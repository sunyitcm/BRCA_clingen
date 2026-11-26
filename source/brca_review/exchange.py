"""Only an allele-matched, version-attributed BRCA Exchange row may supply additional reading material."""
import re
from decimal import Decimal,InvalidOperation
from ._engine.brcalite import core
MISSING=(None,'','-','.','NA','N/A','null')
def numeric(x):
    if isinstance(x,bool) or x in MISSING:return None
    try:n=Decimal(str(x).strip())
    except (TypeError,InvalidOperation,ValueError):return None
    return n if n.is_finite() else None

def detail(doc,linked,release,car,g,c_name):
    rows=doc.get('data') if isinstance(doc,dict) else None
    if not isinstance(rows,list) or not all(isinstance(r,dict) for r in rows):return None
    chrom='17' if g=='BRCA1' else '13'
    hs={b:{h for a in car.get('genomicAlleles',[]) if a.get('referenceGenome')==b for h in a.get('hgvs',[]) if h.startswith(core.NC[b][chrom]+':g.')} for b in core.NC}
    sch={'GRCh38':'38','GRCh37':'37'}
    def brief(i,r):
        rel=r.get('Data_Release') or {}
        if not isinstance(rel,dict):return None
        name,idx=rel.get('name'),rel.get('id')
        if type(name)!=int or type(idx)!=int:return None
        qa=[]
        if r.get('Gene_Symbol')!=g or r.get('HGVS_cDNA')!=c_name or r.get('Genomic_HGVS_38') not in hs['GRCh38']:qa.append('gene_c_or_g38')
        for b,s in sch.items():
            if not hs[b] or all(r.get(k) in MISSING for k in ('Genomic_HGVS_'+s,'Genomic_Coordinate_hg'+s)):qa.append('no_verified_'+b);continue
            for key in ('Genomic_HGVS_'+s,'Genomic_Coordinate_hg'+s):
                val=r.get(key)
                if val in MISSING:continue
                if key.startswith('Genomic_Coordinate'):
                    m=re.fullmatch(r'chr(13|17):g\.(\d+):([ACGT]+)>([ACGT]+)',str(val))
                    if not m:qa.append(key);continue
                    try:val=core.parse('-'.join(m.groups()),b)['hgvs']
                    except ValueError:qa.append(key);continue
                if hs[b] and val not in hs[b]:qa.append(key)
        return dict(row=i,variant_id=r.get('id'),release=name,release_id=idx,release_date=rel.get('date'),issues=qa)
    indexed=[(v,r) for i,r in enumerate(rows) if (v:=brief(i,r))]
    pinned=[(v,r) for v,r in indexed if v['variant_id']==linked and v['release']==release.get('releaseName') and v['release_id']==release.get('releaseId')]
    if len(pinned)!=1:return None
    base,rp=pinned[0]
    if 'gene_c_or_g38' in base['issues']:return None
    major=max(v['release'] for v,r in indexed);latest=[(v,r) for v,r in indexed if v['release']==major]
    # A failed cross-build comparison never appoints this historical row as a clinical data source.
    selected=None
    if not base['issues']:selected=(base,rp)
    elif len(latest)==1 and not latest[0][0]['issues']:selected=latest[0]
    out=dict(ldh_linked_row=base,latest_returned=latest[0][0] if len(latest)==1 else None,method='matched_LDH_release' if selected and selected[0]==base else 'separate_matching_release' if selected else 'none',source_row=None,measures={},splice={},history=[])
    seen=set()
    for v,r in sorted(indexed,key=lambda a:-a[0]['release']):
        n=numeric(r.get('BayesDel_nsfp33a_noAF'))
        if not v['issues'] and n is not None and str(n) not in seen and len(out['history'])<8:
            seen.add(str(n));out['history'].append(dict(release=v['release'],value=str(n),row=v['row'],field='BayesDel_nsfp33a_noAF'))
    if not selected:return out
    obs,row=selected;out['source_row']=obs
    for field,kind in [('Segregation_LR_exLOVD','segregation'),('Sum_family_LR_exLOVD','family_summary_unspecified'),('Posterior_probability_exLOVD','historical_posterior'),('Combined_prior_probablility_exLOVD','historical_prior')]:
        v=row.get(field);n=numeric(v)
        out['measures'][kind]=dict(field=field,raw=v if isinstance(v,(str,int,float)) else None,value=str(n) if n is not None and n>=0 else None,state='not_in_record' if field not in row else 'not_reported' if v in MISSING else 'reported' if n is not None and n>=0 else 'not_numeric',raw_pointer=f"/data/{obs['row']}/{field}")
    keys=('DS_AG_spliceAI','DS_AL_spliceAI','DS_DG_spliceAI','DS_DL_spliceAI')
    if all(numeric(row.get(k)) is not None and 0<=numeric(row[k])<=1 for k in keys):
        out['splice']=dict(max_delta=str(max(numeric(row[k]) for k in keys)),components={k:row[k] for k in keys},mask='not_disclosed',distance=None,reference_alt_scores=None,counts_as_evidence=False,raw_pointer="/data/"+str(obs['row']))
    return out

# Typed exact numerical observations only. Do not parse inequalities, intervals, grading cutoffs or imputed defaults.
NUM=r'(?:[0-9]+(?:\.[0-9]+)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?'
TEXT_TYPE={'segregation':r'(?:(?:co[- ]?)?segregation\s+LR|共分离\s*LR)','combined':r'(?:combined\s+LR|综合\s*LR)','family_history':r'(?:family\s+history\s+LR|家族史\s*LR)'}
def text_lr(lines,typ):
    out=[]
    for line in lines:
        s=line.get('value') or '';pat=TEXT_TYPE[typ]+r'\s*(?:(?:for this variant\s+)?(?:of\s+|=\s*|was\s+|is\s+|(?:均)?为\s*)|[:：]\s*)('+NUM+r')(?![\d.eE%]|\s*[-–~]\s*\d)'
        # A trailing period is punctuation unless followed by a digit.
        pat=pat.replace(r'(?![\d.eE%]|',r'(?![\deE%]|\.\d|')
        if typ=='family_history':
            s=re.sub(r'家族史\s*LR\s*与\s*综合\s*LR\s*均为', '家族史 LR 为', s)
        for m in re.finditer(pat,s,re.I):
            val=m[1]
            if re.search(r'(?:threshold|cutoff|minimum|required|example|阈值|举例)\s*[:：=]?\s*$',s[max(0,m.start()-25):m.start()],re.I) or re.match(r'\s*(?:was|is|as)?\s*(?:the |a )?(?:threshold|cutoff|is required|阈值)',s[m.end():],re.I):continue
            if re.search(r'(?:non[- ]|unverified |no measured |没有|未见)$',s[max(0,m.start()-15):m.start()],re.I):continue
            if numeric(val) is None or numeric(val)<=0:continue
            out.append(dict(value=val,pointer=line.get('pointer'),type=typ))
    return out
