from pathlib import Path
import gzip,json,re,urllib.parse as u
from functools import lru_cache
from . import __version__
from .net import getj,sha
R=Path(__file__).resolve().parents[1]
GENES={'BRCA1':('NM_007294.4','MONDO:0700268','GN092'),'BRCA2':('NM_000059.4','MONDO:0700269','GN097')}
NC={'GRCh38':{'17':'NC_000017.11','13':'NC_000013.11'},'GRCh37':{'17':'NC_000017.10','13':'NC_000013.10'}}
def parse(s,build):
 if build not in NC:raise ValueError('必须明确 GRCh37 或 GRCh38')
 base={'original':s,'assembly':build}
 if '\t' in s:
  fields=s.rstrip('\r\n').split('\t')
  if len(fields)<5:raise ValueError('VCF 至少需要前五列')
  return parse(f'{fields[0]}-{fields[1]}-{fields[3]}-{fields[4]}',build)
 if ':' in s:
  annotated=re.fullmatch(r'(NM_\d+\.\d+)\((BRCA[12])\)(:c\..+)',s)
  if annotated:
   if GENES[annotated[2]][0]!=annotated[1]:raise ValueError('注释基因与规范转录本不符')
   s=annotated[1]+annotated[3]
  m=re.fullmatch(r'(NM_\d+\.\d+):c\.([*+_A-Za-z0-9>\-]+)',s)
  if not m or m[1] not in {v[0] for v in GENES.values()}:raise ValueError('v0.1 仅接收规范指定的版本化 NM:c.；p.、其他版本、复合单倍型进入人工消歧')
  if not re.fullmatch(r'(?:-?\d+|\*\d+)(?:[+-]\d+)?(?:_(?:-?\d+|\*\d+)(?:[+-]\d+)?)?(?:[ACGT]>[ACGT]|delins[ACGT]+|del[ACGT]*|dup[ACGT]*|ins[ACGT]+)',m[2]) or re.search(r'([ACGT])>\1$',s):raise ValueError('不支持或无变化的 c. 表达')
  return dict(base,kind='hgvsc',hgvs=s,gene=next(g for g in GENES if GENES[g][0]==m[1]))
 m=re.fullmatch(r'(?:chr)?(13|17)-(\d+)-([ACGT]+)-([ACGT]+)',s)
 if not m or ',' in s:raise ValueError('格式: chr17-POS-REF-ALT；逐 ALT 输入，未知碱基、SV、复合或其他染色体不支持')
 chrom,pos,ref,alt=m[1],int(m[2]),m[3],m[4]
 if pos<1 or ref==alt or max(len(ref),len(alt))>100:raise ValueError('坐标、REF/ALT或长度超界')
 a,c=ref,alt;start=pos
 while a and c and a[0]==c[0]:a,c,start=a[1:],c[1:],start+1
 while a and c and a[-1]==c[-1]:a,c=a[:-1],c[:-1]
 end=start+len(a)-1; span=str(start) if len(a)==1 else f'{start}_{end}'
 if a and c:expr=f'{start}{a}>{c}' if len(a)==len(c)==1 else f'{span}delins{c}'
 elif a:expr=f'{span}del'
 else:
  if start<=1:raise ValueError('边界插入需人工复核')
  expr=f'{start-1}_{start}ins{c}'
 return dict(base,kind='vcf',chrom=chrom,pos=pos,ref=ref,alt=alt,hgvs=f'{NC[build][chrom]}:g.{expr}',interbase=[start-1,end,a,c])
def requests_for(v):
 b=v['assembly'];host='https://'+('grch37.' if b=='GRCh37' else '')+'rest.ensembl.org';qs='hgvs=1&merged=1&transcript_version=1&numbers=1'+('&mane=1' if b=='GRCh38' else '')+'&variant_class=1'
 z={'car':('https://reg.clinicalgenome.org/allele?'+u.urlencode({'hgvs':v['hgvs']}),None),'info':(host+'/info/data?content-type=application/json',None)}
 if v['kind']=='hgvsc':z['vep']=(host+'/vep/human/hgvs/'+u.quote(v['hgvs'],safe='')+'?'+qs,None)
 else:
  chrom,p,ref,alt=v['chrom'],v['pos'],v['ref'],v['alt']
  z['ref']=(host+f'/sequence/region/human/{chrom}:{p}..{p+len(ref)-1}:1?content-type=application/json',None)
  z['vep']=(host+'/vep/human/region?'+qs,{'variants':[f'{chrom} {p} . {ref} {alt} . . .']})
 return z
def validcar(v,c,reference=None):
 if not isinstance(c,dict) or not re.fullmatch(r'CA\d+',str(c.get('@id','')).split('/')[-1]):return None
 names={g:{h for row in c.get('transcriptAlleles',[]) if row.get('geneSymbol')==g for h in row.get('hgvs',[]) if h.startswith(t[0]+':c.')} for g,t in GENES.items()}
 used=[g for g in names if names[g]]
 if len(used)!=1:return None
 gene=used[0]
 if len(names[gene])!=1:return None
 if v['kind']=='hgvsc':
  if v['gene']!=gene or v['hgvs'] not in names[gene]:return None
 else:
  expected=NC[v['assembly']][v['chrom']]
  ar=[a for a in c.get('genomicAlleles',[]) if a.get('referenceGenome')==v['assembly'] and any(h.startswith(expected+':g.') for h in a.get('hgvs',[]))]
  exact=any(len(a.get('coordinates',[]))==1 and [q.get('start'),q.get('end'),q.get('referenceAllele'),q.get('allele')]==v['interbase'] for a in ar for q in a['coordinates'])
  if not exact:
   if reference is None or len(ar)!=1 or len(ar[0].get('coordinates',[]))!=1:return None
   q=ar[0]['coordinates'][0]
   other=[q.get(k) for k in ('start','end','referenceAllele','allele')]
   if not all(type(n) is int for n in other[:2]) or not all(isinstance(a,str) and re.fullmatch('[ACGT]*',a) for a in other[2:]):return None
   if not reference.equivalent(v['chrom'],v['interbase'],other):return None
 return (c['@id'].split('/')[-1],gene,names[gene])
def ref_matches(v,ref):
 end=v['pos']+len(v['ref'])-1
 return isinstance(ref,dict) and ref.get('seq','').upper()==v['ref'] and ref.get('id')==f"chromosome:{v['assembly']}:{v['chrom']}:{v['pos']}:{end}:1"
def vetted(v,c,vep,ref,reference=None):
 a=validcar(v,c,reference)
 if not a:return None,'CAR 无唯一精确核酸/规范转录本对应'
 if v['kind']=='vcf' and not ref_matches(v,ref):return None,'REF、坐标或序列组装未核对'
 if not isinstance(vep,list) or len(vep)!=1:return None,'VEP 响应未逐输入唯一对账'
 res=vep[0]
 expected=requests_for(v)['vep'][1]; expected=expected['variants'][0] if expected else v['hgvs']
 if res.get('assembly_name')!=v['assembly'] or res.get('input')!=expected:return None,'VEP 组装或原输入不符'
 ca,gene,hh=a
 ts=[t for t in res.get('transcript_consequences',[]) if t.get('gene_symbol')==gene and t.get('transcript_id')==GENES[gene][0] and t.get('hgvsc') in hh]
 if len(ts)!=1:return None,'规范转录本／c. 不一致或多解'
 return (ca,gene,ts[0]),None
@lru_cache(maxsize=1)
def row_index():
 from ... import ruleprep
 return ruleprep.rows()
def lookup(g,h):
 db=row_index();c=h.split(':',1)[-1];found=[]
 for kind in ['T9','T4_splice']:
  for x in db[kind]:
   if x['gene']==g and x['hgvsc']==c:
    source=dict(x,table=kind)
    if kind=='T4_splice':
     source['rna_references']=[r for r in db['ST16'] if r['gene']==g and r['site']==c[:-1] and c[-1] in re.findall(r'\b[ACGT]\b',r['alts'])]
     source['pmids_mentioned']=sorted({p for r in source['rna_references'] for p in r['pmids']})
    found.append(source)
 return found

def ptc_context(g,t):
 terms=set(t.get('consequence_terms',[]))
 if not terms.intersection({'stop_gained','frameshift_variant'}):return []
 match=re.fullmatch(r'c\.(\d+)(?:_(\d+))?(?:[ACGT]>[ACGT]|delins[ACGT]+|del[ACGT]*|dup[ACGT]*|ins[ACGT]+)',t['hgvsc'].partition(':')[2])
 if not match:return []
 start,end=int(match[1]),int(match[2] or match[1])
 exon=t.get('exon','').split('/')[0]
 db=row_index()
 spans=[r for r in db['ST1'] if r['gene']==g and r['exon'].split()[0]==exon and int(r['cdna'][0])<=start<=end<=int(r['cdna'][1])]
 if len(spans)!=1:return []
 events=[r for r in db['T4_event'] if r['gene']==g and re.match(r'E'+re.escape(exon)+r'(?:\D|$)',r['exon']) and r['event']=='PTC']
 return [dict(r,score_counted=False,role='DNA event exon only; verify stop exon separately') for r in spans+events]
def review_prompts(r):
 """Document a manual next step separately from workbook candidates and accepted evidence."""
 f=r.get('facts') or {};t=f.get('transcript') or {};p=t.get('protein_start')
 if (r.get('specification')=='GN092 v1.2' and r.get('status') in ('review_required','draft_reviewed')
     and r.get('gene')=='BRCA1' and f.get('event')=='SNV' and (t.get('hgvsc') or '').startswith('NM_007294.4:')
     and 'missense_variant' in (t.get('consequence_terms') or []) and type(p) is int and 1<=p<=1863
     and not any(lo<=p<=hi for lo,hi in ((2,101),(1391,1424),(1650,1857)))):
  return [{'criterion':'BP1_Strong','state':'review_prompt','score_counted':False,'protein_position':p,
           'specification':'GN092 v1.2 / Figure 1A / Appendix J','reference_url':'https://cspec.genome.network/cspec/ui/svi/doc/GN092?version=1.2',
           'check':'confirm splice prediction and supporting inputs; max delta <= 0.1',
           'spliceai':{'distance':10000,'mask':'unmasked','per_allele_scores':'REF and ALT','max_delta':None,'source':'not_queried_by_pipeline'}}]
 return []
def write_report(folder,v,pieces,links=None,reference=None):
 from .net import writej,utc
 from .report import html_report
 val,err=vetted(v,pieces.get('car'),pieces.get('vep'),pieces.get('ref'),reference)
 identity=validcar(v,pieces.get('car'),reference);gene=identity[1] if identity else v.get('gene')
 out={'schema':'brca-lite.report.v1','engine':'0.1.0','generated_utc':utc(),'status':'review_required' if val else 'blocked_identity_or_annotation','input':v,'gene':gene,'condition':GENES[gene][1] if gene else None,'inheritance':'AD','domain':'germline cancer predisposition','specification':(GENES[gene][2]+' v1.2') if gene else None,'rules_bundle_sha256':sha((R/'rules/bundle.json').read_bytes()),'caid':identity[0] if identity else None,'classification_local':None,'points_accepted':None,'approved_on':None,'published_on':None,'curator':None,'expert_assertions':[],'met':[],'not_met':[],'not_evaluated':[],'candidates':[],'facts':{},'source_files':links or {},'limitations':[]}
 out.update(engine=__version__,run_id=folder.name,classification_draft=None,source_of_sample='unknown',genetic_origin='unknown',somatic_oncogenicity='not_assessed',therapy='not_assessed')
 if pieces.get('ref'):out['facts']['reference_check']=pieces['ref']
 out['engine_manifest_sha256']=sha((R/'rules/engine-0.2.json').read_bytes())
 if err:out['limitations'].append(err)
 if gene:out['not_evaluated']=[k for k,x in getj(R/'rules/catalog.json')[gene].items() if x['enabled']];out['not_applicable_by_spec']=[k for k,x in getj(R/'rules/catalog.json')[gene].items() if not x['enabled']]
 if identity:
  out['representations']={b:[h for a in pieces['car'].get('genomicAlleles',[]) if a.get('referenceGenome')==b for h in a.get('hgvs',[]) if h.startswith(tuple(x+':' for x in NC[b].values()))] for b in ('GRCh37','GRCh38')}
 if val:
  ca,g,t=val;out['facts']['transcript']={k:t.get(k) for k in ('gene_symbol','transcript_id','hgvsc','hgvsp','exon','consequence_terms','protein_start')}
  coords=[q for a in pieces['car'].get('genomicAlleles',[]) if a.get('referenceGenome')==v['assembly'] and any(h.startswith(tuple(n+':g.' for n in NC[v['assembly']].values())) for h in a.get('hgvs',[])) for q in a.get('coordinates',[])]
  if len(coords)==1:
   ref,alt=coords[0].get('referenceAllele'),coords[0].get('allele')
   if isinstance(ref,str) and isinstance(alt,str):
    out['facts']['event']='INS' if not ref else 'DEL' if not alt else 'SNV' if len(ref)==len(alt)==1 else 'MNV' if len(ref)==len(alt) else 'DELINS'
    if v['kind']=='vcf':
     car_interval=[coords[0]['start'],coords[0]['end'],ref,alt]
     out['facts']['genomic_reconciliation']={'input_interbase':v['interbase'],'car_interbase':car_interval,
         'method':'exact_coordinates' if car_interval==v['interbase'] else 'equal_local_alternate_haplotypes'}
  out['facts']['exon_workbook_context']=ptc_context(g,t)
  for x in lookup(g,t['hgvsc']):
   assignment=x.get('code') or x.get('assigned') or ''
   if not re.match(r'^(PS3|BS3|PVS1|BP7)(?:$|\b|_)',assignment) or 'N/A' in assignment:
    out['facts'].setdefault('workbook_rows_without_assignment',[]).append(x)
    continue
   out['candidates'].append(dict(state='requires_review',source=x,code=x.get('code') or x.get('assigned'),strength=x.get('strength'),score_counted=False))
 if pieces.get('info'):out['facts']['ensembl_metadata']=pieces['info']
 if pieces.get('ldh') and identity:
  d=pieces['ldh'].get('data',{})
  if d.get('entId')==identity[0]:out['facts']['ldh_associations_only']=d.get('ld',{})
  else:out['limitations'].append('LDH 等位不一致')
 out['limitations'].append('外部断言、文献、等位频率和证据独立性未完成全量审核。未给本等位签发五分类。')
 out['manual_review_prompts']=review_prompts(out)
 writej(folder/'report.json',out);html_report(folder/'report.html',out);return out
