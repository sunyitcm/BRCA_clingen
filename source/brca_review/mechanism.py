"""Case-specific ENIGMA v1.2 branches. A predictor value, a clinical LR and a code weight are different observations."""
import re
from ._engine.brcalite.core import GENES
NP={'BRCA1':'NP_009225.1','BRCA2':'NP_000050.3'}
AA=r'(?:Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)'
CATS={'frameshift':'移码','stop':'无义／提前终止','canonical_splice':'经典剪接 ±1/2','missense':'错义','inframe':'框内插入、缺失或 delins','silent':'同义','intronic':'非经典内含子','start_loss':'起始丢失','stop_loss':'终止丢失','unknown':'后果待核','ambiguous':'后果不一致／复合'}
LOF={'frameshift','stop','canonical_splice','start_loss'}
def protein_effect(g,text):
    if not isinstance(text,str) or not text.startswith(NP.get(g,'?')+':p.'):return {}
    s=text.split(':p.',1)[1]
    if s.startswith('(') and s.endswith(')'):s=s[1:-1]
    m=re.fullmatch('('+AA+r')([1-9]\d*)(?:'+AA+r')?fs(?:(?:Ter|\*)([1-9]\d*|\?))?',s)
    if m:
        n=int(m[2]);dist=int(m[3]) if m[3] and m[3].isdigit() else None
        return dict(category='frameshift',start=n,stop_aa=n+dist-1 if dist else None,term='frameshift_variant')
    m=re.fullmatch('('+AA+r')([1-9]\d*)('+AA+r'|Ter|\*|=)',s)
    if m:
        cat='stop' if m[3] in ('Ter','*') else 'silent' if m[3] in (m[1],'=') else 'missense'
        return dict(category=cat,start=int(m[2]),stop_aa=int(m[2]) if cat=='stop' else None,term={'stop':'stop_gained','silent':'synonymous_variant','missense':'missense_variant'}[cat])
    m=re.fullmatch(AA+r'([1-9]\d*)(?:_'+AA+r'([1-9]\d*))?(delins(?:'+AA+')+|del|dup|ins(?:'+AA+')+)',s)
    if m:return dict(category='inframe',start=int(m[1]),end=int(m[2] or m[1]),term='inframe_deletion' if m[3]=='del' else 'inframe_insertion' if m[3]=='dup' or m[3].startswith('ins') else 'protein_altering_variant')
    return {}

def classify(t,g):
    from .criteria import DOM,LIMIT
    if g not in GENES or t.get('transcript_id')!=GENES[g][0] or not str(t.get('hgvsc') or '').startswith(GENES[g][0]+':c.'):return dict(category='unknown',domain='unknown',aa=None,stop_aa=None,end=None)
    p=protein_effect(g,t.get('hgvsp'));c=t['hgvsc'].split(':c.',1)[1];ts=set(t.get('consequence_terms') or [])
    m=re.fullmatch(r'(?:\*?\d+|-\d+)[+-]([1-9]\d*)[ACGT]>[ACGT]',c)
    cats=set()
    for term,cat in [('frameshift_variant','frameshift'),('stop_gained','stop'),('start_lost','start_loss'),('stop_lost','stop_loss'),('missense_variant','missense'),('inframe_deletion','inframe'),('inframe_insertion','inframe'),('synonymous_variant','silent')]:
        if term in ts:cats.add(cat)
    if p:cats.add(p['category'])
    cls='unknown'
    if m:cls='canonical_splice' if int(m[1])<=2 else 'intronic'
    elif ts&{'splice_acceptor_variant','splice_donor_variant'}:cls='canonical_splice'
    elif len(cats)==1:cls=cats.pop()
    elif cats:cls='ambiguous'
    # A generic intronic annotation is insufficient on its own to exclude overlap with the native +/-1,2 splice positions.
    start=t.get('protein_start') if type(t.get('protein_start'))==int else p.get('start')
    end=t.get('protein_end') if type(t.get('protein_end'))==int else p.get('end')
    mismatch=type(p.get('start'))==int and type(t.get('protein_start'))==int and p['start']!=t['protein_start']
    if not isinstance(start,int) or not 1<=start<=LIMIT[g] or mismatch:start=None
    if cls!='inframe':end=start
    if cls=='inframe' and (not isinstance(end,int) or not start or not start<=end<=LIMIT[g]):end=None
    domain='unknown'
    if start and end:
        inside=any(l<=start<=end<=h for l,h in DOM[g]);touch=any(start<=h and end>=l for l,h in DOM[g])
        domain='inside' if inside else 'boundary' if touch else 'outside'
    return dict(category=cls,domain=domain,aa=start,end=end,stop_aa=(p.get('stop_aa') if cls in ('frameshift','stop') and not mismatch else None))

def ptc_hint(g,t,index):
    """Project an explicit coding fsTer stop back only across a simple, fully located coding edit.
    This is an exon lookup, not an independent sequence prediction of the stop or of NMD."""
    m=classify(t,g);stop=m['stop_aa'];cat=m['category']
    if not stop or cat not in ('frameshift','stop'):return
    c=t['hgvsc'].split(':c.',1)[1];delta=0;alt_start=None
    if cat=='frameshift':
        mm=re.fullmatch(r'(\d+)(?:_(\d+))?(delins[ACGT]+|del[ACGT]*|ins[ACGT]+|dup[ACGT]*)',c)
        if not mm:return
        lo,hi=int(mm[1]),int(mm[2] or mm[1]);op=mm[3];n=hi-lo+1
        if n<1:return
        if op.startswith('delins'):ins=len(op[6:]);delta=ins-n;alt_start=lo+ins
        elif op.startswith('del'):
            if op[3:] and len(op[3:])!=n:return
            delta=-n;alt_start=lo
        elif op.startswith('dup'):
            if op[3:] and len(op[3:])!=n:return
            delta=n;alt_start=hi+n+1
        else:
            if hi!=lo+1:return
            delta=len(op[3:]);alt_start=hi+delta
        if delta%3==0:return
    else:
        if not re.fullmatch(r'\d+[ACGT]>[ACGT]',c):return
        pos=int(re.match(r'\d+',c)[0])
        if (pos+2)//3!=stop:return
    altcodon=3*stop-2
    if alt_start is not None and altcodon<alt_start:return
    lo,hi=altcodon-delta,altcodon+2-delta
    spans=[]
    for a in index.get('ST1',[]):
        if a.get('gene')!=g:continue
        try:
            l,h=map(int,a['cdna'])
            if l<=lo<=hi<=h:spans.append(a)
        except (TypeError,ValueError,KeyError):pass
    if len(spans)!=1:return
    span=spans[0];exon=span['exon'].split()[0]
    rows=[a for a in index.get('T4_event',[]) if a['gene']==g and a['event']=='PTC' and re.match(r'E'+re.escape(exon)+r'(?:\D|$)',a['exon'])]
    if len(rows)!=1:return
    from .criteria import getj,R
    locks={v.get('sha256'):v.get('url') for v in getj(R/'rules/attachments.lock.json') if v.get('gene')==g}
    return dict(table4_url=locks.get(rows[0].get('sha256')),st1_url=locks.get(span.get('sha256')),gene=g,transcript=t['transcript_id'],first_aa=m['aa'],stop_aa=stop,reference_codon=[lo,hi],exon=span['exon'],span=span,table4=rows[0],role='predicted_stop_exon_not_adoption')

def interpretation(r,k,put):
    """Return branch-level observations; rows are gated by the caller after identity QC."""
    from .criteria import PROTEIN,dec
    f=r.get('facts') or {};g=r.get('selected_gene') or r.get('gene');t=f.get('transcript') or {};m=classify(t,g);cat=m['category'];on=m['domain'];aa=m['aa']
    state='no_evaluable_data';cmp=None;blocked=False
    src=f.get('transcript_sources') or {};txsrc=src.get('consequence_terms') or src.get('hgvsp') or 'vep'
    if cat!='unknown':
        put('规范转录本后果',t['transcript_id']+' · '+CATS[cat]+(' · '+t['hgvsp'] if t.get('hgvsp') else ''),txsrc,'/facts/transcript')
    if aa:
        loc={'inside':'在规范功能域内','outside':'在列出的规范功能域外','boundary':'跨域边界','unknown':'作用范围未核准'}[on]
        show=str(aa) if not m['end'] or m['end']==aa else f"{aa}–{m['end']}"
        put('首个改变的残基' if cat=='frameshift' else '蛋白位置',f'aa {show}；{loc}'+('；只描述起点' if cat=='frameshift' else ''),src.get('protein_start',src.get('hgvsp','vep')),'/facts/transcript/protein_start' if t.get('protein_start') is not None else '/facts/transcript/hgvsp')
    bx=f.get('ldh_brca') or {};score=dec(bx.get('BayesDel_nsfp33a_noAF')); proteinbranch=(cat in ('missense','inframe') and on=='inside')
    if k in ('PP3','BP4') and score is not None:
        ver=bx.get('version_source') or {};label='BayesDel noAF（'+('蛋白分支' if proteinbranch else '背景，不计分')+'）'
        put(label,str(score)+(' · BRCA Exchange release '+str(ver['releaseName']) if isinstance(ver,dict) and ver.get('releaseName') is not None else ''),'ldh_brca','/facts/ldh_brca/BayesDel_nsfp33a_noAF')
    spok=cat in ('silent','missense','inframe','intronic')
    if k in ('PP3','BP4','BP1') and spok:
        dat=f.get('brca_detail') or {};s=dat.get('splice') or {}
        if s.get('max_delta') is not None:
            put('SpliceAI 四项最大 Δ（背景）',str(s['max_delta'])+'；此来源未核准 unmasked、10 kb 及 REF/ALT 明细，不能单独赋码。','brca_detail','/facts/brca_detail/splice')
    if cat in LOF or cat=='stop_loss':
        blocked=True;state='not_applicable'
        cmp=CATS[cat]+'不属于本规范的 '+k+' 蛋白或剪接适用分支。'+('按 PVS1 流程核转录本、提前终止和 NMD；不把 BayesDel 再作一条证据。' if cat in LOF else '按实际机制复核；不能套用错义阈值。')
    elif cat in ('unknown','ambiguous'):
        state='incomplete';cmp='规范转录本后果未核准，不能仅从坐标、DNA 长度或分数推定蛋白分支。'
    elif k in ('PP3','BP4'):
        if proteinbranch:
            cutoff=PROTEIN[g][0 if k=='PP3' else 1];op='≥' if k=='PP3' else '≤'
            if score is not None:
                hits=score>=cutoff if k=='PP3' else score<=cutoff
                state='branch_threshold_met' if hits else 'branch_threshold_not_met';cmp=f'BayesDel noAF {score}，'+('达到' if hits else '未达到')+f' {g} 蛋白分支阈值 {op}{cutoff}。'
            else:state='incomplete';cmp='蛋白分支适用，但未取得该等位的数值型 BayesDel noAF。'
        else:state='branch_not_applicable';cmp=('蛋白域跨度未核准；' if cat=='inframe' and on in ('unknown','boundary') else '不满足域内错义／框内蛋白分支条件；')
        if k=='PP3':cmp+='剪接是另一分支，须用符合设置的 SpliceAI 最大 Δ≥0.2；同一判据只计一次。'
        elif (cat in ('missense','inframe','silent') and on=='outside'):
            cmp+='域外类型应另核 BP1，而不是把低蛋白分计入 BP4。'
        else:cmp+='无剪接影响须核最大 Δ≤0.1；不能用低蛋白分代替剪接预测。'
    elif k=='BP1':
        if cat in ('missense','inframe','silent') and on=='outside':state='incomplete';cmp='域外条件符合；仍须核准该等位的 SpliceAI 最大 Δ≤0.1 及设置。'
        elif cat in ('missense','inframe','silent') and on=='inside':state='does_not_meet';cmp='在规范功能域内，不满足 BP1 的域外条件。'
        elif cat=='intronic':state='not_applicable';blocked=True;cmp='BP1 不用于此内含子变异；按适用性核 BP4 或 PP3 剪接分支。'
    return state,cmp,blocked
