"""One short result plus its rule per applicable criterion. Source-qualified choices are computed on the server."""
import re,math
from .criteria import GENES,PROTEIN,DOM,dec,protein_domain,verified,catalog,lr_policy
from .erepo import PANEL,UUID
from ._engine.brcalite.core import R
from ._engine.brcalite.net import getj
from ._engine.brcalite.score import frequency
from functools import lru_cache
@lru_cache(maxsize=1)
def locks(): return {x['gene']:x for x in getj(R/'rules/attachments.lock.json') if x['name'].startswith('Specifications_Table9')}
RULE={
 'PVS1':'按变异机制、转录本、外显子及 RNA 证据走本基因 PVS1 分支。','PS1':'比较有合格分类的同效应变异；分别核蛋白与剪接分支。','PS3':'采用 ENIGMA Table 9 对该等位的明确功能实验赋码；仅有论文或 RNA 结果不自动给 PS3。','PS4':'按规范核对病例对照及合格的富集统计。','PM2':'本版 Python 仅对 SNV 自动求值；gnomAD v2.1 非癌外显子和 v3.1 非癌的远交人群（outbred）均未见，区域平均深度均 ≥25×。','PM3':'核范可尼贫血表型、另一等位及相位，按 Table 6 汇总。','PM5':'仅按基因组 PTC 和终止所在外显子查 Table 4；核 PVS1 条件。','PP1':'定量共分离 LR，按方向和大小分档；不能用家族史 LR 代替。','PP3':'按域内蛋白影响或合格剪接预测求值。','PP4':'仅用临床校准的综合 LR；Supporting：2.08 ≤ LR < 4.3。','BA1':'合格非创始人群最大 FAF >0.001，覆盖 ≥20×，排除已知致病创始变异。','BS1':'合格非创始人群最大 FAF >0.0001 可用 Strong，>0.00002 至 ≤0.0001 可用 Supporting；另核覆盖和创始效应。','BS2':'按 Table 8 核无相关隐性表型的共现病例，合计独立先证者。','BS3':'采用 ENIGMA Table 9 对该等位明确推荐的正常蛋白功能证据。','BS4':'用定量共分离 LR 衡量反向分离证据，核表型模拟及替代病因。','BP1':'核本基因域外适用变异和 SpliceAI 最大 Δ≤0.1。','BP4':'核本基因蛋白和剪接无影响条件。','BP5':'只用反向的临床校准综合 LR；跨基因共现本身不够。','BP7':'区分同义、内含子及 RNA 实验分支，按本基因条文赋强度。'}
def wb_lock(r,c):
 g=r['gene'];src=c.get('source') or {};lock=locks().get(g,{})
 code,w=c.get('code'),c.get('strength');n=src.get('row');t=r['facts']['transcript']
 e=src.get('checked_assignment') or {}
 return (code in ('PS3','BS3') and w=='Strong' and src.get('table')=='T9' and src.get('sheet')=='Table9_BRCA12VCEP_specs' and
         type(n) is int and n>1 and src.get('sha256')==lock.get('sha256') and src.get('gene')==g and src.get('hgvsc')==t['hgvsc'].split(':',1)[1] and
         src.get('code')==code and src.get('strength')==w and e=={'gene':g,'hgvsc':src['hgvsc'],'code':code,'strength':w,'cells':f'D{n}:E{n}','workbook_sha256':lock['sha256']})
def source(kind,w,refs,units):return dict(kind=kind,strength=w,refs=refs,units=units)
def expert_ok(r,z):
 e=r.get('erepo') or {};x=z.get('expert') or {}
 return (e.get('state')=='matched' and e.get('spec_version')=='1.2' and e.get('review_panel')==PANEL and e.get('caid')==r.get('caid') and
         e.get('condition')==r.get('condition') and e.get('gene')==r.get('gene') and UUID.fullmatch(str(e.get('uuid',''))) and
         all(re.fullmatch('[0-9a-f]{64}',str(e.get(a,''))) for a in ('summary_sha256','detail_sha256')) and not e.get('alerts') and 
         x.get('state')=='Met' and x.get('detail_state')=='Met' and x.get('strength') in z['allowed_strengths'] and x.get('origin')=='summary_and_detail' and x.get('pointer') and not x.get('issues'))
def expert_brief(z):
 x=z.get('expert') or {};txt=' '.join(str(i.get('value') or '') for i in x.get('text',[]));k=z['code']
 if k=='PM2' and 'v4.1' in txt and re.search(r'absent|未见',txt,re.I):
  if re.search('25',txt) and re.search('90',txt):return '专家记录：gnomAD v4.1 未见；>90% 样本深度 ≥25×。' if re.search(r'(?:>\s*90|超过\s*90)',txt) and re.search(r'≥\s*25|>=\s*25',txt) else '专家记录：gnomAD v4.1 未见；覆盖见原文。'
  return '专家记录：gnomAD v4.1 未见。'
 if k in ('PP4','BP5'):
  from .exchange import text_lr
  values=list(dict.fromkeys(row['value'] for row in text_lr(x.get('text') or [],'combined')))
  if len(values)==1:return '专家记录：临床综合 LR '+values[0]+'。'
 return None
def population(r):
 d=(r.get('facts') or {}).get('population_evidence') or {}
 if not isinstance(d,dict) or d.get('caid')!=r.get('caid') or d.get('gene')!=r.get('gene') or d.get('hgvsc')!=(r.get('facts',{}).get('transcript') or {}).get('hgvsc') or not re.fullmatch('[0-9a-f]{64}',str(d.get('file_sha256',''))):return {},{}
 z=d.get('inputs');z=z if isinstance(z,dict) else {};strata=z.get('strata')
 if not isinstance(strata,list) or not strata:return d,{}
 # Imported counts must point to actual named sources; an absent record or null coverage is never an AC of zero.
 if not all(isinstance(a,dict) and isinstance(a.get('source_ref'),str) and re.fullmatch(r'(?:https://[^\s/]+/\S+|LOCAL:[^\s#]+#\S+)',a['source_ref']) for a in strata):return d,{}
 try:
  ans=frequency(dict(z,allele=r['caid'],allele_verified=True,event=r['facts'].get('event')))
  if ans.get('code')=='PM2':
   need={('gnomAD_v2.1_non_cancer','exome'),('gnomAD_v3.1_non_cancer','genome')}
   def valid(a):
    ac=a.get('AC');an=a.get('AN');dp=a.get('mean_depth')
    return (type(ac) is int and ac==0 and type(an) is int and an>0 and type(dp) in (int,float) and math.isfinite(dp) and dp>=25 and
            a.get('qc_pass') is True and a.get('outbred') is True and a.get('absence_verified_full_subset') is True)
   if {(a.get('dataset'),a.get('assay')) for a in strata}!=need or not all(valid(a) for a in strata):return d,{}
  return d,ans
 except (TypeError,KeyError,ValueError,ArithmeticError):return d,{}
def describe(r,z,identity):
 k=z['code'];g=r.get('selected_gene') or r.get('gene');f=r.get('facts') or {};status=r.get('source_status') or {};t=f.get('transcript') or {};ex=f.get('ldh_exome') or {}
 out=dict(result='本次未取得可按这条规范求值的资料。',standard=RULE[k],expert_note=None,option=None,default_selected=False)
 if not identity:return out|{'result':'等位或转录本未核准。'}
 ep=z.get('expert');sources=[];refroot='/facts/';matched=bool(expert_ok(r,z))
 if ep and ep.get('strength') and ep['strength'] not in z['allowed_strengths']:
  z['expert']=None;ep=None;matched=False
 alt=ep and ep.get('state')=='Not Met'
 if ep:
  out['expert_note']=expert_brief(z) or ('专家明确列为不满足。' if alt else '专家记录：'+str(ep.get('strength') or '未披露强度')+'。')
  out['result']='本地未取得可独立求值的资料。'
 if k in ('PS3','BS3'):
  candidates=[c for c in r.get('candidates',[]) if c.get('code')==k];valid=[c for c in candidates if wb_lock(r,c)]
  opposite=any(c.get('code') in ('PS3','BS3') and c.get('code')!=k for c in r.get('candidates',[]))
  # A missing attestation, even when a paper is mentioned, does not create an assignment.
  if len(valid)==len(candidates)==1 and not opposite:
   s=valid[0]['source'];n=s['row'];out['result']=f'Table 9 第 {n} 行明确推荐 {k} Strong。'
   sources=[source('workbook','Strong',[{'source':'ENIGMA Table 9','locator':s['sheet']+f'!D{n}:E{n}','workbook_sha256':s['sha256'],'url':locks()[g]['url'],'pmids':s.get('pmids_mentioned',[])}],[f"workbook:{s['sha256']}:{n}"])]
  elif candidates:out['result']='附件有相关记录，但未核实唯一且一致的原表赋码。'
  if opposite and candidates:return out|{'status':'conflict','result':'同一等位存在方向相反的功能赋码，暂缓计入。'}
 if k in ('PP3','BP4'):
  lo,hi=PROTEIN[g];bx=f.get('ldh_brca') or {};raw=bx.get('BayesDel_nsfp33a_noAF');n=dec(raw);pos=t.get('protein_start')
  out['standard']=(f'{g} 域内错义或框内改变的蛋白分支：BayesDel noAF ≥{lo}。' if k=='PP3' else f'{g} 域内错义或框内改变：BayesDel noAF ≤{hi} 且 SpliceAI 最大 Δ≤0.1。')
  from .mechanism import classify
  mech=classify(t,g);is_miss=mech['category'] in ('missense','inframe');domain={'inside':True,'outside':False}.get(mech['domain']);pos=mech['aa']
  if n is not None and is_miss:
   cutoff=('≥'+str(lo) if k=='PP3' else '≤'+str(hi));msg=('域内' if domain is True else '域外' if domain is False else '域待核')
   out['result']=f'BayesDel noAF（蛋白）{n}（{g} 标准 {cutoff}）；aa {pos}，{msg}。'
   ver=bx.get('version_source') or {};link=bx.get('linked_allele')
   if (k=='PP3' and ((mech['category']=='missense' and f.get('event') in ('SNV','MNV')) or (mech['category']=='inframe' and f.get('event') in ('DEL','INS','DELINS'))) and domain is True and lo<=n<=1 and link and link in (r.get('representations') or {}).get('GRCh38',[]) and
       isinstance(ver,dict) and (ver.get('releaseName') is not None or ver.get('releaseId') is not None) and not set(t.get('consequence_terms',[]))&{'splice_acceptor_variant','splice_donor_variant'}):
    sources=[source('prediction','Supporting',[{'source':'LDH / BRCA Exchange','field':'BayesDel_nsfp33a_noAF','release':ver,'allele':link,'locator':refroot+'ldh_brca/BayesDel_nsfp33a_noAF'},{'source':(f.get('transcript_sources') or {}).get('protein_start','VEP'),'locator':refroot+'transcript/protein_start','transcript':t['transcript_id'],'position':pos}],['predictor:BayesDel:noAF:'+link])]
   if k=='BP4' and domain and n<=hi:out['result']+=' 另需核准剪接分支。'
  elif k=='PP3':out['result']='本次未取得可用的蛋白预测分值或合格的剪接预测。'
 if k=='BP1':
  d=protein_domain(r);pos=t.get('protein_start')
  if 'missense_variant' in (t.get('consequence_terms') or []) and d is not None:
   out['result']=(f'{g} aa {pos}，位于规范功能域外；本次未核实 SpliceAI 最大 Δ 是否 ≤0.1，不能自动计入。' if not d else f'{g} aa {pos}，在规范功能域内，不符合 BP1 的域外条件。')
 if k in ('PM2','BA1','BS1'):
  pop,qualified=population(r);st=status.get('ldh_exome');pairs=(pop.get('inputs') or {}).get('strata') or []
  if pairs:
   flags=[]
   for ds,label in [('gnomAD_v2.1_non_cancer','v2.1 非癌外显子'),('gnomAD_v3.1_non_cancer','v3.1 非癌基因组')]:
    hits=[q for q in pairs if isinstance(q,dict) and q.get('dataset')==ds]
    flags.append(label+'：'+('；'.join(f"AC/AN {q.get('AC','?')}/{q.get('AN','?')}、覆盖 {q.get('mean_depth','?')}×" for q in hits[:2]) if hits else '未取得'))
   out['result']='本地队列：'+'；'.join(flags)+'。'
  elif ex:out['result']='LDH：gnomAD v4.1 外显子 AC/AN '+str(ex.get('AC') if ex.get('AC') is not None else '未返回')+'/'+str(ex.get('AN') if ex.get('AN') is not None else '未返回')+'；未取得规范双队列的覆盖。'
  elif st=='no_link_listed':out['result']='LDH 未列该等位的 gnomAD v4.1 外显子记录；未取得规范双队列的计数及覆盖。'
  elif st in ('http_error','transport_error','unavailable','invalid_json','unresolved_links'):out['result']='LDH 人群频率请求失败或关联未核准；规范双队列资料未取得。'
  else:out['result']='本次没有可用的规范队列计数及覆盖。'
  if qualified.get('code')==k and qualified.get('strength') in z['allowed_strengths']:
   sources=[source('population',qualified['strength'],[{'source':'本地已对账人群数据','locator':refroot+'population_evidence','file_sha256':pop['file_sha256'],'strata_refs':sorted({a['source_ref'] for a in pairs})}],['population:'+r['caid']])]
 if not sources and z['candidates'] and k not in ('PS3','BS3'):out['result']='附件有位点或外显子线索，尚无足以独立采用的判定。'
 if not sources and ep and matched and k not in ('PM2','BA1','BS1','PP3','BP4') and out['expert_note']:out['result']=out['expert_note']
 conflict=bool(sources and ep and (ep['state']!='Met' or ep.get('strength')!=sources[0]['strength'] or ep.get('issues')))
 if conflict:return out|{'status':'conflict','result':out['result']+' 与专家记录不一致，暂缓计入。'}
 if sources:return out|{'option':sources[0],'default_selected':True,'status':'source_ready'}
 if matched:
  e=r['erepo'];ref={'source':'ENIGMA eRepo 个案引用','url':e['url'],'uuid':e['uuid'],'spec_version':e['spec_version'],'record_version':e['record_version'],'detail_sha256':e['detail_sha256'],'locator':ep['pointer'],'pmids':ep.get('pmids',[])}
  return out|{'option':source('expert_citation',ep['strength'],[ref],['expert:'+e['uuid']+':'+k]),'status':'citation_available'}
 if ep and not alt and not matched:
  out['expert_note']=(out.get('expert_note') or '')+' 个案规范版或逐码资料未核全，仅展示。'
 return out|{'status':'expert_not_met' if alt else 'missing'}
