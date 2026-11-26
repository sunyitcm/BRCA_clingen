"""Append-only local human decisions. The frozen scorer is the only classification engine."""
import datetime as dt, hashlib, json, os, re, tempfile, shutil
from pathlib import Path
from copy import deepcopy
from . import __version__
from .criteria import GENES,PROTEIN,VSL,dec,ledger,lr_policy,protein_domain,verified
from ._engine.brcalite.core import R
from ._engine.brcalite.net import writej,getj,utc
from ._engine.brcalite.score import adjudicate,frequency
class ReviewChanged(ValueError):pass
def text(v,maxlen):
 if not isinstance(v,str) or not 1<=len(v.strip())<=maxlen:raise ValueError('审核人、理由、来源或独立证据编号缺失／过长')
 return v.strip()
def canon(s):
 v=s.strip().casefold().rstrip('/')
 m=re.fullmatch(r'(?:pmid\s*[:：]\s*|https?://pubmed\.ncbi\.nlm\.nih\.gov/)([1-9]\d{0,9})',v)
 if m:return 'pmid:'+m[1]
 if v.startswith('https://doi.org/'):return 'doi:'+v.partition('doi.org/')[2]
 return ' '.join(v.split())
def ids(v):
 if not isinstance(v,list) or not 1<=len(v)<=20:raise ValueError('每条需 1–20 个可追溯来源／独立证据单元')
 return [text(s,350) for s in v]
def primary(s,candidate_ids):
 v=canon(s)
 if v.startswith('wb:'):return s in candidate_ids
 if re.search(r'(erepo|clinvar|cspec\.genome\.network|/(?:vcv|rcv|scv)\d+)',v):return False
 return bool(re.fullmatch(r'pmid:[1-9]\d{0,9}|doi:10\.\d{4,9}/\S+|local:[^\s#]+#[^\s]+|https://[^\s/]+/[^\s]+',v))
def splice_ok(v):
 n=dec(v.get('splice_max'))
 if n is None or not 0<=n<=1 or v.get('unmasked_10kb_ref_alt') is not True:raise ValueError('剪接分支需复核 10 kb、unmasked、REF/ALT 和最大 Δ（0–1）')
 return n
def special(r,k,w,o,srcs):
 if k in ('PM2','BS1','BA1'):
  f=dict(o,allele=r['caid'],allele_verified=True,event=r['facts'].get('event'))
  if not isinstance(f.get('strata'),list) or not all(isinstance(z,dict) and z.get('source_ref') in srcs for z in f['strata']):raise ValueError(k+'：需提供对应队列、QC、覆盖、计数和匹配的来源')
  assessed=frequency(f)
  if (assessed.get('code'),assessed.get('strength'))!=(k,w):raise ValueError(k+'：提供的人群数据未通过该码与强度的冻结条件；不能使用外部断言或查无链接代替')
 if k in VSL:
  lr=dec(o.get('lr'))
  if lr is None or lr<=0 or o.get('calibrated') is not True:raise ValueError(k+'：需给出校准 LR，并确认原始资料')
  policy=lr_policy();cfg=policy['criteria'][k]
  if o.get('lr_type')!=cfg['lr_type']:raise ValueError(k+'：家族史等临床综合 LR 与共分离 LR 不可互换')
  limits=policy['directional_thresholds'][cfg['direction']]
  allowed=next((grade for grade in ('Very Strong','Strong','Moderate','Supporting') if (lr>=dec(limits[grade]) if cfg['direction']=='for' else lr<=dec(limits[grade]))),None)
  if w!=allowed:raise ValueError(k+'：所选强度与 LR 的方向或档位不符，不允许降档代算')
  if w=='Very Strong':
   src=o.get('lr_source')
   if not isinstance(src,str) or src not in srcs or not primary(src,set()) or o.get('independence_checked') is not True:raise ValueError(k+'：该档需定位 LR 原始来源，并核查本批独立单元')
   if k in ('PP4','BP5') and o.get('components_checked') is not True:raise ValueError(k+'：请核实综合 LR 的数据类型、组成和权重')
   if k=='BP5' and o.get('not_cross_gene_coobservation') is not True:raise ValueError('BP5：不能仅凭不同癌易感基因的致病变异共现')
   if k in ('PP1','BS4') and o.get('pedigree_context') not in ('single','isolated','multiple'):raise ValueError(k+'：请核对家系背景')
   if k=='BS4' and o.get('phenocopy_reviewed') is not True:raise ValueError('BS4：请核查表型模拟及家系中其他致病变异')
   if k=='PP1':
    es=o.get('effect_source')
    if o.get('effect_kind') not in ('predicted','experimental') or not isinstance(es,str) or es not in srcs or not primary(es,set()):raise ValueError('PP1：该档须另核蛋白或剪接影响预测／实验的原始依据')
    if o.get('pedigree_context') in ('single','isolated') and o.get('alternative_cause_checked') is not True:raise ValueError('PP1：单一家系或隔离人群须评估其他致病变异')
 if k in ('PP3','BP4','BP1'):
  t=r['facts']['transcript'];terms=set(t.get('consequence_terms') or []);miss='missense_variant' in terms;syn='synonymous_variant' in terms;intr='intron_variant' in terms and not terms.intersection({'splice_acceptor_variant','splice_donor_variant'})
  d=protein_domain(r);lo,hi=PROTEIN[r['gene']]
  if k=='PP3' and o.get('branch')=='protein':
   s=dec(o.get('bayes_noaf'))
   if not miss or d is not True or s is None or s<lo or s>1 or o.get('predictor_verified') is not True:raise ValueError('PP3 蛋白分支：需规范转录本域内错义、达本基因阈值的 BayesDel noAF 与资源设置')
  elif k=='PP3' and o.get('branch')=='splice':
   if not (miss or syn or intr) or splice_ok(o)<dec('.2'):raise ValueError('PP3 剪接分支缺适用前提或未达到 0.2')
  elif k=='BP1':
   if not (miss or syn) or d is not False or splice_ok(o)>dec('.1'):raise ValueError('BP1：需域外适用变异及剪接 Δ≤0.1')
  elif k=='BP4':
   s=dec(o.get('bayes_noaf'))
   if not (intr or (d is True and (syn or miss and s is not None and -1<=s<=hi and o.get('predictor_verified') is True))) or splice_ok(o)>dec('.1'):raise ValueError('BP4：需满足该基因适用分支的蛋白及剪接条件')
  else:raise ValueError('需选择可核对的蛋白或剪接分支；其他机制保留人工复核')
 if k=='PM5':
  if not set(r['facts']['transcript'].get('consequence_terms') or []).intersection({'stop_gained','frameshift_variant'}):raise ValueError('PM5_PTC：本位点不属于基因组 PTC 分支')
def policy_hash():
 root=Path(__file__).resolve().parent
 h=hashlib.sha256()
 for n in ('adoption.py','criteria.py','mechanism.py','lineage.py','exchange.py','briefing.py','routes.py','observations.py','selection.py','guidance.py','assets/criteria_zh_v12.json','lr_v12.json','erepo.py','_engine/brcalite/score.py','_engine/brcalite/workbook.py','_engine/rules/catalog.json','_engine/rules/review-20261124.json'):
  h.update(n.encode()+b'\0'+hashlib.sha256((root/n).read_bytes()).digest())
 return h.hexdigest()
def root_for(st,sid):return st.dir/'reviews'/sid
def seal(row):
 return hashlib.sha256(json.dumps({k:v for k,v in row.items() if k!='content_sha256'},sort_keys=True,ensure_ascii=False,allow_nan=False,separators=(',',':')).encode()).hexdigest()
def latest(st,sid):
 base=root_for(st,sid);hits=sorted(x for x in base.glob('r[0-9][0-9][0-9][0-9][0-9][0-9]') if (x/'review.json').is_file())
 prev=None
 for i,path in enumerate(hits,1):
  rev=getj(path/'review.json')
  if (rev.get('content_sha256')!=seal(rev) or rev.get('snapshot')!=sid or rev.get('revision')!=f'r{i:06d}' or rev.get('parent_sha256')!=(prev['content_sha256'] if prev else None)):raise ValueError('本地审核版本链校验失败')
  prev=rev
 return prev
def effective(report,rev):
 if rev and rev.get('schema')=='brca.selection.v2':
  from .selection import present
  return present(report,rev)
 r=deepcopy(report);r.pop('not_applicable_by_spec',None)
 if not rev:
  r['criteria_table']=ledger(r);return r
 rules=rev['decisions'];r['criteria_table']=ledger(r,rules);r.update(review_desk_version=rev['desk_version'],review_policy_sha256=rev['policy_sha256'],is_training=rev['is_training'],met=[x['code'] for x in rules if x['state']=='met'],not_met=[x['code'] for x in rules if x['state']=='not_met'],locally_not_applicable=[x['code'] for x in rules if x['state']=='not_applicable'],not_evaluated=[x['code'] for x in r['criteria_table']['rows'] if x['rule_enabled'] and x['decision'] is None],classification_draft=rev['result']['draft_class'],points_accepted=rev['result']['points'],review_revision=rev['revision'],reviewed_utc=rev['created_utc'],review_result=rev['result'],classification_local=None)
 return r
def evaluate(s,p,old):
 r=s['report'];gene=r.get('selected_gene') or r.get('gene')
 if not verified(r,gene) or r.get('rules_bundle_sha256')!=hashlib.sha256((R/'rules/bundle.json').read_bytes()).hexdigest():raise ValueError('等位、转录本或规则锁未核实，不能采纳')
 if r.get('met') or r.get('not_met') or r.get('points_accepted') is not None:raise ValueError('旧版已有本地计分；请保留原审核，不覆盖旧结果')
 if s['report_sha256']!=p.get('base_hash') or p.get('gene')!=gene or p.get('caid')!=r['caid'] or p.get('assembly')!=s['assembly']:raise ValueError('请求与快照的基因、组装、CA 或原报告哈希不符')
 last=old['revision'] if old else None
 if last!=p.get('previous_revision'):raise ReviewChanged('审核已有新版本，请刷新后查看，旧页面不能覆盖')
 who=text(p.get('reviewer'),80);now=utc();date=now[:10]
 inventory=p.get('decisions')
 if not isinstance(inventory,list) or len(inventory)>28:raise ValueError('判据清单不合法')
 rows={z['code']:z for z in ledger(r)['rows']};clean=[];seen=set();oldrows={z['code']:z for z in (old or {}).get('decisions',[])}
 flags={}
 for f in ('review_complete','independence_reviewed','significant_conflict'):
  val=p.get(f,False)
  if type(val)!=bool:raise ValueError('审核标志不是布尔值')
  flags[f]=val
 for x in inventory:
  if not isinstance(x,dict):raise ValueError('逐项记录格式错误')
  k=x.get('code');st=x.get('state')
  if not isinstance(k,str) or k in seen or not rows.get(k,{}).get('rule_enabled') or st not in ('met','not_met','not_applicable'):raise ValueError('判据停用、重复或审核状态不正确')
  seen.add(k);row=rows[k]
  if x.get('confirmed') is not True:raise ValueError(k+'：须先确认已核查原始证据和适用条件')
  reason=text(x.get('reason'),1800);refs=ids(x.get('source_refs'))
  ck={'code':k,'state':st,'confirmed':True,'reason':reason,'source_refs':refs}
  if st!='met' and not any(primary(v,{a['ref'] for a in row['candidates']}) or v.startswith(row['rule_url']) for v in refs):
   raise ValueError(k+'：不满足／不适用也需要可定位的来源或本基因规范条文')
  if st=='met':
   w=x.get('strength');options=row['allowed_strengths']
   if w not in options:raise ValueError(k+'：强度不在当前规则锁；保留专家原强度但不在本地计分')
   workbook={a['ref'] for a in row['candidates'] if a['strength']==w or not a['strength']}
   if any(ref.startswith('WB:') and ref not in workbook for ref in refs) or not any(primary(ref,workbook) for ref in refs):raise ValueError(k+'：需独立原始证据；整条专家断言或其他位点／判据的附件不能代替')
   units=[canon(z) for z in ids(x.get('independent_units'))];o=x.get('inputs') or {}
   if not isinstance(o,dict):raise ValueError('附加字段须为 JSON 对象')
   special(r,k,w,o,refs);ck.update(strength=w,independent_units=units,inputs=o)
   if k in VSL and w=='Very Strong':ck['strength_policy']=lr_policy()['id']
   if k=='PM5':ck.update(subtype='genomic_PTC',event_exon=text(o.get('event_exon'),15),stop_exon=text(o.get('stop_exon'),15))
  oldx=oldrows.get(k)
  if oldx and all(oldx.get(key)==val for key,val in ck.items()) and all(key in ck for key in oldx if key not in ('reviewer','adjudicated_on')):ck.update(reviewer=oldx['reviewer'],adjudicated_on=oldx['adjudicated_on'])
  else:ck.update(reviewer=who,adjudicated_on=date)
  clean.append(ck)
 change=str(p.get('change_reason') or '').strip()
 if len(change)>1800:raise ValueError('修改说明过长')
 if old and set(oldrows)-seen and not change:raise ValueError('撤回以前的判据须填写此次修改理由')
 if not clean and not old:raise ValueError('尚未选择需要保存的判据')
 enabled={k for k,z in rows.items() if z['rule_enabled']}
 if flags['review_complete'] and (seen!=enabled or not flags['independence_reviewed']):raise ValueError('形成草案须逐项处置全部可用判据并复核独立性；不能把空白当作不满足')
 packet=dict(allele=r['caid'],condition=r['condition'],gene=gene,domain='germline',spec=GENES[gene][2]+'@1.2',origin='manual_adjudication',adjudicated_on=date,case_reviewer=who,criteria=clean,**flags)
 # Only server-validated applications can unlock an official textual strength in the frozen structured list.
 conditional=frozenset(x['code'] for x in clean if x.get('strength_policy')==lr_policy()['id'] and x.get('strength')=='Very Strong')
 result=adjudicate(packet,verified_text_codes=conditional)
 allowed=('evidence inventory / independence review not complete; subtotal only','significant conflict: expert discussion required','stand-alone / significant or unresolved conflict: expert discussion required')
 serious=[v for v in result['issues'] if v not in allowed and not(v=='no accepted evidence' and not any(x['state']=='met' for x in clean))]
 if serious:raise ValueError('未保存；需处理重复、冲突或规则约束：'+'；'.join(serious))
 return dict(schema='brca.local-review.v1',desk_version=__version__,policy_sha256=policy_hash(),is_training=bool(s.get('teaching')),snapshot=s['id'],base_sha256=s['report_sha256'],gene=gene,assembly=s['assembly'],caid=r['caid'],spec=r['specification'],rules_hash=r['rules_bundle_sha256'],revision=f'r{int(last[1:])+1:06d}' if last else 'r000001',parent=last,created_utc=now,reviewer=who,change_reason=change,decisions=clean,flags=flags,result=result)
def save(st,s,p,render):
 old=latest(st,s['id'])
 if p.get('schema')=='brca.selection.v2':
  from .selection import build
  rev=build(s,p,old)
 else:rev=evaluate(s,p,old)
 rev['parent_sha256']=old['content_sha256'] if old else None;base=root_for(st,s['id']);base.mkdir(parents=True,exist_ok=True)
 tmp=Path(tempfile.mkdtemp(prefix='.pending-',dir=base));dst=base/rev['revision']
 try:
  report=effective(s['report'],rev);v=dict(s,report=report,local_review=rev)
  writej(tmp/'report.json',report);rev['saved_report_sha256']=hashlib.sha256((tmp/'report.json').read_bytes()).hexdigest();rev['content_sha256']=seal(rev)
  writej(tmp/'review.json',rev);(tmp/'report.html').write_bytes(render([v]))
  try:tmp.rename(dst)
  except OSError:
   if dst.exists():raise ReviewChanged('另一份审核已先行保存，请刷新')
   raise
 finally:
  if tmp.exists():shutil.rmtree(tmp)
 return v
