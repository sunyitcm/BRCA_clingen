"""A selection carries no clinical signature: attached evidence and explicit manual judgments stay separate; allowed grades and weights are validated on the server."""
import json
from pathlib import Path
from collections import Counter
from copy import deepcopy
from . import __version__
from .criteria import ledger,verified,GENES,R,getj,protein_domain
from ._engine.brcalite.score import P,combine
from ._engine.brcalite.net import utc,sha
from .adoption import ReviewChanged,policy_hash
SCHEMA='brca.selection.v2'
def stats(ds):
 if not ds:return {'points':None,'stand_alone':False}
 stand=any(x['strength']=='Stand Alone' for x in ds)
 return {'points':None if stand else sum((1 if x['code'].startswith('P') else -1)*P[x['strength']] for x in ds),'stand_alone':stand}
def build(s,p,old):
 r=s['report'];g=r.get('selected_gene') or r.get('gene')
 ok=verified(r,g);tab=ledger(r)
 if (not ok and not tab.get('draft_selectable')) or r.get('rules_bundle_sha256')!=sha((R/'rules/bundle.json').read_bytes()):raise ValueError('等位、转录本或本基因规则未核准。')
 if (s['report_sha256'],s['assembly'],r['caid'],g,s['id'])!=(p.get('base_hash'),p.get('assembly'),p.get('caid'),p.get('gene'),p.get('snapshot')):raise ValueError('当前选择与取数快照不符。')
 if r.get('met') or r.get('not_met') or r.get('points_accepted') is not None or old and old.get('schema')!=SCHEMA:raise ValueError('已有旧版审核；原结果只读，请新建查询快照。')
 last=old.get('revision') if old else None
 if last!=p.get('previous_revision'):raise ReviewChanged('已有另一个已保存版本；请刷新后再保存。')
 chosen=p.get('selected');note=p.get('note','')
 if not isinstance(note,str):raise ValueError('备注需为文字，可以留空。')
 if not isinstance(chosen,list) or len(chosen)>19 or not all(isinstance(k,str) for k in chosen) or len(chosen)!=len(set(chosen)):raise ValueError('所选判据无效或重复。')
 manual=p.get('manual',{})
 if not isinstance(manual,dict) or any(not isinstance(k,str) or not isinstance(w,str) or k not in chosen for k,w in manual.items()):raise ValueError('人工强度必须对应本次勾选的判据。')
 table=tab;rows={z['code']:z for z in table['rows']};ds=[];units=set()
 rule=table['guideline_source']
 for k in chosen:
  z=rows.get(k,{});op=(z.get('brief') or {}).get('option')
  if not z or not z.get('rule_enabled'):raise ValueError(k+'：不在本基因启用的规范内。')
  if not ok or z['evidence_view']['blocks_automatic_scoring']:
   w=manual.get(k) or (z['allowed_strengths'][0] if len(z['allowed_strengths'])==1 else None)
   if w is not None and w not in z['allowed_strengths']:raise ValueError(k+'：强度不适用于本基因。')
   ds.append(dict(code=k,state='pending',strength=w,origin='manual',source_refs=[],independent_units=[],stand_alone=False,contribution_points=None,hold='identity_unresolved' if not ok else z['evidence_view'].get('hold_reason') or 'incompatible_observation'))
   continue
  if not op:
   w=manual.get(k)
   if not w or w not in z['allowed_strengths']:raise ValueError(k+'：请明确选择本基因规范允许的人工强度。')
   if k=='PM2' and r.get('facts',{}).get('event') in ('INS','DEL','DELINS'):raise ValueError('PM2 不适用于已知插入、缺失或 delins。')
   if k=='BP1' and protein_domain(r) is True:raise ValueError('BP1 要求位于本基因规范功能域外；本次规范转录本位置在域内。')
   c=z['guideline_zh'];gd=c['grades'][w]
   # A rule locator explains the choice, but is not an independent experiment or patient record.
   ds.append(dict(code=k,state='met',strength=w,origin='manual',source_refs=[],independent_units=[],manual_attestation=True,rule_ref=dict(spec=rule['spec'],rule_key=c['rule_key'],document_url=rule['document_url'],document_sha256=rule['document_sha256'],pointer=gd['pointer'],text_sha256=gd['text_sha256'],conditional_in_text=gd['conditional_in_text'])))
   continue
  if k in manual:raise ValueError(k+'：已有可采纳的来源，不能在同一次选择覆盖其强度。')
  if op['strength'] not in z['allowed_strengths']:raise ValueError(k+'：强度不在本基因规范内。')
  if units.intersection(op['units']):raise ValueError('同一依据不可重复计入。')
  units.update(op['units']);ds.append(dict(code=k,state='met',strength=op['strength'],origin=op['kind'],source_refs=deepcopy(op['refs']),independent_units=op['units']))
 all_decisions=ds;ds=[x for x in ds if x['state']=='met'];pending=[x['code'] for x in all_decisions if x['state']=='pending'];keys={x['code'] for x in ds}
 for group in [('PP3','PVS1'),('PP3','BP4'),('BP1','BP4'),('BP1','PP3'),('BA1','BS1','PM2'),('PS3','BS3'),('PP4','BP5'),('PP1','BS4')]:
  if len(keys.intersection(group))>1:raise ValueError('该组合需消除互斥或重叠：'+'、'.join(k for k in group if k in keys))
 if 'PM3' in keys and keys.intersection(getj(R/'rules/review-20261124.json')[g]['excludes_if_met']):raise ValueError('PM3 与已选人群判据不能同时采用。')
 if 'PM5' in keys and 'PVS1' not in keys:raise ValueError('按此基因规范，PM5_PTC 须核 PVS1。')
 for d in ds:
  d['stand_alone']=d['strength']=='Stand Alone'
  d['contribution_points']=None if d['stand_alone'] else (1 if d['code'].startswith('P') else -1)*P[d['strength']]
 result=stats(ds);local=[x for x in ds if x['origin'] not in ('expert_citation','manual')];man=[x for x in ds if x['origin']=='manual'];quoted=[x for x in ds if x['origin']=='expert_citation'];directions={x['code'][0] for x in ds}
 caveats=[];pattern=None
 if any(rows[d['code']]['brief']['status'] in ('conflict','expert_not_met') for d in man):caveats.append('手工采用的判据与已取得的资料或专家处置不一致；保留两者，合计暂缓。')
 if 'BP7' in keys and any(d['code']=='BP7' and d['strength']=='Supporting' for d in ds) and 'BP4' not in keys:caveats.append('BP7 Supporting 分支须满足 BP4；请先核准并选入 BP4 后合算。')
 seen_papers=set()
 for d in ds:
  papers={str(v) for ref in d['source_refs'] for v in ref.get('pmids',[])}
  if papers&seen_papers:caveats.append('两个判据提及同一研究，不能确认实验或队列独立；合并计算暂缓。');break
  seen_papers.update(papers)
 if len(directions)>1:caveats.append('同时有致病和良性证据，未做独立冲突判定。')
 if keys&{'PP4','BP5'} and keys&{'PP1','BS4','PS4','BS2','PM3'}:caveats.append('综合 LR 可能包含家系或病例数据，组合提示暂缓。')
 if 'PM5' in keys:caveats.append('终止外显子及 PVS1 适用条件需另核。')
 if keys&{'PVS1','PS1','BP7'} and len(keys&{'PVS1','PS1','BP7'})>1:caveats.append('需排查 RNA 证据重复使用。')
 # Pure same-direction pattern, not a signed classification. No implied negatives or invented independence sign-off.
 if ds and not caveats:
  ptn=combine(Counter(x['strength'] for x in ds),next(iter(directions)),ds)
  if ptn!='Uncertain significance':pattern=ptn
 result.update(pending_codes=pending,local_sources=stats(local),quoted_sources=stats(quoted),manual_sources=stats(man),manual_codes=[x["code"] for x in man],local_codes=[x['code'] for x in local],quoted_codes=[x['code'] for x in quoted],selected_combination=pattern,issues=caveats,mode='selected_evidence_only',draft_class=None,local_final=None)
 if caveats:result['points']=None # no validated combined subtotal when there may be double-counting or conflict
 return dict(schema=SCHEMA,desk_version=__version__,policy_sha256=policy_hash(),is_training=bool(s.get('teaching')),snapshot=s['id'],base_sha256=s['report_sha256'],gene=g,assembly=s['assembly'],caid=r['caid'],spec=table['spec'],rules_hash=r['rules_bundle_sha256'],revision=f'r{int(last[1:])+1:06d}' if last else 'r000001',parent=last,created_utc=utc(),note=note,decisions=all_decisions,result=result)
def present(report,rev):
 r=deepcopy(report);r.pop('not_applicable_by_spec',None);ds=rev['decisions'];out=rev['result']
 r['criteria_table']=ledger(r,ds);r.update(review_desk_version=rev['desk_version'],review_policy_sha256=rev['policy_sha256'],review_revision=rev['revision'],reviewed_utc=rev['created_utc'],is_training=rev['is_training'],met=[d['code'] for d in ds if d['state']=='met'],not_met=[],not_evaluated=[z['code'] for z in r['criteria_table']['rows'] if not z['decision'] or z['decision']['state']=='pending'],classification_draft=None,classification_local=None,points_accepted=out['points'],points_accepted_basis='selected_evidence_with_source_breakdown',review_result=out,review_note=rev['note'])
 return r
def population_from_file(path,r):
 if not path:return
 # Explicit machine-readable evidence optional; never scan for files or transmit them.
 b=Path(path).read_bytes()
 if len(b)>12_000_000:raise ValueError('人群证据文件过大。')
 data=json.loads(b);rec=data.get('records') if isinstance(data,dict) and data.get('schema')=='brca.population.v1' else None
 if not isinstance(rec,list):raise ValueError('人群证据文件格式不符。')
 tx=(r.get('facts',{}).get('transcript') or {}).get('hgvsc')
 matches=[q for q in rec if isinstance(q,dict) and q.get('caid')==r.get('caid') and q.get('gene')==r.get('gene') and q.get('hgvsc')==tx]
 if len(matches)>1:raise ValueError('同一等位存在多条人群记录，未采用。')
 if matches:
  q=matches[0];r['facts']['population_evidence']={'gene':r['gene'],'hgvsc':tx,'caid':r['caid'],'file_sha256':sha(b),'inputs':q.get('inputs') or {}}
