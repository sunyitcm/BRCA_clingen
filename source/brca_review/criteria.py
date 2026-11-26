"""One deterministic row per ClinGen main code; used by JSON, HTML and review validation."""
import re,hashlib
from copy import deepcopy
from decimal import Decimal, InvalidOperation
from functools import lru_cache
from pathlib import Path
from ._engine.brcalite.core import GENES,R,parse
from ._engine.brcalite.net import getj
from .erepo import ORDER,parts,blank
from ._engine.brcalite.score import P
CHROM={'BRCA1':'17','BRCA2':'13'}
DOM={'BRCA1':[(2,101),(1391,1424),(1650,1857)],'BRCA2':[(10,40),(2481,3186)]}
LIMIT={'BRCA1':1863,'BRCA2':3418}
PROTEIN={'BRCA1':(Decimal('.28'),Decimal('.15')),'BRCA2':(Decimal('.30'),Decimal('.18'))}
LABEL=dict(zip(ORDER,('功能丧失 / RNA','已知同效应对照','确认新生','异常功能实验','病例对照富集','热点或关键域','对照人群缺失','范可尼贫血与相位','蛋白长度变化','终止外显子（PM5_PTC）','未确认新生','正向共分离','基因错义机制','有害计算预测','正向临床似然比','权威致病断言','独立良性频率','良性人群频率','共现及隐性表型','正常功能实验','反向共分离','域外且无剪接影响','顺反式共现','重复区框内改变','无影响计算预测','反向临床似然比','权威良性断言','同义、内含子及正常 RNA')))
@lru_cache(maxsize=1)
def lr_policy():return getj(Path(__file__).resolve().parent/'lr_v12.json')
VSL=frozenset(lr_policy()['criteria'])
@lru_cache(maxsize=1)
def catalog():return getj(R/'rules/catalog.json')
def scope(g,build,value):
 if g not in GENES or not isinstance(value,str) or len(value)>16384:raise ValueError('请选择基因，并一次输入一个核酸 ALT。')
 value=value.strip()
 if not value or any(c in value for c in '\r\n'):raise ValueError('只接受一个位点；VCF 请逐行、逐 ALT 输入。')
 if re.search(r'[ \t]',value):
  cols=re.split(r'[ \t]+',value)
  if len(cols)<5 or not re.fullmatch(r'(?:chr)?(?:13|17)',cols[0]) or not cols[1].isdigit():raise ValueError('VCF 行应为 CHROM POS ID REF ALT；不能粘贴表头或多行。')
  value='-'.join(cols[i] for i in (0,1,3,4))
 if len(value)>250:raise ValueError('位点过长。')
 v=parse(value,build)
 if v.get('gene',g)!=g or v.get('chrom',CHROM[g])!=CHROM[g]:raise ValueError(f'{g} 入口只接受 chr{CHROM[g]} 或 {GENES[g][0]}:c.；请检查基因与组装。')
 return v
def verified(r,g):
 tx,condition,spec=GENES.get(g,(None,None,None));t=(r.get('facts') or {}).get('transcript') or {}
 return bool(tx and r.get('gene')==g and r.get('condition')==condition and r.get('specification')==spec+' v1.2'
             and r.get('status') in ('review_required','draft_reviewed') and re.fullmatch(r'CA\d+',str(r.get('caid')))
             and t.get('gene_symbol')==g and t.get('transcript_id')==tx and str(t.get('hgvsc') or '').startswith(tx+':c.'))
def protein_domain(r):
 t=r.get('facts',{}).get('transcript') or {};g=r.get('gene');p=t.get('protein_start')
 if g not in GENES or t.get('transcript_id')!=GENES[g][0] or type(p)!=int or not 1<=p<=LIMIT[g]:return None
 return any(a<=p<=b for a,b in DOM[g])
def dec(v):
 try:
  n=Decimal(str(v));return n if n.is_finite() else None
 except (TypeError,ValueError,InvalidOperation):return None
def ledger(r,decisions=()):
 g=r.get('selected_gene') or r.get('gene');meta={'schema':'brca.criteria.v2','gene':g,'spec':GENES[g][2]+' v1.2' if g in GENES else None,'identity_verified':verified(r,g),'rows':[],'alerts':[]}
 if g not in GENES:return meta
 locked=r.get('rules_bundle_sha256')==hashlib.sha256((R/'rules/bundle.json').read_bytes()).hexdigest()
 meta['reviewable']=meta['identity_verified'] and locked and not (r.get('met') or r.get('not_met') or r.get('points_accepted') is not None)
 meta['draft_selectable']=bool(locked and r.get('acquisition')=='live_independent' and r.get('selected_gene')==g and not (r.get('met') or r.get('not_met') or r.get('points_accepted') is not None))
 current=catalog()[g];f=r.get('facts') or {};states=r.get('source_status') or {};expert=r.get('erepo') or blank('not_queried_in_snapshot')
 registry={};review={x['code']:x for x in decisions};ok=meta['identity_verified'];meta['expert_state']=expert['state'];meta['alerts']=list(expert.get('alerts',[]))
 for c in ORDER:
  pol=current[c];enabled=pol['enabled'];z=dict(code=c,label=LABEL[c],rule_enabled=enabled,allowed_strengths=list(pol['strengths']),rule_field=pol['field'],rule_url='https://cspec.genome.network/cspec/ui/svi/doc/'+GENES[g][2]+'?version=1.2',conditional_strengths=['Very Strong'] if c in VSL else [],expert=None,candidates=[],facts=[],gaps=[],local='unassessed' if enabled else 'spec_not_applicable',decision=None,accepted_strength=None,contribution_points=None,stand_alone=False,priority=False)
  z['structured_strengths']=list(pol['strengths']);z['conditional_rule']=None
  if c in VSL and enabled and meta['spec'] in lr_policy()['profiles']:
   policy=lr_policy();cfg=policy['criteria'][c];z['conditional_rule']=dict(policy_id=policy['id'],source=policy['sources'][g],section=cfg['section'],condition=cfg['condition'],lr_type=cfg['lr_type'],direction=cfg['direction'],threshold=policy['directional_thresholds'][cfg['direction']]['Very Strong'])
   if 'Very Strong' not in z['allowed_strengths']:z['allowed_strengths'].insert(0,'Very Strong')
  registry[c]=z
 if not ok:meta['alerts'].append('选择的基因、规范转录本或等位身份尚未核准，禁止本地计分。')
 def fact(c,label,value,pointer,source):
  registry[c]['facts'].append(dict(label=label,value=value,pointer=pointer,source=source))
 if ok:
  if expert.get('state')=='matched' and expert.get('spec_version')!='1.2':meta['alerts'].append('个案专家规范版本未披露或与本地 v1.2 不同，保留两者。')
  if expert.get('classification'):
   def cls(v):return str(v).lower().replace('_',' ').strip()
   if any(rec.get('class') and cls(rec['class'])!=cls(expert['classification']) for rec in r.get('expert_assertions',[])):meta['alerts'].append('ClinVar 专家提交与 eRepo 分类文字不同，请先核疾病、提交方及版次。')
  for er in (expert.get('criteria',[]) if expert.get('state')=='matched' else []):
   c=er.get('code')
   if c in registry:
    registry[c]['expert']=deepcopy(er)
    if er.get('issues'):registry[c]['gaps']+=er['issues'];registry[c]['priority']=True
    if er['state']=='Met':
     registry[c]['priority']=True
     if not registry[c]['rule_enabled'] or er.get('strength') and er['strength'] not in registry[c]['allowed_strengths']:
      registry[c]['gaps'].append('专家采用强度不在本地冻结可选清单；保留原标记，本地不自动计分。')
  for i,cand in enumerate(r.get('candidates') or []):
   src=cand.get('source') or {};base,w=parts(cand.get('code'))
   if base not in registry:meta['alerts'].append('未映射附件原码：'+str(cand.get('code')));continue
   if src.get('gene')!=g or src.get('hgvsc','')!=(f['transcript']['hgvsc'].partition(':')[2]):meta['alerts'].append('附件行身份不符，未采用');continue
   strength=cand.get('strength') or w;key=f"WB:{src.get('sheet')}:{src.get('row')}:{src.get('sha256')}"
   registry[base]['candidates'].append(dict(raw_code=cand.get('code'),strength=strength,sheet=src.get('sheet'),row=src.get('row'),cells=src.get('review_cells') or src.get('cell') or src.get('columns'),sha256=src.get('sha256'),ref=key,pmids=src.get('pmids_mentioned') or [],values=src.get('review_values') or cand.get('verified_cells') or {},pointer='/candidates/'+str(i),workbook_url=(cand.get('workbook') or {}).get('url') or next((a['url'] for a in getj(R/'rules/attachments.lock.json') if a['gene']==g and a['sha256']==src.get('sha256')),None)))
   registry[base]['priority']=True
  d=protein_domain(r);t=f.get('transcript') or {};bx=f.get('ldh_brca') or {};s=bx.get('BayesDel_nsfp33a_noAF');terms=t.get('consequence_terms') or []
  if 'missense_variant' in terms and d is not None:
   if not d:registry['BP1']['priority']=True;registry['BP1']['gaps'].append('规范转录本域外错义；需核 SpliceAI 的设置与最大 Δ≤0.1。')
   if s is not None:
    n=dec(s);c='BP4' if d and n is not None and n<=PROTEIN[g][1] else 'PP3'
    fact(c,'BayesDel noAF（蛋白）',s,'/facts/ldh_brca/BayesDel_nsfp33a_noAF','ldh_brca')
    fact(c,'规范转录本氨基酸位置',t.get('protein_start'),'/facts/transcript/protein_start',(f.get('transcript_sources') or {}).get('protein_start','vep'))
    if d and n is not None and (n>=PROTEIN[g][0] or n<=PROTEIN[g][1]):
     registry[c]['priority']=True;registry[c]['gaps'].append('该分支对照见本行检索结果与标准。')
    for i,cc in enumerate(r.get('candidates') or []):
     src=cc.get('source') or {};values=src.get('review_values') or {}
     if src.get('table')=='T9' and values.get('H') not in (None,''):
      fact(c,'Table 9 剪接预测列（非独立重算）',values['H'],f'/candidates/{i}/source/review_values/H','workbook')
  ex=f.get('ldh_exome') or {}
  for c in ('PM2','BS1','BA1'):
   if ex:
    for k in ('AC','AN','AF','quality'):
     if k in ex:fact(c,'LDH v4.1 exome '+k,ex[k],'/facts/ldh_exome/'+k,'ldh_exome')
   registry[c]['gaps'].append('按本基因冻结 v1.2 核 v2.1 非癌外显子与 v3.1 非癌队列、祖源、适用统计量及深度；v4.1 或无链接不能代替。')
   if c=='PM2' and f.get('event')=='SNV' and (not ex or any(ex.get(k) is None for k in ('AC','AN'))):
    registry[c]['gaps'].append('本次取数没有足以证实目标等位在规范两队列均缺失的计数、覆盖和质控；不能自动采纳 PM2_Supporting。')
   if c=='PM2' and states.get('ldh_exome')=='no_link_listed':fact(c,'LDH 外显子频率关联','no_link_listed','/source_status/ldh_exome','ldh')
  if f.get('exon_workbook_context'):
   for c in ('PVS1','PM5'):
    fact(c,'外显子工作簿线索',f['exon_workbook_context'],'/facts/exon_workbook_context','workbook');registry[c]['priority']=True
  for pr in r.get('manual_review_prompts') or []:
   k,_=parts(pr.get('criterion'))
   if k in registry:registry[k]['priority']=True;registry[k]['gaps'].append(pr.get('check') or '请复核指南条件')
 for k,z in registry.items():
  if k in review:
   x=review[k];z['local']=x['state'];z['decision']=deepcopy(x)
   if x['state']=='met':
    z['accepted_strength']=x['strength'];z['stand_alone']=x['strength']=='Stand Alone'
    if not z['stand_alone']:z['contribution_points']=(1 if k.startswith('P') else -1)*P[x['strength']]
   if x['state']=='met' or (z['expert'] and (x['state']=='met')!=(z['expert']['state']=='Met')):z['priority']=True
   if z['expert'] and ((x['state']=='met')!=(z['expert']['state']=='Met') or x['state']=='met' and z['expert'].get('strength') and x['strength']!=z['expert']['strength']):z['gaps'].append('本地审核与外部判定或强度不同，保留双方依据，复核差异。')
  elif k in (r.get('met') or []) or k in (r.get('not_met') or []):
   z['local']='met' if k in (r.get('met') or []) else 'not_met';z['priority']=True;z['gaps'].append('原报告已有人工处置，强度以原审核记录为准；此历史记录只读。')
  elif z['candidates']:z['local']='candidate'
  elif z['facts'] and z['rule_enabled']:z['local']='facts_only'
  if z['expert'] and not z['gaps'] and z['expert']['state']=='Met' and not z['candidates'] and not z['facts']:z['gaps'].append('只有本条外部断言；本地尚未提供独立来源。')
 for k,v in states.items():
  if v in ('unavailable','http_error','transport_error','unresolved_links','allele_not_verified','unrecognized_entity_url','invalid_json'):meta['alerts'].append(k+'：'+v+'；不能据此排除证据。')
 meta['other_sources']=[]
 if ok:
  for key,label in [('workbook_rows_without_assignment','工作簿有关行未赋码；论文与实验只供复核'),('ldh_associations_only','LDH 关联目录；有链接不等于有合格实验或预测')]:
   if f.get(key):meta['other_sources'].append(dict(label=label,pointer='/facts/'+key,value=f[key]))
 meta['population_prompt']={'code':'PM2','strength':'Supporting','spec':meta['spec'],'status':'local_review_required','message':'两个基因均有 PM2_Supporting。单次无链接或单一版本无变异记录不等于满足；需本基因规范要求的队列和覆盖。'} if ok and f.get('event')=='SNV' else None
 from .briefing import describe
 from . import guidance
 meta["guideline_source"]=guidance.source(g)
 visible=[registry[c] for c in ORDER if current[c]['enabled']]
 for z in visible:
  z['candidates']=[c for c in z['candidates'] if not c['strength'] or c['strength'] in z['allowed_strengths']]
  z['guideline_zh']=guidance.criterion(g,z)
  z['brief']=describe(r,z,ok and locked)
  from .observations import observation
  z['evidence_view']=observation(r,z,ok)
  ev=z['evidence_view']
  if r.get('acquisition')=='live_independent' and (ev['has_observations'] or ev['state']=='does_not_meet') or z.get('decision',{} ) and z['decision']['state']=='pending':z['priority']=True
  if z['brief'].get('option') or z['brief'].get('status')=='conflict':z['priority']=True
  if ok and z['code']=='PM2' and f.get('event')=='SNV' and (states.get('ldh_exome') or f.get('ldh_exome') or f.get('population_evidence')):z['priority']=True
 from .presentation import assessment,reading_guide
 for z in visible:z['assessment']=assessment(z,ok)
 meta['literature']=reading_guide(r,ok)
 meta['rows']=visible;meta['applicable_count']=len(visible)
 meta['all_evidence_missing']=not any(z['brief']['option'] or any(a.get('kind') not in ('searched_no_record','null_field','source_failure','not_reported','not_queried') and a.get('value')!='no_link_listed' for a in z['evidence_view']['observations']) for z in visible)
 return meta
