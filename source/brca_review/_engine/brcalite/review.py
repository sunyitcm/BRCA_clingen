from pathlib import Path
import os
from .core import GENES
from .net import writej,getj,sha,utc
from .score import adjudicate
from .report import html_report
def attach(report,packet,out):
 r=getj(report);p=getj(packet);g=r['gene']
 if r['status']!='review_required' or not r.get('caid') or p.get('allele')!=r['caid'] or p.get('gene')!=g or p.get('condition')!=r['condition'] or p.get('spec')!=GENES[g][2]+'@1.2':raise ValueError('审核单与报告身份或版本不符')
 from .core import R
 if r.get('rules_bundle_sha256')!=sha((R/'rules/bundle.json').read_bytes()):raise ValueError('审核报告的冻结规则包与当前不同')
 q=adjudicate(p);out.mkdir(parents=True,exist_ok=False)
 for name,ref in r['source_files'].items():
  if isinstance(ref,str):r['source_files'][name]=os.path.relpath(report.parent/ref,out)
 if r.get('literature'):
  r['literature']['reading_card']=os.path.relpath(report.parent/r['literature']['reading_card'],out)
 r['parent_run_id']=r.get('run_id');r['run_id']=out.name
 r['review_result']=q;r['points_accepted']=q.get('points');r['classification_local']=None;r['classification_draft']=q['draft_class'];r['reviewed_at_utc']=utc()
 r['previous_report_sha256']=sha(report.read_bytes());r['review_packet_sha256']=sha(packet.read_bytes());r['curator']=p.get('case_reviewer');r['status']='draft_reviewed' if q['draft_class'] else 'review_blocked'
 r['met']=[];r['not_met']=[]
 if q['points'] is not None or q['draft_class']:
  r['met']=[x['code'] for x in p['criteria'] if x.get('state')=='met'];r['not_met']=[x['code'] for x in p['criteria'] if x.get('state')=='not_met']
  r['reviewed_not_applicable']=[x['code'] for x in p['criteria'] if x.get('state')=='not_applicable']
  r['not_evaluated']=[k for k in r['not_evaluated'] if k not in r['met']+r['not_met']+r['reviewed_not_applicable']]
 r['review_lines']=p.get('criteria',[]);r['limitations']+=q['issues']
 r['criterion_dispositions']={}
 for line in r['review_lines']:
  if not isinstance(line,dict) or not isinstance(line.get('code'),str):continue
  code=line['code'];accepted=code in r['met']
  r['criterion_dispositions'][code]={
   'requested_state':line.get('state'),'reviewed_strength':line.get('strength'),
   'accepted':accepted,'contributes_to_subtotal':accepted and r['points_accepted'] is not None,
   'accepted_as_standalone':accepted and code=='BA1','source_refs':line.get('source_refs',[])}
 writej(out/'review.json',p);writej(out/'report.json',r);html_report(out/'report.html',r);return q
