"""Read a single ENIGMA germline assertion. External decisions are never local evidence votes."""
import json, re
from urllib.parse import urlencode
from ._engine.brcalite.core import GENES
from ._engine.brcalite.net import fetch, sha
BASE='https://erepo.clinicalgenome.org/evrepo/api'
PANEL='ENIGMA BRCA1 and BRCA2 VCEP'
UUID=re.compile(r'[0-9a-f]{8}-(?:[0-9a-f]{4}-){3}[0-9a-f]{12}')
ORDER=('PVS1','PS1','PS2','PS3','PS4','PM1','PM2','PM3','PM4','PM5','PM6','PP1','PP2','PP3','PP4','PP5','BA1','BS1','BS2','BS3','BS4','BP1','BP2','BP3','BP4','BP5','BP6','BP7')
def grade(s):
 s=re.sub('[ _-]','',str(s or '').lower().replace('pathogenic','').replace('benign',''))
 return dict(supporting='Supporting',moderate='Moderate',strong='Strong',verystrong='Very Strong',standalone='Stand Alone').get(s)
def parts(raw):
 v=str(raw);m=re.match(r'^(PVS1|[PB][SMPA]\d+)(?=$|[_(\s-])',v,re.I)
 if not m or m[1].upper() not in ORDER:return None,None
 return m[1].upper(),grade(re.sub(r'\(RNA\)|(?:^|_)PTC','',v[m.end():],flags=re.I))
def literature_mentions(notes,known=()):
 # A number without a PMID prefix is not automatically a citation. Retain explicit exclusions.
 out={str(p):{'pmid':str(p),'role':'mentioned','pointers':[]} for p in known if re.fullmatch(r'[1-9][0-9]{0,9}',str(p))}
 idr=r'[1-9][0-9]{0,9}';sep=r'(?:,|，|、|;|；|\band\b)'
 for n in notes:
  txt=n.get('value') or '';ptr=n.get('pointer')
  for m in re.finditer(r'\bPMIDs?\s*[:：]?\s*('+idr+r'(?:\s*'+sep+r'\s*'+idr+r')*)(?![0-9])',txt,re.I|re.ASCII):
   for pid in re.findall(idr,m[1]):
    o=out.setdefault(pid,{'pmid':pid,'role':'mentioned','pointers':[]})
    if ptr and ptr not in o['pointers']:o['pointers'].append(ptr)
    before=txt[max(0,m.start()-30):m.start()];after=txt[m.end():m.end()+65]
    if (re.search(r'排除\s*$|未(?:予)?(?:采用|纳入)\s*$',before) and not re.search(r'(?:未|不|尚未|无法|不能)排除\s*$',before)) or re.match(r'\s*\)?\s*(?:(?:was|were|is|are)\s+(?:therefore\s+)?)?(?:not (?:considered|included|used)|excluded)\b',after,re.I):o['role']='explicitly_excluded'
    elif o['role']=='mentioned' and (re.search(r'(?:专家|本例)纳入\s*$',before) or re.match(r'\s*(?:was|were) (?:explicitly )?included\b',after,re.I)):o['role']='explicitly_included'
 return list(out.values())

def blank(state):return {'state':state,'criteria':[],'alerts':[]}
def select(s,r):
 if not isinstance(s,dict) or not isinstance(s.get('data'),list):return None,'invalid_summary'
 rows=s['data'];target=[]
 for x in rows:
  if not isinstance(x,dict):return None,'invalid_summary'
  if x.get('caId')!=r['caid'] or x.get('gene')!=r['gene'] or x.get('mondoId')!=r['condition'] or x.get('retracted') is not False:continue
  if x.get('ep')!=PANEL or str(x.get('affiliationId'))!='50087' or str(x.get('moi')).lower()!='autosomal dominant inheritance':continue
  if UUID.fullmatch(str(x.get('uuid',''))):target.append(x)
 if len(rows)>=100 or len(target)>1:return None,'ambiguous_or_truncated'
 return (target[0],'matched') if target else (None,'no_matched_assertion')
def outcome_label(x):return {'met':'Met','not met':'Not Met'}.get(str(x or '').strip().lower())
def interpret(x,d,r):
 ans=blank('summary_only');ans.update(uuid=x['uuid'],url=BASE.rsplit('/',1)[0]+'/ui/classification/'+x['uuid'],detail_url=BASE+'/classification/'+x['uuid'],classification=x.get('classification'),published=x.get('publishedDate'),record_version=x.get('docVersion'),gene=x['gene'],condition=x['mondoId'],caid=x['caId'],review_panel=x.get('ep'),spec_version=None)
 items={};al=ans['alerts']
 for prop,status in [('metCodes','Met'),('unMetCodes','Not Met')]:
  if prop in x and not isinstance(x[prop],list):al.append(prop+' 格式未识别');continue
  for i,raw in enumerate(x.get(prop,[])):
   code,w=parts(raw)
   if not code:al.append('未映射专家原码：'+str(raw));continue
   if code in items:al.append('专家摘要中重复或矛盾：'+code);continue
   items[code]=dict(code=code,state=status,strength=w,raw_code=str(raw),origin='summary',text=[],pointer=None,summary_field=f'{prop}/{i}',pmids=[],issues=[])
   if not w and str(raw).strip()!=code:items[code]['issues'].append('原码后缀未完整解析，请读原文')
 ans['criteria']=list(items.values())
 if d is None:return ans
 if not isinstance(d,dict) or d.get('uuid')!=x['uuid'] or (d.get('condition') or {}).get('id')!=r['condition'] or (d.get('variant') or {}).get('id')!='CAR:'+r['caid'] or str((d.get('statementOutcome') or {}).get('label','')).casefold()!=str(x.get('classification')).casefold():
  ans['state']='detail_identity_mismatch';al.append('详情未通过身份或分类对账，只展示已匹配摘要');return ans
 method=d.get('assertionMethod') or {};ans['spec_version']=method.get('version')
 if not isinstance(d.get('evidenceLine'),list):al.append('专家详情没有可解析的 evidenceLine');return ans
 seen=set()
 for n,line in enumerate(d['evidenceLine']):
  for i,y in enumerate(line.get('evidenceItem') or []):
   if not isinstance(y,dict) or 'criterion' not in y:continue
   crit=y['criterion'];code,_=parts(crit.get('label'));status=outcome_label((y.get('statementOutcome') or {}).get('label'))
   if not code or not status:al.append('未映射专家详情判据或状态：'+str(crit.get('label')));continue
   if code in seen:al.append('详情重复：'+code);continue
   seen.add(code);p=f'/evidenceLine/{n}/evidenceItem/{i}'
   z=items.setdefault(code,dict(code=code,state=status,strength=None,raw_code=crit.get('label'),summary_field=None,issues=[]))
   if not z.get('summary_field'):z['issues'].append('摘要未列此项，见详情')
   if z['state']!=status:z['issues'].append('摘要与详情判定不一致')
   explicit=grade((y.get('evidenceStrength') or {}).get('label')) or grade((line.get('evidenceStrength') or {}).get('label'))
   default=grade((crit.get('defaultStrength') or {}).get('label'))
   if z['strength'] and explicit and z['strength']!=explicit:z['issues'].append('摘要与详情采用的强度不一致')
   grade_origin='summary_suffix' if z['strength'] else 'detail_explicit' if explicit else 'criterion_default' if default else 'not_disclosed'
   w=z['strength'] or explicit or default
   notes=[]
   for j,sub in enumerate(y.get('evidenceLine') or []):
    for k,unit in enumerate(sub.get('evidenceItem') or []):
     for h,con in enumerate(unit.get('contribution') or []):
      note=con.get('comments')
      if isinstance(note,str) and note.strip():notes.append({'value':note.strip(),'pointer':f'{p}/evidenceLine/{j}/evidenceItem/{k}/contribution/{h}/comments'})
   z.update(strength=w,strength_origin=grade_origin,origin='summary_and_detail' if z.get('summary_field') else 'detail',pointer=p,text=notes,
            detail_state=status,detail_strength=explicit,criterion_strength=default,
            pmids=sorted(a['pmid'] for a in literature_mentions(notes)),publication_mentions=literature_mentions(notes))
   if not w:z['issues'].append('来源未披露可解析的强度')
 for c,z in items.items():
  if c not in seen:z['issues'].append('仅摘要列出，未在详情找到')
 ans.update(state='matched',criteria=list(items.values()));return ans
def explicit_empty(folder,url):
 p=folder/'http'/'erepo_summary'
 try:
  m=json.loads((p/'receipt.json').read_bytes());a=m['attempts'][-1];name=a['body_file']
  if m['method']!='GET' or m['url']!=url or a.get('status')!=404 or not re.fullmatch(r'failed-\d+\.body',name):return None
  raw=(p/name).read_bytes();o=json.loads(raw)
  if sha(raw)!=a.get('body_sha256') or o.get('data') not in (None,[]):return None
  status=o.get('status') or {};by=o['metadata']['rendered']['by']
  if status.get('code')==404 and status.get('msg')=='No records were found for given query' and by in (BASE+'/summary/srvc','https://erepo.genome.network/evrepo/api/summary/srvc'):return sha(raw)
 except (OSError,KeyError,IndexError,AttributeError,TypeError,ValueError):pass
 return None
def collect(folder,r,reader=fetch):
 if r.get('status')!='review_required' or r.get('gene') not in GENES or not re.fullmatch(r'CA\d+',str(r.get('caid',''))):return blank('blocked_identity')
 url=BASE+'/summary/classifications?'+urlencode(dict(columns='caId',values=r['caid'],matchTypes='exact',pgSize=100,pg=1))
 def request(key,url):
  raw=reader(folder/'http',key,url,attempts=2,timeout=12)
  r.setdefault('source_files',{})[key]='http/'+key+'/receipt.json';r.setdefault('source_status',{})[key]='retrieved' if raw is not None else 'unavailable'
  if raw is None:return None,None
  try:data=json.loads(raw)
  except (ValueError,TypeError,UnicodeDecodeError):
   data={};r['source_status'][key]='invalid_json'
  return data,sha(raw)
 s,h=request('erepo_summary',url);x,state=select(s,r) if s is not None else (None,'summary_unavailable')
 if s is None:
  empty=explicit_empty(folder,url)
  if empty:
   r['source_status']['erepo_summary']='no_records_for_query'
   return dict(blank('explicit_no_records'),summary_url=url,summary_sha256=empty)
 if not x:return dict(blank(state),summary_url=url,summary_sha256=h)
 d,h2=request('erepo_detail',BASE+'/classification/'+x['uuid'])
 try:result=interpret(x,d,r)
 except (ValueError,TypeError,KeyError,AttributeError):
  result=interpret(x,None,r);result['alerts'].append('详情结构异常，请打开原始记录')
 return dict(result,summary_url=url,summary_sha256=h,detail_sha256=h2)
