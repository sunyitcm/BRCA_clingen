"""Verify the rule locators against the reader's official English JSON, offline."""
import argparse,hashlib,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'source'))
from brca_review.guidance import catalog_zh
def get(d,p):
 for k in p.strip('/').split('/'):
  key=k.replace('~1','/').replace('~0','~');d=d[int(key) if isinstance(d,list) else key]
 return d
def check(gene,path):
 x=catalog_zh()['genes'][gene];b=Path(path).read_bytes()
 if hashlib.sha256(b).hexdigest()!=x['document_sha256']:raise ValueError('官方 JSON 的文件哈希与该基因冻结 v1.2 不一致。')
 doc=json.loads(b);rows=doc['ld']['RuleSet'][0]['ld']['CriteriaCode'];enabled={r['entContent']['label'] for r in rows if r['entContent']['applicability']=='Applicable'}
 if set(x['criteria'])!=enabled:raise ValueError('适用判据集合发生变化。')
 for k,c in x['criteria'].items():
  if get(doc,c['instructions']['pointer'].rsplit('/',1)[0])['_uniqueProp']!=c['rule_key']:raise ValueError('条目标识不符：'+k)
  for field in [c['generic_acmg'],c['instructions'],*c['grades'].values()]:
   if hashlib.sha256(get(doc,field['pointer']).encode()).hexdigest()!=field['text_sha256']:raise ValueError('英文条文哈希不符：'+k)
  for grade,v in c['grades'].items():
   descriptor=get(doc,v['pointer'].rsplit('/',1)[0])
   if descriptor['applicability']!='Applicable' or descriptor['strength']!=('Strong' if v['conditional_in_text'] else grade):raise ValueError('原始强度栏不符：'+k)
 return x,doc
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('gene',choices=['BRCA1','BRCA2']);p.add_argument('official_json',type=Path);p.add_argument('--code',help='选填：显示一条原文及本项目中文释义');p.add_argument('--strength',help='选填：显示指定强度的原文；用于 --code');a=p.parse_args()
 try:
  x,doc=check(a.gene,a.official_json)
  print(f"OK · {x['spec']} · {len(x['criteria'])} 条 · {sum(len(c['grades']) for c in x['criteria'].values())} 个强度定位与冻结英文逐一相符")
  if a.code:
   c=x['criteria'][a.code.upper()]
   print(c['rule_key'],c['zh'],sep='\n')
   for w,r in ([(a.strength,c['grades'][a.strength])] if a.strength else c['grades'].items()):
    print(w,r['zh'],r['pointer'],get(doc,r['pointer']),sep='\n')
  elif a.strength:p.error('--strength 需要 --code')
 except (KeyError,ValueError,OSError,IndexError) as e:p.exit(2,str(e)+'\n')
