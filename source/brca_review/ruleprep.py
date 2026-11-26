"""Prepare the locked ENIGMA workbooks on the user's own device. No case data is sent."""
import gzip, hashlib, json, os, re, tempfile
from pathlib import Path
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET
from zipfile import ZipFile
from ._engine.brcalite.workbook import MAIN, REL
BASE=Path(__file__).resolve().parent/'_engine/rules'
_directory=None
_prepared=None
def sha(b):return hashlib.sha256(b).hexdigest()
def load(p):return json.loads(p.read_bytes())
def configure(directory,workbooks=None):
 global _directory,_prepared
 _directory=Path(directory)/'guidelines'/'GN092-GN097-v1.2';_prepared=None
 global _source
 _source=Path(workbooks) if workbooks else None
 from ._engine.brcalite.core import row_index
 row_index.cache_clear()
def parse_rows(source,records):
 data={'T9':[],'T4_splice':[],'T4_event':[],'ST1':[],'ST16':[]}
 for filename in sorted(records):
  a=records[filename];p=source/filename
  with ZipFile(p) as w:
   strings=[''.join(t.text or '' for t in el.iter(MAIN+'t')) for el in ET.fromstring(w.read('xl/sharedStrings.xml'))] if 'xl/sharedStrings.xml' in w.namelist() else []
   rel={x.get('Id'):x.get('Target') for x in ET.fromstring(w.read('xl/_rels/workbook.xml.rels'))}
   for sheet in ET.fromstring(w.read('xl/workbook.xml')).iter(MAIN+'sheet'):
    name=sheet.get('name');st=re.match(r'ST(\d+) ',name);n=int(st[1]) if st else None
    if not (name=='Table9_BRCA12VCEP_specs' or name.startswith('Table 4 -') or n in (1,16)):continue
    t=rel[sheet.get(REL+'id')];member=t.lstrip('/') if t.startswith('/') else 'xl/'+t
    with w.open(member) as inp:
     for _,row in ET.iterparse(inp,events=('end',)):
      if row.tag!=MAIN+'row':continue
      i=int(row.get('r'));v=['']*32
      for cell in row.findall(MAIN+'c'):
       letters=re.sub(r'\d+$','',cell.get('r'));num=0
       for letter in letters:num=num*26+ord(letter)-64
       if num>len(v):continue
       val=cell.findtext(MAIN+'v','');typ=cell.get('t')
       if typ=='s' and val:val=strings[int(val)]
       elif typ=='inlineStr':val=''.join(t.text or '' for t in cell.iter(MAIN+'t'))
       elif typ=='b' and val:val='True' if val=='1' else 'False'
       elif typ in (None,'n') and val and ('.' in val or 'E' in val.upper()):
        try:val=str(float(val))
        except ValueError:pass
       v[num-1]=val.strip()
      row.clear()
      add(data,name,n,i,v,a)
 return data
def add(data,name,n,i,v,a):
 if v[0] not in ('BRCA1','BRCA2'):return
 g=v[0]; base={'gene':g,'sheet':name,'row':i,'sha256':a['sha256']}
 if name=='Table9_BRCA12VCEP_specs' and v[1].startswith('c.'):
  data['T9'].append(dict(base,hgvsc=v[1],hgvsp=v[2],code=v[3],strength=v[4],pmids_mentioned=sorted(set(re.findall(r'\b\d{7,8}\b',' '.join(v[5:14])))),columns='A:N',review_cells=f'F{i}:N{i}'))
 elif name.startswith('Table 4 -'):
  for loc,alt,code,warn in [(4,5,3,2),(15,16,17,18)]:
   if v[loc].startswith('c.') and v[alt]:
    for allele in re.split(r'[,;/ ]+',v[alt]):
     if allele in ('A','C','G','T'):
      data['T4_splice'].append(dict(base,hgvsc=v[loc]+allele,exon=v[1],assigned=v[code],warning=bool(v[warn]),cell=chr(65+code)+str(i)))
  if v[1] and re.search(r'PTC|START|DEL |DUP ',v[10],re.I):
   data['T4_event'].append(dict(base,exon=v[1],event=v[10],code=v[9],qualifier=v[11],cell='J'+str(i)))
 elif n==1 and v[1]:
  data['ST1'].append(dict(base,exon=v[1],cdna=[v[2],v[3]],nmd_bound=v[10],strength=v[11],cell='L'+str(i)))
 elif n==16 and v[1].startswith('c.'):
  data['ST16'].append(dict(base,site=v[1],alts=v[2],pmids=sorted(set(re.findall(r'\b\d{7,8}\b',v[3]))),cell='D'+str(i)))
def atomic(path,body):
 fd,tmp=tempfile.mkstemp(dir=path.parent)
 try:
  with os.fdopen(fd,'wb') as f:f.write(body)
  os.replace(tmp,path)
 finally:
  if os.path.exists(tmp):os.unlink(tmp)
def ensure():
 global _prepared
 if _prepared is not None:return _directory
 if _directory is None:raise ValueError('请先从 start.py 启动并选择本地运行目录。')
 dest=_directory;dest.mkdir(parents=True,exist_ok=True)
 records={}
 for a in load(BASE/'attachments.lock.json'):
  if a['name'].endswith('.xlsx'):records.setdefault(a['name'],a)
 for filename,a in records.items():
  p=dest/filename
  if p.is_file() and sha(p.read_bytes())==a['sha256']:continue
  if _source:
   src=_source/filename
   if not src.is_file():raise ValueError('指定目录缺少冻结工作簿：'+filename)
   b=src.read_bytes()
  else:
   try:
    with urlopen(Request(a['url'],headers={'User-Agent':'BRCA-Review-Desk/0.1'}),timeout=35) as f:b=f.read(15000001)
   except OSError as err:raise ValueError('官方工作簿未取得；可下载三份原表，使用 --workbooks 指定所在目录。') from err
  if sha(b)!=a['sha256']:raise ValueError('来源工作簿哈希不符；未按新规则判读：'+filename)
  atomic(p,b)
 expected=load(BASE/'bundle.json')['row_index_utf8_sha256'];index=dest/'row-index.json.gz'
 if index.exists():
  try:
   raw=gzip.decompress(index.read_bytes())
   if sha(raw)==expected:_prepared=json.loads(raw);return dest
  except (ValueError,OSError,EOFError):pass
 d=parse_rows(dest,records);b=json.dumps(d,ensure_ascii=False,separators=(',',':')).encode('utf-8')
 if sha(b)!=expected:raise ValueError('本地行索引与验证基线不一致；本次停止。')
 atomic(index,gzip.compress(b,mtime=0));_prepared=d
 return dest
def rows():
 ensure()
 return _prepared

