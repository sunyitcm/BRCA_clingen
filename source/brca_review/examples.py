"""Portable, limited public teaching data, separate from private acquisition snapshots."""
import base64,hashlib,json,gzip
from .model import A,version_rows
from .criteria import ledger
def examples():
 p=A/'teaching.json'
 if not p.exists():return []
 bases=json.loads(gzip.decompress((A/"teaching_bases_v2.json.gz").read_bytes()));out=[]
 for src in json.loads(p.read_bytes()):
  r=src['report'];r.pop('not_applicable_by_spec',None);r['criteria_table']=ledger(r)
  raw=base64.b64decode(bases[src['id']],validate=True)
  out.append(dict(id=src['id'],utc=src['utc'],assembly=r['input']['assembly'],original=r['input']['original'],gene=r['gene'],caid=r['caid'],status=r['status'],report=r,report_sha256=hashlib.sha256(raw).hexdigest(),raw_report=base64.b64encode(raw).decode(),versions=version_rows(r,{},None),raw={},blobs={},diff=None,teaching=True))
 return out
