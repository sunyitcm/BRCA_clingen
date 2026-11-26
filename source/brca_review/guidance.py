"""Version-specific, attributed Chinese reading aids. No text retrieved here is variant evidence."""
import json
from pathlib import Path
from functools import lru_cache
@lru_cache(maxsize=1)
def catalog_zh():
 return json.loads((Path(__file__).parent/'assets/criteria_zh_v12.json').read_bytes())
def source(g):
 obj=catalog_zh();doc=obj['genes'][g]
 return {k:doc[k] for k in ('spec','document_sha256','document_url','publication_url')}|{'translation':obj['translation']}
def criterion(g,z):
 c=catalog_zh()['genes'][g]['criteria'][z['code']]
 if c['rule_key']!=z['rule_field'] or set(c['grades'])!=set(z['allowed_strengths']):raise ValueError('本基因中文条文与冻结规范不一致。')
 return c
