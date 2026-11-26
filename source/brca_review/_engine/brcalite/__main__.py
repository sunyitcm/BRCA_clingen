"""Verify the frozen public profile; app entry point: python start.py."""
from .core import R
from .net import getj,sha
def verify():
 for name in ('engine-0.2.json','bundle.json'):
  for path,h in getj(R/'rules'/name)['files'].items():
   if sha((R/path).read_bytes())!=h:raise ValueError('规则文件哈希不符：'+path)
if __name__=='__main__':verify()
