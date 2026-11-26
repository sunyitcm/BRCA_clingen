"""Offline: exercise release boundaries without contacting any data provider."""
import unittest,sys,tempfile,threading,urllib.request as u
from pathlib import Path
from unittest.mock import patch
from http.server import ThreadingHTTPServer
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'source'))
from brca_review.desk import Store,handler
from brca_review import ruleprep
from brca_review._engine.brcalite.net import receipt_headers
from brca_review._engine.brcalite.core import parse,review_prompts
from brca_review.model import compare
class PublicRelease(unittest.TestCase):
 def test_two_entry_network_boundary(self):
  with tempfile.TemporaryDirectory() as tmp:
   st=Store(Path(tmp));srv=ThreadingHTTPServer(('127.0.0.1',0),handler(st));t=threading.Thread(target=srv.serve_forever,daemon=True);t.start()
   try:
    url='http://127.0.0.1:'+str(srv.server_port)
    with patch('brca_review.desk.online',side_effect=AssertionError),patch('brca_review.ruleprep.urlopen',side_effect=AssertionError):
     def get(p):
      with u.urlopen(url+p) as f:return f.read().decode()
     self.assertIn('href="/snapshot"',get('/'));self.assertIn('href="/query"',get('/'))
     self.assertNotIn(st.token,get('/snapshot'));self.assertIn('CA001182',get('/snapshot'))
     self.assertIn(st.token,get('/query?assembly=GRCh37'));self.assertEqual(list(Path(tmp).iterdir()),[])
   finally:srv.shutdown();srv.server_close();t.join()
 def test_integrity_stops_before_index(self):
  with tempfile.TemporaryDirectory() as tmp:
   base=Path(tmp);bad=base/'workbooks';bad.mkdir()
   one=next(a['name'] for a in ruleprep.load(ruleprep.BASE/'attachments.lock.json') if a['name'].endswith('xlsx'))
   (bad/one).write_bytes(b'invalid');ruleprep.configure(base/'runs',bad)
   with self.assertRaises(ValueError):ruleprep.ensure()
   self.assertFalse(list(base.rglob('row-index.json.gz')))
 def test_different_builds_are_not_compared(self):
  a={'caid':'test','assembly':'GRCh37','status':'review_required'}
  self.assertEqual(compare(a,dict(a,assembly='GRCh38'))['rows'],[])
 def test_receipt_filters_headers(self):
  self.assertEqual(receipt_headers({'Content-Type':'application/json','Via':'example','X-Test':'example'}),[('Content-Type','application/json')])
 def test_prompt_is_precise_and_never_counts(self):
  def make(pos,tx='NM_007294.4:c.1A>G',gene='BRCA1'):
   return {'gene':gene,'specification':'GN092 v1.2','status':'review_required','facts':{'event':'SNV','transcript':{'hgvsc':tx,'consequence_terms':['missense_variant'],'protein_start':pos}}}
  for pos in (2,101,1391,1424,1650,1857,0,1864):self.assertEqual(review_prompts(make(pos)),[])
  for pos in (102,400,1390,1425,1649,1858):
   a=review_prompts(make(pos))[0];self.assertEqual(a['criterion'],'BP1_Strong');self.assertFalse(a['score_counted']);self.assertIsNone(a['spliceai']['max_delta'])
  self.assertEqual(review_prompts(make(400,'NM_007294:c.1A>G')),[])
  self.assertEqual(review_prompts(make(400,gene='BRCA2')),[])
 def test_assembly_affects_accession(self):
  self.assertTrue(parse('chr17-41258504-A-C','GRCh37')['hgvs'].startswith('NC_000017.10:'))
  self.assertTrue(parse('chr17-43106487-A-C','GRCh38')['hgvs'].startswith('NC_000017.11:'))
if __name__=='__main__':unittest.main()
