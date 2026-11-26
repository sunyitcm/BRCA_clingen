"""The BRCA2 teaching fixture drives UI contracts; injected human choices are software examples only."""
import copy,hashlib,json,sys,tempfile,unittest,base64
from pathlib import Path
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root/'source'))
from brca_review import guidance,adoption,selection
from brca_review.criteria import ledger
from brca_review.desk import Store,render
from brca_review.examples import examples

def request(s,keys,grades=None,revision=None):
 d=dict(schema=selection.SCHEMA,snapshot=s['id'],base_hash=s['report_sha256'],gene=s['gene'],assembly=s['assembly'],caid=s['caid'],selected=keys,previous_revision=revision,note='')
 if grades is not None:d['manual']=grades
 return d

class ManualSelection(unittest.TestCase):
 def setUp(self):
  temp=tempfile.TemporaryDirectory();self.addCleanup(temp.cleanup);self.store=Store(Path(temp.name));self.s=self.store.find('demo-public-rc8-brca2-38',False)
 def test_every_shown_code_has_official_field_and_permitted_translations(self):
  for s in examples():
   table=ledger(s['report']);doc=table['guideline_source'];self.assertEqual(doc['spec'],s['report']['specification']);self.assertEqual(len(doc['document_sha256']),64);self.assertTrue(doc['document_url'].startswith('https://zenodo.org/records/'))
   self.assertEqual(len(table['rows']),19)
   for z in table['rows']:
    txt=z['guideline_zh'];self.assertEqual(txt['rule_key'],z['rule_field']);self.assertEqual(set(txt['grades']),set(z['allowed_strengths']));self.assertIn('/ld/RuleSet/0/ld/CriteriaCode/',txt['instructions']['pointer']);self.assertGreater(len(txt['zh']),20)
    for k,d in txt['grades'].items():
     self.assertRegex(d['pointer'],r'/strengthDescriptor/\d+/text$');self.assertEqual(len(d['text_sha256']),64);self.assertTrue(d['zh'])
     self.assertEqual(d['conditional_in_text'],k=='Very Strong' and z['code'] in ('PP1','PP4','BS4','BP5'))
 def test_same_revision_stores_and_restores_three_origins_without_inventing_sources(self):
  s=self.s;auto=request(s,['PS3','PP3']);initial=self.store.review(auto);self.assertEqual(initial['report']['points_accepted'],5)
  req=request(s,['PS3','PP3','PM2','PP1'],{'PP1':'Supporting'},'r000001');v=self.store.review(req);out=v['report']['review_result']
  self.assertEqual([out[k]['points'] for k in ('local_sources','quoted_sources','manual_sources')],[5,1,1]);self.assertEqual(out['points'],7);self.assertEqual(out['manual_codes'],['PP1']);self.assertIsNone(v['report']['classification_local'])
  d=next(d for d in v['local_review']['decisions'] if d['origin']=='manual');self.assertEqual((d['contribution_points'],d['source_refs'],d['independent_units']),(1,[],[]));self.assertTrue(d['rule_ref']['pointer'])
  archived=self.store.dir/'reviews'/s['id']/'r000002';body=(archived/'report.json').read_bytes();js=json.loads((archived/'report.html').read_text().split('<script id="boot" type="application/json">')[1].split('</script>')[0])
  self.assertEqual(json.loads(body),js['snapshots'][0]['report']);self.assertEqual(hashlib.sha256(body).hexdigest(),v['local_review']['saved_report_sha256']);self.assertEqual(v['report'],self.store.find(s['id'])['report'])
  rm=self.store.review(request(s,['PS3'],{},'r000002'));self.assertEqual(rm['report']['points_accepted'],4);self.assertEqual(rm['report']['review_result']['manual_codes'],[]);self.assertEqual((archived/'report.json').read_bytes(),body)
 def test_explicit_grades_no_default_no_out_of_scope_no_implicit_negatives(self):
  for p in [request(self.s,['PVS1']),request(self.s,['PVS1'],{'PVS1':'Stand Alone'}),request(self.s,['PP5'],{'PP5':'Supporting'}),request(self.s,['PP1'],{'BS4':'Strong'}),request(self.s,['PP1'],['Strong']),request(self.s,['PS3'],{'PS3':'Strong'}),request(self.s,['BP1'],{'BP1':'Strong'})]:
   with self.assertRaises(ValueError):selection.build(self.s,p,None)
  self.assertEqual(selection.build(self.s,request(self.s,[]),None)['result']['points'],None)
 def test_signed_manual_grades_conditional_grade_and_standalone(self):
  for k,w,num,stand in [('PVS1','Moderate',2,False),('BS2','Supporting',-1,False),('PP1','Very Strong',8,False),('BS4','Very Strong',-8,False),('BA1','Stand Alone',None,True)]:
   v=selection.build(self.s,request(self.s,[k],{k:w}),None);row=v['decisions'][0];out=v['result'];self.assertEqual(row['contribution_points'],num);self.assertEqual(out['stand_alone'],stand);self.assertEqual(out['points'],num)
   self.assertEqual(row['rule_ref']['conditional_in_text'],w=='Very Strong')
 def test_no_combined_conflicting_or_overlap_score(self):
  for keys,grade in [(['PP3','PVS1'],{'PVS1':'Strong'}),(['PM2','BA1'],{'BA1':'Stand Alone'}),(['PS3','BS3'],{'BS3':'Strong'})]:
   with self.assertRaises(ValueError):selection.build(self.s,request(self.s,keys,grade),None)
  for keys,grade in [(['PP4','PP1'],{'PP1':'Strong'}),(['BP7'],{'BP7':'Supporting'}),(['PS3','BS2'],{'BS2':'Supporting'})]:
   out=selection.build(self.s,request(self.s,keys,grade),None)['result'];self.assertTrue(out['issues']);self.assertIsNone(out['points']);self.assertIsNone(out['selected_combination'])
 def test_no_unverified_identity_or_wrong_gene_for_manual(self):
  for k,v in [('caid',None),('gene','BRCA1'),('rules_bundle_sha256','ff')]:
   s=copy.deepcopy(self.s);s['report'][k]=v
   with self.assertRaises(ValueError):selection.build(s,request(s,['PP1'],{'PP1':'Strong'}),None)
 def test_existing_rc2_base_bytes_remain_unchanged_and_new_guidance_distinct(self):
  s=self.s;raw=base64.b64decode(s['raw_report']);self.assertEqual(hashlib.sha256(raw).hexdigest(),s['report_sha256'])
  self.assertNotIn('guideline_zh',json.loads(raw)['criteria_table']['rows'][0]);self.assertIn('guideline_zh',s['report']['criteria_table']['rows'][0])
  # Read the previous contract without a 'manual' key; the same snapshot accepts a subsequent rc3 revision.
  before=self.store.review(request(s,['PS3']));self.assertEqual(before['report']['points_accepted'],4)
  after=self.store.review(request(s,['PS3','PP1'],{'PP1':'Moderate'},'r000001'));self.assertEqual(after['report']['points_accepted'],6);self.assertEqual(after['local_review']['base_sha256'],s['report_sha256'])
