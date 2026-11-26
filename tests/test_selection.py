import sys,unittest,tempfile,copy,json,hashlib,base64
from pathlib import Path
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'source'))
from brca_review import adoption,selection
from brca_review.desk import Store,render
from brca_review.criteria import ledger
from brca_review.examples import examples
from brca_review._engine.brcalite.net import writej
def rows(r):return {z['code']:z for z in ledger(r)['rows']}
def ask(s,selected,note='',prev=None):
 return dict(schema='brca.selection.v2',snapshot=s['id'],gene=s['gene'],caid=s['caid'],assembly=s['assembly'],base_hash=s['report_sha256'],previous_revision=prev,selected=selected,note=note)
class QuickReview(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.st=Store(Path(self.tmp.name))
  self.s=self.st.find('demo-public-rc8-brca2-38',False);self.r=self.s['report']
 def test_only_nineteen_per_gene_and_strengths_restricted(self):
  for s in examples():
   z=rows(s['report']);self.assertEqual(len(z),19)
   self.assertFalse(set(z)&{'PS2','PM1','PM4','PM6','PP2','PP5','BP2','BP3','BP6'})
   self.assertEqual(z['PM2']['allowed_strengths'],['Supporting']);self.assertEqual(z['PS3']['allowed_strengths'],['Strong'])
   if s['gene']=='BRCA2':self.assertEqual({k for k,i in z.items() if i['priority']},{'PS3','PM2','PP3','PP4'})
 def test_default_real_sources_and_neutral_expert_citations(self):
  z=rows(self.r);auto=[k for k,v in z.items() if v['brief']['default_selected']];self.assertEqual(set(auto),{'PS3','PP3'})
  self.assertIn('≥0.30',z['PP3']['brief']['result']);self.assertNotIn('BP4',z['PP3']['brief']['result'])
  self.assertIn('未列',z['PM2']['brief']['result']);self.assertIn('v4.1',z['PM2']['brief']['expert_note'])
  r=self.st.review(ask(self.s,auto))['report'];res=r['review_result'];self.assertEqual((r['points_accepted'],res['selected_combination']),(5,None));self.assertIsNone(r['classification_draft'])
  s=self.st.review(ask(self.s,auto+['PM2','PP4'],prev='r000001'));res=s['report']['review_result']
  self.assertEqual((res['points'],res['local_sources']['points'],res['quoted_sources']['points']),(7,5,2));self.assertEqual(res['selected_combination'],'Likely pathogenic')
  self.assertIsNone(s['report']['classification_local']);self.assertEqual(s['local_review']['note'],'');self.assertNotIn('reviewer',s['local_review'])
  self.assertEqual(s['report'],self.st.find(self.s['id'])['report'])
 def test_sparse_unsupported_overrides_and_removal(self):
  for req in [ask(self.s,['PM2','PM2']),ask(self.s,['PS2']),ask(self.s,['BP1']),dict(ask(self.s,['PS3']),gene='BRCA1')]:
   with self.assertRaises(ValueError):self.st.review(req)
  self.st.review(ask(self.s,['PS3'],'可选备注'));v=self.st.review(ask(self.s,[],prev='r000001'))
  self.assertEqual(v['report']['review_note'],'');self.assertIsNone(v['report']['points_accepted'])
  with self.assertRaises(adoption.ReviewChanged):self.st.review(ask(self.s,['PS3']))
 def test_row_reference_mismatch_or_None_not_a_functional_assignment(self):
  for override in [dict(sha256='0'*64),dict(hgvsc='c.7T>G'),dict(checked_assignment=None),dict(code='None'),dict(checked_assignment={'code':'PS3','strength':'Strong'})]:
   r=copy.deepcopy(self.r);r['candidates'][0]['source'].update(override);z=rows(r)['PS3'];self.assertFalse(z['brief']['default_selected']);self.assertNotEqual((z['brief']['option'] or {}).get('kind'),'workbook')
  r=copy.deepcopy(self.r);r['candidates'].append(copy.deepcopy(r['candidates'][0]));self.assertFalse(rows(r)['PS3']['brief']['default_selected'])
 def test_prediction_dataset_and_identity_gate(self):
  for v in ['0.29','1.1','NaN',None]:
   r=copy.deepcopy(self.r);r['facts']['ldh_brca']['BayesDel_nsfp33a_noAF']=v;self.assertFalse(rows(r)['PP3']['brief']['default_selected'])
  for field,val in [('linked_allele','wrong'),('version_source',{})]:
   r=copy.deepcopy(self.r);r['facts']['ldh_brca'][field]=val;self.assertFalse(rows(r)['PP3']['brief']['default_selected'])
  r=copy.deepcopy(self.r);r['facts']['transcript']['protein_start']=200;self.assertFalse(rows(r)['PP3']['brief']['default_selected'])
  r=copy.deepcopy(self.r);r['rules_bundle_sha256']='wrong';self.assertFalse(any(z['brief']['option'] for z in rows(r).values()))
 def test_expert_unsupported_strength_incomplete_or_disagrees_not_adopted(self):
  for state in ['summary_only','detail_identity_mismatch','unavailable']:
   r=copy.deepcopy(self.r);r['erepo']['state']=state;self.assertIsNone(rows(r)['PM2']['brief']['option'])
  r=copy.deepcopy(self.r);q=next(x for x in r['erepo']['criteria'] if x['code']=='PM2');q['strength']='Strong';z=rows(r)['PM2'];self.assertIsNone(z['expert']);self.assertIsNone(z['brief']['option'])
  r=copy.deepcopy(self.r);q=next(x for x in r['erepo']['criteria'] if x['code']=='PS3');q['state']='Not Met';z=rows(r)['PS3'];self.assertEqual(z['brief']['status'],'conflict');self.assertFalse(z['brief']['default_selected'])
 def test_external_case_is_not_local_frequency_and_proven_population_can_auto_apply(self):
  # Positive and missing-field population tests are synthetic, never part of shipped teaching cases.
  src='LOCAL:synthetic#aggregate';datasets=[('gnomAD_v2.1_non_cancer','exome'),('gnomAD_v3.1_non_cancer','genome')]
  ss=[dict(dataset=d,assay=a,AC=0,AN=12000,mean_depth=25,qc_pass=True,non_founder=True,outbred=True,ancestry='tested outbred controls',source_ref=src,absence_verified_full_subset=True) for d,a in datasets]
  d={'schema':'brca.population.v1','records':[{'gene':self.r['gene'],'caid':self.r['caid'],'hgvsc':self.r['facts']['transcript']['hgvsc'],'inputs':dict(strata=ss,both_subsets_checked=True)}]}
  f=Path(self.tmp.name)/'counts.json';writej(f,d);r=copy.deepcopy(self.r);selection.population_from_file(f,r);z=rows(r)['PM2']['brief'];self.assertEqual(z['option']['kind'],'population');self.assertTrue(z['default_selected'])
  self.assertEqual(z['option']['refs'][0]['file_sha256'],hashlib.sha256(f.read_bytes()).hexdigest())
  for field,val in [('AC',1),('AN',0),('AC','0'),('mean_depth',24.99),('outbred',False),('source_ref','')]:
   q=copy.deepcopy(r);q['facts']['population_evidence']['inputs']['strata'][0][field]=val;self.assertFalse(rows(q)['PM2']['brief']['default_selected'])
  q=copy.deepcopy(r);q['facts']['event']='DEL';self.assertFalse(rows(q)['PM2']['brief']['default_selected'])
  q=copy.deepcopy(r);q['facts']['population_evidence']['inputs']['strata'].pop();self.assertFalse(rows(q)['PM2']['brief']['default_selected'])
  q=copy.deepcopy(self.r);q['caid']='CA00000';selection.population_from_file(f,q);self.assertNotIn('population_evidence',q['facts'])
 def test_clean_exports_no_form_dependencies(self):
  with patch('urllib.request.OpenerDirector.open',side_effect=AssertionError):
   v=self.st.review(ask(self.s,['PS3']));folder=Path(self.tmp.name)/'reviews'/self.s['id']/'r000001';r=json.loads((folder/'report.json').read_bytes());outer=json.loads((folder/'report.html').read_text().split('<script id="boot" type="application/json">')[1].split('</script>')[0]);self.assertEqual(r,outer['snapshots'][0]['report']);self.assertIsNone(outer['token']);self.assertEqual(len(r['criteria_table']['rows']),19);self.assertEqual(r['review_result']['points'],4);self.assertEqual(v['raw_report'],self.s['raw_report'])
