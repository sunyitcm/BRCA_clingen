"""Offline public/constructed examples; no external calls or patient data."""
import copy,hashlib,json,os,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'source'))
from brca_review import adoption,erepo
from brca_review.criteria import scope,ledger,catalog,protein_domain,PROTEIN
from brca_review.desk import Store,render
from brca_review.examples import examples

def item(k='PS3',w='Strong',ref='PMID:33609447',unit='assay-A',inputs=None):
 return dict(code=k,strength=w,state='met',confirmed=True,reason='教学验收：已核对原材料与适用条件',source_refs=[ref],independent_units=[unit],inputs=inputs or {})
def req(s,ds,prev=None):
 return dict(gene=s['gene'],assembly=s['assembly'],caid=s['caid'],snapshot=s['id'],base_hash=s['report_sha256'],previous_revision=prev,reviewer='教学审核',decisions=ds,change_reason='',review_complete=False,independence_reviewed=False,significant_conflict=False)
def strata(absent=True):
 rows=[]
 for ds,assay in [('gnomAD_v2.1_non_cancer','exome'),('gnomAD_v3.1_non_cancer','genome')]:
  rows.append(dict(dataset=ds,assay=assay,AC=0 if absent else 24,AN=10000,ancestry='test-pop',qc_pass=True,non_founder=True,mean_depth=30,absence_verified_full_subset=True,statistic_verified=True,faf=0 if absent else 0.0012,source_ref='LOCAL:two-cohorts#counts'))
 return dict(strata=rows,both_subsets_checked=True,maximum_nonfounder_reviewed=True,known_pathogenic_founder=False)
def four():
 return [item(),item('PP3','Supporting','LOCAL:dbnsfp#noAF','model',dict(branch='protein',bayes_noaf=.34,predictor_verified=True)),
         item('PP4','Supporting','PMID:31853058','families',dict(lr=3.57,calibrated=True,lr_type='combined')),
         item('PM2','Supporting','LOCAL:two-cohorts#counts','population',strata())]
class Review(unittest.TestCase):
 def setUp(self):
  self.dir=tempfile.TemporaryDirectory();self.addCleanup(self.dir.cleanup);self.st=Store(Path(self.dir.name));self.s=self.st.find('demo-public-rc8-brca2-38',False)
 def test_cross_gene_or_multiple_rejected_before_network_or_rule_download(self):
  for g,v in [('BRCA1','chr13-32380106-G-T'),('BRCA2','NM_007294.4:c.181T>G'),('BRCA2','chr13-32380106-G-T,A'),('BRCA2','chr13-123-G-T\nchr13-124-G-T')]:
   with self.subTest(g=g,v=v),patch('brca_review.desk.online',side_effect=AssertionError),patch('brca_review.ruleprep.ensure',side_effect=AssertionError),self.assertRaises(ValueError):self.st.snapshot(v,'GRCh38',None,g)
  for b in ('GRCh37','GRCh38'):
   for g,ch in [('BRCA1','17'),('BRCA2','13')]:self.assertEqual(scope(g,b,'chr'+ch+'-42-A-C')['chrom'],ch)
 def test_19_and_gene_rules_and_evidence_once(self):
  for s in examples():
   r=s['report'];l=r['criteria_table'];rows=l['rows']
   self.assertEqual(len({i['code'] for i in rows}),19);self.assertEqual(sum(x['rule_enabled'] for x in rows),19)
   self.assertTrue(all(r['specification'].split()[0] in x['rule_url'] for x in rows))
   self.assertEqual(len([x for x in rows if x['conditional_strengths']]),4)
  z=self.s['report']['criteria_table']['rows'];self.assertEqual({x['code'] for x in z if x['priority']},{'PS3','PP3','PP4','PM2'})
  for row in z:
   if row['expert']:self.assertEqual(row['expert']['state'],'Met');self.assertFalse(row['decision'])
  self.assertNotEqual(PROTEIN['BRCA1'],PROTEIN['BRCA2'])
 def test_checked_only_is_not_a_decision(self):
  for x in (dict(item(),confirmed=False),dict(item(),reason=''),dict(item(),independent_units=[]),dict(item(),source_refs=[self.s['report']['erepo']['url']])):
   with self.subTest(x=x),self.assertRaises(ValueError):self.st.review(req(self.s,[x]))
  self.assertIsNone(adoption.latest(self.st,self.s['id']))
 def test_one_checkbox_subtotal_reopen_revision_retract(self):
  c=self.s['report']['criteria_table']['rows'];w=next(z for z in c if z['code']=='PS3')['candidates'][0]
  p=req(self.s,[item(ref=w['ref'])]);base=self.s['raw_report'];v=self.st.review(p);self.assertEqual(v['report']['points_accepted'],4);self.assertIsNone(v['report']['classification_draft'])
  self.assertEqual(self.st.find(self.s['id'])['report']['met'],['PS3']);self.assertEqual(v['raw_report'],base)
  with self.assertRaises(adoption.ReviewChanged):self.st.review(p)
  q=req(self.s,[],'r000001')
  with self.assertRaises(ValueError):self.st.review(q)
  q['change_reason']='撤回教学练习';v=self.st.review(q);self.assertIsNone(v['report']['points_accepted'])
  ep=Path(self.dir.name)/'reviews'/self.s['id'];self.assertTrue((ep/'r000001/report.html').is_file())
  self.assertEqual(adoption.latest(self.st,self.s['id'])['parent_sha256'],json.loads((ep/'r000001/review.json').read_text())['content_sha256'])
 def test_record_and_html_json_identical(self):
  v=self.st.review(req(self.s,[item()]));folder=Path(self.dir.name)/'reviews'/self.s['id']/'r000001'
  r=json.loads((folder/'report.json').read_bytes());boot=json.loads((folder/'report.html').read_text().split('<script id="boot" type="application/json">')[1].split('</script>')[0])
  self.assertEqual(r,boot['snapshots'][0]['report']);self.assertEqual(v['report'],r);self.assertIsNone(boot['token'])
 def test_cohort_and_lr_and_prediction_not_assumed(self):
  variants=[item('PM2','Supporting','LOCAL:two-cohorts#counts','pop',{}),item('PP4','Supporting','PMID:31853058','fh',{'lr':3.57,'lr_type':'segregation','calibrated':True}),item('PP4','Strong','PMID:31853058','fh',{'lr':3.57,'lr_type':'combined','calibrated':True}),item('PP3','Supporting','LOCAL:variant#predict','pred',dict(branch='protein',bayes_noaf=.29,predictor_verified=True))]
  q=strata();q['strata'][0]['mean_depth']=24;variants.append(item('PM2','Supporting','LOCAL:two-cohorts#counts','pop',q))
  q=strata();q['strata'][0]['dataset']='v4.1';variants.append(item('PM2','Supporting','LOCAL:two-cohorts#counts','pop',q))
  for x in variants:
   with self.subTest(code=x['code']),self.assertRaises(ValueError):self.st.review(req(self.s,[x]))
 def test_local_scores_only_after_dispositions_and_conflict_gate(self):
  ds=four();p=req(self.s,ds);res=adoption.evaluate(self.s,p,None)['result'];self.assertEqual(res['points'],7);self.assertIsNone(res['draft_class'])
  p.update(review_complete=True,independence_reviewed=True)
  with self.assertRaises(ValueError):adoption.evaluate(self.s,p,None)
  current=self.s['report']['criteria_table']['rows'];existing={z['code'] for z in ds}
  for row in current:
   if row['rule_enabled'] and row['code'] not in existing:ds.append(dict(code=row['code'],confirmed=True,state='not_applicable',source_refs=[row['rule_url']],reason='合成验收；逐项记录适用性，非真实判读'))
  res=adoption.evaluate(self.s,p,None)['result'];self.assertEqual(res['draft_class'],'Likely pathogenic');self.assertIsNone(res['local_final'])
  p['significant_conflict']=True;self.assertIsNone(adoption.evaluate(self.s,p,None)['result']['draft_class'])
 def test_duplicate_units_mutually_exclusive_null_and_unmapped(self):
  p=req(self.s,four()[:2]);p['decisions'][0]['independent_units']=['PMID:33609447'];p['decisions'][1]['independent_units']=['https://pubmed.ncbi.nlm.nih.gov/33609447/']
  with self.assertRaisesRegex(ValueError,'nonindependent'):adoption.evaluate(self.s,p,None)
  for rows in ([item(),item('BS3','Strong',unit='distinct')],[item('PP5','Supporting')],[item('PP4','Very Strong')]):
   with self.assertRaises(ValueError):adoption.evaluate(self.s,req(self.s,rows),None)
  with self.assertRaises(ValueError):adoption.evaluate(self.s,dict(req(self.s,[item()]),gene='BRCA1'),None)
 def test_expert_not_met_is_not_an_accepted_negative(self):
  r=copy.deepcopy(self.s['report']);r['erepo']['criteria']=[dict(code='BS2',state='Not Met',strength='Strong',issues=[])];row=next(x for x in ledger(r)['rows'] if x['code']=='BS2');self.assertFalse(row['priority']);self.assertEqual(row['local'],'unassessed')
 def test_strengths_identity_raw_case_disagreement(self):
  r=self.s['report'];uuid='11111111-2222-3333-4444-555555555555'
  z=dict(gene='BRCA2',mondoId=r['condition'],caId=r['caid'],moi='Autosomal dominant inheritance',ep=erepo.PANEL,affiliationId='50087',uuid=uuid,retracted=False,classification='Benign',metCodes=['BP5_VeryStrong'],unMetCodes=['PS3'],docVersion='1')
  self.assertEqual(erepo.select({'data':[z]},r)[1],'matched')
  self.assertEqual(erepo.select({'data':[z,dict(z)]},r)[1],'ambiguous_or_truncated')
  self.assertIsNone(erepo.select({'data':[dict(z,retracted=True),dict(z,mondoId='different')]},r)[0])
  a=erepo.interpret(z,None,r);self.assertEqual(a['criteria'][0]['strength'],'Very Strong')
  self.assertIsNone(a['criteria'][1]['strength']);self.assertEqual(erepo.interpret(z,{},r)['state'],'detail_identity_mismatch')
 def test_benign_standalone_not_a_negative_number(self):
  x=item('BA1','Stand Alone','LOCAL:two-cohorts#counts','AF',strata(False));res=adoption.evaluate(self.s,req(self.s,[x]),None)['result'];self.assertIsNone(res['points'])
 def test_open_and_save_without_external_api(self):
  with patch('urllib.request.OpenerDirector.open',side_effect=AssertionError):
   self.assertTrue(render(examples()));self.assertEqual(self.st.review(req(self.s,[item()]))['report']['met'],['PS3'])


class EmptyEndpoint(unittest.TestCase):
 def test_official_explicit_no_records_vs_plain_404(self):
  from brca_review._engine.brcalite.net import sha,writej
  r=copy.deepcopy(examples()[0]['report'])
  for explicit in [False,True]:
   with self.subTest(explicit=explicit),tempfile.TemporaryDirectory() as tmp:
    def get(folder,key,url,**kw):
     loc=folder/key;loc.mkdir(parents=True)
     raw=(json.dumps({'status':{'code':404,'msg':'No records were found for given query'},'metadata':{'rendered':{'by':erepo.BASE+'/summary/srvc'}}}).encode() if explicit else b'Not Found')
     (loc/'failed-1.body').write_bytes(raw);writej(loc/'receipt.json',{'url':url,'method':'GET','attempts':[{'status':404,'body_file':'failed-1.body','body_sha256':sha(raw)}]})
    ans=erepo.collect(Path(tmp),r,get)
    self.assertEqual(ans['state'],'explicit_no_records' if explicit else 'summary_unavailable');self.assertFalse(ans['criteria'])
 def test_local_strength_and_points_in_both_formats(self):
  with tempfile.TemporaryDirectory() as tmp:
   st=Store(Path(tmp));base=st.find('demo-public-rc8-brca2-38',False);v=st.review(req(base,[item()]));r=v['report'];x=next(t for t in r['criteria_table']['rows'] if t['code']=='PS3')
   self.assertEqual((x['accepted_strength'],x['contribution_points']),('Strong',4));self.assertTrue(r['is_training']);self.assertEqual(len(r['review_policy_sha256']),64)

class ProfileRegression(unittest.TestCase):
 def test_lr_text_levels_on_both_genes_and_boundaries(self):
  sources=['LOCAL:study#lr','LOCAL:assay#effect']
  for s,g in [(next(v for v in examples() if v['id']=='demo-public-rc8-brca1-37'),'BRCA1'),(next(v for v in examples() if v['id']=='demo-public-rc8-brca2-38'),'BRCA2')]:
   for k,typ,bound,near,expected in [('PP1','segregation','350','349.99',8),('PP4','combined','350','349.99',8),('BS4','segregation','0.00285','0.002851',-8),('BP5','combined','0.00285','0.002851',-8)]:
    with self.subTest(g=g,k=k):
     o=dict(lr=bound,lr_type=typ,calibrated=True,lr_source=sources[0],independence_checked=True,components_checked=True,not_cross_gene_coobservation=True,pedigree_context='single',alternative_cause_checked=True,phenocopy_reviewed=True,effect_kind='experimental',effect_source=sources[1])
     z=item(k,'Very Strong',sources[0],inputs=o);z['source_refs']=sources
     result=adoption.evaluate(s,req(s,[z]),None);self.assertEqual(result['result']['points'],expected);self.assertIsNone(result['result']['draft_class'])
     self.assertIn('V',result['decisions'][0]['strength']);self.assertTrue(result['decisions'][0]['strength_policy'])
     row=next(t for t in adoption.effective(s['report'],result)['criteria_table']['rows'] if t['code']==k);self.assertEqual(row['contribution_points'],expected)
     self.assertIn('Very Strong',row['allowed_strengths']);self.assertNotIn('Very Strong',row['structured_strengths']);self.assertIn(g,adoption.lr_policy()['sources'])
     for bad in [dict(lr=near),dict(lr_type='combined' if typ=='segregation' else 'segregation'),dict(lr_source='LOCAL:elsewhere#x'),dict(independence_checked=False),{'PP1':dict(alternative_cause_checked=False),'PP4':dict(components_checked=False),'BS4':dict(phenocopy_reviewed=False),'BP5':dict(not_cross_gene_coobservation=False)}[k]]:
      with self.assertRaises(ValueError):adoption.evaluate(s,req(s,[dict(z,inputs=dict(o,**bad))]),None)
     if k=='PP1':
      with self.assertRaises(ValueError):adoption.evaluate(s,req(s,[dict(z,inputs=dict(o,effect_source=''))]),None)
     with self.assertRaises(ValueError):adoption.evaluate(s,req(s,[dict(z,strength='Strong')]),None)
 def test_same_predictor_uses_selected_gene_profile_and_both_show_pm2(self):
  a=next(z for z in examples() if z['id']=='demo-public-rc8-brca1-37');b=next(z for z in examples() if z['id']=='demo-public-rc8-brca2-38')
  pp=item('PP3','Supporting','LOCAL:dbnsfp#noaf','model',dict(branch='protein',bayes_noaf='.29',predictor_verified=True))
  self.assertEqual(adoption.evaluate(a,req(a,[pp]),None)['result']['points'],1)
  with self.assertRaises(ValueError):adoption.evaluate(b,req(b,[pp]),None)
  for s in (a,b):
   prompt=s['report']['criteria_table']['population_prompt'];self.assertEqual(prompt['code'],'PM2');self.assertEqual(prompt['spec'],s['report']['specification'])
 def test_old_snapshot_preserves_original_bytes_and_creates_review_separately(self):
  import base64
  from brca_review._engine.brcalite.net import writej
  with tempfile.TemporaryDirectory() as tmp,patch('urllib.request.OpenerDirector.open',side_effect=AssertionError):
   st=Store(Path(tmp));s=copy.deepcopy(next(x for x in examples() if x['id']=='demo-public-rc8-brca2-38'))
   sid='s-20000101T000000Z-00000001';s.update(id=sid,teaching=False);r=s['report']
   for key in ('criteria_table','erepo','selected_gene','is_training'):r.pop(key,None)
   source=(json.dumps(r,ensure_ascii=False)+'\n').encode();s.update(report_sha256=hashlib.sha256(source).hexdigest(),raw_report=base64.b64encode(source).decode())
   folder=Path(tmp)/sid;folder.mkdir();writej(folder/'snapshot.json',s);(folder/'report.json').write_bytes(source);before=(folder/'snapshot.json').read_bytes()
   view=st.find(sid);self.assertEqual(view['report']['criteria_table']['expert_state'],'not_queried_in_snapshot');self.assertEqual(len(view['report']['criteria_table']['rows']),19)
   st.review(req(view,[item()]));self.assertEqual(st.find(sid)['report']['points_accepted'],4);self.assertEqual(st.find(sid,False)['raw_report'],s['raw_report'])
   self.assertEqual((folder/'snapshot.json').read_bytes(),before);self.assertEqual((folder/'report.json').read_bytes(),source)
   with patch('brca_review.adoption.effective',side_effect=AssertionError):self.assertEqual(st.find(sid)['report']['points_accepted'],4)
   sr=Path(tmp)/'reviews'/sid/'r000001'/'report.json';sr.write_bytes(sr.read_bytes()+b' ')
   with self.assertRaisesRegex(ValueError,'哈希'):st.find(sid)
 def test_malformed_erepo_payload_preserves_other_sources_and_does_not_become_absence(self):
  r=copy.deepcopy(examples()[-1]['report'])
  with tempfile.TemporaryDirectory() as td:
   x=erepo.collect(Path(td),r,lambda *args,**kw:b'bad JSON')
   self.assertEqual((x['state'],r['source_status']['erepo_summary']),('invalid_summary','invalid_json'))
   self.assertTrue(r['candidates']);self.assertFalse(x['criteria'])
