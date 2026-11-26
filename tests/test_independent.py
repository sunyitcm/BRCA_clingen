"""Synthetic alleles and service responses; these tests never use a public variant or the network."""
import json,sys,tempfile,unittest
from copy import deepcopy
from pathlib import Path
from unittest.mock import patch
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root/'source'))
from brca_review import criteria, selection, routes
from brca_review._engine.brcalite import core
from brca_review._engine.brcalite.net import writej,sha
from brca_review.desk import Store
from brca_review.examples import examples
from brca_review.model import pack
from brca_review.routes import online,checked_ref,gnomad_record,mapped

def synth(allele,v):
    b=v['assembly'];chrom=v['chrom'];acc=core.NC[b][chrom]
    tc='NM_007294.4:c.10A>G'
    return {'@id':'https://reg.clinicalgenome.org/allele/'+allele,'transcriptAlleles':[{'geneSymbol':'BRCA1','hgvs':[tc],'proteinEffect':{'hgvs':'NP_009225.1:p.Lys4Glu'}}],
       'genomicAlleles':[{'referenceGenome':b,'coordinates':[dict(start=v['pos']-1,end=v['pos'],referenceAllele=v['ref'],allele=v['alt'])],'hgvs':[v['hgvs']]}]}
class Independent(unittest.TestCase):
 def setUp(self):
    temp=tempfile.TemporaryDirectory();self.addCleanup(temp.cleanup);self.dir=Path(temp.name);self.v=core.parse('17-101-A-G','GRCh37')
 def fake(self,answers):
    def fetch(folder,label,url,payload=None,kind='json',**kw):
        target=folder/label;target.mkdir(parents=True);v=answers.get(label);body=v if isinstance(v,bytes) else (json.dumps(v).encode() if v is not None else None)
        if payload is not None:
            data=json.dumps(payload,separators=(',',':')).encode();(target/'request.json').write_bytes(data);req=sha(data)
        else:req=None
        ar=dict(status=200 if body else 500,ok=body is not None)
        if body:
            (target/'body').write_bytes(body);ar.update(body_file='body',body_sha256=sha(body),bytes=len(body))
        writej(target/'receipt.json',dict(label=label,method='GET' if payload is None else 'POST',url=url,at_utc='2000-01-01T00:00:00+00:00',ok=body is not None,request_sha256=req,attempts=[ar]))
        return body
    return fetch
 def run_fetch(self,ans,out='snap'):
    return online(self.dir/out,self.v,True,workbooks=False,fetcher=self.fake(ans))
 def test_all_down_ledger_and_persistent_pending(self):
    st=Store(self.dir/'store')
    def failed(folder,v,details,reference,workbooks):return online(folder,v,details,reference,workbooks=False,fetcher=self.fake({}))
    with patch('brca_review.desk.online',side_effect=failed),patch('brca_review.ruleprep.ensure',side_effect=ValueError('unavailable')),patch('brca_review.desk.collect_portal'):
        s=st.snapshot(self.v['original'],'GRCh37',None,'BRCA1')
    tab=s['report']['criteria_table'];self.assertEqual(len(tab['rows']),19);self.assertTrue(tab['all_evidence_missing'] and tab['draft_selectable']);self.assertFalse(tab['identity_verified'])
    self.assertNotIn('transcript',s['report']['facts']);self.assertFalse(any(z['brief']['default_selected'] for z in tab['rows']))
    p=dict(schema=selection.SCHEMA,selected=['PP1','PM2'],manual={},base_hash=s['report_sha256'],previous_revision=None,gene='BRCA1',assembly='GRCh37',caid=None,snapshot=s['id'])
    saved=st.review(p);r=saved['report'];self.assertEqual(r['met'],[]);self.assertIsNone(r['points_accepted']);self.assertIsNone(r['review_result']['selected_combination']);self.assertEqual(set(r['review_result']['pending_codes']),set(p['selected']));self.assertEqual(set(r['not_evaluated']),{u['code'] for u in tab['rows']})
    self.assertEqual(st.find(s['id'])['report'],r)
 def test_refseq_plus_car_work_even_if_both_veps_and_excel_down(self):
    ans={'car':synth('CA11',self.v),'refseq_ref':b'>NC_000017.10:101-101 synthetic\nA\n'}
    r=self.run_fetch(ans);self.assertTrue(criteria.verified(r,'BRCA1'));self.assertEqual(r['facts']['transcript']['hgvsp'],'NP_009225.1:p.Lys4Glu');self.assertEqual(r['facts']['transcript_sources']['hgvsp'],'car');self.assertEqual(r['source_status']['ldh_med'],'http_error');self.assertEqual(r['source_status']['workbook'],'unavailable')
 def test_mismatch_is_held_but_source_observations_are_preserved(self):
    ans={'car':synth('CA11',self.v),'refseq_ref':b'>NC_000017.10:101-101 synthetic\nC\n'}
    r=self.run_fetch(ans);self.assertFalse(criteria.verified(r,'BRCA1'));self.assertEqual(r['facts']['transcript']['hgvsp'],'NP_009225.1:p.Lys4Glu');self.assertEqual(r['source_status']['ref_check'],'mismatch');self.assertIn('gno_full37',r['source_files']);self.assertIn('ldh_med',r['source_files'])
 def test_graphql_error_null_zero_and_wrong_allele_are_distinct(self):
    q='gno_nc2';fun=lambda d:gnomad_record(q,self.v,d,'BRCA1')
    self.assertEqual(fun({'errors':[{'message':'Variant not found'}],'data':{'v':None}})[1],'no_variant_record')
    self.assertEqual(fun({'errors':[{'message':'rate limited'}],'data':{'v':None}})[1],'graphql_error')
    base=dict(variant_id='17-101-A-G',reference_genome='GRCh37',exome={'ac':0,'an':420,'filters':[]},genome=None)
    rec,st=fun({'data':{'v':base}});self.assertEqual(st,'variant_record');self.assertEqual(rec['cohorts']['exome']['af'],0);self.assertIsNone(rec['cohorts']['genome'])
    self.assertEqual(fun({'data':{'v':dict(base,variant_id='17-101-A-C')}})[1],'allele_not_verified');self.assertEqual(fun({'data':{'v':None}})[1],'allele_not_verified')
 def test_reference_accession_span_and_mapping_validated(self):
    self.assertIsNone(checked_ref(self.v,b'>NC_000017.11:101-101 GRCh38\nA'));self.assertIsNone(checked_ref(self.v,b'>NC_000017.10:102-102\nA'));c=synth('CA11',self.v)
    self.assertEqual(mapped(c,'BRCA1','GRCh37')['hgvs'],self.v['hgvs']);c['genomicAlleles']*=2;self.assertIsNone(mapped(c,'BRCA1','GRCh37'))
 def test_non_cancer_present_vs_full_only_and_zero_score(self):
    r=self.run_fetch({'car':synth('CA11',self.v),'refseq_ref':b'>NC_000017.10:101-101\nA\n'});f=r['facts'];f['ldh_brca']={'BayesDel_nsfp33a_noAF':0}
    row={'state':'variant_record','label':'synthetic','dataset':'gnomad_r2_1_non_cancer','assembly':'GRCh37','cohorts':{'exome':{'ac':1,'an':123,'filters':[],'af':1/123,'populations':[dict(id='sas',ac=1,an=12)],'faf95':{'popmax':None}}}}
    f['gnomad']={'gno_full37':deepcopy(row)};by=lambda:{z['code']:z for z in criteria.ledger(r)['rows']}
    self.assertFalse(by()['PM2']['evidence_view']['blocks_automatic_scoring']);f['gnomad']['gno_nc2']=row
    t=by();self.assertTrue(t['PM2']['evidence_view']['blocks_automatic_scoring']);self.assertIsNone(t['PM2']['brief']['option']);self.assertIn('总 AF',t['BS1']['evidence_view']['comparison']);self.assertTrue(any(a['value']=='0' for a in t['PP3']['evidence_view']['observations']))
    s={'report':r,'report_sha256':'base','assembly':'GRCh37','id':'test'}
    p=dict(selected=['PM2'],base_hash='base',gene='BRCA1',assembly='GRCh37',caid='CA11',snapshot='test',manual={})
    rev=selection.build(s,p,None);self.assertIsNone(rev['result']['points']);self.assertEqual(rev['decisions'][0]['state'],'pending')
    row['cohorts']['exome']['populations'][0]['id']='asj';self.assertFalse(by()['PM2']['evidence_view']['blocks_automatic_scoring'])
 def test_faf_is_compared_only_as_faf_with_ancestry_and_quality(self):
    r=self.run_fetch({'car':synth('CA11',self.v),'refseq_ref':b'>NC_000017.10:101-101\nA\n'});k='gno_nc2'
    rec=dict(state='variant_record',label='test',cohorts={'exome':{'ac':3,'an':800,'filters':[],'af':3/800,'faf95':{'popmax':0.0001,'popmax_population':'nfe'}}});r['facts']['gnomad']={k:rec}
    t={z['code']:z for z in criteria.ledger(r)['rows']};self.assertIn('Supporting 区间',t['BS1']['evidence_view']['comparison']);self.assertIsNone(t['BS1']['brief']['option']);self.assertIn('未超过 0.001',t['BA1']['evidence_view']['comparison'])
 def test_car_loss_gnomad_variant_still_retrieved_without_refseq_fabrication(self):
    v=self.v;row=dict(variant_id='17-101-A-G',reference_genome='GRCh37',exome=dict(ac=2,an=100),genome=None,transcript_consequences=[dict(gene_symbol='BRCA1',transcript_id='ENST00000',hgvsp='p.Lys4Glu')])
    r=self.run_fetch({'gno_full37':{'data':{'v':row}}});self.assertFalse(criteria.verified(r,'BRCA1'));self.assertEqual(r['facts']['gnomad']['gno_full37']['transcripts'][0]['hgvsp'],'p.Lys4Glu');self.assertNotIn('transcript',r['facts']);self.assertEqual(r['source_status']['ldh_med'],'needs_caid')
 def test_malformed_optional_services_do_not_drop_other_sources(self):
    val=synth('CA11',self.v);val['externalRecords']=['schema-changed']
    r=self.run_fetch({'car':val,'refseq_ref':b'>NC_000017.10:101-101\nA\n','info':['bad'],'gno_full37':{'data':['changed']}})
    self.assertEqual(r['source_status']['clinvar'],'invalid_fields');self.assertEqual(r['source_status']['gno_full37'],'graphql_error');self.assertTrue(criteria.verified(r,'BRCA1'));self.assertNotIn('ensembl_metadata',r['facts'])
 def test_versioned_input_marker_is_retained_for_failed_sources(self):
    r=self.run_fetch({});self.assertEqual(r['source_status']['ref_check'],'unavailable');self.assertEqual(r['source_status']['refseq_ref'],'http_error');self.assertEqual(r['source_status']['gno_nc3'],'no_verified_query_coordinate')
    self.assertEqual(r['facts']['gnomad']['gno_nc2']['state'],'http_error')
