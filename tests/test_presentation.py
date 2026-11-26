import re,sys,unittest
from copy import deepcopy
from pathlib import Path
from urllib.parse import urlparse,parse_qs
from xml.etree import ElementTree as ET
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'source'))
from brca_review.criteria import scope,ledger
from brca_review.examples import examples
from brca_review.presentation import reading_guide,view_only
from brca_review.routes import frequency_entity,scv_citations,gnomad_record

class ReviewPresentation(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.base=next(x for x in examples() if x['id']=='demo-public-rc8-brca2-38')
 def get(self):return deepcopy(self.base['report'])
 def test_one_vcf_line_discards_metadata(self):
    for sep in (' ','\t'):
        v=scope('BRCA2','GRCh38',sep.join(('13','32380106','local-id','G','T','.','PASS','private-metadata','GT','0/1'))+'\n')
        self.assertEqual(v['original'],'13-32380106-G-T');self.assertNotIn('private',str(v))
    for val in ('13 1 . A G,T','13 1 . A T\n13 2 . T A','NM_000059.4:p.Asp3073Tyr','17 1 . T G'):
        with self.assertRaises(ValueError):scope('BRCA2','GRCh38',val)
    self.assertEqual(scope('BRCA2','GRCh37','NM_000059.4:c.9217G>T')['kind'],'hgvsc')
 def test_same_rules_distinct_sources(self):
    t=ledger(self.get());rows={z['code']:z for z in t['rows']};self.assertEqual(len(rows),19)
    self.assertEqual((rows['PS3']['assessment']['state'],rows['PP4']['assessment']['state']),('supported','external_met'))
    self.assertIn('workbook',{a['key'] for a in rows['PS3']['assessment']['sources']})
    self.assertEqual(rows['PM2']['assessment']['option_origin'],'expert_citation')
    self.assertTrue(all(z['assessment']['publication'].endswith('GN097?version=1.2') for z in rows.values()))
 def test_literature_no_papers_still_searchable(self):
    r=self.get();r['erepo']={'state':'explicit_no_records'};r['candidates']=[];r['facts'].pop('workbook_rows_without_assignment',None)
    guide=reading_guide(r,True);self.assertEqual(guide['records'],{})
    for link in guide['links'].values():
        term=parse_qs(urlparse(link).query)['term'][0];self.assertIn('BRCA2',term);self.assertIn('c.9217G>T',term);self.assertIn('D3073Y',term)
    t=r['facts']['transcript'];t['hgvsc']='NM_000059.4:c.68+1G>C';t['hgvsp']=None
    self.assertTrue(any('%2B1G%3EC' in u for u in reading_guide(r,True)['links'].values()))
    self.assertFalse(any('term=' in u for u in reading_guide(r,False)['links'].values()))
 def test_scv_reference_is_lead_not_disease_bibliography(self):
    x=ET.fromstring('<VariationArchive><ClassifiedRecord><ClinicalAssertionList><ClinicalAssertion><ClinVarAccession Accession="SCV001" Version="3"/><Classification><GermlineClassification>test</GermlineClassification><Citation><ID Source="PubMed">12345</ID></Citation></Classification><TraitSet><Trait><Citation><ID Source="PubMed">34567</ID></Citation></Trait></TraitSet></ClinicalAssertion></ClinicalAssertionList></ClassifiedRecord></VariationArchive>')
    out=scv_citations(x);self.assertEqual([(q['pmid'],q['scv']) for q in out],[('12345','SCV001.3')])
    r=self.get();r['facts']['clinvar_record']={'germline_citations':out}
    self.assertEqual(ledger(r)['literature']['records']['12345'][0]['origin'],'clinvar_scv')
    self.assertEqual([z['brief']['option'] for z in ledger(r)['rows']],[z['brief']['option'] for z in ledger(self.get())['rows']])
 def test_ldh_true_zero_not_missing_or_wrong_allele(self):
    obj={'locus':{'contig':'chr13','position':32380106,'reference_genome':'GRCh38'},'alleles':['G','T'],'gks_va_freq':{'focusAlleleCount':0,'locusAlleleCount':100,'alleleFrequency':0},'noVariant':True}
    expected={'NC_000013.11:g.32380106G>T'}
    self.assertEqual(frequency_entity('ldh_genome',obj,expected)['AC'],0);self.assertIsNone(frequency_entity('ldh_genome',obj,{'NC_000013.11:g.32380106G>A'}))
    for bad in (None,1):obj['gks_va_freq']['focusAlleleCount']=bad;self.assertIsNone(frequency_entity('ldh_genome',obj,expected))
 def test_service_error_never_equals_absence(self):
    v=scope('BRCA2','GRCh38','13-32380106-G-T');rec,state=gnomad_record('gno_nc3',v,{'data':{'v':None},'errors':[{'message':'Service overloaded'}]},'BRCA2')
    self.assertEqual((state,rec['service_message']),('graphql_error','Service overloaded'))
    r=self.get();r['erepo']={'state':'explicit_no_records'};r['facts']['gnomad']={'gno_nc3':dict(rec,state=state)};row=next(z for z in ledger(r)['rows'] if z['code']=='PM2')
    self.assertIsNone(row['brief'].get('option'));self.assertTrue(any(x['kind']=='source_failure' for x in row['evidence_view']['observations']))
 def test_old_revision_readonly_overlay(self):
    snap=deepcopy(self.base);tab=snap['report']['criteria_table'];tab.pop('literature')
    for z in tab['rows']:z.pop('assessment')
    before=deepcopy(snap);view=view_only(snap);self.assertEqual(snap,before);self.assertEqual(view['report'],before['report']);self.assertEqual(len(view['view_overlay']['assessments']),19)
 def test_excluded_pmid_is_visible_and_cannot_become_a_second_vote(self):
    from brca_review.erepo import literature_mentions
    text='Two studies (PMIDs: 38417439, 39779857). Also PMID 39779848; PMID 39779848 was not considered.'
    refs=literature_mentions([{'value':text,'pointer':'/line'}])
    self.assertEqual({x['pmid']:x['role'] for x in refs},{'38417439':'mentioned','39779857':'mentioned','39779848':'explicitly_excluded'})
    table=ledger(self.get());x=table['literature']['records']['39779848']
    self.assertTrue(any(o.get('use')=='explicitly_excluded' for o in x));self.assertEqual(next(z for z in table['rows'] if z['code']=='PS3')['brief']['option']['strength'],'Strong')
