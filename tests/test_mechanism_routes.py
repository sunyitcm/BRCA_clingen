"""Synthetic mechanism, eligibility and provenance regressions; fixtures do not represent patients."""
import sys,unittest
from copy import deepcopy
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'source'))
from brca_review.mechanism import classify,ptc_hint,interpretation,NP
from brca_review.routes import padded_indel
from brca_review.exchange import detail,text_lr
from brca_review.lineage import clinical_lr
from brca_review.criteria import GENES,ledger
from brca_review.selection import build
from brca_review.examples import examples
def tx(g,c,p=None,terms=(),**extras):
 return dict(transcript_id=GENES[g][0],hgvsc=GENES[g][0]+':c.'+c,hgvsp=NP[g]+':p.'+p if p else None,gene_symbol=g,consequence_terms=list(terms),**extras)
def interp(t,k='PP3',value='-0.10'):
 g=t['gene_symbol'];out=[];r={'gene':g,'facts':{'transcript':t,'ldh_brca':{'BayesDel_nsfp33a_noAF':value,'version_source':{'releaseName':12}}}}
 return interpretation(r,k,lambda *a:out.append(a)),out
def family(r,k='PP1',text=''):
 z={'code':k,'expert':{'text':[{'value':text,'pointer':'/line/1'}]}};out=[]
 return clinical_lr(r,z,lambda *a:out.append(a)),out
class AlleleMechanisms(unittest.TestCase):
 def test_fs_vs_missense_even_if_domain_and_high_score(self):
  for g in NP:
   t=tx(g,'55_56del','(Val19Glyfs*9)')
   m=classify(t,g);self.assertEqual((m['category'],m['aa'],m['stop_aa'],m['domain']),('frameshift',19,27,'inside'))
   for val in ('0.98','-0.12','0'):
    (state,_,held),obs=interp(t,value=val);self.assertTrue(held);self.assertEqual(state,'not_applicable')
    self.assertTrue(any('背景' in x[0] and x[1].startswith(val+' ·') for x in obs))
 def test_inframe_domain_bounds_and_actual_number(self):
  for g in NP:
   t=tx(g,'55_57del','Val19del');m=classify(t,g);self.assertEqual((m['category'],m['domain']),('inframe','inside'))
   (state,_,held),_=interp(t,value='0.4');self.assertEqual((state,held),('branch_threshold_met',False))
  t=tx('BRCA1','298_303del','Gly100_Gln102del');self.assertEqual(classify(t,'BRCA1')['domain'],'boundary')
  self.assertEqual(interp(t,value='0.8')[0][0],'branch_not_applicable')
 def test_intronic_no_amino_acid_and_canonical_sites(self):
  for dist in ['+1','-2','+9','-11']:
   t=tx('BRCA2','51'+dist+'A>G',terms=['intron_variant']);m=classify(t,'BRCA2')
   self.assertIsNone(m['aa']);self.assertEqual(interp(t)[0][2],dist in ('+1','-2'))
   if dist in ('+9','-11'):self.assertIn('剪接',interp(t)[0][1])
 def test_silent_and_hgvs_conflicts_do_not_suggest_protein_score(self):
  t=tx('BRCA1','57G>A','Gly19=');self.assertEqual(interp(t)[0][0],'branch_not_applicable')
  t['consequence_terms']=['stop_gained'];self.assertEqual(classify(t,'BRCA1')['category'],'ambiguous')
  t['transcript_id']='NM_007294.3';self.assertEqual(classify(t,'BRCA1')['category'],'unknown')
 def test_if_protein_positions_disagree_no_false_domain(self):
  t=tx('BRCA1','59G>C','Gly20Ala',protein_start=21)
  self.assertEqual(classify(t,'BRCA1')['domain'],'unknown')
 def test_stop_exon_is_not_first_event_exon(self):
  t=tx('BRCA1','55_56del','Val19GlyfsTer9')
  rows={'ST1':[{'gene':'BRCA1','cdna':[30,64],'exon':'2'},{'gene':'BRCA1','cdna':[65,102],'exon':'3'}],
        'T4_event':[{'gene':'BRCA1','event':'PTC','exon':'E3','row':1}]}
  m=ptc_hint('BRCA1',t,rows);self.assertEqual((m['exon'],m['stop_aa'],m['reference_codon']),('3',27,[81,83]))
  rows['ST1'][1]['cdna']=[82,102];self.assertIsNone(ptc_hint('BRCA1',t,rows))
  t['hgvsp']=NP['BRCA1']+':p.Val19GlyfsTer?';self.assertIsNone(ptc_hint('BRCA1',t,rows))
 def test_target_ref_anchor_must_match(self):
  car={'genomicAlleles':[{'referenceGenome':'GRCh37','hgvs':['NC_000017.10:g.101_102del'],'coordinates':[{'start':100,'end':102,'referenceAllele':'AT','allele':''}]}]}
  good=b'>NC_000017.10:100-102 x\nCAT\n'
  got=padded_indel(car,'BRCA1','GRCh37',lambda *a,**kw:good);self.assertEqual((got['pos'],got['ref'],got['alt']),(100,'CAT','C'))
  for bad in (good.replace(b'CAT',b'CAC'),good.replace(b'.10:',b'.11:'),b'<html>Error</html>',None):
   self.assertIsNone(padded_indel(car,'BRCA1','GRCh37',lambda *a,**kw:bad))
 def test_inapplicable_outranks_failed_or_observed_population_and_manual_click(self):
  s=deepcopy(next(x for x in examples() if x['id']=='demo-public-rc8-brca2-38'));r=s['report'];r['erepo']={'state':'explicit_no_records'};r['candidates']=[];r['facts']['event']='DEL'
  r['facts']['transcript']=tx('BRCA2','55_56del','Val19GlyfsTer9');r['facts']['gnomad']={'gno_nc2':{'label':'X','state':'http_error'}}
  z={a['code']:a for a in ledger(r)['rows']};self.assertEqual(len(z),19)
  for c in ('PM2','PP3','BP4'):
   self.assertEqual(z[c]['assessment']['applicability'],'not_applicable');self.assertIsNone(z[c]['assessment']['locally_meets']);self.assertIsNone(z[c]['brief']['option']);self.assertEqual(z[c]['brief']['result'],z[c]['evidence_view']['comparison'])
  p=dict(snapshot=s['id'],base_hash=s['report_sha256'],gene=s['gene'],assembly=s['assembly'],caid=s['caid'],previous_revision=None,selected=['PM2','PP3','BP4'],manual={'PM2':'Supporting'})
  saved=build(s,p,None);self.assertIsNone(saved['result']['points']);self.assertTrue(all(d['state']=='pending' and d['hold']=='spec_inapplicable' for d in saved['decisions']))
 def test_unknown_population_is_not_numeric_zero_or_negative(self):
  r=deepcopy(next(x for x in examples() if x['id']=='demo-public-rc8-brca2-38')['report']);r['erepo']={'state':'no_records'};r['facts']['event']='SNV';r['facts'].pop('population_evidence',None)
  for status in ('no_variant_record','http_error','no_verified_query_coordinate'):
   r['facts']['gnomad']={'gno_nc3':{'label':'X','state':status}}
   z=next(a for a in ledger(r)['rows'] if a['code']=='PM2')
   self.assertEqual(z['assessment']['applicability'],'applicable');self.assertIsNone(z['assessment']['locally_meets']);self.assertIsNone(z['brief']['option'])
class Provenance(unittest.TestCase):
 def test_typed_lr_not_generic_lr_threshold_or_history_for_pp1(self):
  lines=[{'value':s,'pointer':'/data'} for s in ('A combined LR for this variant of 3.57 (Family History LR=3.57).','PP1_Strong LR ≥18.7; co-segregation LR >=2.08; segregation LR of 2.08 is the threshold.','家族史 LR 与综合 LR 均为 3.57。')]
  self.assertEqual(text_lr(lines,'segregation'),[])
  for k in ('combined','family_history'):self.assertEqual([v['value'] for v in text_lr(lines,k)],['3.57','3.57'])
  self.assertEqual(text_lr([{'value':'Co-segregation LR = 5.2 (six relatives).'}],'segregation')[0]['value'],'5.2')
 def test_quarantine_one_historical_row_no_merge_or_probability_conversion(self):
  g='BRCA2';c='NM_000059.4:c.57G>A';h38='NC_000013.11:g.100G>A';h37='NC_000013.10:g.90G>A'
  car={'genomicAlleles':[{'referenceGenome':'GRCh38','hgvs':[h38]},{'referenceGenome':'GRCh37','hgvs':[h37]}]}
  a={'id':19,'Gene_Symbol':g,'HGVS_cDNA':c,'Genomic_HGVS_38':h38,'Genomic_HGVS_37':'NC_000017.10:g.90G>A','Data_Release':{'name':10,'id':11},'Segregation_LR_exLOVD':'305'}
  b=dict(a,id=20,Genomic_HGVS_37=h37,Data_Release={'name':11,'id':12},Segregation_LR_exLOVD='-',Posterior_probability_exLOVD='0.991')
  v=detail({'data':[a,b]},19,{'releaseName':10,'releaseId':11},car,g,c)
  self.assertEqual(v['method'],'separate_matching_release');self.assertEqual(v['ldh_linked_row']['issues'],['Genomic_HGVS_37']);self.assertIsNone(v['measures']['segregation']['value'])
  self.assertEqual(v['measures']['historical_posterior']['value'],'0.991')
  (lr,_,cmp),_=family({'facts':{'brca_detail':v}});self.assertEqual(lr['reported_values'],[]);self.assertIn('数值为空',cmp)
  self.assertIsNone(detail({'data':[a,b]},19,{'releaseName':9,'releaseId':11},car,g,c))
 def test_family_and_segregation_cannot_be_counted_twice(self):
  r={'erepo':{'state':'matched'}};s='Combined LR for this variant of 3.57 (based on family history LR=3.57).'
  self.assertEqual(family(r,'PP1',s)[0][0]['reported_values'],[]);(lr,_,cmp),obs=family(r,'PP4',s)
  self.assertEqual(lr['reported_values'],['3.57']);self.assertIn('Supporting',cmp);self.assertTrue(any(x[-1]=='expert_component' for x in obs))
