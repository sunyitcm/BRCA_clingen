"""Score curator-accepted, source-linked independent lines, on the germline profile only."""
from decimal import Decimal,InvalidOperation
from collections import Counter
from datetime import date
import re
from .core import GENES,R
from .net import getj
P={'Supporting':1,'Moderate':2,'Strong':4,'Very Strong':8}
def bin_points(p):
 return 'Pathogenic' if p>=10 else 'Likely pathogenic' if p>=6 else 'Benign' if p<=-7 else 'Likely benign' if p<=-2 else 'Uncertain significance'
def combine(g,d,e):
 v,s,m,u=[g[k] for k in ('Very Strong','Strong','Moderate','Supporting')]
 if d=='P':
  if (v>=1 and (s>=1 or m>=1 or u>=2)) or (s>=1 and (m>=3 or (m>=2 and u>=2) or (m>=1 and u>=4))) or (s>=2 and (m>=1 or u>=2)) or s>=3:return 'Pathogenic'
  if (v>=1 and u>=1) or s>=2 or (s>=1 and (m>=1 or u>=2)) or m>=3 or (m>=2 and u>=2) or (m>=1 and u>=4):return 'Likely pathogenic'
 else:
  if g['Stand Alone'] or (v>=1 and (s+m+u>=1)) or s>=2 or (s and (m>=2 or (m and u) or u>=3)):return 'Benign'
  if v or (s and (u or m)) or (m and u) or u>=2 or any(x.get('multi_type_confirmed') is True and x['strength']=='Strong' for x in e):return 'Likely benign'
 return 'Uncertain significance'

def _filled(s):return isinstance(s,str) and bool(s.strip())
def _ids(xs):return isinstance(xs,list) and bool(xs) and all(_filled(x) for x in xs)
def has_provenance(x,packet):
 try:date.fromisoformat(x.get('adjudicated_on',packet.get('adjudicated_on')))
 except (TypeError,ValueError):return False
 return x.get('origin',packet.get('origin'))=='manual_adjudication'
def adjudicate(packet,*,verified_text_codes=frozenset()):
 g=packet.get('gene'); v=GENES.get(g);codes=packet.get('criteria',[])
 base={'local_final':None,'draft_class':None,'points':None,'mode':None,'issues':[]}
 if not v or not isinstance(codes,list) or packet.get('spec')!=v[2]+'@1.2' or packet.get('condition')!=v[1] or not _filled(packet.get('allele')) or packet.get('domain')!='germline':return base|{'issues':['scope/variant/specification mismatch']}
 if not re.fullmatch(r'CA\d+|TEST_ONLY_[A-Za-z0-9_]+',packet['allele']):return base|{'issues':['canonical allele or explicitly synthetic identifier required']}
 allowed=getj(R/'rules/catalog.json')[g];gates=[];seen=set();all_codes=set();units=set();accepted=[];counts={k:Counter() for k in ('P','B')}
 for x in codes:
  if not isinstance(x,dict):gates.append('malformed evidence');continue
  if not has_provenance(x,packet):gates.append('missing manual_adjudication origin or ISO review date')
  state=x.get('state');key=x.get('code')
  if not _filled(key) or key in all_codes:gates.append('invalid or repeated criterion, including negative entries');continue
  all_codes.add(key)
  if state in ('not_met','not_applicable'):
   if x.get('code') not in allowed or not _filled(x.get('reviewer')) or not _filled(x.get('reason')) or not _ids(x.get('source_refs')):gates.append('unsupported negative decision / missing audit trail')
   continue
  if state!='met':gates.append('unresolved or invalid evidence state: '+str(state));continue
  k=x.get('code','');w=x.get('strength');desc=allowed.get(k,{}); ids=x.get('independent_units');refs=x.get('source_refs')
  # Internal review adapter may authorize these four v1.2 prose levels only after checking the source, LR and additional criteria.
  text_level=w=='Very Strong' and k in ('PP1','PP4','BS4','BP5') and k in verified_text_codes
  if k in seen or not desc.get('enabled') or (w not in desc.get('strengths',[]) and not text_level):gates.append('duplicate, disabled or unmapped code/strength: '+str(k));continue
  if (not isinstance(ids,list) or not ids or not all(isinstance(z,str) and z.strip() for z in ids) or not isinstance(refs,list) or not refs or not all(isinstance(z,str) and z.strip() for z in refs) or not _filled(x.get('reviewer')) or not _filled(x.get('reason'))):gates.append('missing basis/source/reviewer/reason: '+k);continue
  if units.intersection(ids) or len(set(ids))!=len(ids):gates.append('nonindependent evidence: '+k)
  seen.add(k);units.update(ids);accepted.append(x);counts[k[0]][w]+=1
 pm3_policy=getj(R/'rules/review-20261124.json')[g]
 if 'PM3' in seen and seen.intersection(pm3_policy['excludes_if_met']):gates.append('PM3 population exclusion: '+pm3_policy['field'])
 for group in [('PP3','PVS1'),('PP3','BP4'),('BA1','BS1','PM2'),('BP1','BP4'),('BP1','PP3'),('PS3','BS3'),('PP4','BP5'),('PP1','BS4')]:
  if len(seen.intersection(group))>1:gates.append('requires dependency adjudication: '+','.join(group))
 if any(x['code']=='PM5' and (not x.get('event_exon') or not x.get('stop_exon') or x.get('event_exon')!=x.get('stop_exon')) for x in accepted):gates.append('PM5_PTC: exon absent or cross-exon source discrepancy; hold')
 if 'PM5' in seen and ('PVS1' not in seen or any(x.get('subtype')!='genomic_PTC' for x in accepted if x['code']=='PM5')):gates.append('PM5 repurposed only as genomic PTC in applicable PVS1 context')
 if not accepted:gates.append('no accepted evidence')
 if gates:return base|{'issues':gates}
 stand='BA1' in seen;pos=sum(counts['P'].values());neg=sum(counts['B'].values())
 subtotal=None if stand else sum((1 if x['code'][0]=='P' else -1)*P[x['strength']] for x in accepted)
 base['points']=subtotal
 if packet.get('significant_conflict') is True:return base|{'mode':'discuss','issues':['significant conflict: expert discussion required']}
 if packet.get('review_complete') is not True or packet.get('independence_reviewed') is not True or not _filled(packet.get('case_reviewer')):return base|{'issues':['evidence inventory / independence review not complete; subtotal only']}
 if str(packet['allele']).startswith('CA'):
  missing=sorted(k for k,rule in allowed.items() if rule['enabled'] and k not in all_codes)
  if missing:return base|{'issues':['applicable criteria lack a review disposition: '+', '.join(missing)]}
 if (pos and neg) or packet.get('uninformative_prediction_reviewed') is True:
  if stand or min(pos,neg)>1 or packet.get('significant_conflict') is not False:return base|{'mode':'discuss','issues':['stand-alone / significant or unresolved conflict: expert discussion required']}
  return base|{'mode':'ENIGMA_conflict_points','draft_class':bin_points(subtotal)}
 return base|{'mode':'ENIGMA_Table3','draft_class':combine(counts['P' if pos else 'B'],'P' if pos else 'B',accepted)}

def frequency(f):
 """Candidate branch for positively established site-level facts; never certifies data by existence of a field."""
 from math import isfinite
 def res(status,code=None,strength=None,why=''):return dict(status=status,code=code,strength=strength,reason=why,score_counted=False)
 if f.get('allele_verified') is not True or not f.get('allele'):return res('insufficient',why='identity unverified')
 if f.get('event') not in ('SNV','INS','DEL','DELINS','MNV'):return res('insufficient',why='unsupported event')
 rows=f.get('strata')
 if not isinstance(rows,list) or not rows:return res('insufficient',why='no strata')
 passing=[];abs_ok=set();assessed=set();observed=False
 for z in rows:
  if not isinstance(z,dict):return res('insufficient',why='invalid stratum')
  label=(z.get('dataset'),z.get('assay'))
  if label not in (('gnomAD_v2.1_non_cancer','exome'),('gnomAD_v3.1_non_cancer','genome')):continue
  if not (z.get('qc_pass') is True and z.get('non_founder') is True and z.get('source_ref') and z.get('ancestry')):continue
  dp=z.get('mean_depth');ac=z.get('AC');an=z.get('AN')
  if type(dp) not in (float,int) or not isfinite(dp) or type(ac)!=int or type(an)!=int or not (an>0 and 0<=ac<=an):continue
  assessed.add(label)
  observed=observed or ac>0
  if dp>=20 and z.get('faf') is not None and z.get('statistic_verified') is True:
   try:x=Decimal(str(z['faf']))
   except InvalidOperation:return res('insufficient',why='invalid FAF')
   if not x.is_finite() or not (0<=x<=1) or (ac==0 and x!=0):return res('insufficient',why='invalid FAF/counts')
   passing.append(x)
  if ac==0 and dp>=25 and z.get('absence_verified_full_subset') is True:abs_ok.add(label)
 if passing and f.get('known_pathogenic_founder') is False and f.get('maximum_nonfounder_reviewed') is True:
  x=max(passing)
  if x>Decimal('.001'):return res('candidate','BA1','Stand Alone')
  if x>Decimal('.0001'):return res('candidate','BS1','Strong')
  if x>Decimal('.00002'):return res('candidate','BS1','Supporting')
 if (len(abs_ok)==2 and f.get('event')=='SNV' and
     f.get('both_subsets_checked') is True and not any(z.get('AC',0)>0 for z in rows if z.get('dataset') in ('gnomAD_v2.1_non_cancer','gnomAD_v3.1_non_cancer') and isinstance(z.get('AC'),int))):return res('candidate','PM2','Supporting')
 if not assessed:return res('insufficient',why='no eligible source strata with QC, depth, non-founder and allele counts')
 freq_assessed=(f.get('known_pathogenic_founder') is True or
                bool(passing) and f.get('known_pathogenic_founder') is False and f.get('maximum_nonfounder_reviewed') is True)
 absence_assessed=observed or f.get('event') in ('INS','DEL','DELINS')
 if not freq_assessed or not absence_assessed:return res('insufficient',why='one or more population branches lack eligible inputs or coverage')
 return res('no_eligible_candidate',why='reviewed FAF/founder status and observed counts/event do not support an eligible code')
