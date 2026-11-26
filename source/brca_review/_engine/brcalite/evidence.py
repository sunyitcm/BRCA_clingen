import json,re,urllib.parse as u
from pathlib import Path
from xml.etree import ElementTree as ET
from .core import validcar, requests_for,write_report,ref_matches,parse,NC
from .net import fetch,getj,writej,sha
from .report import html_report
def enrich_ldh(name,obj,expect_g,variant=None,reference=None):
 d=obj.get('data',obj)
 if name=='ldh_exome':
  g=d.get('gks_va_freq',{});l=d.get('locus',{});a=d.get('alleles',[])
  if len(a)!=2 or l.get('reference_genome')!='GRCh38':return None
  try:linked=parse(f"{l.get('contig')}-{l.get('position')}-{a[0]}-{a[1]}",'GRCh38')
  except ValueError:return None
  if linked['hgvs'] not in expect_g:
   if not (reference and variant and variant['kind']=='vcf' and variant['assembly']=='GRCh38' and variant['chrom']==linked['chrom'] and ref_matches(linked,reference.check(linked)) and reference.equivalent(linked['chrom'],linked['interbase'],variant['interbase'])):return None
  return {'source':g.get('derivedFrom'),'AC':g.get('focusAlleleCount'),'AN':g.get('locusAlleleCount'),'AF':g.get('alleleFrequency'),'quality':g.get('qualityMeasures'),'role':'context_only; not v1.2 dual non-cancer cohorts'}
 if d.get('entType')!={'ldh_brca':'BrcaExchangeRecord','ldh_consequence':'AlleleMolecularConsequenceStatement'}.get(name):return None
 if d.get('entId') not in expect_g:
  if not (reference and variant and variant['kind']=='vcf' and variant['assembly']=='GRCh38'):return None
  chrom=variant['chrom'];accession=NC['GRCh38'][chrom]
  deletion=re.fullmatch(re.escape(accession)+r':g\.(\d+)(?:_(\d+))?del',str(d.get('entId')))
  if not deletion:return None
  start,end=int(deletion[1])-1,int(deletion[2] or deletion[1])
  if not 0<=start<end or end-start>100 or abs(start-variant['interbase'][0])>10000:return None
  left=[start,end,reference.sequence(chrom,start,end),'']
  if not reference.equivalent(chrom,left,variant['interbase']):return None
 ec=d.get('entContent',{})
 return {'entity':d.get('entType'),'linked_allele':d.get('entId'),'version_source':d.get('entContentSrc'),'modified':d.get('modified'),'revision':d.get('rev'),'BayesDel_nsfp33a_noAF':ec.get('BayesDel_nsfp33a_noAF'),'role':'source fact; no PS3, no automatic PP3'}
def external(xml,vid):
 root=ET.fromstring(xml);v=root.find('.//VariationArchive')
 if v is None or v.get('VariationID')!=str(vid):return []
 out=[]
 for s in v.findall('./ClassifiedRecord/ClinicalAssertionList/ClinicalAssertion'):
  acc=s.find('ClinVarAccession');cl=s.find('Classification')
  if acc is None or cl is None or cl.findtext('ReviewStatus','').lower() not in ('reviewed by expert panel','practice guideline') or not cl.findtext('GermlineClassification'):continue
  out.append({'vcv':v.get('Accession')+'.'+v.get('Version'),'scv':acc.get('Accession')+'.'+acc.get('Version'),'submitter':acc.get('SubmitterName'),'condition':s.findtext('./TraitSet/Trait/Name/ElementValue'),'class':cl.findtext('GermlineClassification'),'review':cl.findtext('ReviewStatus'),'evaluated':cl.get('DateLastEvaluated'),'source_updated':acc.get('DateUpdated'),'vcep_version':None,'counts_as_criterion':False})
 return out
def cvurl(vid):return 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid=true&retmode=xml&id='+str(vid)
def online(folder,v,details,reference=None,previous=None):
 from .archive import replay_fetch
 folder.mkdir(parents=True,exist_ok=False);pieces={};links={}; warnings=[];states={}
 prior=None
 if previous is not None:
  prior=getj(previous/'report.json')
  if any(prior['input'].get(k)!=v.get(k) for k in ('assembly','kind','hgvs','chrom','pos','ref','alt')):
   raise ValueError('旧运行与本次输入或组装不符')
  prevref=prior['facts'].get('reference_check',{})
  if v['kind']=='vcf' and (bool(reference)!=(prevref.get('source')=='local_fasta') or
     reference and reference.provenance['manifest_sha256']!=prevref.get('manifest_sha256')):
   raise ValueError('旧运行与本次参考来源或清单不符')
 def visit(k,url,data=None,kind='json'):
  if previous is None:
   raw=fetch(folder/'http',k,url,data,kind)
  else:
   raw=replay_fetch(previous,folder/'http',k,url,data,kind)
  links[k]='http/'+k+'/receipt.json'
  states[k]='retrieved'
  if raw is None:
   receipt=getj(folder/links[k])
   if not receipt['attempts']:
    states[k]=receipt.get('state','not_requested');warnings.append(k+' 暂缓，见请求记录');return None
   attempt=receipt['attempts'][-1]
   # A generic 404 (endpoint or path) alone cannot establish allele absence.
   states[k]='http_error' if attempt.get('status') else 'transport_error'
   if k=='car' and attempt.get('status') in (200,400,404) and attempt.get('body_file'):
    try:problem=getj(folder/'http'/k/attempt['body_file'])
    except (ValueError,UnicodeError):problem={}
    if isinstance(problem,dict) and problem.get('errorType')=='NotFound':states[k]='registry_not_found'
   warnings.append(k+' 获取失败，见请求记录');return None
  return json.loads(raw) if kind=='json' else raw
 req=requests_for(v)
 if v['kind']=='vcf':
  if reference:
   try:
    pieces['ref']=reference.check(v)
    states['ref']='local_verified_file'
   except ValueError as err:
    pieces['ref']=None;states['ref']='local_region_unavailable';warnings.append(str(err))
  else:
   pieces['ref']=visit('ref',*req['ref'])
  if not ref_matches(v,pieces['ref']):
   states['ref_check']='mismatch' if pieces['ref'] else 'unavailable'
   warnings.append('原始 REF 无法核准，已停止后续请求')
  else:states['ref_check']='matched'
 if states.get('ref_check','matched')=='matched':
  for k in ('car','info','vep'):
   pieces[k]=visit(k,*req[k])
 auth=validcar(v,pieces.get('car'),reference)
 if auth:
  ca,g,_=auth; base='https://ldh.clinicalgenome.org/ldh/Variant/id/'+ca
  pieces['ldh']=visit('ldh',base+'?detail=low')
 r=write_report(folder,v,pieces,links,reference)
 if details and auth and r['status']=='review_required':
  car=pieces['car'];known_g={s for a in car.get('genomicAlleles',[]) if a.get('referenceGenome')=='GRCh38' for s in a.get('hgvs',[])}
  med=visit('ldh_med',base+'?detail=med') or {};d=med.get('data',{});rels=d.get('ld',{}) if d.get('entId')==ca else {}
  for k,typ in [('ldh_exome','GnomADExomesV4.1'),('ldh_consequence','AlleleMolecularConsequenceStatement'),('ldh_brca','BrcaExchangeRecord')]:
   arr=rels.get(typ,[])
   if not isinstance(arr,list) or len(arr)!=1:
    states[k]='no_link_listed' if d.get('entId')==ca and arr==[] else 'unresolved_links'
    continue
   url=arr[0].get('dataIri') or arr[0].get('ldhIri') or ''
   path_pattern=(rf'/ldh/dss/cg/ns/gnomad_v4_1/set/exome/id/{ca}/data' if k=='ldh_exome' else rf'/ldh/{typ}/id/[0-9]+')
   parts=u.urlsplit(url)
   if parts.scheme!='https' or parts.netloc not in ('ldh.genome.network','ldh.clinicalgenome.org') or not re.fullmatch(path_pattern,parts.path) or parts.query or parts.fragment:
    states[k]='unrecognized_entity_url';continue
   obs=visit(k,url)
   if obs:
    facts=enrich_ldh(k,obs,known_g,v,reference)
    if facts:r['facts'][k]=facts
    else:states[k]='allele_not_verified';warnings.append(k+' linked allele not verified; omitted')
  ids={x['variationId'] for x in car.get('externalRecords',{}).get('ClinVarVariations',[]) if str(x.get('variationId','')).isdigit()}
  if len(ids)==1:
   vid=next(iter(ids));xml=visit('clinvar',cvurl(vid),kind='xml')
   if xml:
    r['expert_assertions']=external(xml,vid)
    cv=ET.fromstring(xml).find('.//VariationArchive')
    if cv is not None and cv.get('VariationID')==str(vid) and re.fullmatch('VCV[0-9]+',cv.get('Accession','')) and cv.get('Version','').isdigit():
     r['facts']['clinvar_record']={'accession':cv.get('Accession')+'.'+cv.get('Version'),'updated':cv.get('DateLastUpdated'),'role':'allele_record_not_an_evidence_code'}
  else:states['clinvar']='no_unique_variation_id'
 r['source_status']=states
 r['acquisition']='replayed_http' if previous else 'live_http'
 if prior is not None:
  r['previous_http_report']={'sha256':sha((previous/'report.json').read_bytes()),'engine':prior.get('engine'),'run_id':prior.get('run_id')}
 r['details_requested']=details
 r['source_files']=links;r['limitations'].extend(warnings)
 writej(folder/'report.json',r);html_report(folder/'report.html',r);return r
