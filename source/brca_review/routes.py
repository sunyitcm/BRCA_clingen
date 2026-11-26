"""Independent read-only acquisitions. Unverified observations never grant a score."""
import json
import re
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import urlencode, urlsplit
from xml.etree import ElementTree as ET
from ._engine.brcalite import core, evidence
from ._engine.brcalite.net import fetch, getj, writej

DATASETS = {
    'gno_full37': ('gnomad_r2_1', 'GRCh37', 'v2.1 全队列'),
    'gno_full38': ('gnomad_r4', 'GRCh38', 'v4 全队列'),
    'gno_nc2': ('gnomad_r2_1_non_cancer', 'GRCh37', 'v2.1 非癌外显子'),
    'gno_nc3': ('gnomad_r3_non_cancer', 'GRCh38', 'v3.1 非癌基因组'),
}
PROTEIN={'BRCA1':'NP_009225.1','BRCA2':'NP_000050.3'}
# Ancestries selected for comparison; do not treat every aggregate or founder stratum as outbred.
OUTBRED={'afr','amr','eas','nfe','sas'}
API='https://gnomad.broadinstitute.org/api'

def mapped(car,g,b):
    chrom='17' if g=='BRCA1' else '13'
    nc=core.NC[b][chrom]
    rows=[a for a in car.get('genomicAlleles',[]) if a.get('referenceGenome')==b and any(h.startswith(nc+':g.') for h in a.get('hgvs',[]))]
    if len(rows)!=1 or len(rows[0].get('coordinates',[]))!=1:return
    q=rows[0]['coordinates'][0];start,end=q.get('start'),q.get('end');ref,alt=q.get('referenceAllele'),q.get('allele')
    if type(start)!=int or start<0 or end!=start+1 or not all(isinstance(s,str) and re.fullmatch('[ACGT]',s) for s in (ref,alt)) or ref==alt:return
    k=core.parse(f'{chrom}-{end}-{ref}-{alt}',b)
    return k if k['hgvs'] in rows[0]['hgvs'] else None

def padded_indel(car,g,b,visit):
    """Obtain the target-build VCF anchor from the versioned reference, never from the other assembly."""
    chrom='17' if g=='BRCA1' else '13';acc=core.NC[b][chrom]
    hits=[a for a in car.get('genomicAlleles',[]) if a.get('referenceGenome')==b and any(h.startswith(acc+':g.') for h in a.get('hgvs',[]))]
    if len(hits)!=1 or len(hits[0].get('coordinates',[]))!=1:return
    q=hits[0]['coordinates'][0];st,en,ref,alt=(q.get(k) for k in ('start','end','referenceAllele','allele'))
    if type(st)!=int or type(en)!=int or st<=0 or not all(isinstance(v,str) and re.fullmatch('[ACGT]{0,99}',v) for v in (ref,alt)) or en!=st+len(ref) or bool(ref)==bool(alt):return
    ask=dict(assembly=b,chrom=chrom,pos=st,ref='N'*(len(ref)+1))
    raw=visit('anchor_'+b[-2:],'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+urlencode(dict(db='nuccore',id=acc,seq_start=st,seq_stop=en,strand=1,rettype='fasta',retmode='text')),kind='text')
    hit=checked_ref(ask,raw);seq=(hit or {}).get('seq','')
    if not seq or not re.fullmatch('[ACGT]+',seq) or seq[1:]!=ref:return
    site=core.parse(f"{chrom}-{st}-{seq}-{seq[0]+alt}",b)
    return site if site['interbase']==[st,en,ref,alt] and site['hgvs'] in hits[0]['hgvs'] else None

def checked_ref(v,body):
    if not isinstance(body,bytes):return
    l=body.decode('ascii','replace').strip().splitlines()
    end=v['pos']+len(v['ref'])-1;acc=core.NC[v['assembly']][v['chrom']]
    if not l:return
    head=l[0].split(' ',1)[0]
    allowed={f">{acc}:{v['pos']}-{end}"}
    if end==v['pos']:allowed.add(f">{acc}:{end}")
    seq=''.join(l[1:]).upper()
    if head not in allowed or len(seq)!=len(v['ref']) or not re.fullmatch('[ACGTN]+',seq):return
    return dict(id=f"chromosome:{v['assembly']}:{v['chrom']}:{v['pos']}:{end}:1",seq=seq,source='NCBI RefSeq EFetch',accession=acc,header=l[0])

def car_tx(car,auth):
    _,g,names=auth;h=next(iter(names));acc=core.GENES[g][0];np=PROTEIN[g]
    rows=[t for t in car.get('transcriptAlleles',[]) if t.get('geneSymbol')==g and h in t.get('hgvs',[])]
    ps={p for t in rows for p in [(t.get('proteinEffect') or {}).get('hgvs')] if isinstance(p,str) and p.startswith(np+':p.')}
    p=next(iter(ps)) if len(ps)==1 else None
    from .mechanism import protein_effect
    effect=protein_effect(g,p)
    return dict(gene_symbol=g,transcript_id=acc,hgvsc=h,hgvsp=p,protein_start=effect.get('start'),protein_end=effect.get('end'),exon=None,consequence_terms=[effect['term']] if effect.get('term') else [])

def vep_tx(v,data,g):
    if not isinstance(data,list) or len(data)!=1 or not isinstance(data[0],dict):return
    d=data[0];expected=core.requests_for(v)['vep'][1]
    if d.get('input')!=(expected['variants'][0] if expected else v['hgvs']) or d.get('assembly_name')!=v['assembly']:return
    if not isinstance(d.get('transcript_consequences'),list):return
    ts=[t for t in d['transcript_consequences'] if isinstance(t,dict) and t.get('gene_symbol')==g and t.get('transcript_id')==core.GENES[g][0] and (t.get('hgvsc') or '').startswith(core.GENES[g][0]+':c.')]
    if len(ts)==1:return {k:ts[0].get(k) for k in ('gene_symbol','transcript_id','hgvsc','hgvsp','exon','consequence_terms','protein_start','protein_end')}

def variant_id(v):return f"{v['chrom']}-{v['pos']}-{v['ref']}-{v['alt']}"

def query_args(k,v):
    ds,b,_=DATASETS[k];vid=variant_id(v)
    detail='ac an filters populations{id ac an} faf95{popmax popmax_population}'
    assay=('exome{'+detail+'}' if k=='gno_nc2' else 'genome{'+detail+'}' if k=='gno_nc3' else 'exome{ac an filters} genome{ac an filters} transcript_consequences{gene_symbol transcript_id transcript_version hgvsc hgvsp consequence_terms}')
    return {'query':'{v:variant(variantId:'+json.dumps(vid)+',dataset:'+ds+'){variant_id reference_genome '+assay+'}}'}

def gnomad_record(k,v,doc,g):
    ds,b,label=DATASETS[k];vid=variant_id(v)
    rec=dict(dataset=ds,assembly=b,variant_id=vid,label=label,cohorts={},transcripts=[],record_url='https://gnomad.broadinstitute.org/variant/'+vid+'?dataset='+ds)
    if not isinstance(doc,dict):return rec,'not_available'
    errors=doc.get('errors') or [];data=doc.get('data') or {}
    if not isinstance(data,dict) or not isinstance(errors,list):return rec,'graphql_error'
    row=data.get('v')
    if errors:
        if 'v' in data and row is None and all(str(e.get('message','')).strip().lower()=='variant not found' for e in errors if isinstance(e,dict)) and all(isinstance(e,dict) for e in errors):
            return rec,'no_variant_record'
        if any(isinstance(e,dict) and str(e.get('message','')).strip().lower()=='service overloaded' for e in errors):rec['service_message']='Service overloaded'
        return rec,'graphql_error'
    if not isinstance(row,dict) or (row.get('variant_id'),row.get('reference_genome'))!=(vid,b):
        return rec,'allele_not_verified'
    for key in ('exome','genome'):
        if key not in row:continue
        grp=row[key]
        if grp is None:rec['cohorts'][key]=None;continue
        if not isinstance(grp,dict):return rec,'invalid_fields'
        ac,an=grp.get('ac'),grp.get('an')
        if type(ac)!=int or type(an)!=int or not 0<=ac<=an:return rec,'invalid_fields'
        if grp.get('populations') is not None and (not isinstance(grp['populations'],list) or not all(isinstance(p,dict) for p in grp['populations'])):return rec,'invalid_fields'
        if grp.get('faf95') is not None and not isinstance(grp['faf95'],dict):return rec,'invalid_fields'
        if grp.get('filters') is not None and not isinstance(grp['filters'],list):return rec,'invalid_fields'
        grp=dict(grp,af=ac/an if an else None)
        rec['cohorts'][key]=grp
    if row.get('transcript_consequences') is not None and not isinstance(row['transcript_consequences'],list):return rec,'invalid_fields'
    rec['transcripts']=[x for x in (row.get('transcript_consequences') or []) if isinstance(x,dict) and x.get('gene_symbol')==g]
    return rec,'variant_record'

def frequency_entity(k,obj,known,variant=None,reference=None):
    # Lookup is scoped to an exact CA, dataset collection and allele; absence in a link directory is never zero.
    d=obj.get('data',obj);out=evidence.enrich_ldh('ldh_exome',obj,known,variant,reference)
    if not out:return None
    a,n,p=out.get('AC'),out.get('AN'),out.get('AF')
    dataset=out.get('source') or {}
    if not isinstance(dataset,dict) or dataset.get('version') not in (None,'4.1'):return None
    if out.get('quality') is not None and not isinstance(out['quality'],dict):return None
    if type(a)!=int or type(n)!=int or not 0<=a<=n or n<=0:return None
    if p is not None and (type(p) not in (int,float) or not 0<=p<=1):return None
    if 'noVariant' in d and (type(d['noVariant']) is not bool or (d['noVariant'] and a!=0)):return None
    out.update(_build='GRCh38',assay=k[4:])
    if 'noVariant' in d:out['noVariant']=d['noVariant']
    return out

def scv_citations(cv):
    out=[];seen=set()
    for index,case in enumerate(cv.findall('./ClassifiedRecord/ClinicalAssertionList/ClinicalAssertion')):
        acc=case.find('ClinVarAccession');cl=case.find('Classification')
        if acc is None or cl is None or cl.find('GermlineClassification') is None:continue
        scv=acc.get('Accession','');version=acc.get('Version','')
        if not re.fullmatch('SCV[0-9]+',scv) or not version.isdigit():continue
        for path in ('Classification/Citation/ID','ObservedIn/ObservedData/Citation/ID','Citation/ID'):
            for ident in case.findall(path):
                pmid=(ident.text or '').strip()
                if ident.get('Source')!='PubMed' or not re.fullmatch(r'[1-9][0-9]{0,9}',pmid) or (scv,pmid,path) in seen:continue
                seen.add((scv,pmid,path));out.append(dict(pmid=pmid,scv=scv+'.'+version,locator='ClinicalAssertion['+str(index+1)+']/'+path))
    return out

def online(folder,v,details,reference=None,*,workbooks=True,fetcher=None):
    folder.mkdir(parents=True,exist_ok=False); states={};links={};ev={};pieces={};mapped_requests={};gn={};fetcher=fetcher or fetch
    b=v['assembly'];g=v.get('gene') or ('BRCA1' if v['chrom']=='17' else 'BRCA2');other='GRCh38' if b=='GRCh37' else 'GRCh37';req=core.requests_for(v)
    def visit(k,url,data=None,kind='json'):
        raw=fetcher(folder/'http',k,url,data,kind,attempts=2,timeout=12)
        links[k]='http/'+k+'/receipt.json'
        if raw is not None:
            states[k]='retrieved'
            return json.loads(raw) if kind=='json' else raw
        m=getj(folder/links[k]);a=(m.get('attempts') or [{}])[-1];status=a.get('status')
        states[k]=('http_error' if status and status>=400 else 'invalid_response' if status else m.get('state','transport_error'))
        if k=='car' and a.get('body_file') and status in (200,400,404):
            try:
                if getj(folder/'http'/k/a['body_file']).get('errorType')=='NotFound':states[k]='registry_not_found'
            except (ValueError,AttributeError):pass
        return None
    def independent(key,call,*args):
        try:return call(*args)
        except (TypeError,AttributeError,KeyError,IndexError,ValueError,ET.ParseError):
            states[key]='invalid_fields'
            return None
    def getref():
        if v['kind']!='vcf':return
        if reference:
            try:rr=reference.check(v);states['ref']='local_verified_file'
            except ValueError:rr=None;states['ref']='local_region_unavailable'
        else:
            rr=visit('ref',*req['ref'])
            if not isinstance(rr,dict) or not isinstance(rr.get('seq'),str) or not rr.get('id'):
                acc=core.NC[b][v['chrom']];args=dict(db='nuccore',id=acc,seq_start=v['pos'],seq_stop=v['pos']+len(v['ref'])-1,strand=1,rettype='fasta',retmode='text')
                rr=checked_ref(v,visit('refseq_ref','https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+urlencode(args),kind='text'))
        pieces['ref']=rr
        states['ref_check']='matched' if core.ref_matches(v,rr) else 'mismatch' if rr else 'unavailable'
    def gnrun(k,vv):
        raw=visit(k,API,query_args(k,vv))
        row,s=gnomad_record(k,vv,raw,g)
        if raw is not None:states[k]=s
        row['state']=states[k];gn[k]=row
    def put(k):pieces[k]=visit(k,*req[k])
    def ldh(ca,car):
        base='https://ldh.clinicalgenome.org/ldh/Variant/id/'+ca
        raw=visit('ldh_med',base+'?detail=med')
        if not isinstance(raw,dict) or not isinstance(raw.get('data'),dict) or raw['data'].get('entId')!=ca:
            for k in ('ldh_exome','ldh_genome','ldh_brca','ldh_consequence'):states[k]='directory_unavailable'
            return
        pieces['ldh']=raw;rels=raw['data'].get('ld');known={h for a in car.get('genomicAlleles',[]) if a.get('referenceGenome')=='GRCh38' for h in a.get('hgvs',[])}
        for k,typ in [('ldh_exome','GnomADExomesV4.1'),('ldh_genome','GnomADGenomesV4.1'),('ldh_brca','BrcaExchangeRecord'),('ldh_consequence','AlleleMolecularConsequenceStatement')]:
            arr=rels.get(typ,[]) if isinstance(rels,dict) else None
            if not isinstance(arr,list) or len(arr)!=1:
                states[k]='no_link_listed' if arr==[] else 'unresolved_links';continue
            url=(arr[0].get('dataIri') or arr[0].get('ldhIri') or '') if isinstance(arr[0],dict) else ''
            path=rf'/ldh/dss/cg/ns/gnomad_v4_1/set/{k[4:]}/id/{ca}/data' if k in ('ldh_exome','ldh_genome') else rf'/ldh/{typ}/id/[0-9]+'
            u=urlsplit(url)
            if u.scheme!='https' or u.netloc not in ('ldh.genome.network','ldh.clinicalgenome.org') or not re.fullmatch(path,u.path) or u.query or u.fragment:states[k]='unrecognized_entity_url';continue
            obj=visit(k,url)
            if obj is not None:
                try:
                    if k in ('ldh_exome','ldh_genome'):
                        value=frequency_entity(k,obj,known,v,reference)
                    else:value=evidence.enrich_ldh(k,obj,known,v,reference)
                except (TypeError,AttributeError,KeyError,ValueError):value=None
                if value:
                    ev[k]=value
                    if k=='ldh_brca':
                        from .exchange import detail
                        content=obj.get('data',obj).get('entContent') or {};vid=content.get('variantId');release=value.get('version_source') or {}
                        h=next(iter(auth[2]))
                        if type(vid)==int and vid>0 and isinstance(release,dict) and release.get('releaseId') is not None:
                            u='https://brcaexchange.org/backend/data/variant/?'+urlencode({'variant_id':vid})
                            try:
                                rr=visit('brca_detail',u)
                                if rr is not None:
                                    candidate=detail(rr,vid,release,car,g,h)
                                    if candidate:ev['brca_detail']=candidate
                                    else:states['brca_detail']='allele_or_release_unverified'
                            except (ValueError,TypeError,KeyError,AttributeError,OSError):
                                states['brca_detail']='invalid_fields'
                        else:states['brca_detail']='no_verified_record_id'
                    elif k=='ldh_consequence':
                        rows=(obj.get('data',obj).get('entContent') or {}).get('preferredTranscripts',[])
                        hits=[t for t in rows if isinstance(t,dict) and t.get('id')==core.GENES[g][0] and isinstance(t.get('molecularConsequence'),str)]
                        if len(hits)==1:value['verified_refseq_consequence']=hits[0]['molecularConsequence']
                else:states[k]='allele_not_verified'
    def clinvar(car):
        ids={str(x['variationId']) for x in (car.get('externalRecords') or {}).get('ClinVarVariations',[]) if str(x.get('variationId','')).isdigit()}
        if len(ids)!=1:states['clinvar']='no_unique_variation_id';return
        vid=ids.pop();xml=visit('clinvar',evidence.cvurl(vid),kind='xml')
        if xml:
            cv=ET.fromstring(xml).find('.//VariationArchive')
            if cv is None or cv.get('VariationID')!=vid:states['clinvar']='allele_not_verified';return
            if re.fullmatch('VCV[0-9]+',cv.get('Accession','')) and cv.get('Version','').isdigit():
                ev['clinvar_record']={'accession':cv.get('Accession')+'.'+cv.get('Version'),'updated':cv.get('DateLastUpdated'),'germline_citations':scv_citations(cv)}
                agg=cv.find('./ClassifiedRecord/Classifications/GermlineClassification')
                if agg is not None:
                    ev['clinvar_record']['germline_aggregate']={k:agg.findtext(k) for k in ('Description','ReviewStatus','Explanation') if agg.findtext(k)}
                ev['expert_assertions']=evidence.external(xml,vid)
    jobs=[]
    with ThreadPoolExecutor(max_workers=4) as ex:
        jobs.append(ex.submit(getref));native=ex.submit(put,'vep');jobs.append(ex.submit(put,'info'));ca_job=ex.submit(put,'car')
        def dispatch(vv,mode):
            full='gno_full37' if vv['assembly']=='GRCh37' else 'gno_full38';sub='gno_nc2' if vv['assembly']=='GRCh37' else 'gno_nc3'
            for k in ([full,sub] if mode=='native' else [sub]):jobs.append(ex.submit(gnrun,k,vv))
        if v['kind']=='vcf':dispatch(v,'native')
        ca_job.result();car=pieces.get('car') or {}
        try:
            if not isinstance(car,dict):raise TypeError
            if not all(isinstance(car.get(k,[]),list) and all(isinstance(row,dict) for row in car.get(k,[])) for k in ('genomicAlleles','transcriptAlleles')):raise TypeError
            for row in car.get('genomicAlleles',[])+car.get('transcriptAlleles',[]):
                if not isinstance(row.get('hgvs',[]),list) or not all(isinstance(h,str) for h in row.get('hgvs',[])):raise TypeError
                if not isinstance(row.get('coordinates',[]),list) or not all(isinstance(q,dict) for q in row.get('coordinates',[])):raise TypeError
                if row.get('proteinEffect') is not None and not isinstance(row['proteinEffect'],dict):raise TypeError
            auth=core.validcar(v,car,reference)
        except (TypeError,AttributeError,KeyError,ValueError):
            car={};pieces['car']=None;auth=None;states['car']='invalid_fields'
        if auth and auth[1]==g:
            for build in (b,other):
                mv=mapped(car,g,build)
                if mv is None and (build!=b or v['kind']=='hgvsc'):mv=padded_indel(car,g,build,visit)
                if not mv:continue
                mapped_requests[build]=mv
                if build!=b or v['kind']!='vcf':dispatch(mv,'native' if build==b else 'cross')
            if details:jobs.extend([ex.submit(independent,'ldh_med',ldh,auth[0],car),ex.submit(independent,'clinvar',clinvar,car)])
        native.result()
        mt=mapped_requests.get('GRCh38')
        if b=='GRCh37' and mt and not vep_tx(v,pieces.get('vep'),g):
            mr=core.requests_for(mt)
            def mapped_vep(k):pieces['mapped_'+k]=visit('mapped_'+k,*mr[k])
            jobs.extend(ex.submit(mapped_vep,k) for k in ('vep','info'))
        for j in jobs:j.result()
    for target in ('gno_full37' if b=='GRCh37' else 'gno_full38','gno_nc2','gno_nc3'):
        states.setdefault(target,'no_verified_query_coordinate')
    auth=auth if auth and auth[1]==g else None
    # Only reference-consistent native reports may contribute evidence. Linked assertions still have saved receipts.
    identity=bool(auth and (v['kind']=='hgvsc' or states.get('ref_check')=='matched'))
    # Initialise the original report envelope without invoking the workbook or the VEP dependency gate.
    for meta in ('info','mapped_info'):
        if meta in pieces and pieces[meta] is not None and (not isinstance(pieces[meta],dict) or not isinstance(pieces[meta].get('releases'),list)):
            pieces[meta]=None;states[meta]='invalid_fields'
    seed={k:pieces[k] for k in ('car','info','ref','ldh') if isinstance(pieces.get(k),dict)}
    r=core.write_report(folder,v,seed,links,reference);f=r['facts']
    r['limitations']=[s for s in r['limitations'] if s not in ('CAR 无唯一精确核酸/规范转录本对应','REF、坐标或序列组装未核对','VEP 响应未逐输入唯一对账')]
    r['selected_gene']=g;r['acquisition']='live_independent';r['details_requested']=details
    if not auth:
        states.setdefault('ldh_med','needs_caid');states.setdefault('clinvar','needs_variation_id')
    tcar=car_tx(car,auth) if auth else None;tn=vep_tx(v,pieces.get('vep'),g)
    tm=vep_tx(mapped_requests['GRCh38'],pieces.get('mapped_vep'),g) if 'mapped_vep' in pieces else None
    annotations={k:t for k,t in [('car',tcar),('vep',tn),('mapped_vep',tm)] if t}
    f['source_transcripts']=annotations
    if not isinstance(f.get('ensembl_metadata'),dict):f.pop('ensembl_metadata',None)
    if isinstance(pieces.get('mapped_info'),dict):f['mapped_ensembl_metadata']=pieces['mapped_info']
    for k in ('vep','mapped_vep'):
        if k in pieces and states.get(k)=='retrieved' and k not in annotations:states[k]='no_matched_transcript'
    def same_protein(x,y):
        if not x.get('hgvsp') or not y.get('hgvsp'):return True
        strs=[t['hgvsp'].replace('(','').replace(')','').replace('*','Ter') for t in (x,y)]
        if strs[0]==strs[1]:return True
        # A synonym reported as e.g. Gly12Gly vs Gly12= is equal only on the *same* accession and codon.
        normal=r'([^:]+:p\.)([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=)'
        mm=[re.fullmatch(normal,s) for s in strs]
        return bool(all(mm) and all(m[4] in (m[2],'=') for m in mm) and mm[0].groups()[:3]==mm[1].groups()[:3])
    wrong=tcar and any(t['hgvsc']!=tcar['hgvsc'] or not same_protein(t,tcar) for t in (tn,tm) if t)
    if wrong:identity=False;states['transcript_check']='disagreement'
    source='vep' if tn else 'mapped_vep' if tm else 'car' if tcar else None
    if source:
        t=dict(annotations[source]);f['transcript']=t;f['transcript_sources']={k:source for k,val in t.items() if val is not None}
        if not t.get('hgvsp') and tcar and tcar['hgvsc']==t.get('hgvsc') and tcar.get('hgvsp'):
            t['hgvsp']=tcar['hgvsp'];f['transcript_sources']['hgvsp']='car'
    else:t={}
    f['gnomad']=gn
    if auth:
        f['allele_identity']={'caid':auth[0],'gene':g,'method':'CAR exact input or locally verified equivalent haplotype','reference_status':states.get('ref_check','transcript_input')}
    if identity and t:
        # Only fill omitted fields on the same transcript and exactly reconciled c.; preserve field-level provenance.
        if tcar and tcar['hgvsc']==t.get('hgvsc') and same_protein(t,tcar):
            for key in ('protein_start','protein_end','consequence_terms'):
                if not t.get(key) and tcar.get(key):
                    t[key]=tcar[key];f['transcript_sources'][key]='car'
        named=(ev.get('ldh_consequence') or {}).get('verified_refseq_consequence')
        if not t.get('consequence_terms') and named:
            t['consequence_terms']=[named];f['transcript_sources']['consequence_terms']='ldh_consequence'
    if not identity:
        if ev:f['linked_observations_unverified']=dict(ev)
        r['limitations'].append('等位、REF 或规范转录本尚未核准；分别保留取得的注释，不计分。')
    else:
        r['status']='review_required';r['expert_assertions']=ev.pop('expert_assertions',[])
        f.update(ev)
        if v['kind']=='vcf':
            coord=[q for a in car.get('genomicAlleles',[]) if a.get('referenceGenome')==b and any(h.startswith(core.NC[b][v['chrom']]+':g.') for h in a.get('hgvs',[])) for q in a.get('coordinates',[])]
            if len(coord)==1:
                q=coord[0];iv=[q.get(k) for k in ('start','end','referenceAllele','allele')]
                f['genomic_reconciliation']={'input_interbase':v['interbase'],'car_interbase':iv,'method':'exact_coordinates' if v['interbase']==iv else 'equal_local_alternate_haplotypes'}
        chrom='17' if g=='BRCA1' else '13';coords=[q for a in car.get('genomicAlleles',[]) if a.get('referenceGenome')==b and any(h.startswith(core.NC[b][chrom]+':g.') for h in a.get('hgvs',[])) for q in a.get('coordinates',[])]
        if len(coords)==1:
            a,c=coords[0]['referenceAllele'],coords[0]['allele'];f['event']='INS' if not a else 'DEL' if not c else 'SNV' if len(a)==len(c)==1 else 'MNV' if len(a)==len(c) else 'DELINS'
        if workbooks:
            from .mechanism import ptc_hint
            hint=ptc_hint(g,t,core.row_index())
            if hint:f['stop_exon_review']=hint
            f['exon_workbook_context']=core.ptc_context(g,dict(t,exon=t.get('exon') or ''));matches=core.lookup(g,t['hgvsc']);states['workbook']='matched_rows' if matches else 'matched_stop_exon_only' if hint else 'no_matching_row'
            for x in matches:
                code=x.get('code') or x.get('assigned') or ''
                if re.match(r'^(PS3|BS3|PVS1|BP7)(?:$|\b|_)',code) and 'N/A' not in code:r['candidates'].append(dict(state='requires_review',source=x,code=code,strength=x.get('strength'),score_counted=False))
                else:f.setdefault('workbook_rows_without_assignment',[]).append(x)
    if not workbooks:states['workbook']='unavailable'
    else:states.setdefault('workbook','identity_unresolved')
    r['manual_review_prompts']=core.review_prompts(r);r['source_files']=links;r['source_status']=states
    writej(folder/'report.json',r);return r
