"""Two-stage retrieval: citation discovery and position-addressable JATS body evidence."""
import hashlib,json,re,unicodedata,urllib.parse as u
from xml.etree import ElementTree as ET
from pathlib import Path
from .net import fetch,writej,sha
EPMC='https://www.ebi.ac.uk/europepmc/webservices/rest/'
AMINO=dict(zip('Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val Ter'.split(),'ARNDCQEGHILKMFPSTWYV*'))
def norm(s):return ' '.join(unicodedata.normalize('NFKC',s).replace('→','>').replace('−','-').split())
def text(n):return norm(''.join(n.itertext())) if n is not None else ''
def patterns(aliases):
 if not aliases or any(not a or len(a)>120 or '"' in a or '\n' in a for a in aliases):raise ValueError('provide explicit validated aliases')
 return [(a,re.compile(r'(?<![A-Za-z0-9])'+re.escape(norm(a))+r'(?![A-Za-z0-9])',re.I)) for a in aliases]
def scan(raw,gene,aliases,pmid=None,pmcid=None,*,context=True):
 rx=patterns(aliases);root=ET.fromstring(raw)
 for n in root.iter():n.tag=n.tag.rsplit('}',1)[-1]
 ids={n.get('pub-id-type'):n.text for n in root.findall('./front/article-meta/article-id')}
 if pmid and ids.get('pmid')!=str(pmid):raise ValueError('PMID mismatch or unavailable in article')
 if pmcid and str(ids.get('pmcid') or ids.get('pmc') or '').removeprefix('PMC')!=str(pmcid).removeprefix('PMC'):raise ValueError('PMCID mismatch')
 pmid=ids.get('pmid');pmcid=ids.get('pmcid') or ids.get('pmc');pmcid=('PMC'+str(pmcid).removeprefix('PMC')) if pmcid else None
 body=root.find('body'); lic=root.find('./front/article-meta/permissions/license')
 url_license=lic.get('{http://www.w3.org/1999/xlink}href','') if lic is not None else ''
 licensed=bool(re.match(r'^https?://creativecommons\.org/(licenses/by/|publicdomain/(zero|mark)/)',url_license))
 authors=root.findall('./front/article-meta/contrib-group/contrib/name/surname')
 obj={'pmid':pmid,'pmcid':pmcid,'doi':ids.get('doi'),'title':text(root.find('./front/article-meta/title-group/article-title')),'first_author':text(authors[0]) if authors else None,'xml_sha256':sha(raw),'license_url':url_license,'pubmed_url':('https://pubmed.ncbi.nlm.nih.gov/'+pmid+'/') if pmid else None,'hits':[],'supplement_links':[],'state':'fulltext_no_hit','distinct_patients':None}
 if body is None or not text(body):obj['state']='body_unavailable';return obj
 parent={x:n for n in body.iter() for x in n};article_text=text(body);gene_rx=re.compile(r'\b'+re.escape(gene)+r'\b',re.I)
 for i,node in enumerate(body.iter(),1):
  if node.tag not in ('p','tr'):continue
  chain=[node]
  while chain[-1] in parent:chain.append(parent[chain[-1]])
  if any(a.tag in ('ref-list','supplementary-material','sub-article') for a in chain):continue
  if node.tag=='p' and any(a.tag in ('tr','td','th') for a in chain):continue
  if node.tag=='p' and any(n.tag in ('table','ref-list','supplementary-material') for n in node.iter()):continue
  para=(' | '.join(text(a) for a in node if a.tag in ('td','th')) if node.tag=='tr' else text(node))
  hits=[a for a,pat in rx if pat.search(para)]
  if not hits:continue
  sec=next((a for a in chain if a.tag=='sec'),None);tab=next((a for a in chain if a.tag in ('table-wrap','fig')),None)
  sec_titles=[text(a.find('title')) for a in reversed(chain) if a.tag=='sec']
  anchor=node.get('id') or (tab.get('id') if tab is not None else None) or (sec.get('id') if sec is not None else '')
  hit={'aliases':hits,'section':' / '.join(sec_titles),'type':node.tag,'body_node_ordinal':i,'xml_id':anchor,'paragraph_sha256':sha(para.encode()),'unit_label':text(tab.find('label')) if tab is not None else None,'gene_context':'same_unit' if gene_rx.search(para) else ('other_body_unit' if gene_rx.search(article_text) else 'unverified'),'rna_protein_or_dna_match':'manual','curation_state':'candidate','review_url':('https://pmc.ncbi.nlm.nih.gov/articles/'+pmcid+'/'+('#'+u.quote(anchor) if anchor else '')) if pmcid else obj['pubmed_url']}
  if licensed:hit['paragraph']=para
  else:hit['excerpt_state']='link_only_license_not_allowlisted'
  if context:
   flags=[]
   if re.search(r'\b'+('BRCA2' if gene=='BRCA1' else 'BRCA1')+r'\b',para,re.I):flags.append('同一单元出现另一 BRCA 基因，需核对等位归属')
   if re.search(r'\b(?:positive|negative) control\b',para,re.I):flags.append('包含实验对照用语')
   if re.search(r'\b(?:not (?:detected|observed|found)|no carriers?)\b',para,re.I):flags.append('包含否定或未检出用语')
   if hit['gene_context']!='same_unit':flags.append('本单元没有目标基因名称')
   hit.update(review_flags=flags,score_counted=False,match_type='核酸名称' if any(a.startswith('c.') for a in hits) else '仅别名，需核对核酸')
   if licensed and node.tag=='tr' and tab is not None:
    hit['table_caption']=text(tab.find('caption'))
    hit['table_headers']=[[text(c) for c in row if c.tag in ('td','th')] for row in tab.findall('.//thead/tr')]
    if not hit['table_headers']:hit['review_flags'].append('未取得明确表头，单位和分母需查原表')
  obj['hits'].append(hit)
 for n in root.iter('supplementary-material'):
  for a in n.iter():
   if a.tag not in ('media','graphic','supplementary-material','inline-supplementary-material'):continue
   h=a.get('{http://www.w3.org/1999/xlink}href')
   if h and h not in obj['supplement_links']:obj['supplement_links'].append(h)
 if obj['hits']:obj['state']='body_or_table_hits'
 return obj
def _legacy_discover(out,gene,aliases,pmids,limit,previous):
 if gene not in ('BRCA1','BRCA2') or not 1<=limit<=25:raise ValueError('gene/limit invalid')
 patterns(aliases)
 # Only public allele expressions may become remote search terms.
 allele_pattern=r'(?:[cpg]\.\(?)?(?:[-*]?\d[0-9A-Za-z_*+>\-]*|[A-Za-z]{1,3}\d+[A-Za-z*][A-Za-z*0-9]*)\)?'
 if any(not re.fullmatch(allele_pattern,a) for a in aliases):raise ValueError('联网别名须为核酸或蛋白变异表达')
 if any(not re.fullmatch(r'\d{1,10}',x) for x in pmids):raise ValueError('invalid PMID')
 out.mkdir(parents=True,exist_ok=False);hits=[];queue=[];cands={};caps={}
 def acquire(label,url,kind='json'):
  if previous is not None:
   from .archive import replay_fetch
   return replay_fetch(previous,out/'http',label,url,kind=kind)
  return fetch(out/'http',label,url,kind=kind,timeout=60 if label=='epmc_body' else 25)
 def j(label,url):
  r=acquire(label,url)
  return json.loads(r) if r else {}
 eutils='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?'
 expr=gene+' AND ('+' OR '.join('"'+x+'"' for x in aliases)+')'
 pub=j('pubmed',eutils+u.urlencode({'db':'pubmed','term':expr,'retmax':limit,'retmode':'json'})).get('esearchresult',{})
 caps['pubmed']={'total':pub.get('count'),'download_cap':limit,'status':'ok' if pub else 'failed'}
 def meta(label,query):
  return j(label,EPMC+'search?'+u.urlencode({'query':query,'format':'json','resultType':'core','pageSize':limit}))
 query=gene+' AND ('+' OR '.join('BODY:"'+a+'"' for a in aliases)+')'
 ep=meta('epmc_body',query);rows=ep.get('resultList',{}).get('result',[])
 caps['epmc_body']={'total':ep.get('hitCount'),'next_cursor':ep.get('nextCursorMark'),'download_cap':limit,'status':'ok' if 'hitCount' in ep else 'failed'}
 for it in rows:
  if it.get('pmid'):cands[str(it['pmid'])]=it
  else:queue.append({'id':it.get('id'),'state':'no_pmid_review','pmcid':it.get('pmcid')})
 requested=list(dict.fromkeys(pmids+list(cands)+pub.get('idlist',[])))
 for i,p in enumerate(requested[:limit]):
  it=cands.get(p)
  if not it:
   m=meta('map_'+p,'EXT_ID:'+p+' AND SRC:MED');it=next((x for x in m.get('resultList',{}).get('result',[]) if x.get('pmid')==p),{})
  pc=it.get('pmcid');base={'pmid':p,'pmcid':pc,'pubmed_url':'https://pubmed.ncbi.nlm.nih.gov/'+p+'/'}
  if not it:
   queue.append(base|{'state':'metadata_unavailable'});continue
  if pc and not re.fullmatch(r'PMC[0-9]+',pc):
   queue.append(base|{'state':'invalid_pmcid'});continue
  if not pc or it.get('isOpenAccess')!='Y':queue.append(base|{'state':'metadata_only_review_fulltext_rights'});continue
  raw=acquire('full_'+p,EPMC+pc+'/fullTextXML',kind='xml')
  if not raw:queue.append(base|{'state':'fulltext_fetch_failed'});continue
  try:rec=scan(raw,gene,aliases,p,pc,context=False)
  except (ValueError,ET.ParseError) as e:queue.append(base|{'state':'id_or_xml_problem','reason':str(e)});continue
  if it.get('doi') and rec.get('doi') and it['doi'].casefold()!=rec['doi'].casefold():
   queue.append(base|{'state':'doi_mismatch'});continue
  hits.append(rec)
 queue.extend({'pmid':p,'state':'not_processed_cap'} for p in requested[limit:])
 result={'queries':{'pubmed':expr,'europepmc':query},'coverage':caps,'fulltexts':hits,'review_queue':queue,'limit':limit,'completeness':'bounded retrieval; supplements linked, not searched automatically'}
 result.update(gene=gene,aliases=aliases,run_id=out.name,acquisition='replayed_http' if previous else 'live_http')
 from .report import reading_html
 writej(out/'literature.json',result);reading_html(out/'reading.html',result);return result

def discover(out,gene,aliases,pmids,limit,previous=None,*,citations=None,discovery_limit=None,pool=20):
 from .net import getj
 from .paper_queue import collect
 if previous and getj(previous/'literature.json').get('schema')!='brca-lite.reading-tasks.v1':
  return _legacy_discover(out,gene,aliases,list(dict.fromkeys([*pmids,*(citations or {})])),limit,previous)
 if discovery_limit is None:discovery_limit=0 if pmids or citations is not None else limit
 return collect(out,gene,aliases,pmids,limit,previous,citations or {},discovery_limit,pool)

def report_literature(report,out,pmids,limit,previous=None,discovery_limit=0,pool=20):
 """Start a reading queue from reconciled allele expressions and cited table rows."""
 tx=report['facts']['transcript']
 cdna=tx['hgvsc'].split(':',1)[1]
 aliases=[cdna]
 protein=(tx.get('hgvsp') or '').partition(':p.')[2].strip('()')
 if protein and re.fullmatch(r'[A-Z][a-z]{2}\d+[A-Za-z0-9*]+',protein):
  aliases.append(protein)
  missense=re.fullmatch(r'([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})',protein)
  if missense and missense[1] in AMINO and missense[3] in AMINO:
   aliases.append(AMINO[missense[1]]+missense[2]+AMINO[missense[3]])
 from .reading_sources import assignments
 citations=assignments(report)
 seeds=list(dict.fromkeys([*pmids,*citations]))
 result=discover(out,report['gene'],aliases,pmids,limit,previous,citations=citations,discovery_limit=discovery_limit,pool=pool)
 return {'reading_card':'papers/reading.html','pmids_seeded':seeds,'aliases':aliases,
         'passages':sum(len(x['hits']) for x in result['fulltexts']),
         'papers_scanned':sum(x.get('state')!='body_unavailable' for x in result['fulltexts']),'coverage':result['coverage'],
         'reading_tasks':result.get('reading_tasks',[]),'budgets':result.get('budgets',{}),'score_counted':False}
