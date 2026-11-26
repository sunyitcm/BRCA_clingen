"""Fetch cited studies first and keep optional discovery on an independent budget."""
import json
import re
import urllib.parse as u
from xml.etree import ElementTree as ET

from .archive import replay_fetch
from .net import fetch, writej

EPMC='https://www.ebi.ac.uk/europepmc/webservices/rest/'
PUBMED='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?'


def collect(out,gene,aliases,pmids,user_limit,previous,citations,discover_limit,pool):
    from .literature import patterns,scan
    from .report import reading_html
    patterns(aliases)
    if gene not in ('BRCA1','BRCA2') or not 1<=user_limit<=25 or not 0<=discover_limit<=25 or not 1<=pool<=100:
        raise ValueError('invalid literature search scope/budget')
    allele=r'(?:[cpg]\.\(?)?(?:[-*]?\d[0-9A-Za-z_*+>\-]*|[A-Za-z]{1,3}\d+[A-Za-z*][A-Za-z*0-9]*)\)?'
    if any(not re.fullmatch(allele,a) for a in aliases) or any(not re.fullmatch(r'\d{1,10}',p) for p in [*pmids,*citations]):
        raise ValueError('联网查询只接受核酸、蛋白变异表达或 PMID')
    out.mkdir(parents=True,exist_ok=False)
    def acquire(label,url,kind='json'):
        if previous is not None:
            return replay_fetch(previous,out/'http',label,url,kind=kind)
        return fetch(out/'http',label,url,kind=kind,timeout=60 if label=='epmc_body' else 25)
    def query(label,url):
        raw=acquire(label,url)
        value=json.loads(raw) if raw else None
        return value if isinstance(value,dict) else {}
    def epmc(label,expr,size):
        return query(label,EPMC+'search?'+u.urlencode({'query':expr,'format':'json','resultType':'core','pageSize':size}))
    def records(response):
        container=response.get('resultList')
        rows=container.get('result') if isinstance(container,dict) else None
        if (isinstance(rows,list) and all(isinstance(row,dict) for row in rows)
            and type(response.get('hitCount')) is int and response['hitCount']>=len(rows)):
            return rows
        return None
    tasks={}
    def task(pmid,source):
        item=tasks.setdefault(pmid,{'pmid':pmid,'sources':[],'state':'pending','score_counted':False,'human_review':'not_reviewed'})
        if source not in item['sources']:item['sources'].append(source)
        return item

    for p,origins in citations.items():
        for origin in origins:task(p,origin)
    users=list(dict.fromkeys(pmids))
    for p in users:task(p,{'origin':'user_seed','role':'人工指定'})
    coverage={};queries={};metadata={};queue=[];fulltexts=[]
    required=set(citations)|set(users)
    selected=list(dict.fromkeys([*citations,*users[:user_limit]]))
    for p in users[user_limit:]:
        if p not in citations:tasks[p]['state']='not_processed_user_budget'
    if discover_limit:
        expr=gene+' AND ('+' OR '.join('"'+x+'"' for x in aliases)+')'
        queries['pubmed']=expr
        pub=query('pubmed',PUBMED+u.urlencode({'db':'pubmed','term':expr,'retmax':pool,'retstart':0,'sort':'relevance','retmode':'json'})).get('esearchresult',{})
        ok=isinstance(pub,dict) and str(pub.get('count','')).isdigit() and isinstance(pub.get('idlist'),list)
        pubids=[p for p in pub['idlist'] if isinstance(p,str) and re.fullmatch(r'\d{1,10}',p)] if ok else []
        coverage['pubmed']={'status':'ok' if ok else 'service_or_fields_failed','total':pub.get('count') if ok else None,'idlist':pubids,'pool_per_source':pool,'sort':'relevance','retstart':0,
                            **{k:pub.get(k) for k in ('querytranslation','warninglist','errorlist') if ok}}
        expr=gene+' AND ('+' OR '.join('BODY:"'+a+'"' for a in aliases)+')'
        queries['europepmc']=expr
        ep=epmc('epmc_body',expr,pool)
        rows=records(ep)
        ok=rows is not None
        epids=[]
        for row in rows if ok else []:
            p=str(row.get('pmid') or '')
            if re.fullmatch(r'\d{1,10}',p):
                epids.append(p)
                if row.get('source')=='MED' and str(row.get('id'))==p:metadata[p]=row
            else:queue.append({'id':row.get('id'),'state':'no_pmid_review'})
        coverage['epmc_body']={'status':'ok' if ok else 'service_or_fields_failed','total':ep.get('hitCount') if ok else None,
                              'idlist':epids,'pool_per_source':pool,'sort':'service_default_relevance','next_cursor':ep.get('nextCursorMark')}
        for origin,ids in [('discovery_epmc',epids),('discovery_pubmed',pubids)]:
            for rank,p in enumerate(ids,1):task(p,{'origin':origin,'rank':rank,'role':'关键词发现'})
        discovered=[p for p in dict.fromkeys(epids+pubids) if p not in required]
        selected+=discovered[:discover_limit]
        for p in discovered[discover_limit:]:tasks[p]['state']='not_processed_discovery_budget'
    # Exact MED lookups in bounded groups. Missing in a truncated/malformed reply is unknown.
    need=[p for p in selected if p not in metadata]
    for i in range(0,len(need),20):
        group=need[i:i+20]
        expr='SRC:MED AND ('+' OR '.join('EXT_ID:'+p for p in group)+')'
        data=epmc('ids_'+str(i//20+1),expr,len(group))
        rows=records(data)
        valid=rows is not None
        for row in rows if valid else []:
            p=str(row.get('pmid',''))
            if p in group and row.get('source')=='MED':metadata[p]=row
        for p in group:
            if p not in metadata:
                tasks[p]['state']='no_record' if valid and len(rows)==data['hitCount'] and all(r.get('pmid') in group and r.get('source')=='MED' for r in rows) else 'metadata_service_or_fields_failed'
    for p in selected:
        info=metadata.get(p)
        if not info:continue
        todo=tasks[p];pc=info.get('pmcid')
        todo.update(pmcid=pc,doi=info.get('doi'),title=info.get('title'))
        if not info.get('title') or info.get('isOpenAccess') not in ('Y','N'):
            todo['state']='metadata_fields_missing';continue
        if info['isOpenAccess']=='N' or (not pc and info.get('inPMC')=='N'):
            todo['state']='no_accessible_fulltext_here';continue
        if not re.fullmatch(r'PMC\d+',str(pc)):
            todo['state']='metadata_fields_missing';continue
        raw=acquire('full_'+p,EPMC+pc+'/fullTextXML','xml')
        if raw is None:
            todo['state']='fulltext_service_failure';continue
        try:
            rec=scan(raw,gene,aliases,p,pc)
        except (ValueError,ET.ParseError) as error:
            todo.update(state='id_or_xml_problem',reason=str(error));continue
        if info.get('doi') and (rec.get('doi') or '').casefold()!=info['doi'].casefold():
            todo['state']='doi_mismatch_or_missing';continue
        fulltexts.append(rec)
        todo.update(state=rec['state'],hits=len(rec['hits']),supplements='linked_not_scanned' if rec['supplement_links'] else 'not_listed')
    for item in tasks.values():
        if item['state']!='body_or_table_hits':queue.append({'pmid':item['pmid'],'state':item['state']})
    result={'schema':'brca-lite.reading-tasks.v1','gene':gene,'aliases':aliases,
            'queries':queries,'coverage':coverage,'reading_tasks':list(tasks.values()),'fulltexts':fulltexts,'review_queue':queue,
            'limit':user_limit,'budgets':{'user':user_limit,'guideline':'all','discovery':discover_limit,'pool_per_source':pool},
            'acquisition':'replayed_http' if previous else 'live_http','run_id':out.name}
    writej(out/'literature.json',result);reading_html(out/'reading.html',result)
    return result
