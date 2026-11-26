"""Read one selected PMID. Article retrieval is independent of clinical scoring."""
import datetime as dt, json, re, secrets
from pathlib import Path
from urllib.parse import quote
from xml.etree import ElementTree as ET
from .model import A, read, sha
from ._engine.brcalite.literature import AMINO, scan
from ._engine.brcalite.paper_queue import collect
from ._engine.brcalite.reading_sources import assignments
from ._engine.brcalite.net import utc, writej

EV=r'e-[0-9]{8}T[0-9]{6}Z-[a-f0-9]{10}'
PID=r'[1-9][0-9]{0,9}'
def eid():return 'e-'+dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+secrets.token_hex(5)
def links(pmid,pc=None,doi=None,supp=None):
    out={'PubMed':'https://pubmed.ncbi.nlm.nih.gov/'+pmid+'/','Europe PMC 题录':'https://europepmc.org/article/MED/'+pmid}
    if re.fullmatch(r'PMC[0-9]+',pc or ''):out['PMC 全文']='https://pmc.ncbi.nlm.nih.gov/articles/'+pc+'/'
    if isinstance(doi,str) and re.fullmatch(r'10\.[0-9]{4,9}/[^\s<>]+',doi):
        out['期刊 / DOI']='https://doi.org/'+quote(doi,safe='/():._-')
    if supp and supp.get('url') and supp.get('source')=='GEO GSE117159, Findlay et al., Supplementary Table 1':
        out['GEO 补充资料']='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117159'
    return out
def allele_terms(r):
    tx=(r.get('facts') or {}).get('transcript') or {};hg=tx.get('hgvsc') or ''
    if ':c.' not in hg or r.get('gene') not in ('BRCA1','BRCA2'):raise ValueError('先核对 BRCA 等位及版本化转录本。')
    a=[hg.split(':',1)[1]];prot=(tx.get('hgvsp') or '').partition(':p.')[2].strip('()')
    if re.fullmatch(r'[A-Z][a-z]{2}\d+[A-Za-z0-9*]+',prot):
        a.append(prot);m=re.fullmatch(r'([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})',prot)
        if m and m[1] in AMINO and m[3] in AMINO:a.append(AMINO[m[1]]+m[2]+AMINO[m[3]])
    return a
def projected(receipt):
    r=read(receipt);aa=r.get('attempts',[])
    return {'url':r['url'],'utc':r['finished_utc'],'ok':r.get('ok',False),'state':r.get('state'),'retry_at_utc':r.get('retry_at_utc') or (aa[-1].get('retry_at_utc') if aa else None),'attempts':[{'status':a.get('status'),'error':a.get('error'),'retry_at_utc':a.get('retry_at_utc'),'seconds':a.get('elapsed_s'),'body_sha256':a.get('body_sha256')} for a in aa]}
def stored_xml(root,event,p):
    if not re.fullmatch(EV,event or '') or not re.fullmatch(PID,p or ''):raise ValueError('文献索引格式错误')
    return root/'reading_events'/event/'query/http'/('full_'+p)/'body'
def verified_cache(root,p,gene,aliases):
    path=root/'article_cache'/(p+'.json')
    if not path.exists():return None
    info=read(path);data=stored_xml(root,info['xml_event'],p).read_bytes()
    if info['pmid']!=p or sha(data)!=info['sha256']:raise ValueError('已存正文哈希不符')
    rec=scan(data,gene,aliases,p,info['pmcid'])
    if rec['state']=='body_unavailable' or (rec.get('doi') or '').casefold()!=(info.get('doi') or '').casefold():raise ValueError('已存文章身份不符')
    return (info,rec)
def lookup(root,s,pmid,refresh=False,collect_fn=None):
    if not isinstance(pmid,str) or not re.fullmatch(PID,pmid) or type(refresh)!=bool:raise ValueError('请输入有效 PMID')
    r=s['report']
    if r.get('status') not in ('review_required','draft_reviewed') or not r.get('caid'):raise ValueError('等位身份未核准，不能归属文献。')
    from .presentation import origins as source_roles
    aliases=allele_terms(r);gene=r['gene'];origins=source_roles(r).get(pmid,[]);sup=next((x['supplement'] for x in origins if x.get('supplement')),None)
    ident=eid();d=root/'reading_events'/ident;d.mkdir(parents=True,exist_ok=False)
    base={'id':ident,'snapshot':s['id'],'pmid':pmid,'gene':gene,'allele':s['caid'],'assembly':s['assembly'],'aliases':aliases,'origins':origins or [{'origin':'user_seed','role':'人工指定，待核相关性'}],
          'created_utc':utc(),'source_utc':None,'score_counted':False,'network_state':'not_requested','receipts':{},'article':None,'xml_event':None,'supplement':sup}
    cache=None
    try:cache=verified_cache(root,pmid,gene,aliases)
    except (OSError,ValueError,TypeError,KeyError,ET.ParseError) as ex:base['cache_warning']='已有缓存未通过完整性或身份检查；不展示其原文。'
    rec=None;task={}
    if cache and not refresh:
        info,rec=cache;base.update(state='verified_cache',source_utc=info['utc'],xml_event=info['xml_event'])
    else:
        result=(collect_fn or collect)(d/'query',gene,aliases,[pmid],1,None,{},0,1)
        task=next((t for t in result.get('reading_tasks',[]) if t['pmid']==pmid),{})
        state=task.get('state','metadata_service_or_fields_failed');base.update(state=state,network_state=state)
        receipts={f.parent.name:projected(f) for f in sorted((d/'query/http').glob('*/receipt.json'))};base['receipts']=receipts
        rec=next((v for v in result.get('fulltexts',[]) if v.get('pmid')==pmid and v.get('state') in ('fulltext_no_hit','body_or_table_hits')),None)
        if rec:
            raw=stored_xml(root,ident,pmid).read_bytes()
            if sha(raw)!=rec['xml_sha256']:raise ValueError('已获取正文哈希不符')
            rr=receipts.get('full_'+pmid,{})
            if not rr.get('ok') or not rr['attempts'] or rr['attempts'][-1].get('body_sha256')!=rec['xml_sha256']:raise ValueError('正文与来源回执不符')
            base.update(xml_event=ident,source_utc=rr['utc'])
            if re.match(r'^https?://creativecommons\.org/(licenses/by/|publicdomain/(zero|mark)/)',rec.get('license_url') or ''):
                cp=root/'article_cache';cp.mkdir(exist_ok=True)
                tmp=cp/(ident+'.tmp');writej(tmp,{'pmid':pmid,'pmcid':rec['pmcid'],'doi':rec.get('doi'),'utc':rr['utc'],'xml_event':ident,'sha256':rec['xml_sha256']});tmp.replace(cp/(pmid+'.json'))
        elif cache and state in ('fulltext_service_failure','metadata_service_or_fields_failed') and (not task.get('pmcid') or task['pmcid']==cache[0]['pmcid']) and (not task.get('doi') or task['doi'].casefold()==(cache[0].get('doi') or '').casefold()):
            info,rec=cache;base.update(state='failed_live_with_cache',xml_event=info['xml_event'],source_utc=info['utc'])
    base['article']=rec
    pc=(rec or task).get('pmcid');doi=(rec or task).get('doi')
    base.update(links=links(pmid,pc,doi,sup),title=(rec or task).get('title'),finished_utc=utc())
    writej(d/'event.json',base)
    return base
