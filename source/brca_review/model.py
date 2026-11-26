"""Presentation adapter: never assigns an ACMG code or changes a frozen report."""
from __future__ import annotations
import base64, hashlib, json, re
from pathlib import Path
from ._engine.brcalite.core import GENES, R
from ._engine.brcalite.net import fetch, writej
HERE=Path(__file__).resolve().parent
A=HERE/'assets'
BUILDS=('GRCh37','GRCh38')
def read(p): return json.loads(Path(p).read_bytes())
def sha(b): return hashlib.sha256(b).hexdigest()
def b64(b): return base64.b64encode(b).decode('ascii')
def coord_for(b): return 'chr17-'+('41258504' if b=='GRCh37' else '43106487')+'-A-C'
def collect_portal(folder,report):
    gene=report.get('gene')
    if gene in GENES: fetch(folder/'http','cspec_current','https://cspec.genome.network/cspec/ui/svi/doc/'+GENES[gene][2],kind='html',attempts=1,timeout=12)
def portal_version(body,gene):
    if not body or gene not in GENES: return None
    text=body.decode('utf-8','replace')
    if 'id="lhs-'+GENES[gene][2]+'"' not in text: return None
    m=re.findall(r'Interpretation Guidelines for '+gene+r' Version\s+([0-9]+(?:\.[0-9]+)+)\s*<',text)
    return m[0] if len(set(m))==1 else None
def capture(folder,share,blobs):
    rows={}
    for f in sorted((folder/'http').glob('*/receipt.json')):
        m=read(f);att=(m.get('attempts') or [{}])[-1];key=f.parent.name
        fn=att.get('body_file');raw=(f.parent/fn).read_bytes() if fn else None
        digest=sha(raw) if raw is not None else None
        if digest!=att.get('body_sha256'): raise ValueError('响应校验失败：'+key)
        reqfile=f.parent/'request.json';rb=reqfile.read_bytes() if reqfile.is_file() else None
        if (sha(rb) if rb is not None else None) != m.get('request_sha256'): raise ValueError('请求体校验失败：'+key)
        include=raw is not None and (not share or key!='cspec_current')
        if include: blobs[digest]=b64(raw)
        rows[key]={'sha256':digest,'size':len(raw) if raw is not None else None,'body':bool(include),'method':m['method'],'url':m['url'],
        'request_sha256':m.get('request_sha256'),'request':None if rb is None else rb.decode('utf-8'),'utc':m.get('finished_utc'),'ok':bool(m['ok']),
        'status':att.get('status'),'error':att.get('error'),'attempts':len(m.get('attempts',[])),
        'content_type':next((v for k,v in att.get('response_headers',[]) if k.lower()=='content-type'),None),
        'note':('分发样例仅留指南网页标题核验、网址和原始哈希；本机联网运行另存整页。' if raw is not None and not include else None)}
    return rows
LABELS={'erepo_summary':'eRepo 专家摘要','erepo_detail':'eRepo 专家逐条证据','cspec_current':'CSpec 当前发布页','car':'CAR 等位身份','ref':'Ensembl 参考碱基','info':'VEP 主机数据 release','vep':'VEP 变异注释','ldh':'LDH 等位目录','ldh_med':'LDH 实体目录','ldh_consequence':'LDH 分子后果','ldh_exome':'LDH · gnomAD 外显子','ldh_genome':'LDH · gnomAD 基因组','workbook':'ENIGMA 规范附件','population_file':'本地人群证据','gno_nc2':'gnomAD v2.1 非癌外显子','gno_nc3':'gnomAD v3.1 非癌基因组','gno_full37':'gnomAD v2.1 全队列','gno_full38':'gnomAD v4 全队列','ldh_brca':'LDH · BRCA Exchange','clinvar':'ClinVar VCV / SCV'}
def version_rows(r,receipts,current):
    f=r['facts'];t=f.get('transcript') or {};co=f.get('ldh_consequence') or {};ex=f.get('ldh_exome') or {};bx=f.get('ldh_brca') or {};ns=GENES.get(r.get('gene'),('',None,None))[2]; assertn=r.get('expert_assertions') or []
    cs='https://cspec.genome.network/cspec/ui/svi/doc/'+ns if ns else None
    def rv(key,label,version,scope,meaning,keys,observed=None,source_date=None,link=None):
        return dict(key=key,label=label,version=version,scope=scope,meaning=meaning,raw=keys,fields=observed,
             source_date=source_date,url=link,utc=next((receipts[k]['utc'] for k in keys if k in receipts),None))
    build=r['input']['assembly']; local=r.get('specification')
    a=[rv('rules','本地 CSpec 规则',local,'冻结规范 · 胚系 BRCA1/2','用于解释候选判据；导入新规则须复核。',[],{'rule_lock':r.get('rules_bundle_sha256')},link=(cs+'?version=1.2' if cs else None)),
       rv('cspec_current','CSpec 官方发布页',current,'本次联网检查；离线打开不会更新','页面版本与冻结规范可不同；未可靠读出时显示未知。',['cspec_current'],link=cs)]
    releases=f.get('ensembl_metadata',{}).get('releases') or []
    a.append(rv('vep','Ensembl / VEP',' / '.join(map(str,releases)) or None,build+' · 本次主机','这是主机数据 release。未查询的插件及其版本不能从此数推断。',['info','vep','ref'],{k:t.get(k) for k in ('hgvsc','hgvsp','consequence_terms')} if t and (not f.get('source_transcripts') or 'vep' in f['source_transcripts']) else None))
    if f.get('source_transcripts') and 'vep' in f['source_transcripts']:a[-1]['fields']=f['source_transcripts']['vep']
    a.append(rv('car','ClinGen CAR',None,'跨组装等位身份','CA 是稳定标识；此请求未披露数据库发布号。',['car'],r.get('caid')))
    if f.get('source_transcripts',{}).get('car'):a[-1]['fields']={'caid':r.get('caid'),'transcript':f['source_transcripts']['car']}
    if 'mapped_vep' in receipts:
        a.append(rv('mapped_vep','Ensembl VEP · CAR 映射的 GRCh38 等位',' / '.join(map(str,(f.get('mapped_ensembl_metadata') or {}).get('releases',[]))) or None,'GRCh38 · 另行查询','同一规范转录本交叉核对；不是原 GRCh37 主机的输出。',['mapped_info','mapped_vep'],f.get('source_transcripts',{}).get('mapped_vep')))
    if 'refseq_ref' in receipts:
        rr=f.get('reference_check') or {}
        a.append(rv('refseq_ref','NCBI RefSeq · 原组装 REF',rr.get('accession'),build,'核对版本化染色体序列的指定位置与正链 REF。',['refseq_ref'],rr or None))
    for k,row in (f.get('gnomad') or {}).items():
        cohorts={a:({q:v.get(q) for q in ('ac','an','af','filters','faf95') if q in v} if isinstance(v,dict) else None) for a,v in row.get('cohorts',{}).items()}
        a.append(rv(k,'gnomAD · '+row['label'],row['dataset'],row['assembly'],'dataset 是本次接口选择；未披露的数据小版本不能用 UTC 补齐。未见变异记录不等于 AC=0。',[k],{'variant_id':row['variant_id'],'state':row['state'],'cohorts':cohorts},link=row['record_url']))
    sources=co.get('version_source') or []; sv=[str(s.get('source',''))+' '+str(s['ver']) for s in sources if isinstance(s,dict) and s.get('ver') is not None]
    a.append(rv('ldh_co','LDH · 分子后果',' / '.join(sv) or None,' / '.join(str(s.get('genomeBuild','')) for s in sources if isinstance(s,dict)) or '组装未披露','版本仅属于这一等位的这一实体；不能外推给所有 LDH 记录。',['ldh_consequence'],{'linked_allele':co.get('linked_allele'),'revision':co.get('revision')} if co else None,co.get('modified')))
    locus=ex.get('source') or {}; a.append(rv('exome','LDH · gnomAD exome',locus.get('label'),ex.get('_build') or '本次实体未确认组装','AC/AN/AF 属于这一数据集。未取得规范要求的队列，不能据此自动给 PM2。',['ldh_exome'],{k:ex.get(k) for k in ('AC','AN','AF','quality')} if ex else None))
    gg=f.get('ldh_genome') or {}
    if gg or 'ldh_genome' in receipts:
        a.append(rv('genome','LDH · gnomAD genome',(gg.get('source') or {}).get('label'),gg.get('_build') or '未核准','与外显子分开。这里的 AC=0 只表示此数据集给出的计数；不能替代规范队列和覆盖。',['ldh_genome'],{k:gg[k] for k in ('AC','AN','AF','noVariant','quality') if k in gg} or None))
    bs=bx.get('version_source') or {}
    a.append(rv('exchange','LDH · BRCA Exchange',('releaseName '+str(bs['releaseName'])+' / ID '+str(bs.get('releaseId')) if 'releaseName' in bs else None),'等位实体 · '+str(bx.get('linked_allele') or '未确认'),'数据发布日与 LDH 实体修改日不同。预测分数不等于功能实验。',['ldh_brca'],{'source_release_date':bs.get('releaseDate'),'BayesDel_noAF':bx.get('BayesDel_nsfp33a_noAF')} if bx else None,bx.get('modified')))
    detail=f.get('brca_detail') or {}; row=detail.get('source_row') or {}; linked=detail.get('ldh_linked_row') or {}
    if detail or 'brca_detail' in receipts:
        msg='逐行核基因、转录本和两组装。数据库中的先验、后验及家族合计不能冒充共分离 LR。'
        if linked.get('issues'):
            msg+=' LDH 对应历史行跨组装字段不一致，已隔离：'+', '.join(linked['issues'])+'。'+('另列通过独立核验的版次；不混合两行。' if row else '未使用历史行的 LR。')
        a.append(rv('brca_detail','BRCA Exchange · 完整历史记录',('release '+str(row['release'])+' / ID '+str(row['release_id']) if row else None),'只属于此次返回并核准的行',msg,['brca_detail'],{k:detail.get(k) for k in ('method','source_row','ldh_linked_row','latest_returned','measures','history')} if detail else None,row.get('release_date')))
    for buildname,short in [('GRCh37','37'),('GRCh38','38')]:
        if 'anchor_'+short in receipts:a.append(rv('anchor_'+short,'NCBI RefSeq · '+buildname+' VCF 左锚',None,buildname,'按 CAR 此组装间碱基坐标核原参考及锚；不能替代另一组装或凭固定偏移换算。',['anchor_'+short]))
    a.append(rv('clinvar','ClinVar · 专家记录','；'.join(sorted({str(x['vcv'])+' / '+str(x['scv']) for x in assertn})) or None,'聚合 VCV 与提交 SCV 的记录版本','保留提交方、疾病与原评估日；历史断言不计为独立判据。',['clinvar'],[{k:x.get(k) for k in ('scv','class','review','condition','submitter','evaluated')} for x in assertn] or None,' / '.join(sorted({x['evaluated'] for x in assertn if x.get('evaluated')})) or None))
    record=f.get('clinvar_record') or {}
    if record.get('accession'):
        row=a[-1];row['label']='ClinVar · 等位记录／专家提交';row['version']=record['accession']+('；' + row['version'] if row['version'] else '')
        row['source_date']=record.get('updated');row['fields']={'vcv':record['accession'],'expert_submissions':row['fields'] or []}
        if record.get('germline_aggregate'):row['fields']['aggregate']=record['germline_aggregate']
        row['meaning']='VCV 是此等位的记录版次；专家 SCV 另列。没有合格专家提交时仍有 VCV。两者都不自动计分。'
    ep=r.get('erepo') or {}
    if ep:
        a.append(rv('erepo','ClinGen eRepo · ENIGMA',ep.get('record_version'),r.get('condition'),'个案断言；按 CA、基因、疾病与专家组匹配。可明确选择引用，单列为专家来源；不充当本地独立取数。',['erepo_summary','erepo_detail'],{'state':ep.get('state'),'classification':ep.get('classification'),'spec_version':ep.get('spec_version'),'codes':[(y.get('code'),y.get('state'),y.get('strength')) for y in ep.get('criteria',[])]},ep.get('published'),ep.get('url')))
    return a
def pack(folder,ident,share=False):
    r=read(folder/'report.json')
    blobs={};raw=capture(folder,share,blobs)
    pd=folder/'http/cspec_current/body';cur=portal_version(pd.read_bytes() if pd.is_file() else None,r.get('gene'))
    rows=version_rows(r,raw,cur)
    for v in rows:
        if v['key']=='exome' and r['facts'].get('ldh_exome') and (folder/'http/ldh_exome/body').exists():
            rg=read(folder/'http/ldh_exome/body').get('locus',{}).get('reference_genome');v['scope']=str(rg) if rg else '未披露'
        if v['key']=='rules':
            v['utc']=None
    from ._engine.brcalite.reading_sources import assignments
    roles=assignments(r);locks=read(R/'rules/attachments.lock.json')
    for x in r.get('candidates',[]):
        s=x['source'];tname=('Specifications_Table9' if s.get('table')=='T9' else 'Specifications_Table4')
        at=next((l for l in locks if l['gene']==r['gene'] and l['name'].startswith(tname)),{})
        x['workbook']={k:at.get(k) for k in ('url','sha256','utc','name')}
        if s.get('review_values'):x['verified_cells']=s['review_values']
        x['pmid_roles']={p:'；'.join(z['role'] for z in entries if z.get('sheet')==s.get('sheet') and z.get('row')==s.get('row')) for p,entries in roles.items() if p in s.get('pmids_mentioned',[])}
    def ts(s): return max((v.get('utc') or '' for v in s.values()),default='') or r.get('generated_utc')
    return dict(id=ident,utc=ts(raw),assembly=r['input']['assembly'],original=r['input']['original'],gene=r.get('gene'),caid=r.get('caid'),status=r['status'],
        report=r,versions=rows,raw=raw,blobs=blobs,report_sha256=sha((folder/'report.json').read_bytes()),raw_report=b64((folder/'report.json').read_bytes()),diff=None)
def brief(s): return {k:s[k] for k in ('id','utc','assembly','original','gene','caid','status')} | {'gene':s.get('report',{}).get('selected_gene') or s.get('gene')}
def compare(old,new):
    if not old: return {'message':'首次查询或未选基线。','rows':[]}
    if (not new['caid'] or old['caid']!=new['caid'] or old['assembly']!=new['assembly'] or any(s['status']!='review_required' for s in (old,new))):
        return {'message':'组装或已核准等位不同，或身份未通过；不做证据比较。','rows':[]}
    oldv={r['key']:r for r in old['versions']};res=[]
    for v in new['versions']:
        prev=oldv.get(v['key']);k=v['key']
        if not prev: continue
        vals=lambda a:{n:a[n] for n in ('version','scope','fields','source_date')}
        flags=[];changed=[]
        if vals(prev)!=vals(v): changed=[n for n in vals(v) if vals(prev)[n]!=vals(v)[n]]
        for n in v['raw']:
            if n=='ref' and n not in new['raw'] and n not in old['raw']: continue
            l,h=old['raw'].get(n),new['raw'].get(n)
            if not l or not h or not l['ok'] or not h['ok']: flags.append('不可核对');continue
            if any(l[m]!=h[m] for m in ('method','url','request_sha256')): flags.append('请求不同');continue
            flags.append('原文相同' if l['sha256']==h['sha256'] else '原文变化')
        incomplete=('不可核对' in flags or (k=='cspec_current' and not v['version']))
        res.append(dict(label=v['label'],state=('信息不全' if incomplete else '可见字段变化' if changed else '可见字段相同'),raw=' / '.join(dict.fromkeys(flags)) or '本地锁',changes={n:{'old':prev[n],'new':v[n]} for n in changed}))
    return {'message':'对照 '+old['id']+' · '+(old.get('utc') or ('公开教学示例，不含采集时间' if old.get('teaching') else '取数时间未提供'))+'。同 CA、同组装；仅比较列出的字段，原文变化不等于致病性变化。','rows':res}

LABELS.update(brca_detail='BRCA Exchange 分版完整记录',anchor_37='NCBI GRCh37 参考锚',anchor_38='NCBI GRCh38 参考锚',refseq_ref='NCBI RefSeq 原组装参考',mapped_vep='VEP（映射 GRCh38）',mapped_info='映射主机版本',gno_full37='gnomAD v2.1 全队列',gno_full38='gnomAD v4 全队列',gno_nc2='gnomAD v2.1 非癌外显子',gno_nc3='gnomAD v3.1 非癌基因组',ref_check='原组装 REF 核对',workbook='冻结工作簿')
