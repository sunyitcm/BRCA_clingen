"""Case-level observations and human reading aids; never award points here."""
import re
from urllib.parse import urlencode
from ._engine.brcalite.reading_sources import assignments

STATE={
 'does_not_meet':('not_met','本次不符合'), 'not_applicable':('not_applicable','本位点不适用'),
 'branch_not_applicable':('partial','这一分支不适用'), 'branch_threshold_met':('partial','分支达标；其余条件待核'),
 'branch_threshold_not_met':('partial','此分支未达标'), 'searched_no_row':('no_record','已查该表，未见对应行'),
 'incomplete':('partial','资料不足以确认满足'), 'observed':('partial','有观察值，条件待核'),
 'lookup_failed':('unavailable','来源未查成，无法判断'), 'lookup_incomplete':('unavailable','取数不完整，无法判断'),
 'expert_not_met':('external_not_met','专家标为不满足'), 'source_conflict':('conflict','来源数值冲突，暂缓计分'), 'identity_unresolved':('unverified','等位待核，不可计分')}
def assessment(z,ok):
    v=z['evidence_view'];b=z['brief'];op=b.get('option') or {};k=op.get('kind');refs=[]
    def src(key,label,pointer=None,url=None,kind='observation'):
        if key and not any((x['key'],x.get('url'))==(key,url) for x in refs):refs.append(dict(key=key,label=label,pointer=pointer,url=url,kind=kind))
    for c in z['candidates']:
        src('workbook',str(c.get('sheet') or 'ENIGMA 表')+' 第 '+str(c.get('row'))+' 行',c.get('pointer'),c.get('workbook_url'))
    if b.get('status')=='conflict':state,title='conflict','来源冲突，暂缓计分'
    elif not ok:state,title='unverified','等位待核，不可计分'
    elif v['state'] in ('does_not_meet','not_applicable') or v['blocks_automatic_scoring']:state,title=STATE.get(v['state'],('pending','待核，不计分'))
    elif k=='expert_citation':state,title='external_met','专家采用 · 本地未独立证实'
    elif k:state,title='supported',{'workbook':'符合 · 附件明确赋码','prediction':'符合 · 已达到所用规则','population':'符合 · 规范队列已核准'}.get(k,'符合 · 有可追溯依据')
    elif b.get('status')=='expert_not_met':state,title='external_not_met','专家标为不满足'
    else:state,title=STATE.get(v['state'],('not_assessed','未取得可评资料'))
    for x in v['observations']+z['facts']:
        key=x.get('source')
        if key=='workbook' and any(s['key']==key for s in refs):continue
        src(key,key,x.get('pointer'),kind=x.get('kind','observation'))
    if z.get('expert'):src('erepo_detail','eRepo 个案判据',None)
    # The case comparator and local adoption are independent: a manual decision cannot rewrite the retrieval result.
    return dict(state=state,label=title,
        applicability=v.get('applicability') or ('not_applicable' if state=='not_applicable' else 'undetermined' if state in ('unverified','not_assessed','unavailable','partial') else 'applicable'),
        locally_meets=True if state=='supported' else False if state=='not_met' else None,
        clinical_lr=v.get('clinical_lr'),hold_reason=v.get('hold_reason'),comparison=v.get('comparison'),standard=v.get('standard'),sources=refs,
        publication=z['rule_url'],selected=(z.get('decision') or {}).get('state'),
        option_strength=op.get('strength'),option_origin=k)

def origins(r):
    result={str(p):list(rows) for p,rows in assignments(r).items() if re.fullmatch(r'[1-9]\d{0,9}',str(p))}
    er=r.get('erepo') or {}
    if er.get('state')=='matched':
        for row in er.get('criteria') or []:
            from .erepo import literature_mentions
            for ref in literature_mentions(row.get('text') or [],row.get('pmids') or ()):
                role=ref['role'];p=ref['pmid'];msg={'explicitly_excluded':'专家原文明确排除本研究，不作为该条的支持；仍保留供复核。','explicitly_included':'专家原文明确纳入本研究；属于该判据所用文献，不另加分。','mentioned':'个案说明提及；具体采用与否须看专家原文。'}[role]
                result.setdefault(p,[]).append(dict(origin='erepo',criterion=row.get('code'),use=role,role='eRepo '+str(row.get('code') or '')+'：'+msg,pointers=ref['pointers']))
    for q in ((r.get('facts') or {}).get('clinvar_record') or {}).get('germline_citations') or []:
        result.setdefault(q['pmid'],[]).append(dict(origin='clinvar_scv',role='ClinVar '+q['scv']+' 提交引用；尚未确认实验或队列是否适用',locator=q['locator']))
    return result

def reading_guide(r,verified):
    from .papers import allele_terms
    g=r.get('selected_gene') or r.get('gene');aliases=[]
    if verified:
        try:aliases=[a for a in allele_terms(r) if re.fullmatch(r'[A-Za-z0-9+.*_=>?\-]{1,160}',a)]
        except ValueError:pass
    p='https://pubmed.ncbi.nlm.nih.gov/';m='https://pmc.ncbi.nlm.nih.gov/search/'
    result=dict(allele_verified=bool(verified),aliases=aliases,links={},records=origins(r) if verified else {},note='人工检索链接不会自动发起检索。核对核酸 ALT、转录本、实验对照和队列重叠；论文数量不决定判据强度。')
    if not aliases or g not in ('BRCA1','BRCA2'):
        result['links']={'PubMed 入口（无已核准别名）':p,'PMC 全文入口':m};return result
    ti=lambda a:'"'+a+'"[Title/Abstract]'
    term=ti(g)+' AND ('+' OR '.join(ti(a) for a in aliases)+')'
    full=g+' AND ('+' OR '.join('"'+a+'"' for a in aliases)+')'
    result['links']={'PubMed · 按此等位':p+'?'+urlencode({'term':term}),'PubMed · 限功能或 RNA 线索':p+'?'+urlencode({'term':term+' AND (functional[Title/Abstract] OR assay[Title/Abstract] OR splicing[Title/Abstract] OR RNA[Title/Abstract] OR "saturation genome editing"[Title/Abstract])'}),'PMC · 搜开放全文':m+'?'+urlencode({'term':full})}
    return result

def view_only(s):
    """Compatibility overlay for old sealed revisions; never change their report or hash."""
    r=s.get('report') or {};t=r.get('criteria_table') or {}
    if not t or ('literature' in t and all('assessment' in z for z in t.get('rows',[]))):return s
    overlay=dict(literature=reading_guide(r,t.get('identity_verified',False)),assessments={})
    for z in t.get('rows',[]):
        if z.get('evidence_view') and z.get('brief'):
            overlay['assessments'][z['code']]=assessment(z,t.get('identity_verified',False))
    return dict(s,view_overlay=overlay)
