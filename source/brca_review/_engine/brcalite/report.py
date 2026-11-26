import json,html
import re
from .net import writej
READING_STATES={
 'body_or_table_hits':'已取得正文，有段落或表行命中','fulltext_no_hit':'已取得正文，未命中别名',
 'fulltext_service_failure':'全文服务请求失败','no_accessible_fulltext_here':'本接口无可获取的开放全文',
 'no_record':'指定题录无记录','metadata_fields_missing':'题录缺核验字段',
 'metadata_service_or_fields_failed':'题录服务或字段核验失败','id_or_xml_problem':'标识或 XML 核验失败',
 'doi_mismatch_or_missing':'全文 DOI 核验失败','body_unavailable':'未取得有效正文',
 'not_processed_discovery_budget':'超过发现名额，未处理','not_processed_user_budget':'超过指定文献名额，未处理'}
def html_report(path,r):
 e=lambda s:html.escape('—' if s is None else str(s))
 def tbl(headers,rows):
  return '<table><tr>'+''.join('<th>'+e(s)+'</th>' for s in headers)+'</tr>'+''.join('<tr>'+''.join('<td>'+e(s)+'</td>' for s in row)+'</tr>' for row in rows)+'</table>'
 def candidate_row(item):
  source=item['source'];code=str(item['code'])
  match=re.match(r'[A-Z]+[0-9]+',code);base=match[0] if match else code
  review=r.get('criterion_dispositions',{}).get(base,{})
  # Existing reviewed reports also carry per-code met and accepted subtotals.
  accepted=base in r.get('met',[])
  contributes=review.get('contributes_to_subtotal',accepted and r.get('points_accepted') is not None)
  disposition=('已接受 '+str(review.get('reviewed_strength') or '') if accepted else
               '未满足' if base in r.get('not_met',[]) else '不适用' if base in r.get('reviewed_not_applicable',[]) else
               '审核阻断' if review else '待审核')
  contribution='是（按判据计一次）' if contributes else '独立判据' if accepted and base=='BA1' else '否'
  return (code+' / '+str(item.get('strength') or '见原表'),source['sheet']+' '+str(source['row']),
          str(source.get('review_cells') or source.get('cell'))+' '+', '.join(source.get('pmids_mentioned',[])),disposition,contribution,r.get('classification_draft') or '尚无草案')
 f=r['facts'];t=f.get('transcript',{});html0='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>BRCA 变异证据审核单</title><style>body{max-width:1000px;margin:32px auto;padding:0 20px;font:16px/1.6 system-ui,sans-serif;color:#202a38}h2{margin:28px 0 6px}table{border-collapse:collapse;width:100%;margin:8px 0}td,th{padding:8px 10px;border:1px solid #bbc4cd;text-align:left;overflow-wrap:anywhere}th{background:#eef3f5}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f2f5f6;padding:15px}.tag{color:#784c00}a{color:#154876}</style>'''
 h=[html0,'<h1>BRCA 胚系变异证据审核单</h1><p class="tag">本地研究草案 · 未签发 · 未提交 ClinGen</p>',
 tbl(['变异身份','本次对象'],[('输入',r['input']['original']),('本次注释组装',r['input']['assembly']),('GRCh37 g.',', '.join(r.get('representations',{}).get('GRCh37',[]))),('GRCh38 g.',', '.join(r.get('representations',{}).get('GRCh38',[]))),('基因／转录本',t.get('hgvsc')),('蛋白',t.get('hgvsp')),('CA',r.get('caid')),('疾病 / 遗传方式',str(r.get('condition'))+' / AD'),('本地规范',r.get('specification')),('本地判读状态',{'review_required':'待专家审核','blocked_identity_or_annotation':'身份或转录本未通过','draft_reviewed':'已复核的计算草案','review_blocked':'审核尚未完成或存在冲突'}.get(r['status'],r['status'])),('草案分级',r.get('classification_draft') or '尚无草案分类'),('已接受证据积分',r.get('points_accepted')),('签发状态','未签发')]),
 '<h2>外部既有断言</h2>',tbl(['来源','版本及提交方','原分类','原疾病'],[(x.get('vcv'),x.get('scv','')+' · '+str(x.get('submitter')),x.get('class'),x.get('condition')) for x in r['expert_assertions']]) if r['expert_assertions'] else '<p>此次未导入合格外部断言。</p>','<h2>待审核判据</h2>',
 tbl(['原表建议','工作表及物理行','回看位置／关联 PMID','当前判据处置','判据贡献小计','本次草案'],[candidate_row(x) for x in r['candidates']]),
 '<p>已满足：'+e(', '.join(r['met']) or '暂无审核结果')+'；经审核不满足：'+e(', '.join(r['not_met']) or '暂无审核结果')+'</p><p>尚未评估：'+e(', '.join(r['not_evaluated']))+'</p>']
 if r.get('literature'):
  lit=r['literature']
  tasks=lit.get('reading_tasks',[])
  if tasks:
   for origin,label in [('guideline','规范引用'),('user_seed','人工指定')]:
    selected=[t for t in tasks if any(s['origin']==origin for s in t['sources'])]
    scanned=sum(t['state'] in ('body_or_table_hits','fulltext_no_hit') for t in selected)
    h+=['<p>'+label+'已扫描全文 '+str(scanned)+' / '+str(len(selected))+'</p>']
   h+=[tbl(['PMID','来源及用途','获取状态','人工复核'],[(t['pmid'],' / '.join(s['role'] for s in t['sources']),READING_STATES.get(t['state'],t['state']),'待复核') for t in tasks])]
  for name,src in lit.get('coverage',{}).items():
   if src.get('status')!='ok':h+=['<p>'+e(name)+' · '+e(src.get('status'))+'</p>']
  h+=['<h2>本次全文检索</h2><p><a href="'+e(lit['reading_card'])+'">打开全文阅读卡</a> · '+e(lit['papers_scanned'])+' 篇已扫描 · '+e(lit['passages'])+' 处段落或表行命中</p>']
 else:
  h+=['<p>本次运行未检索全文。</p>']
 for cand in r['candidates']:
  row=cand['source']
  if 'review_values' in row:
   h+=['<h3>Table 9 第 '+e(row['row'])+' 行的实验说明</h3>',tbl(['单元格','原表内容'],[(k+str(row['row']),v) for k,v in row['review_values'].items()])]
 if f.get('workbook_rows_without_assignment'):
  h+=['<h2>原表未赋码或列为不适用</h2>',tbl(['工作表','物理行','原值'],[(row['sheet'],row['row'],row.get('code') or row.get('assigned')) for row in f['workbook_rows_without_assignment']])]
 if 'frequency_review' in f:
  freq=f['frequency_review']['result']
  h+=['<h2>人群频率审核</h2>',tbl(['状态','候选判据','强度','说明'],[(freq['status'],freq['code'],freq['strength'],freq['reason'])])]
 if r.get('caid'):
  from .core import GENES
  g=r.get('gene');ca=r.get('caid');ns=GENES[g][2] if g in GENES else ''
  h+=['<p>核对原站：<a href="https://cspec.genome.network/cspec/ui/svi/doc/'+e(ns)+'?version=1.2">CSpec v1.2</a> · <a href="https://reg.clinicalgenome.org/allele/'+e(ca)+'">CAR 等位</a> · <a href="https://ldh.clinicalgenome.org/ldh/Variant/id/'+e(ca)+'?detail=low">LDH 关联</a></p>']
 if r.get('review_lines'):
  h+=['<h2>审核记录</h2>',tbl(['代码','强度','状态','理由／来源'],[(x.get('code'),x.get('strength'),x.get('state'),str(x.get('reason',''))+' '+str(x.get('source_refs',[]))) for x in r['review_lines'] if isinstance(x,dict)])]
 if r.get('source_status'):
  h+=['<h2>本次来源状态</h2>',tbl(['来源','取数或核对结果'],r['source_status'].items())]
 if r.get('facts',{}).get('exon_workbook_context'):
  h+=['<h2>外显子工作表复核</h2>',tbl(['工作表','单元格','外显子','表内信息'],[(x['sheet'],x['cell'],x['exon'],x.get('strength') or x.get('event')) for x in r['facts']['exon_workbook_context']])]
 for hit in r.get('literature_examples',[]):
  h.extend(['<h2>全文定位</h2><p>',e(hit['section'])+' · <a href="https://pubmed.ncbi.nlm.nih.gov/'+e(hit['pmid'])+'/">PMID '+e(hit['pmid'])+'</a>；'+e(hit['summary_zh'])+' <a href="https://pmc.ncbi.nlm.nih.gov/articles/'+e(hit['pmcid'])+'/#'+e(hit['xml_id'])+'">复读原文</a></p>'])
 vrows=[]
 ens=f.get('ensembl_metadata',{}); ex=f.get('ldh_exome',{}); co=f.get('ldh_consequence',{}); bx=f.get('ldh_brca',{})
 if ens:vrows.append(('Ensembl 数据发布',', '.join(map(str,ens.get('releases',[]))),'指主机数据 release，不代表每个插件版本'))
 for src in (co.get('version_source') or []) if isinstance(co.get('version_source'),list) else []:
  vrows.append(('LDH 已抓取的分子后果',str(src.get('source'))+' '+str(src.get('ver')),'本等位实体的来源版本'))
 if ex:vrows.append(('gnomAD / LDH',str((ex.get('source') or {}).get('label'))+' exome','AC '+str(ex.get('AC'))+' / AN '+str(ex.get('AN'))+'；不能替代规范指定的双队列'))
 if bx:
  bs=bx.get('version_source') or {}
  if isinstance(bs,dict):vrows.append(('BRCA Exchange / LDH','releaseName '+str(bs.get('releaseName'))+'；ID '+str(bs.get('releaseId')),'BayesDel noAF '+str(bx.get('BayesDel_nsfp33a_noAF'))+'；预测不是功能实验证据'))
 if vrows:h.extend(['<h2>不同来源各是哪一版？</h2>',tbl(['来源','版本','本次怎样使用'],vrows)])
 if r['candidates'] and r.get('gene'):
  from .core import R
  from .net import getj
  urls={z['name']:z['url'] for z in getj(R/'rules/attachments.lock.json') if z['gene']==r.get('gene')}
  h+=['<p>回到官方表格查原行：']
  for name,uri in urls.items():
   if name.startswith('Specifications_Table') and any(x['source']['table']==('T9' if 'Table9' in name else 'T4_splice') for x in r['candidates']):
    h+=[' <a href="'+e(uri)+'">'+e(name.split('_')[1])+'</a> ']
  h+=['</p>']
  pm=sorted({i for x in r['candidates'] for i in x['source'].get('pmids_mentioned',[]) if str(i).isdigit()})
  if pm:h+=['<p>本行提到的文献（可能含排除项）：'+ ' · '.join('<a href="https://pubmed.ncbi.nlm.nih.gov/'+e(i)+'/">'+e(i)+'</a>' for i in pm)+'</p>']
 h.extend(['<details><summary>查看原始字段与请求记录</summary><pre>',e(json.dumps({'facts':f,'receipts':r['source_files']},ensure_ascii=False,indent=2)),'</pre></details><h2>人工任务</h2><ul>',*['<li>'+e(s)+'</li>' for s in r['limitations']],'</ul><p>完整结构：<a href="report.json">report.json</a></p></html>'])
 path.write_text(''.join(h),encoding='utf-8')

def reading_html(path,result):
 """Render locally licensed passages and identifiers; every hit remains a reading task."""
 esc=lambda v:html.escape(str(v or ''))
 def url(s,label):
  import urllib.parse
  u=urllib.parse.urlsplit(s or '')
  if u.scheme in ('http','https') and u.hostname in ('pubmed.ncbi.nlm.nih.gov','pmc.ncbi.nlm.nih.gov','www.ncbi.nlm.nih.gov','creativecommons.org','doi.org'):
   return '<a href="'+esc(s)+'">'+esc(label)+'</a>'
  return esc(label)
 a=['<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>BRCA 全文复核卡</title><style>body{max-width:900px;margin:36px auto;padding:0 18px;font:16px/1.65 system-ui,sans-serif;color:#22323e}article{margin:22px 0;border-top:2px solid #287369}blockquote{border-left:3px solid #82a79e;padding-left:1em;margin-left:0}small{overflow-wrap:anywhere;color:#485467}a{color:#14576d}details{margin:.7em 0}</style><h1>全文复核卡</h1><p>以下是原文段落或表行的定位。摘要和文末参考文献不算正文命中；复核研究设计与重复病例后才能判断证据方向。</p>']
 for t in result.get('reading_tasks',[]):
  a+=['<p>'+url('https://pubmed.ncbi.nlm.nih.gov/'+t['pmid']+'/','PMID '+t['pmid'])+' · '+esc(READING_STATES.get(t['state'],t['state']))+' · '+esc(' / '.join(s['role'] for s in t['sources']))+'</p>']
  for s in t['sources']:
   if s.get('sheet'):a+=['<p>'+esc(s['sheet'])+' · '+esc(s.get('cells'))+'</p>']
   if sup:=s.get('supplement'):
    a+=['<p>'+url('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117159','GEO GSE117159')+' · '+esc(sup['source'])+' · '+esc(sup['sheet'])+' 行 '+esc(sup['row'])+' · '+esc(sup['identity'])+' '+esc(sup['fields'])+'</p><small>定位时工作簿 SHA-256 '+esc(sup['sha256'])+' · 待人工复核原实验</small>']
 for f in result.get('fulltexts',[]):
  a+=['<article><h2>'+esc(f.get('title'))+'</h2><p>'+esc(f.get('first_author'))+' 等 · '+url(f.get('pubmed_url'),'PMID '+str(f.get('pmid')))+' · '+esc(f.get('pmcid'))+' · '+url('https://doi.org/'+str(f.get('doi')),'DOI '+str(f.get('doi')))+'</p><p>再利用条款：'+url(f.get('license_url'),f.get('license_url') or '未确认')+'</p>']
  if f.get('pmid')=='30040829' and f.get('pmcid')=='PMC6057642' and f.get('doi')=='10.1371/journal.pone.0201086' and 'c.181T>G' in result.get('aliases',[]):
   a+=['<p>复核 Kowalik 的 Results 和 Table 1 时，请核对同一转诊队列的 24 人，勿将段落与表行累计为两组。</p>']
  for h in sorted(f.get('hits',[]),key=lambda h:('results' not in str(h.get('section')).lower(),h.get('body_node_ordinal',0))):
   a+=['<h3>'+esc(h.get('section') or '未标章节')+' · '+esc(h.get('unit_label') or ('正文段' if h.get('type')=='p' else '表行'))+'</h3><p>'+url(h.get('review_url'),'跳到原文位置')+' · 命中：'+esc(' / '.join(h['aliases']))+'</p>']
   if h.get('review_flags') or h.get('match_type'):a+=['<p>'+esc(' · '.join([h.get('match_type',''),*h.get('review_flags',[])]))+'</p>']
   if h.get('table_headers'):a+=['<p>'+esc(h.get('table_caption'))+'</p><p>表头 '+esc(' / '.join(' | '.join(row) for row in h['table_headers']))+'</p>']
   if h.get('paragraph'):a+=['<blockquote>'+esc(h['paragraph'])+'</blockquote>']
   else:a+=['<p>仅列定位；请由有权访问的审阅者阅读全文。</p>']
   a+=['<small>正文节点 '+esc(h.get('body_node_ordinal'))+' · ID '+esc(h.get('xml_id'))+' · 段 SHA-256 '+esc(h.get('paragraph_sha256'))+'</small>']
  a+=['<p>补充材料：'+esc(len(f.get('supplement_links',[])))+' 个链接，尚未自动检索表内内容。</p></article>']
 if not result.get('fulltexts'):a+=['<p>此次未取得可核对的全文，详见人工队列。不能据此排除相关证据。</p>']
 a+=['<details><summary>检索范围与人工队列</summary><pre style="white-space:pre-wrap">'+esc(json.dumps({k:result.get(k) for k in ('queries','limit','coverage','review_queue')},ensure_ascii=False,indent=2))+'</pre></details><a href="literature.json">结构化记录</a></html>']
 path.write_text(''.join(a),encoding='utf-8')
