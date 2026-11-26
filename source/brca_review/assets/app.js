'use strict';
const boot=JSON.parse(document.getElementById('boot').textContent), $=s=>document.getElementById(s), docs=new Map(boot.snapshots.map(x=>[x.id,x])), menu=new Map(boot.index.map(x=>[x.id,x]));
const esc=x=>String(x??'—').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c])), jstr=x=>JSON.stringify(x,null,2);
const safeLink=(u,label)=> typeof u==='string'&&u.startsWith('https://')?'<a target="_blank" rel="noopener noreferrer" href="'+esc(u)+'">'+esc(label)+'</a>':esc(label);
const statusZh={matched_stop_exon_only:'未见逐等位行；已定位终止外显子规则行',allele_or_release_unverified:'等位或发布行未核准',retrieved:'已取得',variant_record:'已收录，身份一致',no_variant_record:'此数据集明确未见记录',no_link_listed:'目录已查，未列对应实体',no_matching_row:'已查工作簿，未见匹配行',matched_rows:'有匹配行',matched:'已核对',mismatch:'坐标或 REF 不符',no_matched_transcript:'未取得匹配的规范转录本',transport_error:'连接或读取失败',http_error:'HTTP 请求失败',graphql_error:'服务返回查询错误',unavailable:'未取得',needs_caid:'未取得精确 CAid，未查',needs_variation_id:'无可核对的 Variation ID，未查',no_verified_query_coordinate:'无已核准的该组装坐标，未查',identity_unresolved:'身份未核准',directory_unavailable:'上游实体目录未取得',no_unique_variation_id:'未取得唯一 Variation ID',allele_not_verified:'返回等位未能对账',no_records_for_query:'查询成功，未见匹配记录'};
const names=boot.labels, live=Boolean(boot.token)&&location.protocol==='http:', tabs=['evidence','papers','versions','raw','usage'];
let gene=boot.gene||'BRCA1',snap=null, tab='evidence', build=boot.assembly, rawKey='vep', busy=false, selected=boot.selected, seq=0;
const jsonView=x=>'<pre class="cell">'+esc(jstr(x))+'</pre>', chip=(s,w=false)=>'<span class="chip'+(w?' warn':'')+'">'+esc(s)+'</span>';
const table=(fields,rows)=>'<div class="scroll"><table><thead><tr>'+fields.map(f=>'<th>'+esc(f)+'</th>').join('')+'</tr></thead><tbody>'+rows.map(r=>'<tr>'+r.map(d=>'<td>'+d+'</td>').join('')+'</tr>').join('')+'</tbody></table></div>';
const gotoRaw=k=>'<button type="button" class="link" data-raw="'+esc(k)+'">'+esc(names[k]||k)+'</button>';
const rawLinks=ks=>ks.filter(k=>snap&&snap.raw[k]).map(gotoRaw).join(' · ');
const facts=(pairs)=>'<div class="grid">'+pairs.map(p=>'<div class="fact"><span>'+esc(p[0])+'</span><b>'+esc(p[1])+'</b></div>').join('')+'</div>';
const download=(bytes,fn)=>{const u=URL.createObjectURL(new Blob([bytes],{type:'application/octet-stream'}));const a=document.createElement('a');a.href=u;a.download=fn;a.click();setTimeout(()=>URL.revokeObjectURL(u),8000)};
const bytes=s=>Uint8Array.from(atob(s),c=>c.charCodeAt(0));
const textOf=b=>new TextDecoder().decode(bytes(b));
function flip(t){tab=t;tabs.forEach(k=>$(k+'Pane').hidden=k!==t);document.querySelectorAll('[data-tab]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.tab===t)));if(t==='raw')rawRender();if(t==='papers')paperTab()}
function byDate(b){return [...menu.values()].filter(d=>d.assembly===b&&(d.report?.selected_gene||d.gene)===gene).sort((a,b)=>(b.utc||'').localeCompare(a.utc||''))}
function lists(){document.querySelectorAll('[data-gene]').forEach(x=>{x.setAttribute('aria-pressed',String(x.dataset.gene===gene));x.disabled=busy||saveBusy});const chr=gene==='BRCA1'?'17':'13',tx=gene==='BRCA1'?'NM_007294.4':'NM_000059.4';$('vlabel').textContent=gene+' · chr'+chr+' · 一次一个 ALT';$('variant').placeholder='chr'+chr+'-POS-REF-ALT 或 '+tx+':c.…';
$('inputGuide').innerHTML='<b>优先 VCF 位点：</b>选原始组装，输入 chr'+chr+'-POS-REF-ALT，或粘贴一条 VCF 数据行（只取 CHROM、POS、REF、ALT）。也支持 <b>'+esc(tx)+':c.…</b>。每次一个 ALT；单独的 p. 不能用来确定核酸等位。';document.querySelectorAll('[data-build]').forEach(b=>{b.setAttribute('aria-pressed',String(b.dataset.build===build));b.disabled=busy||saveBusy});const h=$('history');h.innerHTML='';for(const d of byDate(build)){h.add(new Option((d.id.startsWith('demo-')?'样例':'本机')+' · '+(d.caid||'未定等位')+' · '+d.original+' · '+(d.utc?d.utc.replace('T',' ').slice(0,19)+' UTC':d.id.startsWith('demo-')?'公开教学':'未提供取数时间'),d.id))}if(selected)h.value=selected;h.disabled=busy||saveBusy||!h.options.length}
async function select(id){if(typeof saveBusy!=='undefined'&&saveBusy)return;const stamp=++seq;selected=id;lists();$('msg').textContent='';if(!id){snap=null;$('variant').value='';paint();return}if(!docs.has(id)){if(location.protocol==='file:'){snap=null;paint();return}try{const resp=await fetch('/api/snapshot/'+encodeURIComponent(id));if(!resp.ok)throw Error('快照不可读取');docs.set(id,await resp.json())}catch(e){if(stamp===seq){snap=null;paint();$('msg').textContent=String(e)}return}}if(stamp!==seq)return;paperMessage='';snap=docs.get(id);$('variant').value=snap.original;rawKey='current_report';paint()}
function paint(){const r=snap?.report;const ident=$('identity');$('mode').textContent=live?'查询 · 手动联网':'快照 · 只读';$('query').disabled=!live||busy;$('qhint').textContent=live?'只在点击时发送该公开可查询的核酸位点；先查原组装；CAR 有唯一对应时可另查另一组装，分别标明来源。':'本文件不会联网更新；需要查新位点或复查，请看「使用与安装」。';
 if(!snap){ident.innerHTML='<p>此文件未附该组装的快照。可以选回已有快照；在本机模式下也可用原生坐标新查。</p>';$('evidencePane').innerHTML='<h2>暂无该组装的证据</h2>'; $('versionsPane').innerHTML='<h2>尚无可显示的版本</h2>';rawRender();paperTab();return}
 const f=r.facts||{},fallback=Object.values(f.gnomad||{}).flatMap(x=>x.transcripts||[]).filter(t=>t.hgvsp&&t.gene_symbol===gene),t=f.transcript||{},rep=r.representations||{}, ok=!!r.criteria_table?.identity_verified;
 ident.innerHTML='<div class="row between"><div><div class="eyebrow">当前等位　/　'+esc(snap.assembly)+'</div><h2>'+esc(t.hgvsc||r.input.original)+'</h2></div>'+chip(ok?'等位与注释已对账':'身份或注释待核',!ok)+'</div><p>'+esc(r.gene)+'　·　'+esc(t.hgvsp||(!t.hgvsc&&fallback.length?null:'未返回蛋白变化'))+'　·　'+esc(r.caid)+'　·　'+esc(r.condition)+'　/　'+esc(r.inheritance)+'</p><details><summary>核对两套组装的表示与原输入</summary>'+table(['项目','原值'],[['原始输入',esc(r.input.original)],['GRCh37',esc((rep.GRCh37||[]).join(' / ')||'未确认')],['GRCh38',esc((rep.GRCh38||[]).join(' / ')||'未确认')],['本次注释主机',esc(snap.raw.info?.url||'未取到')]])+'</details>'+decisionHeader(r)+'<small>生成此快照的引擎：'+esc(r.review_desk_version||r.engine||'未披露')+' · '+(snap.utc?'快照：'+esc(snap.utc)+' UTC':'公开教学示例，不含采集时间')+'　·　'+esc(snap.id)+'　·　打开时间不替代来源版本。</small>';
 if(r.acquisition==='live_independent'){
  const tx=f.transcript_sources||{},k=tx.hgvsp||tx.hgvsc;
  if(k)ident.innerHTML+='<p class="smallnote">本页规范转录本／蛋白：'+esc(names[k]||k)+(f.source_transcripts?.car?.hgvsp&&f.source_transcripts.car.hgvsp===t.hgvsp&&k!=='car'?'；CAR 蛋白描述一致':'')+'。'+(snap.raw[k]?gotoRaw(k):'')+'</p>';
  else if(fallback.length){const unique=[...new Set(fallback.map(v=>[v.transcript_id,v.hgvsc,v.hgvsp].filter(Boolean).join(' · ')))];ident.innerHTML+='<details open><summary>gnomAD 返回的转录本注释（未核准规范 RefSeq）</summary><p>'+unique.slice(0,8).map(esc).join('<br>')+'</p></details>';}
  if(r.source_status?.ref_check==='mismatch')ident.innerHTML+='<p class="note">输入 REF 或参考坐标不匹配；查询结果保留供对照，不计分。</p>';
 }
 if(!ok&&r.source_status?.ref_check==='unavailable'){const reason=snap.raw?.ref?.status||snap.raw?.ref?.http_status;ident.innerHTML+='<p class="note"><b>参考碱基未核实，暂不计分。</b> '+esc(snap.assembly)+' 参考服务本次无可用响应'+(reason?'（'+esc(reason)+'）':'')+'。保留原输入，可稍后重查；要换组装，请先取得同一等位在目标组装的坐标。</p>';}
 evidence();versions();paperTab();if(tab==='raw')rawRender()}
/* REVIEW_LIB */
function shownFields(v,r){if(v.key!=='erepo'||!v.fields)return v.fields;const rows=r.criteria_table?.rows||[];return {...v.fields,codes:rows.filter(z=>z.expert&&(!z.expert.strength||z.allowed_strengths.includes(z.expert.strength))).map(z=>[z.code,z.expert.state,z.expert.strength])};}
function versions(){const r=snap.report;let h='<div class="eyebrow">02　/　PROVENANCE</div><h2>这个位点用了各库的哪一版？</h2><p class="muted tight">下面是快照采集时各接口实际提供的版本。用户打开页面不会自动覆盖它。规则版本、主机 release、数据集版和断言版分别保留。</p><p class="note tight">联网复查指今天能从所选接口读到的内容。gnomAD 官网队列、LDH 实体和专家发布记录分别列出；同一来源的不同子集不合并，未披露的版本不补写。</p>';
 const cur=snap.versions.find(v=>v.key==='cspec_current')?.version,locked=(r.specification||'').split(' v')[1];
 if(!cur||cur!==locked)h+='<p class="note callout">'+esc(!cur?'官方现行规范版本未核实。':'官方规范与本地冻结规则不同；请人工审核后再更新规则。')+'</p>';
 h+=table(['来源 / 适用范围','冻结于此快照的版本','解释与原始数据'],snap.versions.filter(v=>v.key!=='erepo'||r.erepo?.state==='matched').map(v=>['<b>'+esc(v.label)+'</b><br><small>'+esc(v.scope)+'</small>',esc(v.version??'未披露或未取得')+(v.source_date?'<br><small>来源内日期：'+esc(v.source_date)+'</small>':'')+(v.utc?'<br><small>取数 UTC：'+esc(v.utc.substring(0,19).replace('T',' '))+'</small>':''),esc(v.meaning)+'<br>'+rawLinks(v.raw)+(v.url?' '+safeLink(v.url,'查看原站'):'')+(v.fields!=null?'<details><summary>保留字段</summary>'+jsonView(shownFields(v,r))+'</details>':'')]));
 const ex=r.facts.ldh_exome||{};if(ex.AN!=null)h+='<div class="note tight">本等位的 gnomAD exome：AC '+esc(ex.AC)+' / AN '+esc(ex.AN)+'；需按原规范核对人群、队列和覆盖。</div>';
 h+='<details><summary>本次各接口的取数状态（含未取得或未列出）</summary>'+table(['接口','状态'],Object.entries(r.source_status||{}).map(([k,v])=>[esc(names[k]||k),esc(statusZh[v]||v)]))+'</details><h3 class="sp">与上次选定快照的比较</h3>';
 const d=snap.diff;h+='<p>'+esc(d?.message||(live?'从本页重查同一等位后在此显示差异。':'目前是冻结快照；下载本机包后可重查并比较。'))+'</p>';
 if(d?.rows?.length)h+=table(['来源','字段比较','响应字节 / 变化'],d.rows.map(x=>[esc(x.label),esc(x.state),esc(x.raw)+(Object.keys(x.changes).length?'<details><summary>前值 / 后值</summary>'+jsonView(x.changes)+'</details>':'')]));
 $('versionsPane').innerHTML=h}
function rawRender(){const el=$('rawPane');if(!snap){el.innerHTML='<h2>暂无原始数据</h2>';return}
 const opts=['current_report',...Object.keys(snap.raw),'report'];if(!opts.includes(rawKey))rawKey=opts[0];el.innerHTML='<div class="eyebrow">03　/　RAW DATA</div><h2>原始响应与请求凭据</h2><p class="muted tight">原始字节可下载。美化 JSON 仅改变阅读排版；报告 JSON 是本地派生文件。</p><div class="row bar"><label for="rawsel">来源</label><select id="rawsel"></select><button id="rawDown">下载原始字节</button><button id="pretty" aria-pressed="false">美化 JSON</button></div><div id="rawmeta" class="sp"></div><pre id="rawbody" class="cell rawpre"></pre>';
 opts.forEach(k=>$('rawsel').add(new Option(k==='current_report'?'本次审核报告 JSON（与页面一致）':k==='report'?(snap.teaching?'公开教学基线 JSON（已去采集时间）':'取数时原报告 JSON（冻结）'):names[k]||k,k)));$('rawsel').value=rawKey;$('rawsel').onchange=e=>{rawKey=e.target.value;rawRender()};
 if(rawKey==='current_report'){
  const plain=jstr(snap.report)+'\n';$('rawbody').textContent=plain;$('pretty').disabled=true;$('rawDown').textContent='下载本次审核 JSON';$('rawDown').disabled=false;
  $('rawmeta').textContent='本页与导出读取同一份 criteria_table；原始查询字节另见“取数时原报告”。';
  $('rawDown').onclick=()=>download(new TextEncoder().encode(plain),snap.id+'-'+(snap.local_review?.revision||'view')+'-report.json');return
 }
 const m=rawKey==='report'?{sha256:snap.report_sha256,ok:true,note:snap.teaching?'公开教学摘录；此 SHA-256 对应发布包内的示例字节。':'本地生成的结构化报告；不是外部响应。'}:snap.raw[rawKey],b=rawKey==='report'?snap.raw_report:snap.blobs[m.sha256];
 const plain=b?textOf(b):'',mtype=m.content_type||'';$('rawDown').disabled=!b;let obj=null;try{obj=JSON.parse(plain)}catch(e){}$('pretty').disabled=obj===null;$('rawbody').textContent=b?plain:(m.note||'未取得响应字节；请查看取数状态。');
 $('rawmeta').innerHTML='<p class="break tight">'+(m.method?esc(m.method)+' · '+safeLink(m.method==='GET'?m.url:null,m.url):'本地文件')+'</p><p class="tight">'+(m.method?'HTTP '+esc(m.status)+' · '+esc(m.ok?'成功':'失败')+' · '+esc(m.utc)+' UTC · ':'')+esc(m.size??(b?bytes(b).length:'—'))+' bytes</p><p class="break tight"><b>SHA-256</b> <code>'+esc(m.sha256)+'</code></p>'+(m.request?'<details><summary>原始请求体及哈希</summary><pre class="cell">'+esc(m.request)+'</pre><code class="break">'+esc(m.request_sha256)+'</code></details>':'')+(m.note?'<p>'+esc(m.note)+'</p>':'')+(m.error?'<p>'+esc(m.error)+'</p>':'');
 $('pretty').onclick=e=>{const on=e.target.getAttribute('aria-pressed')==='true';e.target.setAttribute('aria-pressed',String(!on));$('rawbody').textContent=on?plain:jstr(obj)};
 $('rawDown').onclick=()=>{const ext=rawKey==='report'||obj!==null?'json':mtype.includes('xml')?'xml':'txt';download(bytes(b),snap.id+'-'+rawKey+'.'+ext)}}
document.querySelectorAll('[data-tab]').forEach(b=>b.onclick=()=>flip(b.dataset.tab));document.querySelectorAll('[data-build]').forEach(b=>b.onclick=()=>{if(busy||saveBusy)return;build=b.dataset.build;select(byDate(build)[0]?.id||null)});
document.addEventListener('click',e=>{const a=e.target.closest('[data-raw]');if(a){rawKey=a.dataset.raw;flip('raw');$('rawPane').scrollIntoView({block:'start'})}});
$('history').onchange=e=>select(e.target.value);
$('query').onclick=async()=>{if(!live||busy||saveBusy)return;let variant=$('variant').value.trim();if(/[\r\n]/.test(variant)){$('msg').textContent='请一次粘贴一行、一项 ALT。';return}if(/[ \t]/.test(variant)){const c=variant.split(/[ \t]+/);if(c.length<5||!/^(chr)?(13|17)$/.test(c[0])||!/^\d+$/.test(c[1])||!c.slice(3,5).every(x=>/^[ACGT]+$/.test(x))){$('msg').textContent='VCF 行需 CHROM POS ID REF ALT；多 ALT 请先拆分。';return}variant=[c[0],c[1],c[3],c[4]].join('-');$('variant').value=variant}if(!variant){$('msg').textContent='请填写原始组装的坐标或带版本 c.。';return}busy=true;lists();$('query').disabled=true;$('variant').disabled=true;$('msg').textContent='正在查询公开服务并留存原文，通常需要数十秒；可继续查看下方原快照。';try{const res=await fetch('/api/query',{method:'POST',headers:{'Content-Type':'application/json','X-Review-Token':boot.token},body:JSON.stringify({gene,assembly:build,variant,previous:snap?.id||null})});const s=await res.json();if(!res.ok)throw Error(s.error||'请求未完成');docs.set(s.id,s);menu.set(s.id,s);await select(s.id);$('msg').textContent=s.status==='review_required'?'已另存新快照；「各库版本与差异」列有本次取数 UTC。选有同等位、同组装基线时可见差异。':'已保存本次取数，身份或注释未通过；不可据此判读。'}catch(e){$('msg').textContent='查询未完成：'+e.message}finally{busy=false;lists();$('variant').disabled=false;$('query').disabled=!live}};
// Each paper action is user-initiated; opening an exported HTML never issues network requests.
const paperStates={body_or_table_hits:'已取得正文，有段落或表行命中',fulltext_no_hit:'已取得正文，所列别名未命中',fulltext_service_failure:'全文请求失败，可转人工',metadata_service_or_fields_failed:'题录服务失败或返回字段不足',no_accessible_fulltext_here:'此接口未提供可获取的开放全文',metadata_fields_missing:'题录字段不足，不能可靠取全文',no_record:'本次精确题录查询无记录',id_or_xml_problem:'文章身份或 XML 校验未通过',doi_mismatch_or_missing:'DOI 不一致或缺失，未采纳',body_unavailable:'文件中未找到可读正文',failed_live_with_cache:'本次联网失败；下列内容来自已核验旧缓存',verified_cache:'已核验的历史缓存；本次未联网',pending:'尚未取全文'};
let paperBusy=false,paperMessage='';
function lit(){return snap?.report?.criteria_table?.literature||snap?.view_overlay?.literature||{}}
function groupPapers(){
 if(!snap)return {};const list={};
 const orig=lit().records||snap.reading_assignments||{};
 for(const [p,roles] of Object.entries(orig))if(/^[1-9]\d{0,9}$/.test(p)){
   const supp=roles.find(o=>o.supplement)?.supplement;
   list[p]={pmid:p,state:'pending',origins:roles,supplement:supp,links:{PubMed:'https://pubmed.ncbi.nlm.nih.gov/'+p+'/','Europe PMC 题录':'https://europepmc.org/article/MED/'+p,...(supp?.source?.startsWith('GEO GSE117159,')?{'GEO 补充资料':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117159'}:{})}};
 }
 const demo=boot.paper_demo||{},sample=demo.record;
 if(snap.teaching===true&&sample&&snap.caid===demo.caid&&snap.report?.gene===demo.gene&&snap.report?.facts?.transcript?.hgvsc===demo.hgvsc&&!Object.hasOwn(snap.reading_events||{},sample.pmid)){
   const p=sample.pmid,old=list[p]||{};list[p]={...old,...sample,links:{...old.links,...sample.links},origins:[...(old.origins||[]),...sample.origins]};
 }
 for(const [p,event] of Object.entries(snap.reading_events||{})){const original=list[p]||{};const all=[...(original.origins||[]),...(event.origins||[])];const seen=new Set();list[p]={...original,...event,links:{PubMed:'https://pubmed.ncbi.nlm.nih.gov/'+encodeURIComponent(p)+'/','Europe PMC 题录':'https://europepmc.org/article/MED/'+encodeURIComponent(p),...original.links,...event.links},origins:all.filter(o=>{const id=[o.origin,o.role,o.criterion,o.row,o.locator].join('/');if(seen.has(id))return false;seen.add(id);return true})};}
 return list;
}
function paperPriority(e){
 const src=e.origins||[],excl=src.some(v=>v.use==='explicitly_excluded'||/在本行排除/.test(v.role||'')),used=src.some(v=>v.use==='explicitly_included'||/表内采用的实验/.test(v.role||''));
 if(excl)return used?{rank:1,label:'核查来源',why:'有来源采用，也有来源排除；先核对各自针对的实验和判据。'}:{rank:3,label:'核对排除项',why:'有来源明确排除这项实验或研究；先核排除范围与理由。'};
 if(e.supplement?.sheet&&e.supplement?.row)return {rank:0,label:'优先读',why:'已有作者补表的工作表和行号，可直接核对。'};
 const c=(snap?.report?.facts?.transcript?.hgvsc||'').split(':')[1]||'',hits=e.article?.hits||[];
 if(c&&hits.some(v=>v.aliases?.includes(c)&&(v.type==='tr'||/\bresults?\b|结果/i.test(v.section||''))&&!(v.review_flags||[]).some(x=>x.includes('另一 BRCA'))))return {rank:0,label:'优先读',why:'在表行或结果段找到当前 c. 名称；仍须确认基因、ALT 和观察对象。'};
 if(used||hits.length||src.some(v=>v.origin==='guideline'||v.origin==='erepo'))return {rank:1,label:'核查来源',why:'有指南、专家或正文线索，先确认它是否提供本等位的直接观察。'};
 return {rank:2,label:'阅读线索',why:'先核原文、补表和核酸等位，再决定是否有可用数据。'};
}
function orderedPapers(){
 const rank=x=>Math.min(6,...(x.origins||[]).map(o=>o.origin==='guideline'?0:o.origin==='erepo'&&o.use!=='explicitly_excluded'?1:o.origin==='erepo'?2:o.origin==='clinvar_scv'?4:3));
 return Object.values(groupPapers()).sort((a,b)=>paperPriority(a).rank-paperPriority(b).rank||rank(a)-rank(b)||String(a.pmid).localeCompare(String(b.pmid)));
}
function readingStrip(compact){
 const g=lit(),all=orderedPapers(),n=all.length,limit=compact?5:0;
 let h='<section class="readstrip"><h3>文献与人工复核</h3>';
 if(compact){h+='<p>'+(n?'本快照列出 '+n+' 个 PMID：'+all.slice(0,limit).map(v=>safeLink('https://pubmed.ncbi.nlm.nih.gov/'+v.pmid+'/','PMID '+v.pmid)).join('、')+(n>limit?' 等':''):'本快照尚无关联 PMID，可按此位点人工查找。')+' <button type="button" class="link" data-read-guide>查看文献与检索说明</button></p>';}
 h+='<p class="readinglinks">'+Object.entries(g.links||{'PubMed':'https://pubmed.ncbi.nlm.nih.gov/','PMC':'https://pmc.ncbi.nlm.nih.gov/'}).map(([name,url])=>safeLink(url,name)).join(' · ')+'</p>';
 if(!compact)h+='<p>检索别名：'+esc((g.aliases||[]).join('、')||'未核准')+'。PubMed 查题名、摘要等题录字段；位点只在正文或补表时，请用 PMC 或期刊原文。功能／RNA 筛选只是可选的窄检索。</p><p>'+esc(g.note||'请核对同一核酸 ALT、实验对照及独立队列。')+' 即便蛋白变化相同，也要回到核酸核对。</p>';
 return h+'</section>';
}
document.addEventListener('click',e=>{if(e.target.closest('[data-read-guide]')){flip('papers');$('papersPane').scrollIntoView({block:'start'})}});
// Read-only navigation derived from this paper event. Never infers evidence strength.
function readWhere(e){
 const orig=e.origins||[],s=e.supplement,a=e.article||{},cdna=(snap?.report?.facts?.transcript?.hgvsc||'').split(':')[1]||'';
 const dna=v=>Boolean(cdna&&v.aliases?.some(x=>x===cdna));
 const results=v=>/\bresults?\b|结果/i.test(v.section||'');
 const order=v=>(dna(v)?0:4)+(v.type==='tr'?0:results(v)?1:2);
 const hits=(a.hits||[]).slice().sort((x,y)=>order(x)-order(y));
 const leads=[],contexts=[];
 if(orig.some(v=>v.use==='explicitly_excluded'||/在本行排除/.test(v.role||'')))
   contexts.push('有来源在对应判据中明确排除这项研究或实验。先查排除原因及更正；不能凭文献命中将其重新计入。');
 if(s?.sheet&&s?.row)leads.push({label:s.sheet+' 第 '+s.row+' 行（作者补表）',url:s.url,reason:'已有该等位的补表定位。先核转录本、REF/ALT、列标题与实验对照。'});
 if(hits.length){
   const first=hits[0];const second=hits.find(v=>v!==first&&results(v)&&v.type!=='tr'&&!(v.review_url&&v.review_url===first.review_url));
   [first,...(second?[second]:[])].forEach(v=>{
     const where=[v.unit_label,v.section].filter(Boolean).join(' · ')||(v.type==='tr'?'正文表行':'正文段落');
     const why=dna(v)?'命中当前 c. 名称；仍需核对基因、分母与等位。':'只命中其他别名；先确认是否对应当前核酸 ALT。';
     leads.push({label:where,url:v.review_url,reason:(v.type==='tr'?'对照表头看该行；':results(v)?'查看作者如何描述实验或病例；':'核查这处是否只是背景引用；')+why});
   });
 }else if(a.state==='fulltext_no_hit'||e.state==='fulltext_no_hit'){
   contexts.push('已核验正文 XML，这组别名没有命中。请查作者补表、旧命名或勘误；正文未命中不等于文献没有数据。');
 }else if(!e.article&&(e.state==='pending'||!e.state)){
   contexts.push('尚未取得正文定位。先看推荐此文的来源；在本机可点「取正文（优先缓存）」。');
 }else if(!e.article){
   contexts.push('当前没有可定位的正文。请从本卡的 PubMed、PMC 或期刊链接查看原文与补表。');
 }
 const codes=orig.map(v=>String(v.criterion||''));
 if(codes.some(c=>/^(PS3|BS3)\b/.test(c)))contexts.push('复核功能证据：变异的具体实验、阳性／阴性对照和校准；RNA 与蛋白实验分开看。');
 if(codes.some(c=>/^(PP1|BS4|PP4|BP5|PS4)\b/.test(c)))contexts.push('复核临床数据：LR 类型、原始尺度、家系或队列来源；组成项不能和已包含它的综合值重复使用。');
 const priority=paperPriority(e);
 return {hits,html:'<aside class="note callout" data-reading-priority="'+esc(priority.rank)+'"><strong>先读哪里 · '+esc(priority.label)+'</strong> <small>阅读顺序，不计分</small><p class="tight">'+esc(priority.why)+'</p>'+
 (leads.length?'<ol>'+leads.map(q=>'<li>'+safeLink(q.url,q.label)+'：'+esc(q.reason)+'</li>').join('')+'</ol>':'')+
 contexts.map(q=>'<p class="tight">'+esc(q)+'</p>').join('')+'</aside>'};
}
function paperTab(){
 const pane=$('papersPane');if(!pane)return;
 if(!snap){pane.innerHTML='<h2>先选择有等位核对结果的快照</h2>';return}
 const entries=orderedPapers(),active=live&&['review_required','draft_reviewed'].includes(snap.status);
 let h='<div class="eyebrow">READING　/　PUBLICATIONS</div><h2>全文与人工阅读</h2><p>先看指南指定的文章，再按需取一篇全文。题录、可定位的正文、补充表是不同的材料；未取得正文不赋分，也不阻断变异查询。</p>';
 h+=readingStrip(false);
 if(active)h+='<div class="controls"><label for="paperInput">指定 PMID　</label><input id="paperInput" inputmode="numeric" maxlength="10" placeholder="如 30040829"> <button type="button" data-paper-manual="yes" '+(paperBusy?'disabled':'')+'>查这篇（优先缓存）</button></div>';
 else h+='<p class="muted tight">这是只读页面。下列官方链接可手工打开；如要重新查询，请启动本机 Python。</p>';
 h+='<p class="status" aria-live="polite" id="paperNote">'+esc(paperMessage)+'</p>';
 if(!entries.length)h+='<p>未在本快照来源中提取到关联 PMID。可用上方预填检索式；找到文章后可输入 PMID 查看题录、正文状态与定位。</p>';
 for(const e of entries){const p=String(e.pmid),a=e.article,hrefs=e.links||{PubMed:'https://pubmed.ncbi.nlm.nih.gov/'+encodeURIComponent(p)+'/'};
 h+='<article class="detail"><div class="row between"><h3>PMID '+esc(p)+'</h3>'+chip(e.state==='teaching_locators'?'公开教学定位 · 未调用本机缓存':(paperStates[e.state]||e.state),!['body_or_table_hits','verified_cache','fulltext_no_hit','teaching_locators'].includes(e.state))+'</div>'+(e.title&&!a?'<p><b>'+esc(e.title)+'</b></p>':'')+'<p>'+esc([...new Set((e.origins||[]).map(x=>x.role))].join('；'))+'</p><div class="row tight">'+Object.entries(hrefs).map(([k,v])=>safeLink(v,k)).join('　')+'</div>';
 const guide=readWhere(e);h+=guide.html;
 if(active)h+='<p class="row sp"><button type="button" data-paper="'+esc(p)+'" '+(paperBusy?'disabled':'')+'>取正文（优先缓存）</button><button type="button" data-paper="'+esc(p)+'" data-refresh="yes" '+(paperBusy?'disabled':'')+'>重新联网</button></p>';
 if(e.supplement){const z=e.supplement;h+='<p class="note tight">独立补表定位：'+esc(z.source)+' · '+esc(z.sheet)+' 第 '+esc(z.row)+' 行。'+esc(z.identity)+'。请人工核对；整篇正文未取得也可查这条线索。</p>'}
 if(a){h+='<p><b>'+esc(a.title)+'</b> · '+esc(a.first_author)+' 等</p><p class="tight">PMCID '+esc(a.pmcid)+'；'+esc((a.hits||[]).length)+' 处定位；全文来源 UTC：'+esc(e.source_utc||(e.state==='teaching_locators'?'公开教学不附采集时间':'未提供'))+'。</p>';
 const hits=guide.hits;
 if(p==='30040829'&&snap.caid==='CA001182'&&hits.some(x=>x.unit_label==='Table 1')&&hits.some(x=>x.type==='p'&&/Results/i.test(x.section||'')))h+='<p class="note callout">Results 和 Table 1 的 24 人来自同一转诊队列。两处定位只供交叉核对。</p>';
 for(const v of hits)h+='<details><summary>'+esc(v.section||'正文章节未标明')+' · '+esc(v.unit_label||'正文段')+' · '+esc((v.aliases||[]).join(' / '))+'</summary><p>'+safeLink(v.review_url,'跳到原文位置')+'</p>'+(v.table_headers?.length?'<p class="tight">表头：'+esc(v.table_headers.map(a=>a.join(' | ')).join('; '))+'</p>':'')+(v.paragraph?'<blockquote>'+esc(v.paragraph)+'</blockquote>':'<p>仅列定位；请在原站阅读。</p>')+((v.review_flags||[]).length?'<p class="note tight">复读提示：'+esc(v.review_flags.join('；'))+'</p>':'')+'<p class="tight">段落或表行 SHA-256：<code class="break">'+esc(v.paragraph_sha256)+'</code></p></details>';
 h+='<p class="tight break">XML SHA-256：<code>'+esc(a.xml_sha256)+'</code>　'+safeLink((a.license_url||'').replace(/^http:/,'https:'),'原文许可')+'</p>'}
 if(e.cache_warning)h+='<p class="note callout">'+esc(e.cache_warning)+'</p>';
 if(e.receipts)h+='<details><summary>本次请求、尝试次数与错误记录</summary>'+jsonView(e.receipts)+'</details>';
 if(e.id){h+='<p class="tight">阅读记录 '+esc(e.id)+' · '+esc(e.finished_utc)+'；本地目录 reading_events/'+esc(e.id)+'/report.html 可离线复读。</p>';if(live&&e.xml_event)h+='<a download href="/api/paper-xml/'+encodeURIComponent(e.id)+'">下载已核准的 JATS 原始字节</a>'}
 h+='</article>'}
 pane.innerHTML=h;
}
async function askPaper(pmid,refresh){
 if(!live||paperBusy||!snap)return;
 if(!/^[1-9][0-9]{0,9}$/.test(pmid)){paperMessage='请输入 1–10 位正整数 PMID。';paperTab();return}
 const s=snap;paperBusy=true;paperMessage='正在核对 PMID '+pmid+'。可继续查看判据与版本。';paperTab();
 try{const res=await fetch('/api/paper',{method:'POST',headers:{'Content-Type':'application/json','X-Review-Token':boot.token},body:JSON.stringify({snapshot:s.id,pmid,refresh})});const z=await res.json();if(!res.ok)throw Error(z.error||'获取失败');s.reading_events??={};s.reading_events[pmid]=z;paperMessage=snap?.id===s.id?'已另存阅读记录；本地判据和积分未改变。':''}
 catch(e){s.reading_events??={};if(!s.reading_events[pmid])s.reading_events[pmid]={pmid,state:'fulltext_service_failure',links:{PubMed:'https://pubmed.ncbi.nlm.nih.gov/'+pmid+'/'},origins:[{role:'请人工核对；本次页面请求未完成留档。'}]};paperMessage=snap?.id===s.id?('查询未完成：'+e.message+'。可打开该 PMID 的原站链接。'):''}
 finally{paperBusy=false;paperTab()}
}
document.addEventListener('click',e=>{const b=e.target.closest('[data-paper],[data-paper-manual]');if(b){const p=b.dataset.paper||$('paperInput')?.value.trim()||'';askPaper(p,b.dataset.refresh==='yes')}});

document.querySelectorAll('[data-gene]').forEach(b=>b.onclick=()=>{if(busy||saveBusy)return;gene=b.dataset.gene;select(byDate(build)[0]?.id||null)});
const requestedGene=new URLSearchParams(location.search).get('gene');if(['BRCA1','BRCA2'].includes(requestedGene)){gene=requestedGene;selected=byDate(build)[0]?.id||null}
const requestedBuild=new URLSearchParams(location.search).get('assembly');if(['GRCh37','GRCh38'].includes(requestedBuild)){build=requestedBuild;selected=byDate(build)[0]?.id||null}
const wantedSnapshot=new URLSearchParams(location.search).get('snapshot');if(wantedSnapshot&&byDate(build).some(s=>s.id===wantedSnapshot))selected=wantedSnapshot;
function followCriterionLink(){
 const key=location.hash.slice(1);
 if(!/^criterion-[A-Z]{2,3}[1-9][0-9]?$/.test(key))return false;
 const row=document.getElementById(key);
 if(!row||!row.hasAttribute('data-code'))return false;
 flip('evidence');
 for(let d=row.closest('details');d;d=d.parentElement?.closest('details'))d.open=true;
 row.focus({preventScroll:true});row.scrollIntoView({block:'center'});return true;
}
window.addEventListener('hashchange',followCriterionLink);
select(selected||byDate(build)[0]?.id||null).then(()=>{if(!followCriterionLink()&&new URLSearchParams(location.search).get('tab')==='papers')flip('papers')});
