const unsaved=new Map();let saveBusy=false;
const means={workbook:'附件明确推荐',prediction:'本地阈值已满足',population:'规范队列已核准',expert_citation:'引用专家',manual:'手工补入'};
const nsign=n=>n==null?'—':(n>0?'+':'')+String(n);
const summaryOf=obj=>obj?.stand_alone?'BA1，单独判据':nsign(obj?.points);
function decisionHeader(r){
 const e=r.erepo||{},L=r.criteria_table||{},out=r.review_result||{},rev=snap.local_review,codes=(L.rows||[]).filter(z=>z.expert?.state==='Met'&&z.allowed_strengths.includes(z.expert.strength)).map(z=>z.code+' '+z.expert.strength);
 const matched=e.state==='matched'&&Boolean(e.classification);
 const extra=matched?'<details class="expert-extra" id="expertExtra"><summary>eRepo 专家补充：<b>'+esc(e.classification)+'</b></summary><p>'+esc(codes.join('；')||'逐码未核全')+'</p><p>'+esc(e.review_panel||'ENIGMA')+' · 规范 '+esc(e.spec_version||'未披露')+' · '+safeLink(e.url,'打开此专家个案')+'。在下方明确勾选的专家引用才参与本地计算。</p></details>':'';
 let right='<b>尚未保存</b><div>已有确切来源的项，会在下方预选。</div>';
 if(!L.reviewable&&!rev&&!L.identity_verified)right='<b>身份未核准，不能计分</b><div>可勾选并保存待核清单，计分前需核准等位。</div>';
 if(rev?.schema==='brca.selection.v2'){
  right='<b>'+(!rev.decisions.length?'未采纳判据':rev.decisions.every(d=>d.state==='pending')?'已保存待核，尚未计分':out.issues?.length?'组合待核':summaryOf(out)+(out.stand_alone?'':' 分'))+'</b><div>本地来源 '+esc(summaryOf(out.local_sources))+'；专家引用 '+esc(summaryOf(out.quoted_sources))+'；手工补入 '+esc(summaryOf(out.manual_sources))+'</div>';
  if(rev.decisions.length){const by=new Map((L.rows||[]).map(z=>[z.code,z]));right+='<ul class="savedlines">'+rev.decisions.map(d=>'<li><strong>'+esc(d.code)+(d.strength?' '+esc(d.strength):'')+'</strong> · '+esc(means[d.origin]||'历史审核')+' · '+(d.state==='pending'?'待核，不计分':d.strength==='Stand Alone'?'独立良性':esc(nsign(d.contribution_points??by.get(d.code)?.contribution_points)))+'</li>').join('')+'</ul>';}
  if(out.selected_combination)right+='<div>所选组合：'+esc(out.selected_combination)+'（研究提示）</div>';
 }else if(rev||r.legacy_read_only||r.points_accepted!=null)right='<b>历史审核：'+esc(r.classification_draft||'未形成草案')+'</b><div>原积分 '+esc(nsign(r.points_accepted))+'；保留原记录。</div>';
 return '<div class="grid brief-head"><div class="fact"><span>已保存的所选证据'+(snap.teaching?' · 教学':'')+'</span>'+right+'</div></div>'+extra;
}
function buf(){const key=snap.id+'/'+(snap.local_review?.revision||'base');
 if(!unsaved.has(key)){const rev=snap.local_review,rows=snap.report.criteria_table?.rows||[];unsaved.set(key,{keys:new Set(rev?rev.decisions.filter(d=>d.state==='met'||d.state==='pending').map(d=>d.code):rows.filter(z=>z.brief?.default_selected).map(z=>z.code)),manual:Object.fromEntries((rev?.decisions||[]).filter(d=>d.origin==='manual').map(d=>[d.code,d.strength])),note:rev?.note||''})}
 return unsaved.get(key);}
const active=()=>live&&!!(snap?.report?.criteria_table?.reviewable||snap?.report?.criteria_table?.draft_selectable)&&!snap.report.legacy_read_only&&(!snap.local_review||snap.local_review.schema==='brca.selection.v2');
function ruleParts(z){
 const g=snap.report.criteria_table.gene,doc=boot.rule_guide?.[g];
 return {rule:z.guideline_zh||doc?.criteria?.[z.code],source:snap.report.criteria_table.guideline_source||doc};
}
const officialByHash=new Map();
const lookupPointer=(root,loc)=>loc.replace(/^\//,'').split('/').reduce((a,k)=>a?.[k.replace(/~1/g,'/').replace(/~0/g,'~')],root);
function englishHtml(z){
 const {rule,source}=ruleParts(z),doc=officialByHash.get(source?.document_sha256);if(!rule||!source)return '';
 if(!doc)return '<button type="button" class="link" data-english="'+esc(z.code)+'">对照此条英文（选择已下载的规范 JSON）</button>';
 const fields=[['VCEP 适用说明',rule.instructions],...Object.entries(rule.grades).map(([g,r])=>[g,r])],seen=new Set();
 return fields.map(([label,h])=>{if(seen.has(h.pointer))return '';seen.add(h.pointer);const en=lookupPointer(doc,h.pointer);return '<div lang="en"><b>'+esc(label)+'</b><p class="english">'+esc(en)+'</p></div>';}).join('');
}
function manualPick(z,w,can){
 const levels=z.allowed_strengths||[],only=levels.length===1,selected=only?levels[0]:w||'',{rule}=ruleParts(z),grade=rule?.grades?.[selected];
 return '<div class="manualpick" '+(!can||!buf().keys.has(z.code)?'hidden':'')+'><label><b>手工补入</b> · 强度 '+(only?'<strong>'+esc(selected)+'</strong>':'<select data-grade="'+esc(z.code)+'" aria-label="'+esc(z.code)+' 人工强度"><option value="">请选择</option>'+levels.map(g=>'<option '+(g===selected?'selected ':'')+'value="'+esc(g)+'">'+esc(g)+(rule?.grades?.[g]?.conditional_in_text?'（正文条件档）':'')+'</option>').join('')+'</select>')+'</label><p class="gradehint">'+esc(grade?.zh||'按此位点的证据选择；不同档位见下方条文。')+'</p></div>';
}
function lrLinkout(s){
 const r=s.report||{},build=s.assembly,g=r.gene,L=r.criteria_table||{};
 const cfg=({GRCh38:{BRCA1:['hg38','17','NC_000017.11'],BRCA2:['hg38','13','NC_000013.11']},
              GRCh37:{BRCA1:['hg19','17','NC_000017.10'],BRCA2:['hg19','13','NC_000013.10']}})[build]?.[g];
 const base='https://hgdownload.soe.ucsc.edu/hubs/enigma/',help=safeLink(base+'enigma.html','LR 轨道说明及修订记录');
 const locs=r.representations?.[build];let target='';
 if(L.identity_verified===true&&/^CA\d+$/.test(r.caid||'')&&cfg&&Array.isArray(locs)&&locs.length===1&&typeof locs[0]==='string'){
   const [db,chr,acc]=cfg,prefix=acc+':g.';const m=locs[0].startsWith(prefix)&&/^(\d+)(?:_(\d+))?(?:[ACGT]>[ACGT]|delins[ACGT]+|del[ACGT]*|dup[ACGT]*|ins[ACGT]+)$/.exec(locs[0].slice(prefix.length));
   if(m){
     const lo=Number(m[1]),hi=Number(m[2]||m[1]);
     if(Number.isSafeInteger(lo)&&Number.isSafeInteger(hi)&&lo>0&&hi>=lo&&hi-lo<=100){
       const p=new URLSearchParams({db,position:'chr'+chr+':'+Math.max(1,lo-10)+'-'+(hi+10),hubUrl:base+'hub.txt',hub_5570236_enigma:'show',hub_5570236_BRCAmla:'full'});
       target=safeLink('https://genome.ucsc.edu/cgi-bin/hgTracks?'+p,'UCSC 本位点：逐条展开（'+build+'）');
     }
   }
 }
 return '<p class="case-sources"><b>UCSC 轨道：</b>'+(target?target+' · ': '当前没有唯一核准的本组装区间；')+help+'。点同一核酸 ALT 的条目查看 LR 分量。</p>';
}
function parsonsSource(code){
 const combined=code==='PP4'||code==='BP5',seg=code==='PP1'||code==='BS4';
 if(!combined&&!seg)return '';
 const xls=safeLink('https://pmc-oa-opendata.s3.amazonaws.com/PMC6772163.1/HUMU-40-1557-s001.xlsx','下载 Parsons 2019 原始 Excel');
 const paper=safeLink('https://pmc.ncbi.nlm.nih.gov/articles/PMC6772163/','论文及补充资料');
 const which=combined?'K 列 · 综合 LR':'F 列 · 共分离 LR';
 const note=combined?'按综合值的方向与规范阈值核 PP4／BP5；不要将作者表与 UCSC 对同一资料的汇总重复计分。':'只用独立的共分离证据；若该分量已计入采用的综合 LR，不可再次计分。';
 return '<div class="case-sources" data-parsons-column="'+(combined?'K':'F')+'"><b>作者原表（手工查）：</b>'+xls+' · '+paper+'。<strong>'+which+'</strong>。<details><summary>怎样查这张表</summary><p>打开 <code>SuppT1 All Data</code>，在 B 列按 c. 查找，核对 A 列基因。旧表用 BRCA1 NM_007294.3、BRCA2 NM_000059.3，须先确认与当前转录本是同一核酸等位。用 <code>SuppT2 Data Sources</code> 追查数据来源。</p><p>'+note+' L/M 列是先验／后验概率。链接不代表本位点已命中；空白或未收录不计为反向证据。</p></details></div>';
}
function rowHtml(z){
 const b=z.brief||{},op=b.option,x=z.expert,can=active(),rv=z.decision;
 const {rule,source}=ruleParts(z),manual=!op,w=can?buf().manual[z.code]:rv?.origin==='manual'?rv.strength:'';
 const hard=z.evidence_view?.hold_reason==='spec_inapplicable';
 const checked=can?buf().keys.has(z.code):['met','pending'].includes(rv?.state);
 let grade=rv?.strength||op?.strength||w||(x?.strength&&z.allowed_strengths.includes(x.strength)?x.strength:null);
 if(!grade&&checked&&z.allowed_strengths.length===1)grade=z.allowed_strengths[0];
 const marker=rv?(rv.state==='pending'?'已保存待核 · 不计分':'已保存 · '+(means[rv.origin]||'历史审核')):op?means[op.kind]:b.status==='conflict'?'资料冲突':b.status==='expert_not_met'?'专家未采用':'';
 const box='<label class="pickone"><input type="checkbox" data-adopt="'+esc(z.code)+'" '+(checked?'checked ':'')+(!can?'disabled ':'')+'aria-label="计入 '+esc(z.code)+'"> '+(hard?'留待核':'计入')+'</label>';
 const head='<strong>'+esc(z.code)+(grade?' <span class="rowweight">'+esc(grade)+'</span>':'')+'</strong><span class="rowlabel">'+esc(z.label)+'</span>';
 let src='';
 if(rule&&source){
  src+=(!manual||hard?'<p><b>条文中文释义：</b>'+esc(rule.zh)+'</p>':'')+'<dl class="zhgrades">'+z.allowed_strengths.map(g=>'<div><dt>'+esc(g)+'</dt><dd>'+esc(rule.grades[g]?.zh)+'</dd></div>').join('')+'</dl>';
  src+='<p>'+safeLink(source.publication_url,'ClinGen '+source.spec)+' · '+safeLink(source.document_url,'下载官方英文 JSON')+' · 项目中文释义，非官方译本。</p>';
  src+='<details><summary>英文原文与 JSON 定位</summary><div data-original="'+esc(z.code)+'">'+englishHtml(z)+'</div><p>字段标识：<code>'+esc(rule.rule_key)+'</code></p>'+jsonView({document_sha256:source.document_sha256,instructions:rule.instructions,grades:Object.fromEntries(Object.entries(rule.grades).map(([k,v])=>[k,{pointer:v.pointer,text_sha256:v.text_sha256,conditional_in_text:v.conditional_in_text}]))})+'</details>';
 }else src='<p><b>本基因标准：</b>'+esc(b.standard||'见规范条文')+' '+safeLink(z.rule_url,'ClinGen CSpec')+'</p>';
 if(['PP1','BS4'].includes(z.code)){
   const guide=snap.report.criteria_table?.literature;
   if(guide?.allele_verified&&guide.aliases?.length){
    const term=snap.report.criteria_table.gene+' AND ('+guide.aliases.map(a=>'"'+a+'"').join(' OR ')+') AND (segregation OR cosegregation OR pedigree OR "likelihood ratio")';
    src+='<p>'+safeLink('https://pmc.ncbi.nlm.nih.gov/search/?term='+encodeURIComponent(term),'PMC 全文：按此等位找共分离、家系及 LR')+' · '+safeLink('https://bioexp.net/cosegregation/','ClinGen 条文所列 CAL-Leiden 方法')+'。这是人工复核入口，不代表本文已报告该等位 LR。</p>';}
 }
 if(z.conditional_rule){const q=z.conditional_rule;src+='<p>正文条件档 Very Strong：LR '+(q.direction==='for'?'≥':'≤')+esc(q.threshold)+'，'+esc(q.condition)+'</p>';}
 if(op){src+='<p>采用方式：'+esc(means[op.kind])+'。</p>';}
 for(const c of z.candidates||[]){src+='<p>'+safeLink(c.workbook_url,'ENIGMA 工作簿')+' · '+esc(c.sheet)+' 第 '+esc(c.row)+' 行'+(c.pmids?.length?' · '+c.pmids.map(p=>safeLink('https://pubmed.ncbi.nlm.nih.gov/'+p+'/','PMID '+p)).join('、'):'')+'</p>';}
 const stop=snap.report.facts?.stop_exon_review;
 if(stop&&['PVS1','PM5'].includes(z.code)){src+='<p>'+safeLink(stop.table4_url,'Table 4 · 官方工作簿')+' · '+safeLink(stop.st1_url,'ST1 · 官方附件')+'</p><details><summary>终止落点、附件行及 SHA-256</summary>'+jsonView(stop)+'</details>';}
 const links=[...new Set((z.facts||[]).map(f=>f.source).concat((z.evidence_view?.observations||[]).map(f=>f.source)).concat(x?['erepo_detail','erepo_summary']:[]))];if(links.length)src+='<p>'+rawLinks(links)+'</p>';
 if(x){src+='<p>'+safeLink(snap.report.erepo?.url,'专家个案')+' · '+(x.state==='Not Met'?'明确不满足':'采用 '+esc(x.strength||'未注明'))+'</p>';
  if(x.text?.length)src+='<details><summary>专家原文或教学摘录</summary>'+x.text.map(n=>'<p class="subblock">'+esc(n.value)+'</p>').join('')+'</details>';}
 if(z.evidence_view?.clinical_lr){src+='<details><summary>临床 LR 类型、原记录定位与分档</summary>'+jsonView(z.evidence_view.clinical_lr)+'</details>';}
 if((z.facts||[]).length||z.candidates?.length||rv||op)src+='<details><summary>原始字段、行号与哈希</summary>'+jsonView({facts:z.facts,attachment:z.candidates,source:op?.refs,saved:rv})+'</details>';
 const view=z.evidence_view||{},modern=snap.report.acquisition==='live_independent',values=view.observations||[];
 const display=v=>{const val=typeof v.value==='object'?'<details><summary>查看返回字段</summary>'+jsonView(v.value)+'</details>':esc(v.value);return '<div class="measure"><span>'+esc(v.label)+'：</span>'+val+'</div>';};
 const vv=values.filter(v=>v.kind!=='expert_assertion'),isContext=v=>z.code==='PM2'?hard||['gno_full37','gno_full38','ldh_exome','ldh_genome'].includes(v.source):['BS1','BA1'].includes(z.code)&&!v.label.includes('FAF')&&v.source!=='population_file';
 const context=vv.filter(isContext),groups=vv.filter(v=>!isContext(v)).map(display).join('')+(context.length?'<details class="popcontext"><summary>'+(hard?'背景人群计数与查询状态，不用于 PM2（':'查看其他人群计数（')+context.length+'）</summary>'+context.map(display).join('')+'</details>':'');
 let actual=values.length?groups:'';
 let outcome=view.comparison?'<p class="casecomparison"><b>对照结果：</b>'+esc(view.comparison)+'</p>':'';
 if(modern&&!values.length&&!view.comparison&&['missing',undefined].includes(b.status))actual='<div class="muted">'+(view.state==='lookup_incomplete'?'未取得可评资料。':'暂无可评资料。')+'</div>';
 if(view.standard&&(values.length||view.comparison)){if(hard||view.clinical_lr)src='<p><b>'+(view.clinical_lr?'规范 LR 阈值，不是本位点实测值：':'完整适用条件：')+'</b>'+esc(view.standard)+'</p>'+src;else outcome+='<div class="casenote">标准：'+esc(view.standard)+'</div>';}
 const mini=modern||outcome?'':z.code==='PM2'?'<div class="casenote">本版标准：SNV；v2.1 非癌外显子和 v3.1 非癌均未见，区域平均深度均 ≥25×。</div>':z.code==='PP4'&&op?.strength==='Supporting'?'<div class="casenote">标准：2.08 ≤ 临床综合 LR < 4.3。</div>':'';
 if(view.blocks_automatic_scoring&&modern)outcome+='<div class="casenote">可勾选保存为「待核」；本次不计入积分。</div>';
 const note=(b.expert_note&&b.result!==b.expert_note&&(!op||op.kind==='expert_citation'||z.code==='PM2'))?'<div class="casenote">'+esc(b.expert_note)+'</div>':'';
 const assess=z.assessment||snap.view_overlay?.assessments?.[z.code]||{};
 const tag=assess.label?'<div class="assessment" data-state="'+esc(assess.state)+'"><b>本次判定：'+esc(assess.label)+'</b></div>':'';
 const via=(assess.sources||[]).map(o=>{
  if(o.key==='workbook'&&o.url)return safeLink(o.url,o.label);
  if(snap.raw?.[o.key])return gotoRaw(o.key);
  if(o.key==='erepo_detail'&&snap.report.erepo?.url)return safeLink(snap.report.erepo.url,'eRepo 专家个案');
  const item=(snap.report.facts?.gnomad||{})[o.key];if(item?.record_url)return safeLink(item.record_url,names[o.key]||o.label);
  return '<button type="button" class="link" data-raw="current_report">'+esc(names[o.key]||o.label)+(snap.teaching?'（教学摘录）':'（报告字段）')+'</button>';
 });
 const assoc=via.length?'<div class="case-sources"><b>来源／查询记录：</b>'+[...new Set(via)].join(' · ')+'</div>':'';
 const directRule=parsonsSource(z.code)+(['PP4','BP5'].includes(z.code)?lrLinkout(snap):'')+'<span class="case-sources">'+safeLink(z.rule_url,'本基因规范')+'</span>';
 return '<div class="crit brief-row" tabindex="-1" id="criterion-'+esc(z.code)+'" data-code="'+esc(z.code)+'">'+box+'<div><div class="crithead">'+head+(marker?'<small class="tag">'+esc(marker)+'</small>':'')+'</div>'+tag+(hard?outcome:'')+'<div class="finding">'+(rv?.origin==='manual'?'<b>自动检索：</b>':'')+(actual||esc(b.result||'尚无可用取数。'))+'</div>'+(hard?'':outcome)+mini+note+assoc+directRule+(manual&&rule&&!hard?'<p class="zhlead"><b>指南中文：</b>'+esc(rule.zh)+'</p>':'')+(manual&&can&&!hard?manualPick(z,w,can):'')+'<details class="sourced"><summary>逐档条文、英文与出处</summary>'+src+'</details></div></div>';
}
function evidence(){
 const r=snap.report,L=r.criteria_table;if(!L){$('evidencePane').innerHTML='<p>未获得本基因判据表。</p>';return}
 const rows=L.rows.filter(z=>z.rule_enabled),rel=rows.filter(z=>z.priority),rest=rows.filter(z=>!z.priority);let h='<h2>本次证据</h2><p class="muted">'+esc(L.spec||'规范未核准')+' · 本基因 '+rows.length+' 条可用主判据；本例优先显示 '+rel.length+' 条。</p>';
 h+='<details class="scorelegend"><summary>怎样读 BayesDel、临床 LR 和证据积分</summary><p>BayesDel noAF 是预测器分数，可为负；不是患病概率。共分离 LR 评估家系分离（PP1／BS4）；临床综合 LR 可含家族史等组成（PP4／BP5）。资料相同不能重复使用。本地 +1、+2、+4、+8 等表示被采纳判据的强度，不是从预测分数或数据库条目直接换算的后验概率。</p></details>';
 if(active()){h+='<p class="smallnote">有自行核实的位点证据，可展开其余判据手工计入。没有资料的项不会预选；指南条文只说明使用条件。</p>';h+='<div class="briefsave" id="reviewform"><div class="row"><button type="button" data-save class="primary">校验、保存并计算</button> <span id="pending"></span></div><label class="single-note">本次备注（选填）<input type="text" data-note placeholder="可留空直接保存" value="'+esc(buf().note)+'"></label><p id="saveNote" role="status" aria-live="polite"></p></div>'}
 else if(r.legacy_read_only)h+='<p class="note">此为旧版已保存审核，只读；原始报告保留。要使用新流程可另存查询快照。</p>';
 if(!live)h+='<p class="smallnote">这是只读快照；复选框反映已经保存的选择。要勾选或补录，请在本机启动后的「查询」页打开同一快照。</p>';
 if(L.alerts.length)h+='<details><summary>取数提示（'+L.alerts.length+'）</summary><p>'+L.alerts.map(esc).join('；')+'</p></details>';
 h+=rel.length?rel.map(rowHtml).join(''):'<p>本次未检出可优先展示的证据。</p>';
 if(rest.length)h+='<details id="allRows" '+(!L.identity_verified||L.all_evidence_missing||!rel.length?'open':'')+'><summary>'+(!rel.length?'本基因全部 ':'展开其余 ')+rest.length+' 条可用判据</summary>'+rest.map(rowHtml).join('')+'</details>';
 h+=readingStrip(true);
 if(L.other_sources?.length){h+='<details><summary>其他已返回的线索（没有自动赋码）</summary>';for(const a of L.other_sources){
  if(a.pointer==='/facts/ldh_associations_only'){h+='<p>LDH 实体目录列有以下关联。目录命中不代表已读过实体或取得判据；下列链接可手工复核。</p><ul>';for(const [key,list] of Object.entries(a.value||{})){if(!Array.isArray(list)||!list.length)continue;const urls=[...new Set(list.map(v=>v?.dataIri||v?.ldhIri).filter(u=>typeof u==='string'&&u.startsWith('https://')))];h+='<li>'+esc(key)+'（'+list.length+'） '+urls.slice(0,2).map((u,i)=>safeLink(u,'实体 '+(i+1))).join(' · ')+'</li>';}h+='</ul>'+rawLinks(['ldh_med']);}
  else if(Array.isArray(a.value)){h+='<p>'+esc(a.label)+'</p><ul>'+a.value.map(v=>'<li>'+esc(v.sheet||v.table)+' 第 '+esc(v.row)+' 行'+(v.pmids_mentioned?.length?' · '+v.pmids_mentioned.map(id=>safeLink('https://pubmed.ncbi.nlm.nih.gov/'+encodeURIComponent(id)+'/','PMID '+id)).join('、'):'')+'</li>').join('')+'</ul>';}
 }h+='<details><summary>所有原始关联字段</summary>'+jsonView(L.other_sources)+'</details></details>';}
 const cv=r.facts?.clinvar_record,agg=cv?.germline_aggregate;
 if(agg)h+='<details><summary>ClinVar 已有提交 · '+esc(cv.accession)+'（不计分）</summary><p>'+esc(agg.Description)+' · '+esc(agg.Explanation||agg.ReviewStatus)+'</p>'+rawLinks(['clinvar'])+'</details>';
 if(r.expert_assertions?.length)h+='<details><summary>ClinVar 专家提交</summary>'+table(['分类','提交版本'],r.expert_assertions.map(x=>[esc(x.class),esc(x.scv)]))+'</details>';
 h+='<p class="sp"><button type="button" data-ledger-download>下载本页报告 JSON</button></p>';
 if(!live)h+='<p class="smallnote">快照只读。'+(location.protocol==='http:'&&['127.0.0.1','localhost'].includes(location.hostname)?'<a href="/query?gene='+encodeURIComponent(L.gene)+'&amp;assembly='+encodeURIComponent(snap.assembly)+'&amp;snapshot='+encodeURIComponent(snap.id)+'">在本机查询页采纳并保存</a>':'启动本机程序后，在「查询」里选同一示例即可直接保存。')+'</p>';
 if(snap.local_review?.schema==='brca.selection.v2'){const v=snap.local_review;h+='<p class="smallnote">已保存 '+esc(v.revision)+' · '+esc(v.created_utc)+' UTC'+(v.note?' · 备注：'+esc(v.note):'')+'。未勾选仍为未采纳。</p>';if(v.result.issues.length)h+='<p class="note">'+v.result.issues.map(esc).join('；')+'</p>';}
 $('evidencePane').innerHTML=h;pending();
}
function pending(){if($('pending')&&snap){const b=buf(),keys=b.keys,rows=snap.report.criteria_table.rows,missing=rows.filter(z=>keys.has(z.code)&&snap.report.criteria_table.identity_verified&&!z.evidence_view?.blocks_automatic_scoring&&!z.brief?.option&&!b.manual[z.code]);$('pending').textContent='本次选择 '+keys.size+' 项'+(keys.size?'：'+[...keys].join('、'):'')+'。'+(missing.length?'待选强度：'+missing.map(z=>z.code).join('、')+'。':'');}}
document.addEventListener('input',e=>{if(e.target.matches('[data-note]')&&snap&&!saveBusy)buf().note=e.target.value;});
document.addEventListener('change',e=>{
 if(!active()||saveBusy||!e.target.matches('[data-adopt],[data-grade]'))return;
 const b=buf(),code=e.target.dataset.adopt||e.target.dataset.grade,z=snap.report.criteria_table.rows.find(n=>n.code===code),row=e.target.closest('[data-code]');if(!z||!row)return;
 if(e.target.matches('[data-adopt]')){
  if(e.target.checked){b.keys.add(code);if(!z.brief?.option&&z.allowed_strengths.length===1)b.manual[code]=z.allowed_strengths[0];}
  else {b.keys.delete(code);delete b.manual[code];const sel=row.querySelector('[data-grade]');if(sel)sel.value='';}
 }else if(b.keys.has(code)){if(e.target.value)b.manual[code]=e.target.value;else delete b.manual[code];}
 const pick=row.querySelector('.manualpick');if(pick){pick.hidden=!b.keys.has(code);const {rule}=ruleParts(z);pick.querySelector('.gradehint').textContent=rule?.grades?.[b.manual[code]]?.zh||'按此位点的证据选择；不同档位见下方条文。';}
 pending();
});
function chooseOriginal(code){
 const z=snap.report.criteria_table.rows.find(r=>r.code===code),src=ruleParts(z).source,box=document.querySelector('[data-original="'+code+'"]');
 if(!src||!box)return;const f=document.createElement('input');f.type='file';f.accept='.json,application/json';
 f.onchange=async()=>{if(!f.files.length)return;try{
  const file=f.files[0];if(file.size>3000000)throw Error('文件过大。');if(!globalThis.crypto?.subtle)throw Error('本浏览器不能核验哈希，请在本机查询页打开。');
  const data=await file.arrayBuffer(),h=[...new Uint8Array(await crypto.subtle.digest('SHA-256',data))].map(v=>v.toString(16).padStart(2,'0')).join('');
  if(h!==src.document_sha256)throw Error('文件与本基因 v1.2 的 SHA-256 不同，未使用。');
  officialByHash.set(h,JSON.parse(new TextDecoder().decode(data)));
  if(snap.report.criteria_table.guideline_source?.document_sha256===h||boot.rule_guide?.[snap.report.criteria_table.gene]?.document_sha256===h)
   for(const row of snap.report.criteria_table.rows){const target=document.querySelector('[data-original="'+row.code+'"]');if(target)target.innerHTML=englishHtml(row);}
 }catch(err){box.textContent=err.message+' ';const retry=document.createElement('button');retry.type='button';retry.dataset.english=code;retry.textContent='重新选取';box.append(retry);}};f.click();
}
document.addEventListener('click',e=>{const c=e.target.closest('[data-english]');if(c)chooseOriginal(c.dataset.english);});
async function saveLocal(){
 if(!active()||saveBusy||busy)return;const s=snap,rv=s.local_review,b=buf();
 const missing=s.report.criteria_table.rows.find(z=>b.keys.has(z.code)&&s.report.criteria_table.identity_verified&&!z.evidence_view?.blocks_automatic_scoring&&!z.brief?.option&&!b.manual[z.code]);
 if(missing){const row=document.querySelector('[data-code="'+missing.code+'"]');if(row?.closest('#allRows'))$('allRows').open=true;row?.querySelector('select')?.focus();$('saveNote').textContent=missing.code+'：请选择人工强度，再保存。';return;}
 const body={schema:'brca.selection.v2',snapshot:s.id,base_hash:s.report_sha256,caid:s.caid,gene:s.report.criteria_table.gene,assembly:s.assembly,previous_revision:rv?.revision||null,selected:[...b.keys],manual:Object.fromEntries(Object.entries(b.manual).filter(([k,v])=>b.keys.has(k)&&typeof v==='string'&&v.length)),note:b.note};
 saveBusy=true;lists();const pane=$('evidencePane');pane.inert=true;$('saveNote').textContent='正在保存…';
 try {const reply=await fetch('/api/review',{method:'POST',headers:{'Content-Type':'application/json','X-Review-Token':boot.token},body:JSON.stringify(body)}),v=await reply.json();if(!reply.ok)throw Error(v.error||'保存失败');
  unsaved.delete(s.id+'/'+(rv?.revision||'base'));docs.set(v.id,v);snap=v;paint();if($('saveNote'))$('saveNote').textContent='已保存 '+v.local_review.revision+'。'+(v.local_review.result.issues.length?'组合待核。':'结果见上方合计。');}
 catch(e){if($('saveNote'))$('saveNote').textContent='未保存：'+e.message;}finally{saveBusy=false;pane.inert=false;lists();}}
document.addEventListener('click',e=>{if(e.target.closest('[data-save]'))saveLocal();if(e.target.closest('[data-ledger-download]')&&snap)download(new TextEncoder().encode(jstr(snap.report)+'\n'),snap.id+'-'+(snap.local_review?.revision||'view')+'-report.json');});
