"""Describe what the source actually measured, including nonqualifying observations."""
import re
from .routes import OUTBRED
from .criteria import protein_domain,PROTEIN,dec

BAD={'http_error','transport_error','graphql_error','unavailable','host_cooldown','invalid_response','invalid_fields','allele_not_verified','unresolved_links','directory_unavailable'}
def numbers(d):
    return f"AC {d['ac']:,} / AN {d['an']:,}"+(f"；AF {d['af']:.6g}" if d.get('af') is not None else '')
def observation(r,z,ok):
    k=z['code'];f=r.get('facts') or {};g=r.get('selected_gene') or r.get('gene')
    br=z['brief'];ep=z.get('expert');gn=f.get('gnomad') or {};states=r.get('source_status') or {};items=[];gaps=[];comparison=None;hold=False;lr=None;hold_reason=None
    external=f.get('linked_observations_unverified') or {};state='no_evaluable_data';relevant=False
    def put(label,value,source,ptr,kind='observed'):
        nonlocal relevant
        if value is not None:items.append(dict(label=label,value=value,source=source,pointer=ptr,kind=kind));relevant=True
    if k in ('PM2','BS1','BA1'):
        seen=[];unusable=False;eligible_faf=[];noncancer_record=False;statname=''
        for key in ('gno_nc2','gno_nc3','gno_full37','gno_full38'):
            q=gn.get(key)
            if not q:continue
            label=q['label'];s=q.get('state');ptr='/facts/gnomad/'+key;target=key in ('gno_nc2','gno_nc3')
            if s=='no_variant_record':
                put(label,'该版变异接口明确未见记录；没有给出 AN 或覆盖。',key,ptr,'searched_no_record')
                continue
            if s=='no_verified_query_coordinate':
                put(label,'未查询：缺少此组装已核准、带 REF 锚的 VCF 等位；收录和计数未知。',key,ptr+'/state','not_queried');continue
            if s!='variant_record':
                unusable=True
                if s in BAD:
                    detail='（'+q['service_message']+'）' if q.get('service_message') else ''
                    put(label,'查询未完成'+detail+'；收录情况和 AC 未知。',key,ptr+'/state','source_failure')
                continue
            noncancer_record|=target
            for kind,a in q.get('cohorts',{}).items():
                suf='' if target else (' 外显子' if kind=='exome' else ' 基因组')
                if a is None:
                    put(label+suf,'未返回该队列记录；AC 未知。',key,ptr+'/cohorts/'+kind,'null_field');continue
                base=ptr+'/cohorts/'+kind
                msg=numbers(a)
                if a.get('filters'):msg+='；过滤标记 '+','.join(map(str,a['filters']))
                put(label+suf,msg,key,base)
                if not target:continue
                pops=a.get('populations') or []
                vals=[(i,x) for i,x in enumerate(pops) if x.get('id') in OUTBRED and type(x.get('ac'))==int and type(x.get('an'))==int and 0<x['ac']<=x['an']]
                if vals and a.get('filters')==[]:
                    seen.extend([x for _,x in vals])
                    if k=='PM2':
                        for i,x in vals:
                            put(label+' · '+x['id'].upper(),f"AC {x['ac']:,} / AN {x['an']:,}",key,base+'/populations/'+str(i))
                ff=a.get('faf95') or {};n=dec(ff.get('popmax'));p=ff.get('popmax_population')
                if k in ('BA1','BS1') and n is not None and 0<=n<=1:
                    put(label+' · FAF95',str(n)+'（'+str(p or '祖源未披露')+'）',key,base+'/faf95')
                    if p in OUTBRED and a.get('filters')==[]:eligible_faf.append(n)
        for source,label in [('ldh_exome','外显子'),('ldh_genome','基因组')]:
            data=f.get(source) or external.get(source) or {}
            if type(data.get('AC'))==int and type(data.get('AN'))==int:
                note='；实体明确 noVariant=true' if data.get('noVariant') is True else ''
                qc=data.get('quality') or {};depth=qc.get('meanDepth');cov20=qc.get('fractionCoverage20x')
                if type(depth) in (int,float) and depth>=0:note+=f'；该实体平均深度 {depth:.3g}×'
                if type(cov20) in (int,float) and 0<=cov20<=1:note+=f'；≥20× 覆盖比例 {100*cov20:.2f}%'
                put('LDH · v4.1 '+label+'（背景）',f"AC {data['AC']:,} / AN {data['AN']:,}"+(f"；AF {data['AF']:.6g}" if type(data.get('AF')) in (int,float) else '')+note,source,'/facts/'+('' if f.get(source) else 'linked_observations_unverified/')+source)
        # Display provided local strata as observations regardless of the result; imports have their own source pointers.
        from .briefing import population
        pop,qualified=population(r)
        for i,x in enumerate((pop.get('inputs') or {}).get('strata') or []):
            if not isinstance(x,dict):continue
            vals=[a+' '+str(x[n]) for n,a in (('AC','AC'),('AN','AN'),('faf','FAF'),('mean_depth','区域平均深度 ×')) if x.get(n) is not None]
            if (k=='PM2' and (x.get('dataset'),x.get('assay')) in (('gnomAD_v2.1_non_cancer','exome'),('gnomAD_v3.1_non_cancer','genome')) and x.get('outbred') is True and x.get('qc_pass') is True and x.get('source_ref') and type(x.get('AC'))==int and type(x.get('AN'))==int and 0<x['AC']<=x['AN']):seen.append(x)
            if vals:put(str(x.get('dataset','本地队列'))+' · '+str(x.get('ancestry','')),'；'.join(vals),'population_file',f'/facts/population_evidence/inputs/strata/{i}')
        if k=='PM2' and (f.get('event') in ('DEL','INS','DELINS') or f.get('event')=='MNV' and 'delins' in str((f.get('transcript') or {}).get('hgvsc'))):
            state='not_applicable';hold=True;hold_reason='spec_inapplicable'
            comparison='规范排除插入、缺失及 delins。本例 DNA 事件为 '+f['event']+'；PM2 不适用，即使两队列覆盖充分且 AC=0 也不计分。其余人群数值仅供核查。'
        elif k=='PM2' and seen:
            state='does_not_meet';comparison='规定队列的远交人群中已观察到 ALT，不满足“两队列均未见”。单次观察不提供致病或良性分。'
            hold=True;hold_reason='incompatible_observation'
        elif k=='PM2' and gn and qualified.get('code')!='PM2':
            state='incomplete';comparison='现有结果未证实规范两队列均缺失且区域平均深度均 ≥25×。'
        elif k in ('BA1','BS1') and gn and qualified.get('code')!=k:
            state='incomplete'
            if eligible_faf:
                n=max(eligible_faf);w='超过 0.001' if k=='BA1' and n>dec('.001') else '未超过 0.001' if k=='BA1' else '超过 Strong 阈值 0.0001' if n>dec('.0001') else '落入 Supporting 区间 (0.00002, 0.0001]' if n>dec('.00002') else '未超过 Supporting 下限 0.00002'
                comparison=f'已返回的合格祖源 FAF95 中，最大值 {n}，{w}。仍须核队列、覆盖 ≥20× 及创始变异例外。'
            else:comparison='没有取得可与本条阈值比较的非创始人群 FAF；总 AF 不能代替 FAF。'
        if state!='not_applicable' and unusable and not any(x['kind'] not in ('source_failure','not_queried') for x in items):state='lookup_failed';comparison='人群接口未查成；不能判断是否收录。'
    if k in ('PP3','BP4','BP1') and ok:
        from .mechanism import interpretation
        state,comparison,hold=interpretation(r,k,put)
        if hold:hold_reason='spec_inapplicable'
    if k in ('PVS1','PM5') and ok:
        from .mechanism import classify,CATS
        t=f.get('transcript') or {};mech=classify(t,g);pt=f.get('stop_exon_review') or {}
        if mech['category'] in ('frameshift','stop','canonical_splice','start_loss'):
            src=(f.get('transcript_sources') or {}).get('hgvsp') or (f.get('transcript_sources') or {}).get('consequence_terms','vep')
            put('转录本后果',CATS[mech['category']]+' · '+str(t.get('hgvsp') or t.get('hgvsc')),src,'/facts/transcript')
            if mech['stop_aa']:put('推导的突变蛋白终止位置','第 '+str(mech['stop_aa'])+' 个密码子（据 p. 描述推导，不是 NMD 实验）',src,'/facts/transcript/hgvsp','derived')
            if pt:
                a=pt['table4'];b=pt['span'];lo,hi=pt['reference_codon']
                put('终止密码子映回参考转录本',f"c.{lo}–{hi}；外显子 {pt['exon']}",'workbook','/facts/stop_exon_review','derived')
                put('冻结规范行',f"Table 4 第 {a['row']} 行（{a['cell']}） · "+str(a['code'] if k=='PVS1' else a['qualifier'])+f"；ST1 第 {b['row']} 行。",'workbook','/facts/stop_exon_review/table4','workbook_context')
                comparison='提前终止密码子所在外显子已定位。核对 Table 4 完整行、可挽救转录本和 NMD 条件；此为分支候选，尚未采纳。' if k=='PVS1' else '根据终止所在外显子核 PM5_PTC，先满足适当的 PVS1；单凭发生移码的位置不能赋码。'
            elif mech['category']!='canonical_splice':comparison='已有提前终止／起始丢失线索。终止所在外显子、转录本挽救和 NMD 等未全部核准；不自动给强度。'
            else:comparison='经典剪接位点走 PVS1 剪接流程；不得用预测剪接造成的 PTC 直接套 PM5_PTC。'
            state='incomplete'
    for item in z.get('facts',[]):
        if any(a['pointer']==item['pointer'] for a in items) or (k in ('PM2','BA1','BS1') and item['source']=='ldh_exome'):continue
        put(item['label'],item['value'],item['source'],item['pointer'])
    for c in z.get('candidates',[]):
        put('指南附件',f"{c['sheet']} 第 {c['row']} 行："+str(c.get('raw_code') or '')+' '+str(c.get('strength') or ''),'workbook',c['pointer'])
    if ep:
        put('eRepo 专家个案',('明确不满足' if ep['state']=='Not Met' else '采用 '+str(ep.get('strength') or '未注明强度') if ep['state']=='Met' else ep['state']),'erepo_detail',ep.get('pointer'),'expert_assertion')
        if ep['state']=='Not Met' and state=='no_evaluable_data':state='expert_not_met'
    if k in ('PP4','BP5','PP1','BS4') and ok:
        from .lineage import clinical_lr
        lr,br['standard'],comparison=clinical_lr(r,z,put)
        state='observed' if lr['observations'] else 'expert_not_met' if ep and ep.get('state')=='Not Met' else 'no_evaluable_data'
        if lr.get('conflicting'):state='source_conflict';hold=True;hold_reason='source_conflict'
    if not items and k in ('PS3','BS3') and ok:
        if states.get('workbook') in ('no_matching_row','matched_stop_exon_only') or states.get('workbook')=='matched_rows' and not any((a.get('source') or {}).get('table')=='T9' for a in r.get('candidates',[])) and not any(a.get('table')=='T9' for a in f.get('workbook_rows_without_assignment',[])):
            state='searched_no_row';comparison='冻结 Table 9 未见该等位行；不能据此推断功能实验阴性。'
        elif states.get('workbook')=='unavailable':state='lookup_failed';comparison='冻结工作簿未取得，未能查核。'
    if ok and comparison and (hold or k in ('PP1','BS4','PP4','BP5','PVS1','PM5') and not br.get('option')):
        br['result']=comparison
    if hold:
        br.update(option=None,default_selected=False,status=hold_reason or 'incompatible_observation')
    elif br.get('option') and state not in ('does_not_meet','not_applicable'):
        state='source_available'
    elif state=='no_evaluable_data' and items and lr is None:state='observed'
    if not ok and items:state='identity_unresolved';comparison='已保存该查询的返回值；等位身份未完成核准，不计分。';hold=True
    if not ok and not items and states and any(s in BAD for s in states.values()) and not comparison:state='lookup_incomplete'
    return dict(state=state,applicability=('not_applicable' if state=='not_applicable' else 'applicable' if ok and k=='PM2' and f.get('event')=='SNV' else None),clinical_lr=lr,hold_reason=hold_reason,observations=items,standard=br['standard'],comparison=comparison,blocks_automatic_scoring=hold,has_observations=relevant)
