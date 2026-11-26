"""Source-attributed clinical LR observations, kept outside the score engine."""
from .exchange import text_lr,numeric
from .criteria import lr_policy
def clinical_lr(r,z,put):
    k=z['code'];fam=k in ('PP1','BS4');typ='segregation' if fam else 'combined'
    display='定量共分离 LR' if fam else '临床综合 LR'
    cfg=lr_policy()['criteria'][k];bounds=lr_policy()['directional_thresholds'][cfg['direction']]
    std=('仅比较'+display+'；规范分档（不是实测值）：')+'；'.join(w+(' ≥' if cfg['direction']=='for' else ' ≤')+bounds[w] for w in ('Supporting','Moderate','Strong','Very Strong'))+'。'
    f=r.get('facts') or {};detail=f.get('brca_detail') or {};row=detail.get('source_row') or {}
    m=(detail.get('measures') or {}).get('segregation') if fam else None
    expert=z.get('expert') or {};data=text_lr(expert.get('text') or [],typ) if (r.get('erepo') or {}).get('state')=='matched' else []
    measured=[];read=[]
    for it in data:
        put('专家个案引用 · '+display,it['value'],'erepo_detail',it['pointer'],'expert_measure')
        measured.append(('专家引用',it['value']));read.append(dict(source='erepo_detail',**it))
    if not fam and expert:
        for it in text_lr(expert.get('text') or [],'family_history'):
            put('专家个案 · 家族史 LR（综合值的组成，不作 PP1）',it['value'],'erepo_detail',it['pointer'],'expert_component')
    if m:
        if m.get('value') is not None:
            put('BRCA Exchange · 已报告的共分离 LR',m['value']+' · release '+str(row.get('release')),'brca_detail','/facts/brca_detail/measures/segregation','external_measure')
            measured.append(('BRCA Exchange（历史 exLOVD 字段，方法待核）',m['value']));read.append(dict(source='brca_detail',value=m['value'],type=typ,pointer=m['raw_pointer']))
        else:
            msg='原字段为占位符 '+repr(m['raw']) if m.get('state')=='not_reported' else '返回行没有此字段' if m.get('state')=='not_in_record' else '该字段无可解析数值'
            put('BRCA Exchange · Segregation_LR_exLOVD',msg+'；没有可用 LR。','brca_detail','/facts/brca_detail/measures/segregation','not_reported')
    none='补证：查该核酸等位已发表的家系共分离分析，或用有亲缘、携带状态和表型的家系按附录 I 建模。家族史 LR、综合 LR 和病例篇数不能代替。' if fam else '补证：查经校准且组成可追溯的临床综合 LR；单一未校准组件、先验概率、后验概率均不能代入。'
    if not measured:
        if fam and not m and r.get('source_status',{}).get('brca_detail') in ('http_error','transport_error','invalid_fields','allele_or_release_unverified'):
            none='该次 BRCA Exchange 完整记录未核准；'+none
        cmp='本次没有取得此等位的'+display+'；数值为空，不能当作 0 或 1。'+none
    else:
        vals={numeric(v) for _,v in measured};cmp=''
        if len(vals)==1:
            n=vals.pop();w=next((w for w in ('Very Strong','Strong','Moderate','Supporting') if n>0 and (n>=numeric(bounds[w]) if cfg['direction']=='for' else n<=numeric(bounds[w]))),None)
            cmp='本位点'+display+' = '+measured[0][1]+'；'+('数值落入 '+w+' 档' if w else '不直接提供本方向强度')+'。须核定量方法、家系／队列和是否已用于其他判据；不在本地重算此 LR。'
        else:cmp='来源报告的 LR 不一致；保留各值，不能挑选或相乘。'
    return dict(type=typ,conflicting=len({numeric(v) for _,v in measured})>1,observations=read,reported_values=[a['value'] for a in read],thresholds=bounds),std,cmp

