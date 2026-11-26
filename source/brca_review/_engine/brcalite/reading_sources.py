"""Reading assignments with allele- and workbook-specific provenance; no scoring."""
FINDLAY = {
    'url':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117159/suppl/GSE117159_Supplementary_Table_1.xlsx',
    'source':'GEO GSE117159, Findlay et al., Supplementary Table 1',
    'sha256':'570d84bffef01a465b57476e3a84d62550bb0c1c5d1948795b451ee981f1da36',
    'sheet':'Sheet1','header_row':4,'row':778,
    'identity':'B778:E778 = 17, 41258504 (hg19), A, C. F778 NM_007294.3; J778 c.181T>G.',
    'fields':'P778 function.score.mean; Q778 func.class; 按相同基因组等位回看原实验。',
    'status':'located_for_human_review','score_counted':False,
}
ROLES_388 = {'30209399':('表内采用的实验','F388/K388'),
             '32546644':('该 cDNA 实验在本行排除','F388/L388'),
             '32741062':('剪接／RNA 线索','G388')}


def assignments(report):
    result={}
    rows=[c['source'] for c in report.get('candidates',[])]
    rows+=report.get('facts',{}).get('workbook_rows_without_assignment',[])
    for row in rows:
        for pmid in row.get('pmids_mentioned',[]):
            role='原表提及，角色待核'
            cells=row.get('review_cells') or row.get('cell')
            entry={'origin':'guideline','criterion':row.get('code') or row.get('assigned'),
                   'sheet':row['sheet'],'row':row['row'],'sha256':row['sha256']}
            if (report.get('caid')=='CA001182' and report['gene']=='BRCA1' and
                row.get('hgvsc')=='c.181T>G' and row['row']==388 and row.get('table')=='T9' and
                row['sha256']=='b023f69c4d3131d474d90badc08012f130fa9bbdb47a9ece88ea7e3581da6bc8' and pmid in ROLES_388):
                role,cells=ROLES_388[pmid]
                if pmid=='30209399':entry['supplement']=FINDLAY
            elif row.get('table')=='T4_splice':
                role='ST16 剪接／RNA 线索'
                entry['st16_cells']=[r['cell'] for r in row.get('rna_references',[]) if pmid in r['pmids']]
            entry.update(role=role,cells=cells)
            result.setdefault(pmid,[]).append(entry)
    return result
