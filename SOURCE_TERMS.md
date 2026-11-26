# 软件与数据来源

BRCA Review Desk 0.2.0rc8 汇集公开来源的必要教学字段、定位与校验指纹。不含患者输入、作者工作簿、全文 XML、整库镜像或真实运行目录。

| 来源 | 用途与分发方式 |
|---|---|
| [ClinGen ENIGMA BRCA1 GN092](https://cspec.genome.network/cspec/ui/svi/doc/GN092?version=1.2) / [BRCA2 GN097](https://cspec.genome.network/cspec/ui/svi/doc/GN097?version=1.2) | 锁定 v1.2。随包为本项目中文释义（非官方译本）、字段定位与校验信息。读者可从 [GN092 存档](https://zenodo.org/records/21434315) / [GN097 存档](https://zenodo.org/records/21434343) 获取官方 JSON，从规范发布页获取完整工作簿。 |
| [Parsons 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6772163/)、[Li](https://pubmed.ncbi.nlm.nih.gov/31853058/)、[UCSC ENIGMA](https://hgdownload.soe.ucsc.edu/hubs/enigma/enigma.html) | 逐判据链接和读表指引。原表、整条轨道均不随包；后续更新可能改变外站内容。 |
| [ClinGen CAR](https://reg.clinicalgenome.org/)、[LDH](https://ldh.clinicalgenome.org/)、[eRepo](https://erepo.clinicalgenome.org/evrepo/ui/)、[Ensembl REST](https://rest.ensembl.org/)、[GRCh37 REST](https://grch37.rest.ensembl.org/)、[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)、[gnomAD](https://gnomad.broadinstitute.org/) | 数据归各来源管理；公开可查不等于任意批量再分发。教学用已登记的公开等位；历史页只列 4 个最小 Ensembl 响应（主机版本与参考碱基）及指纹。其余实际响应仅在使用者主动取数后落入其本机目录。 |
| [PubMed](https://pubmed.ncbi.nlm.nih.gov/)、[Europe PMC](https://europepmc.org/)、[Kowalik 等 2018](https://doi.org/10.1371/journal.pone.0201086)（[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)） | 阅读顺序教学只分发章节、表名、核酸别名、定位和指纹，不附论文全文或原始 XML。其他文章、功能实验、预测工具及补表遵循各自许可。[SpliceAI Lookup](https://spliceailookup.broadinstitute.org/) 只作人工查证入口。 |

**代码许可。** 本候选包未附经权利人确认的开源代码许可证，不声明 MIT、Apache 或其他授权。拟公开仓库时由维护者确认自有源码的权属和拟授予的许可。第三方数据权利与代码许可分开；下载地址、内容哈希或注明来源都不能替代原许可。
