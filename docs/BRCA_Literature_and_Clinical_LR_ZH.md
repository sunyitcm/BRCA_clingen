# BRCA1/2 临床 LR 与文献｜从判据行查原表

适用规范：[BRCA1 GN092 v1.2](https://cspec.genome.network/cspec/ui/svi/doc/GN092?version=1.2) 与 [BRCA2 GN097 v1.2](https://cspec.genome.network/cspec/ui/svi/doc/GN097?version=1.2)。本章随 0.2.0rc8 分发，区分原始研究、汇总轨道、专家引用和本地实测；不产生临床分类。返回 [首页](../README.md) 或 [API／脚本说明](../METHODS_API_ZH.md)。

## 1. 从 PP4/BP5 或 PP1/BS4 判据行，查工作簿与轨道

**从判据进原表：** “判据与出处” → PP4/BP5（综合 LR）或 PP1/BS4（仅共分离）→ “下载 Parsons 2019 原始 Excel”。未在优先展示区时，展开其余判据即可。工作簿保持作者原样，由使用者下载；下载后可断网用 Excel 查找。打开 **SuppT1 All Data**，在 B 列查 c. 并核 A 列基因，再分别读 K 或 F 列。用第二张表追数据来源。仅打开链接不会标记命中或改变本地分数。

**打开实例：** [UCSC Genome Browser 中的 BRCA2 c.9217G>T（GRCh38）](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr13%3A32380096-32380116&hubUrl=https%3A%2F%2Fhgdownload.soe.ucsc.edu%2Fhubs%2Fenigma%2Fhub.txt&hub_5570236_enigma=show&hub_5570236_BRCAmla=full)。默认请求只将指定 LR 子轨道设为 full（逐条显示），并打开 ENIGMA 父组；要查看分量仍需点击同一核酸等位的条目。另见 [UCSC ENIGMA 汇总轨道及更新说明](https://hgdownload.soe.ucsc.edu/hubs/enigma/enigma.html)。轨道名为 `BRCAmla`，文件名为 `BRCAmfa.bb`。[hg38](https://hgdownload.soe.ucsc.edu/hubs/enigma/hg38/BRCAmfa.bb) 和 [hg19](https://hgdownload.soe.ucsc.edu/hubs/enigma/hg19/BRCAmfa.bb) 表示同一批数据的不同组装，不能相加。它把多篇研究整理为来源分量和综合 LR，是派生汇总，不是新的独立队列，也不是规范 Table 7。Table 7 是 **PP4/BP5 使用示例**。

| 需要的资料 | 从哪里读 | 怎样用 |
|---|---|---|
| 多因素综合、或单独共分离分量 | [Parsons 2019，PMID 31131967](https://pmc.ncbi.nlm.nih.gov/articles/PMC6772163/)；[作者原始 Excel](https://pmc-oa-opendata.s3.amazonaws.com/PMC6772163.1/HUMU-40-1557-s001.xlsx) | 涵盖 1,395 个 BRCA1/2 变异，主要是错义/内含子，734 个有可进行后验分类的数据。`SuppT1 All Data` F=共分离，G=病理，H=共现，I=家族史，J=病例对照，K=综合 LR；L/M 是先验/后验概率。再用 `SuppT2 Data Sources` 找各研究中心、家系/肿瘤数。原表参考转录本为 NM_007294.3 / NM_000059.3，不能盲拼版本或蛋白名。 |
| 个人及家族肿瘤史 LR | [Li 2020，PMID 31853058](https://pubmed.ncbi.nlm.nih.gov/31853058/)；[BRCA1 补表 3](https://pmc-oa-opendata.s3.amazonaws.com/PMC7118020.1/41436_2019_729_MOESM2_ESM.xls)、[BRCA2 补表 4](https://pmc-oa-opendata.s3.amazonaws.com/PMC7118020.1/41436_2019_729_MOESM3_ESM.xls) | 文章报告对 4,644 名携带 2,383 种 VUS/VLP 的受检先证者进行分析。读位点行的 LOG(LR) 和 N_Probands；这是经校准的病史证据，不是亲属间共分离。 |
| 病例–对照成分 | [Zanti 2025，PMID 40413188](https://pubmed.ncbi.nlm.nih.gov/40413188/)，尤其是作者的 Supplementary Data 4 | 来自 BRIDGES、CARRIERS、UK Biobank。ccLR、队列 OR、家系 LR 不是同一量。轨道在 2026 年 11 月修订：作者标为 `No evidence*` 的小携带者数 ccLR 仍可见，但不再乘入综合 LR。不要把同一病例-对照数据再作为 PS4 与综合 LR 累加。 |
| 新的家系共分离 | ClinGen 条文附录 I；[CAL-Leiden 方法，PMID 41882377](https://pubmed.ncbi.nlm.nih.gov/41882377/) | 找不到发表数据、自己又有信息足够的家系时，可用亲缘关系、可用的携带/非携带状态、癌种、发病年龄和末次随访年龄，按外显率与招募条件作定量分析。没有这些资料则留空；单个肿瘤基因检测或“多位亲属患癌”不能代替。 |

## 2. 规范怎样读 LR？

LR 比较“该变异致病”与“该变异中性”两种假设下观察到这些临床数据的可能性。**PP1/BS4** 评价定量共分离；**PP4/BP5** 评价经临床校准、来源和组成可核的综合 LR。综合值可以只有一个有信息的分量。指南也允许满足条件的未发表数据；并不要求各实验室重复建队列。发表即给分、按文献数量加分、把 LOG(LR)/后验概率当线性 LR 都不成立。

| 方向 | Supporting | Moderate | Strong | Very Strong |
|---|---:|---:|---:|---:|
| 支持致病（PP1 或 PP4，按对应 LR 类型） | ≥2.08 | ≥4.3 | ≥18.7 | ≥350 |
| 支持良性（BS4 或 BP5，按对应 LR 类型） | ≤0.48 | ≤0.23 | ≤0.05 | ≤0.00285 |

取满足条件的最高档，不累加各档。0.48 < LR < 2.08 不支持任一方向；这与“没取到 LR”不同。对 PP1_VeryStrong，规范另要求预测或实验证据支持对蛋白或剪接的影响；若结果来自单一家系或隔离人群，要排查其他病因。共分离若已含在综合 LR 中，不能再以同一批成员计 PP1/BS4。数据有重叠、方向冲突或同一变异只在外部专家摘要中出现时，保留来源逐项核查，不挑最大值。

## 3. 两个公开例子：看哪一行，能回答什么？

| 等位 | 精确定位与观察 | 可用结论 |
|---|---|---|
| BRCA2 NM_000059.4:c.9217G>T | Li [补表 4](https://pmc-oa-opendata.s3.amazonaws.com/PMC7118020.1/41436_2019_729_MOESM3_ESM.xls)，Sheet1 行 2009：D 列 LOG(LR)=0.5526850819587708，E 列 N_Probands=1。`10**0.5526850819587708=3.57013864…`。UCSC 行的 `liLRs`、`familyHistoryCombinedLR` 和 `LLR` 都是 3.5701，`segregationCombinedLR` 为空。[eRepo 专家个案](https://erepo.clinicalgenome.org/evrepo/ui/classification/237f2389-fa0c-4474-9f7a-50a8f9a3542d) 写为 3.57。 | 与专家的 PP4_Supporting 一致。这几个数字是同一条病史证据的不同表示，不是三次支持，不提供 PP1。UCSC 的 `LLR` 虽带两个 L，文件定义是**线性综合 LR**。 |
| BRCA2 c.3G>A | Parsons 原始表 `SuppT1 All Data` 行 745：F745 共分离 67.6211；K745 综合 157.4935，包含了前者。第二张表可追来源家系。 | 这是确有已发表共分离数值的例子，独立复核时不用自己重建作者队列。若采用综合值，不能再把 67.62 作为一条独立证据。不是对这个变异作新的临床分类。 |

**本次字段覆盖清点（一个组装为分母）：**13,182 条轨道记录，5,873 条的综合 LR 栏非空，564 条的共分离栏非空。非空不等于已满足某档强度。这些数字只描述下列 SHA-256 对应的轨道版本，不是所有可能 BRCA 变异的覆盖率。

**原始查询，可在浏览器打开：**[BRCA2 示例所在 GRCh38 小区间](https://api.genome.ucsc.edu/getData/track?hubUrl=https%3A%2F%2Fhgdownload.soe.ucsc.edu%2Fhubs%2Fenigma%2Fhub.txt&genome=hg38&track=BRCAmla&chrom=chr13&start=32380104&end=32380108)。这是区间 GET，返回第一条不一定是目标；必须核组装、核酸等位和转录本。规范版号与轨道更新独立，下游至少保存来源、表行/字段、组装、返回尺度、UTC、文件 SHA-256，以及哪些分量已经进入综合值。

核查快照（SHA-256）：

- hg38 BRCAmfa.bb：`833ca06d1dce0a0e412ee2e7fc3f345caf4037760caa310923cf58d2bb63182c`
- hg19 BRCAmfa.bb：`f8120dd7b06e69c88a059990df2d432abf6867cfdf625c929be40228c28baf56`
- Li BRCA2 补表 4：`0846fd4b5db3e900d40ff561e447ccc49c76654bc2b81af8fa9a187f5ec8b1c3`
- Parsons 补表工作簿：`76d25b0e1934d46ef5818158c61c8646c07ec14add50c31b71e0bc96505da20a`

## 4. 文献卡为什么建议先读 Table 1 或 Results？

公开教学例中，BRCA1 c.181T>G 的 [Kowalik，PMID 30040829](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0201086) 确实取得开放全文并在正文／表行中找到 6 处。优先读 Table 1 对应该变异的一行，然后读 Results 说明计数和招募方法。两处涉及同一批携带者，不产生两条独立的证据。方法与引言也可能命中名称，但不是该等位的新观察。查到题录仍可能无法下载正文，此时可使用卡片上的原文链接。

本版本机：“全文与人工阅读” → 目标 PMID → “取正文（优先缓存）”。打开阅读卡后表行排在前面，同类再优先 Results，有时能跳到原表行；只有快照里确实有定位，才显示命中。具体教学位点的同队列提醒是已核实写入的提示；通用扫描不会替每篇文献判断家系或队列是否独立。本版静态教学页显示公开示例定位，不会自动取全文。

**阅读卡的用法：** 先看首行的工作优先级，再按下列原则复核。

| 首屏提示 | 在什么情况下显示 | 人工先核哪里 |
|---|---|---|
| 优先读：有可对账的位点数据 | 同一核酸 ALT；指南/专家实际采用，或研究的原始表行可定位 | 展示证据类型、PMID、表/工作表及行号、观察值、单位或尺度，以及采用/排除理由。 |
| 值得复核：命中原文 | 同一等位出现在 Results、家系图、表格或实验；结论或独立性未核 | 看成员与基因型、分母、对照、原文单位；标注是否只是另一个结果的转述。 |
| 线索：仅引用、题录或未取到全文 | 只有 ClinVar 引用、摘要、基因/蛋白别名、关键词，或全文超时 | 标注为什么尚不能评价；给 PubMed、PMC、期刊和补表入口。一次重试后可改走手工。 |
| 已排除 | 专家或原作者明确排除该项实验或旧数据被更正 | 保留原文理由和更正来源；不要因为引文多或有“功能实验”字样重新计入。 |

读 LR 优先看等位级补表、独立先证者/家系、计算模型和数据来源；读 PS3/BS3 优先看指南 Table 9 及相应的原始实验/对照。两种任务的“必读论文”未必相同。
