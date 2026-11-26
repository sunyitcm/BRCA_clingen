# BRCA Review Desk · BRCA1/2 胚系变异证据审核台

**0.2.0rc8 · 候选发布版**　｜　Python 3.10+ 标准库　｜　GRCh37 / GRCh38　｜　ClinGen ENIGMA VCEP v1.2

从一个核酸等位出发，在同一页核对 **等位 → 判据 → 观察值与原始出处 → 所选证据与组合**。适用对象是复核遗传性肿瘤易感变异的遗传学与临床生信人员。输出是有出处的研究级审核记录，不自动签发临床分类，也不按体细胞致癌性或用药指南计分。

## 1. 先看，或在本机查询

| 目的 | 入口 | 是否联网 |
|---|---|---|
| 无安装试读 | 解压**完整 ZIP**，双击 [START_HERE.html](START_HERE.html)。选「打开 BRCA2 教学位点的 PP4，手查 Parsons 2019 原表」，可直接到该行。 | 打开不联网；点作者 Excel、论文或 UCSC 时才访问外站。 |
| 查一个位点、取正文、保存勾选 | 在含 `start.py` 的目录运行下面的命令；在打开的浏览器选择「查询」。 | 查位点或取未缓存的正文时，Python 才请求对应公共服务。勾选、保存和计算在本机完成。 |

Mac / Linux：

```bash
python3 start.py
```

Windows PowerShell：`py -3 start.py`；已激活 Miniforge 的三个平台均可用 `python start.py`。无需安装 Perl、Docker、VEP cache 或填写 API key。终端保持打开，结束按 Ctrl+C；端口被占可用 `python3 start.py --port 0` 并打开它显示的网址。`127.0.0.1` 只指运行程序的这台电脑。GitHub Pages 可展示 `docs/` 下的只读教学页，不能代替本机保存或查询。

## 2. 明确组装与核酸等位

**在线查询提示（2026-11-26 04:01:14 UTC）**　本次从 Mac 请求 Ensembl GRCh37 VEP 的[公开测试位点](https://grch37.rest.ensembl.org/vep/human/hgvs/NC_000017.10:g.7577107A%3ET)，收到 HTTP 500 和 HTML 错误页。当前建议优先使用已核实的 **GRCh38 输入**。手头只有 GRCh37 坐标时，请保留原始组装，核实同一等位的 GRCh38 坐标与 REF 后再查询。

先选 BRCA1（chr17）或 BRCA2（chr13），再选原始组装。优先输入 **chr-POS-REF-ALT** 或单行 VCF；也支持带版本的 HGVS c.。每次一个 ALT。以下三种写法是同一公开教学等位：

| GRCh38 | GRCh37 | 规范转录本 |
|---|---|---|
| `chr13-32380106-G-T` | `chr13-32954243-G-T` | `NM_000059.4:c.9217G>T` |

本例应核至 BRCA2、CA10579818、p.Asp3073Tyr。不要只切换组装而保留旧坐标，也不要只凭蛋白名称认定是同一核酸等位。CAR、VEP、LDH、gnomAD 等分别取证，原组装与映射组装分开标记。GRCh37 官方 REST 有时会返回 500 或 HTML；这属于服务失败，不能写成未收录、AC=0 或原组装注释成功。详见[脚本与接口](METHODS_API_ZH.md#grch37)。

## 3. 判据与出处：先看 PP4，再回作者表

两基因各列冻结规范启用的 **19 条主判据**。已有线索的行在前，其余在「展开其余可用判据」。从首页走直达链接，不必展开全表，也不必先勾选「计入」。

| 判据行 | 点击哪里 | 打开后看哪里 |
|---|---|---|
| **PP4 / BP5** | **作者原表（手工查）→ 下载 Parsons 2019 原始 Excel**；旁边有论文入口。 | `SuppT1 All Data`：A 列核基因、B 列查核酸 c.、**K 列取综合 LR**。根据方向及阈值复核。 |
| **PP1 / BS4** | 同一作者原表 | 同一张表的 **F 列是共分离 LR**。若已包含在采用的综合值里，不能再次计分。 |
| **PP4 / BP5** | **UCSC 本位点：逐条展开**，另有方法与修订说明 | 打开原组装的 ENIGMA 组与 LR 子轨道（`full`）。点核酸 ALT 相同的条目看来源分量。身份或位置不唯一时不给位点链接；仍可手查原表。 |

[Parsons 2019 原始 Excel](https://pmc-oa-opendata.s3.amazonaws.com/PMC6772163.1/HUMU-40-1557-s001.xlsx) · [原论文／补充资料，PMID 31131967](https://pmc.ncbi.nlm.nih.gov/articles/PMC6772163/) · [UCSC 轨道方法及更新](https://hgdownload.soe.ucsc.edu/hubs/enigma/enigma.html)。

**两个不同的读表示例：** BRCA2 `c.3G>A` 在 Parsons 第一张表第 745 行：F≈67.62、K≈157.49；这是分量与综合值。BRCA2 `c.9217G>T` 的教学例中，3.57 来自 [Li 的 BRCA2 补表 4](https://pmc-oa-opendata.s3.amazonaws.com/PMC7118020.1/41436_2019_729_MOESM3_ESM.xls)的家族史数据，不应归给 Parsons 或 PP1。下载链接是查证入口，不代表程序已检索该表、找到了本等位或自动采用了分数。作者表原用 NM_007294.3 / NM_000059.3，检索时核回同一核酸等位。

其他常用行：功能实验证据 **PS3/BS3** 先查规范 Table 9 及采用／排除理由；**PVS1/PM5_PTC** 核转录本、机制、Table 4 与 ST1；**PM2** 核指定人群、覆盖与事件适用范围。点每行「逐档条文、英文与出处」核对条文和工作表。冻结规范：[BRCA1 GN092](https://cspec.genome.network/cspec/ui/svi/doc/GN092?version=1.2) · [BRCA2 GN097](https://cspec.genome.network/cspec/ui/svi/doc/GN097?version=1.2)。

## 4. 读原文，再决定是否计入

「全文与人工阅读」按工作优先级提示先看哪一表行或 Results。可从 [BRCA1 离线试读](docs/snapshot.html?gene=BRCA1&assembly=GRCh38&snapshot=demo-public-rc8-brca1-38&tab=papers) 看 PMID 30040829 的 Table 1 和 Results 为什么属于同一队列。在线版点「取正文（优先缓存）」；最多重试一次或转 PubMed、PMC、期刊和补表手查。阅读优先级、论文篇数和证据强度是三件事。

在本机「查询」页勾选有依据的判据，按需选择强度，可留一条选填备注，然后点 **校验、保存并计算**。附件／本地规则、引用专家、人工补入的贡献分开列。固定 BRCA2 教学习题：先保存预选的 PS3、PP3 得 +5，再明确勾 PM2 和 PP4 的专家引用得 +7（其中专家引用 +2）。今天联网查回的字段不保证与教学题相同。

## 5. 回看版本、审核记录和校验值

「各库版本与差异」分列组装、转录本、来源版本和队列／实体版本；「原始数据」可核查请求、响应与 SHA-256。查询 UTC 不替代数据库版本。分发的公开示例不含制作者的操作时间，但保留原数据源发布的日期；你自己运行产生的 UTC、原始字节和审核修订仍存于本机。前期试用版的教学审核继续用原副本回看；公开教学摘录采用新 ID 与新哈希，不会套用旧勾选。真实位点的既有快照不随包分发，也不因解压被改写。

| 进一步阅读 | 内容 |
|---|---|
| [解读员手册](docs/Reviewer_Guide_ZH.md) | 各判据如何显示、勾选和保存；阴性、不适用与未评价。 |
| [临床 LR 与文献](docs/BRCA_Literature_and_Clinical_LR_ZH.md) | 原表、UCSC、新研究之间如何辨别独立性；作者表逐列及检索实例。 |
| [脚本、API 与 Python 判读](METHODS_API_ZH.md) | 各来源实际 GET/POST、curl 例子、失败和缓存的处理。 |
| [本版说明](RELEASE_NOTES.md) · [验收与复现](docs/Acceptance_ZH.md) · [数据来源与权利](SOURCE_TERMS.md) | 安装包完整性、离线验证、公开范围及第三方来源。 |
