# BRCA 审核台：脚本、API 与 Python 判读方法

供维护者、临床生信人员和希望核查某次取数的解读员使用。[返回解读员 README](README.md)。描述的是 **BRCA Review Desk 0.2.0rc8 / BRCA-lite 0.2.1rc1 冻结内核**。本版在 PP4/BP5 生成逐条展开的 UCSC 链接，并在 PP4/BP5（K 列）、PP1/BS4（F 列）直链 Parsons 2019 的工作簿；诊疗证据规则与上版相同；UCSC 和作者原表用于人工对照。所有例子均为公开等位。

## 一、GET、POST、curl、Python 和 HTML 各负责什么

**GET / POST 是 HTTP 请求方法；REST 是接口设计方式；curl 是可以发送两种请求的命令行客户端。** 应用启动后由 Python 标准库 `urllib.request` 发请求，不会启动 curl、Perl 或本地 VEP。浏览器的“联网查询”先把当前等位 POST 给自己的 Python 服务；Python 再按来源选择 GET、POST 或读取本机文件。HTML 是显示层，JSON、XML、FASTA、XLSX 是内容格式，都不能只靠文件扩展名判断成功。

| 源／步骤（实际入口） | 方法、输入、输出 | 必须核的字段／失败时的含义 |
|---|---|---|
| Ensembl 主机版本、坐标 REF | GET `/info/data`、`/sequence/region/human/{region}` → JSON。可配置锁定的本地 FASTA 清单；远程 REF 不可用可走 NCBI EFetch GET（FASTA）。 | 主机 `releases`；参考组装、区间、正链、实际 REF。只有匹配才通过 VCF 的 REF 核对。 |
| VEP：VCF 型输入 | **POST** `/vep/human/region`，JSON `variants` 数组（本界面每次一 ALT）；见下面原样示范。 | 核 `assembly_name`、`input`、唯一规范转录本、`hgvsc`；`hgvsp` 和外显子只在可得时显示。未返回该转录本与 HTTP 错误分别登记。 |
| VEP：版本化 HGVSc 输入 | **GET** `/vep/human/hgvs/{经URL编码的NM:c.}` → JSON。两路径均启用 `hgvs=1&merged=1&transcript_version=1&numbers=1&variant_class=1`，GRCh38 另启用 `mane=1`。 | 与规范 NM 及 CAR 的 c./p. 对账。当前默认不请求 REVEL、CADD、SpliceAI 等插件，不能以 VEP 返回 200 证明这些分数存在。 |
| ClinGen Allele Registry（CAR） | **GET** `https://reg.clinicalgenome.org/allele?hgvs=...`；VCF 先构造版本化 NC:g.；输入 c. 时查该 c.。 | 取 `@id` 的 CAid、`genomicAlleles`、`transcriptAlleles` 和外部 ID。用基因、规范转录本和等位核准；重复区需等价单倍型核对。只查询，不创建新等位。 |
| ClinGen LDH；必要时 BRCA Exchange | **GET** `https://ldh.clinicalgenome.org/ldh/Variant/id/{CAid}?detail=med`；再 GET 目录中该实体的受信 URL。若有核准记录 ID 及发布信息，再 GET BRCA Exchange 的 `/backend/data/variant/?variant_id=...`。 | 核 CA、记录类型、实体版本、GRCh38 表示或记录 release。目录本身不等于分数。用到的实体包括 gnomAD v4.1 外显子／基因组、BRCA Exchange、分子后果。“没列链接”、字段占位、明确为 0 分开保存。LDH 不需要先有 VEP 蛋白名。 |
| gnomAD Browser | **POST GraphQL** 到 `https://gnomad.broadinstitute.org/api`；提交 `variantId` 和明确的 `dataset`。本机只读，不是 REST GET。 | 全队列 `gnomad_r2_1`（GRCh37）或 `gnomad_r4`（GRCh38）；历史非癌外显子 `gnomad_r2_1_non_cancer` 和非癌基因组 `gnomad_r3_non_cancer`。核组装、等位、测序类型、AC/AN、过滤和适用人群；覆盖/FAF 单独核。GraphQL 的 HTTP 200 仍可同时有 `errors`；`data.v=null` 不是 AC=0。 |
| ClinVar；ClinGen eRepo | ClinVar：CAR 有唯一 Variation ID 才 GET NCBI `efetch.fcgi?db=clinvar&rettype=vcv&is_variationid=true&retmode=xml&id=...`。eRepo：GET `/evrepo/api/summary/classifications`（CA 精确匹配），再 GET `/evrepo/api/classification/{uuid}`。 | ClinVar 核 `VariationArchive`、VCV/SCV 版本、胚系断言、引用的 PMID。eRepo 再核基因、MONDO、显性遗传、ENIGMA VCEP、未撤回及逐码详情；个案匹配才展示专家补充。聚合星级、题录、其他疾病的断言都不自动计分。 |
| CSpec 与三份 Excel | 未锁版本的当前发布页 **GET HTML** 只作版本对照；实际代码／组合使用冻结文件。官方工作簿按固定 URL GET 或读用户提供的同版文件，初次核 SHA-256，再以标准库 `zipfile` + XML 建索引。 | 核基因、规范版、文件哈希、工作表、物理行和字段；不能验证时不按未知新表授分。Table 4／9、ST1／16 是现版提取的部分；作者 Li 的旧式二进制 `.xls` 不是这套 XLSX 解析器的输入。 |
| 文献与链接 | 单篇按钮：Europe PMC `search` **GET** 以 `SRC:MED AND EXT_ID:PMID` 的等价表达取题录，再按 PMCID **GET** `/{PMCID}/fullTextXML`。另有 PubMed／PMC 预填网页链接；底层批量模块可选 ESearch／正文发现，但按钮不会自动把整套候选论文都抓一遍。 | 核 PMID、PMCID、已知 DOI、可读正文、原 XML 哈希。题录、仅摘要、无开放全文、XML 成功但无别名命中、只提供补表链接、服务失败各有不同含义。 |
| UCSC `BRCAmla`、Li／Parsons 作者原表 | **目前为人工或独立核查**。本版依据已核准的 CA、原组装和唯一 NC:g. 为 PP4/BP5 构造 Genome Browser 链接，ENIGMA 父组 show、BRCAmla 子轨道 full；不能可靠定位就不生成该位点链接，仍可手查作者原表和方法。另可 GET UCSC `/getData/track` JSON。判据行链接到 PMC Open Data 中的作者 Excel；点击原表只发文件下载 GET，不在工作簿 URL 附加当前位点，也不导入评分。作者表为下载文件；**未写入本版常规联网流水线**。 | 同一查询区间可出现多个等位；作者稿与汇总轨道可能复用队列。记录有效数字和原始尺度，区分“精确收录”和“可用于本判据”。链接能打开不证明程序已读过它。 |

## 二、从输入到已保存的判据：按代码核对

源码在 `source/brca_review/`。文件名按发行包给出，前缀 `_engine/brcalite/` 指内置冻结引擎。

| 顺序 | 对应脚本 | 可验证的行为 |
|---|---|---|
| 1. 输入、选择基因和组装 | `start.py`、`desk.py`、`criteria.py`、`_engine/brcalite/core.py` | 一次一 ALT；VCF 行只取 CHROM/POS/REF/ALT，别把整份样本表粘入。选基因、染色体、规范 NM 和组装须一致。不会只凭 p. 反推一个核酸等位。 |
| 2. 独立取数和对账 | `routes.py`、`_engine/brcalite/net.py`、`erepo.py`、`exchange.py` | 原组装 VEP、CAR、REF、版本可同时启动，有本组装 VCF 坐标时 gnomAD 可直接查。LDH 在核准 CA 后查询；ClinVar 看唯一 Variation ID。一个上游失败，其他来源的成功和失败均留档。跨组装资料单独标来源，不能依赖坐标相近。 |
| 3. 事实整理及逐码求值 | `ruleprep.py`、`mechanism.py`、`observations.py`、`briefing.py`、`lineage.py`、`criteria.py` | 从机制、转录本、冻结工作簿和保留的分量生成 `criteria_table.rows`。例如：Table 9 的明确推荐与 PMIDs 分开；BayesDel noAF 走有机制和域限制的规则；PP1 与 PP4 读不同 LR。仅有氨基酸位置不能通过所有条件；未知、阈值未达和本规则不适用不同。 |
| 4. 保存与组合 | `selection.py`、`adoption.py`、`_engine/brcalite/score.py` | 本机 POST `/api/review` 校验快照哈希、身份、允许强度、重复及互斥。Supporting/Moderate/Strong/Very Strong 的分值绝对值分别为 1、2、4、8，致病方向为正、良性方向为负；BA1 独立显示。组合提示按冻结的**强度组合**评估，不能只凭总分转换分类。来源分列；混合方向或可能重叠可暂缓合计，绝不把所有候选都默认计入。 |
| 5. 展示与复读 | `model.py`、`presentation.py`、`papers.py`、`_engine/brcalite/{paper_queue,literature}.py`、`assets/{app,criteria}.js` | 归档 `report.json` 和快照；审核另存修订。浏览器渲染同一报告的逐码行。“先读哪里”从本阅读事件已有的来源、工作表行、章节和命中别名生成；没有调用大模型，也不写入评分。 |

**谁决定什么？** 已锁定表行明确给出同等位代码/强度、或本地规则所需条件完整，软件可预选；只有专家已发布个案时，可明确选择“引用专家”；新研究的队列独立性、实验校准、家系定量输入和相互矛盾的资料，需要人工核。算法只拦截已编码的重复／互斥；没有报警不是已证实样本独立。未确认项可以留在待核清单。

## 三、手工复核：可复制的 Mac/Linux curl

以下命令在独立工作目录运行，只发公开教学等位。Windows 用户可直接用本机网页；以下是 Bash 写法，不要把反斜杠续行原样粘到 PowerShell。命令模拟核心请求，不替代应用的等位对账、归档和评分。

**也可以运行随包示范。** 在项目根目录执行 `bash scripts/demo_vep_curl.sh`（Mac/Linux，可选）。它只查询下面这个公开 BRCA2 等位，在系统的临时目录新建独立结果目录，打印路径并调用 [check_public_vep.py](scripts/check_public_vep.py) 验证回显、组装和转录本。已有两份响应时，可在响应目录用 Python 单独运行校验脚本，不会联网。正式审核台仍通过 `python3 start.py` 启动，不依赖这两个脚本。

**VCF 坐标 → GRCh38 VEP（POST）；版本化 c. → 同一服务（GET）。** 两种操作都是只读。

```bash
curl -sS --max-time 30 -D vep-vcf.headers -o vep-vcf.json -w 'VCF: HTTP %{http_code}, %{content_type}\n' \
  -H 'Accept: application/json' -H 'Content-Type: application/json' \
  --data-binary '{"variants":["13 32380106 . G T . . ."]}' \
  'https://rest.ensembl.org/vep/human/region?hgvs=1&merged=1&transcript_version=1&numbers=1&mane=1&variant_class=1'
curl -sS --max-time 30 -D vep-c.headers -o vep-c.json -w 'HGVSc: HTTP %{http_code}, %{content_type}\n' \
  -H 'Accept: application/json' \
  'https://rest.ensembl.org/vep/human/hgvs/NM_000059.4%3Ac.9217G%3ET?hgvs=1&merged=1&transcript_version=1&numbers=1&mane=1&variant_class=1'
```

先看**各自的 HTTP 状态和类型**：若为 500/HTML，文件仍可能名为 `.json`，将其当作错误原文；不要当作空结果。若为 200 再核内容（不能只看第一条转录本）。下面这个小脚本仅用于此教学输入：

```python
# 保存为 check_public_vep.py；运行 python3 check_public_vep.py
import json
from pathlib import Path

for name, original in [('vep-vcf.json', '13 32380106 . G T . . .'),
                       ('vep-c.json', 'NM_000059.4:c.9217G>T')]:
    try:
        items = json.loads(Path(name).read_bytes())
    except (OSError, ValueError) as exc:
        raise SystemExit(f'{name}: 文件不可读或不是 JSON，请看 HTTP 状态和响应头') from exc
    if not (isinstance(items, list) and len(items) == 1 and isinstance(items[0], dict)):
        raise SystemExit(f'{name}: 未得到唯一输入结果')
    rec = items[0]
    if (rec.get('assembly_name'), rec.get('input')) != ('GRCh38', original):
        raise SystemExit(f'{name}: 组装或输入回显不符')
    hits = [t for t in (rec.get('transcript_consequences') or []) if isinstance(t, dict)
            and t.get('gene_symbol') == 'BRCA2' and t.get('transcript_id') == 'NM_000059.4'
            and t.get('hgvsc') == 'NM_000059.4:c.9217G>T']
    if len(hits) != 1:
        raise SystemExit(f'{name}: 没有唯一匹配的规范转录本')
    print(name, hits[0]['hgvsc'], hits[0].get('hgvsp') or '未返回蛋白描述')
```

**CAR 和 LDH 用 GET，可直接在浏览器看 JSON。** [CAR：按版本化 c. 查询](https://reg.clinicalgenome.org/allele?hgvs=NM_000059.4%3Ac.9217G%3ET) → 核到 CA10579818 → [LDH：按该 CA 取实体目录](https://ldh.clinicalgenome.org/ldh/Variant/id/CA10579818?detail=med)。不要把页面上第一条 gnomAD 或分子后果实体直接当作目标，还要核实体 ID、来源、组装及版本。同行若需重放本机一次请求，优先复制那次回执的 URL 和 POST 请求体。

**UCSC 手工核查：** [GRCh38 本位点 API 原始 JSON](https://api.genome.ucsc.edu/getData/track?hubUrl=https%3A%2F%2Fhgdownload.soe.ucsc.edu%2Fhubs%2Fenigma%2Fhub.txt&genome=hg38&track=BRCAmla&chrom=chr13&start=32380105&end=32380106)。这里的 API 区间是 0 起算、右端不含；Browser 显示和 VCF POS 是 1 起算。此区间实测返回两条：`c.9219_9228del`（1.1329）与 `c.9217G>T`（3.5701）。匹配名称／规范化核酸 ALT 后再取数；不要按返回顺序或仅按 p. 选取。UCSC 的 `LLR` 字段在**这条轨道中是线性合并 LR**；与 Li 作者表的 `LOG(LR)` 不同。

<a id="grch37"></a>
## 四、两组装的实测边界

以下是同一公开 BRCA2 等位在一次历史探测中的响应；用于说明不同返回状态，不代表各服务当前可用率。公开版不附探测时钟或现场日志：

| 实测路径 | HTTP / 实际内容 | 使用决定 |
|---|---|---|
| GRCh38 `rest.ensembl.org`：VCF POST、HGVSc GET、`info/data` GET | 各 200 JSON；两种输入均核得 `NM_000059.4:c.9217G>T`、`NP_000050.3:p.Asp3073Tyr`；主机 release 116 | 可对该等位复核。服务将来可更新；这不是所有来源冻结在 116。 |
| GRCh37 `grch37.rest.ensembl.org`：原生 VCF POST、`info/data` GET | 两个请求均为 500、text/html | 记服务错误。保留 GRCh37 入口，保留为原组装；勿写“无变异”或“端口已更改”。 |
| 单独探测 `grch37-archive.rest.ensembl.org`，同一 GRCh37 VCF | 200 JSON；归档 release 99，六条转录本后果中包括 `NM_000059.3:c.9217G>T`／`NP_000050.2:p.Asp3073Tyr`，没有所需的 NM .4 | 可当有明确版本的历史线索。不是程序默认路由，不能静默替换或把 .3 改写成 .4。 |
| CAR 的 HGVSc GET、LDH 目录 GET、UCSC 区间 GET | 三项均为 200 JSON；前两项对上 CA10579818；UCSC 精确匹配到线性合并 LR 3.5701 | 对上述查询字段有效；不能据此推断 gnomAD、Europe PMC 或其他位点均正常。 |

同样的 500 出现在页面上时，先区分**外站回执**与**本机保存失败**（例如磁盘写入错误），不要改组装来消除报错。对原生 GRCh37，如 CAR 唯一对账且能构造 GRCh38 表示而原 VEP 缺指定转录本，可独立取 GRCh38 VEP 并标 `mapped_vep`；这不是原生 37 成功，也不是静默 liftover。原始 VCF 的 REF 无法核准或规范 c./p. 冲突时，保留各源但不计分。没有唯一注册 CA 时，LDH 不能按此路径取得；不等于 VEP 没有注释或等位是新发现。

## 五、LR 和文献：复核要留下什么

- **哪里找到，哪里没找到。** 分开记 CSpec 版本、作者／PMID、DOI、勘误、原文件名／SHA-256、工作表名／行列、HGVS 或组装与 REF/ALT、是否精确匹配、字段名／原始尺度、取数 UTC。UCSC 的两个组装是同一资源的投影；“全表未见”只可用于已完整核过的版本，附件下载失败记“未核实”。
- **不要比较错列。** Parsons 原表 `SuppT1 All Data` 的 F=共分离、I=家族史、K=综合 LR、L/M=先验／后验概率；`SuppT2 Data Sources` 追来源。旧 NM 版本须用同一核酸等位消歧。Li 表 3/4 的 D=LOG(LR)，E=先证者数；对该来源按 `10**D` 还原，保留原值。未知底数不猜。中性区间 0.48 < LR < 2.08、没有值、作者明确排除、另一等位有值，是四种状态。
- **怎样对旧表、轨道、新文章作决定。** 先核同一等位，再核每个分量的研究和原始人群；标“沿用、替换、已含、独立、待核”。新发布不等于新队列；历史值和新值都保留。有实质差别或跨阈值时不选最大值、不自行相乘。组合已包含的共分离或病例–对照证据不能再按 PP1/BS4 或 PS4 重计。表格用于校准阈值时，只作方法依据。详细实例和直达链接见 README。
- **先读哪里。** 后端从开放 JATS 的正文段和表行用当前核酸／蛋白别名查找，记录章节、表头（可得时）、单元定位、命中别名及哈希。图片/PDF 和普通外链补表**不在**这次 XML 扫描范围。只有许可符合项目条件时才内嵌原段；否则给定位与原站链接。更新后的阅读卡优先已知补表行，或当前 c. 命中的表行／Results；只有蛋白别名时明确提示核回 ALT。PMID 30040829 的队列说明是已核实的具体示例，不能把其他论文的段落自动去重成独立患者。
- **缓存、限流与可及性。** `papers.py` 先按 PMID 找可复用缓存，核原文 SHA、PMID/PMCID、已知 DOI，再用当前等位别名重新扫描。点击“重新联网”失败而有可校验旧缓存时，明确标注旧来源 UTC。位点数据/文献的公共 HTTP 读请求通常至多两次（首发 + 一次临时故障重试）；429/5xx/传输错误按原因处理，有 Retry-After 先尊重，过长不在现场等待。CSpec 发布页探测只尝试一次；工作簿下载是独立路径。免登录不等于无限额、稳定 SLA 或取得第三方全文再分发权。

## 六、如何复验、交接和升级

在含 `start.py` 的目录运行：

```bash
python3 verify_files.py
python3 start.py --self-test
python3 -m unittest discover -s tests
```

`--self-test` 是断网教学自检，不是测试外站是否在线。首次联网查询前可单独运行 `python3 start.py --prepare-rules`（只准备／核验三份工作簿，不发位点）；已有指定文件可用 `--workbooks` 传文件夹。可用 `--data-dir` 指定仓库以外的运行目录。`--ref37`、`--ref38` 接受**对应 FASTA 的项目清单 JSON**，不是直接传一个 .fa 文件；原生组装和参考文件必须相配。`--population-file` 为本地、可核准的双队列资料，不是上传人群数据的接口。

现场从“原始数据”核请求状态与链接；本机 `s-…/http/<来源>/receipt.json` 记实际方法、URL、尝试、响应头和原文 SHA，POST 另有 `request.json`；失败可能有 `failed-1.body`。原查询保留 `report.json`、`snapshot.json`；文献在 `reading_events/` 留存，与评分修订分开。不要把运行目录提交到 GitHub。轮到升级时，先核新规范、作者勘误／轨道更新日期及是否许可分发，再用命中、错 ALT、无信息、失败、重复证据等公开对照重放；保留原版快照。项目自检不是临床性能或跨人群有效性验证。

**官方接口与方法入口：** [VEP POST](https://rest.ensembl.org/documentation/info/vep_region_post) · [VEP HGVS GET](https://rest.ensembl.org/documentation/info/vep_hgvs_get) · [GRCh37 归档接口](https://grch37-archive.rest.ensembl.org/) · [CAR](https://reg.clinicalgenome.org/) · [LDH](https://ldh.clinicalgenome.org/) · [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25497/) · [Europe PMC API](https://europepmc.org/RestfulWebService) · [UCSC 轨道方法／修订说明](https://hgdownload.soe.ucsc.edu/hubs/enigma/enigma.html) · [GN092](https://cspec.genome.network/cspec/ui/svi/doc/GN092?version=1.2) · [GN097](https://cspec.genome.network/cspec/ui/svi/doc/GN097?version=1.2)。

专项阅读：[先读哪里、临床 LR 与作者原始表](docs/BRCA_Literature_and_Clinical_LR_ZH.md)。
