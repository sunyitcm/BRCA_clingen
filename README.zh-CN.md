# BRCA_clingen：BRCA1/BRCA2 胚系变异注释与证据判读方法

[English](README.en.md) · [双语首页](.github/README.md) · [rc8 随包解读员说明](README.md)

这是一套独立的 Python 工具，供分子遗传学家和临床生信人员按 **ClinGen ENIGMA BRCA1/2 VCEP** 的基因专属 ACMG/AMP 规范复核遗传性肿瘤易感变异。Ensembl VEP 回答“这个等位在指定转录本上造成什么后果”；ClinGen 规范回答“哪些观察能用作哪一档证据”。报告把这两个问题分开，也让读者能回到原始来源。程序不签发临床结论，不作体细胞用药或致癌性分级。

## 用一个公开位点试读

1. **先打开离线示例。** 在 GitHub 仓库的 **Releases → v0.2.0-rc.8 → Assets** 下载 `brca-review-desk-0.2.0-rc.8.zip` 和 `SHA256SUMS.txt`。完整解压产品 ZIP 后，双击 `START_HERE.html`。直接在 GitHub 网页点 HTML 文件不会运行本机审核台。若维护者尚未公开 Release，首页的方法说明仍可独立阅读。
2. **再做一次本机查询。** 在含 `start.py` 的目录运行 `python3 start.py`（Windows：`py -3 start.py`；Python 3.10+）。浏览器选 **BRCA2 → GRCh38**，输入 `chr13-32380106-G-T`。同一等位的 GRCh37 输入是 `chr13-32954243-G-T`，也支持带版本的 `NM_000059.4:c.9217G>T`。应识别为 CA10579818、p.Asp3073Tyr。一次一个 ALT；坐标与组装一起填写。单独蛋白名不能确定核酸等位。
3. **对照一条证据。** 打开示例的 **PP4**，可从[专家个案](https://erepo.clinicalgenome.org/evrepo/ui/classification/237f2389-fa0c-4474-9f7a-50a8f9a3542d)读到临床 LR 3.57。它来自 [Li 等的家族史资料（PMID 31853058）](https://pubmed.ncbi.nlm.nih.gov/31853058/)。同一行也给出 Parsons 2019 作者原表的手查入口；这个链接不表示 3.57 来自 Parsons，此数也不是 PP1 的共分离 LR。专家已采用与本机已选择分别记录。

## 如何从注释走到逐项判读

| 问题 | 取什么、怎样判断 | 如何核查 |
|---|---|---|
| 是否同一等位？ | 优先用原始组装的 VCF 坐标或一行 VCF，也可用带版本的 HGVS c.。VEP 注释转录本；ClinGen Allele Registry（CAR）核对核酸等位。ClinGen Linked Data Hub（LDH）可按 CAid 独立查询，不需要先得到蛋白名称。 | 看组装、REF/ALT、转录本版本和 CAid。原组装失败和另一组装成功分别记录。 |
| 数据能支持哪档证据？ | 频率看 gnomAD 的发布版、指定队列、AC、AN 和覆盖；功能实验先看规范 Table 9 对该等位的推荐与排除理由。ClinVar 提交和符合条件的 eRepo 专家个案分别列来源。 | 在「判据与出处」查看观察值、规则阈值及报告原始字段。未查到、请求失败、测量值确实未达标各有不同含义。不能把空响应当作 AC=0。 |
| 下一篇论文先读哪里？ | 阅读卡优先提示正文段落、表行或补表，附 PMID。没有可读全文时可到 PubMed、PMC 或期刊手查。PP4/BP5 还可对照 [UCSC ENIGMA 临床 LR 轨道说明](https://hgdownload.soe.ucsc.edu/hubs/enigma/enigma.html)和作者原表。 | 查同一 ALT、样本／队列有无重叠；论文数量不决定证据强度。综合 LR 和已经计入的分量不能重复使用。 |
| 本机结论如何形成？ | 每个基因可展开冻结规范启用的 **19 条主判据**。勾选要采用的来源及允许的强度，点「校验、保存并计算」。报告保留未采用项与选择记录。 | 可核对已接受小计与组合草案；仅看总分不能代替组合规则或遗传学专家签发。 |

## 版本、运行方式与再核查

BRCA Review Desk 0.2.0rc8 固定使用 [BRCA1 GN092 v1.2](https://cspec.genome.network/cspec/ui/svi/doc/GN092?version=1.2) 和 [BRCA2 GN097 v1.2](https://cspec.genome.network/cspec/ui/svi/doc/GN097?version=1.2)。离线示例是固定快照；主动查询则读取当时各公共服务的响应。打开「各库版本与差异」和「原始数据」，对照可取得的来源版本、组装、转录本、查询 UTC 和响应 SHA-256。时间戳不代替数据库版本。GRCh37 服务如返回 HTTP 500 或 HTML，相应来源记为服务失败，不能据此判定等位未收录。

本机只需 Python 标准库，不需 Perl、Docker、VEP cache 或 API key。主动查询时会向相应公共服务发送该核酸等位；只试读离线教学不必联网。界面以中文为主。

[脚本、REST 和 curl 示例](METHODS_API_ZH.md) · [解读员手册](docs/Reviewer_Guide_ZH.md) · [临床 LR 与文献](docs/BRCA_Literature_and_Clinical_LR_ZH.md) · [验收与复现](docs/Acceptance_ZH.md) · [第三方来源与许可说明](SOURCE_TERMS.md)。本项目并非 ClinGen 或 ENIGMA 官方软件。
