# BRCA Review Desk 0.2.0rc8

**预发布 · 建议标签：`v0.2.0-rc.8`**  
适用于 BRCA1/2 相关癌症易感性的胚系变异研究审核。Python 3.10+，只用标准库。一个包支持 GRCh37 和 GRCh38。

### 从哪里开始

解压完整资产 `brca-review-desk-0.2.0-rc.8.zip`，双击 **START_HERE.html**。选择 BRCA2 的 PP4 教学入口，直接到作者原表链接；或在终端运行 `python3 start.py` 进入本机查询（Windows：`py -3 start.py`）。

### 本版的查证路径

- PP4/BP5 行提供 Parsons 2019 原始 Excel 的直接下载、论文入口、综合 LR 的 K 列提示；PP1/BS4 提示同表 F 列的共分离 LR。手工打开来源不会勾选或计分。
- PP4/BP5 位点链接按已核准的组装和等位打开 UCSC，将 ENIGMA 父组设为 show、LR 子轨道设为 full。查看分量仍需点目标 ALT。保留轨道方法与修订记录链接。
- 保留 BRCA1 的旧版离线样例、两套组装、文献阅读优先级和可核对的表行／段落定位。静态页没有本机写入或后台抓取功能。
- 公共教学摘录使用独立 ID；不随包携带采集时钟、工作目录或审核历史。各来源自身的发布日期、科学标识与指纹保留。真实查询与人工采纳的行为及冻结的判读规则保持不变。

### 校验、复现与来源

完整 ZIP 的 SHA-256 见同一 Release 的 `SHA256SUMS.txt`。解压后可运行 `python3 verify_files.py` 与 `python3 start.py --self-test`。完整复核见 [验收说明](docs/Acceptance_ZH.md)。GitHub 的 “Source code” 是平台按标签生成的另一个归档，哈希不与手工上传的资产通用。

[首页及使用说明](README.md) · [脚本与 API](METHODS_API_ZH.md) · [数据与权利说明](SOURCE_TERMS.md)。公开原表按原站链接取得；本包未镜像第三方工作簿或数据库。网络故障、查无记录、资料不足和不适用各自保留，不据此自动判定。
