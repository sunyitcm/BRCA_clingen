# rc8：包完整性和解读员验收

这套验收检查分发包、软件行为和证据出处，不代表对临床诊断性能的验证。

## 1. 下载与静态试用

把版本 ZIP 和同一 Release 的 `SHA256SUMS.txt` 放到同一目录。macOS 用 `shasum -a 256 -c SHA256SUMS.txt`，Linux 用 `sha256sum -c SHA256SUMS.txt`。Windows PowerShell 用 `Get-FileHash .\brca-review-desk-0.2.0-rc.8.zip -Algorithm SHA256` 与校验文件逐字比较。校验值用于检查字节是否一致，不是发行者签名。

解压后在含 `start.py` 的目录执行（Windows 把 `python3` 换成 `py -3`）：

```bash
python3 verify_files.py
python3 start.py --self-test
python3 -m unittest discover -s tests
```

| 手工验收 | 怎样做 | 应看到什么 |
|---|---|---|
| PP4 手查原表 | 断网打开 START_HERE.html → BRCA2 的 PP4 入口；再联网点「下载 Parsons 2019 原始 Excel」 | 先定位到 PP4 行，K 列提示与论文链接直接可见。另查 BP5 也是 K 列；PP1/BS4 为 F 列。工作表 SuppT1 All Data；BRCA2 c.3G>A 在 745 行。下载不是命中，不勾选不计分。 |
| UCSC | 对同一 BRCA2 c.9217G>T，在两个输入组装的 PP4 行打开 UCSC | 分别指向 hg19 或 hg38，ENIGMA 组 show、BRCAmla full。外站若要求验证，在浏览器完成；点同一核酸 ALT 查看分量。 |
| 阅读优先级 | 首页点 BRCA1 的「试读：Table 1 与 Results」 | 同一 PMID 有 6 处公开教学定位，先读该表与 Results；不把同队列按两条证据计入。 |
| 保存与重开 | 运行本机程序，选固定 BRCA2 教学快照。先保存默认 PS3、PP3，再加入专家 PM2、PP4 并保存 | 先 +5，再 +7（专家引用 +2）。重开可见保存值。纯静态教学页只读。 |

## 2. 公共教学样例与真实查询

发布样例是从公开信息整理的教学摘录。样例、阅读卡和旧离线页不附制作者的查询 UTC；示例报告哈希对应**本公开版本的字节**。真实数据库发布日期、论文年份、规范发布号和内容指纹保留。使用者在本机新查仍记录真实 UTC、响应和哈希。早期内部教学记录不自动映射为新示例的已采纳证据。

## 3. 维护者复建 GitHub Release

从与目标版本匹配的源码目录运行 `python3 scripts/build_release.py --output-dir dist`。脚本先逐项核对 `MANIFEST.sha256`，只收录清单和清单自身；不会遍历打包工作目录。输出单一带版本顶层目录的 ZIP 和单独的 `SHA256SUMS.txt`。文件排序、权限和 ZIP 头时间固定，不携带文件系统用户名、扩展属性或实际修改时间。相同源码在同一 Python/zlib 工具链下应得到相同字节。清单更新应逐项审查，不能把运行目录加入清单。

建议在已确认权属及许可的干净仓库中，从对应提交创建标签 `v0.2.0-rc.8`，在 **Releases → Draft a new release** 勾选 **pre-release**。使用 [发布说明](../RELEASE_NOTES.md)，上传上述两件资产。GitHub 会另外生成源码归档。包内不带仓库历史、账号、远程地址或发布凭据。[GitHub Releases 说明](https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases)
