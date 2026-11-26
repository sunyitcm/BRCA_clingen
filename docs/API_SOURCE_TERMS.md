# API 与数据使用说明 / API and data use

BRCA Review Desk 以 Python 向公开服务发出只读请求。对每条资料，请区分提供服务的机构、原始数据来源、适用的数据条款及本项目代码许可。免登录是一种访问方式，不是对全部结果的再分发授权。

| 对象 | 本项目如何使用、如何注明来源 |
|---|---|
| **Ensembl VEP 软件及 REST 服务** | 通过 Ensembl REST 取得注释；本版不分发 VEP 的 Perl 软件或 cache。[VEP release/116 的代码许可](https://github.com/Ensembl/ensembl-vep/blob/release/116/LICENSE)是 **Apache License 2.0**。访问服务参见 [EMBL-EBI 使用条款](https://www.ebi.ac.uk/about/terms-of-use/)；[Ensembl 数据说明](https://www.ensembl.org/info/about/legal/disclaimer.html)区分 Ensembl 项目生成的数据和可能有独立限制的第三方数据／软件。请引用 Ensembl VEP 和实际使用的数据源，记录可取得的发布版；不要把整个 JSON 或自有客户端标成 Apache 或 MIT。 |
| **ClinGen ENIGMA 规范及专家整理内容** | 本版固定 [BRCA1 GN092 v1.2](https://zenodo.org/records/21434315) 和 [BRCA2 GN097 v1.2](https://zenodo.org/records/21434343)。这两个标准的存档都明确援引 **CC0 1.0**。[ClinGen 条款和引用方式](https://www.clinicalgenome.org/docs/terms-of-use/)说明了其对公开专家整理内容的政策，并请求在可行时注明 ClinGen 和访问日期。引用具体标准时写明 ENIGMA VCEP、基因、版本及存档链接；引用某条专家判读时保留对应个案链接。中文释义由本项目整理，并非官方译本。 |
| **ClinGen CAR、LDH 与外部证据** | [CAR](https://reg.clinicalgenome.org/) 用于核对等位身份，可用 CAid 指回记录；[LDH](https://ldh.clinicalgenome.org/ldh/ui/about) 提供外部来源的链接和部分字段。引用 LDH 条目时同时标明提供原数据的数据库、可获得的版本及具体记录。gnomAD、BRCA Exchange、预测资源、论文和作者补表等分别按原来源条件使用。ClinGen 对自有整理内容的 CC0 声明不能替代第三方许可；提供原表链接不等于获得镜像整表的许可。 |

**本项目代码。** 本候选版没有附加由权利人选定的开源代码许可证。调用这些 API 不要求给本项目自动加 MIT，也不授予对第三方材料的转许可。详见[随包来源说明](../SOURCE_TERMS.md)。

**本机与联网。** 打开离线示例不需要向 API 提交位点。主动查询会把该次核酸等位或其标识发往相应服务，服务端可能看到网络出口 IP。公开示例与使用者在自己设备上产生的查询、响应和审查记录分开。引用数据库可记实际访问的日历日期，并同时注明版本；实验室内部审计仍可保留精确 UTC。

**English.** This Python client queries the Ensembl and ClinGen services. The Ensembl VEP software uses Apache 2.0; EMBL-EBI service terms and the originating databases’ data conditions apply separately. The two archived BRCA specifications above explicitly state CC0 1.0. ClinGen requests attribution and, where appropriate, the access date. A CAR identifier or an LDH link does not relicense external records, prediction scores or authors’ supplements. This release candidate does not declare an open-source license for the project’s own code. Offline examples require no API submission; a live query sends that allele or identifier to the relevant service.
