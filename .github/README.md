# BRCA_clingen | BRCA1/BRCA2 胚系变异注释与证据判读 · Germline variant annotation and interpretation

[中文：使用与判读方法](../README.zh-CN.md) · [English: quick start and methods](../README.en.md)

**BRCA Review Desk** 是给分子遗传学家和临床生信人员使用的本机审核工具。这套 **BRCA1/BRCA2 变异注释方法**用 Ensembl VEP 核对转录本后果，再按 ClinGen ENIGMA 的基因专属 ACMG/AMP 规范逐条核对证据。输入一个 VCF 位点或带版本的 HGVS c.，可以看到哪些判据有依据、观察值和阈值是什么、哪些仍缺资料。适用于遗传性肿瘤易感性的**胚系**判读。

**BRCA Review Desk** is a local Python tool for germline BRCA1 and BRCA2 variant annotation and evidence review in hereditary cancer predisposition. It pairs Ensembl Variant Effect Predictor (VEP) annotations with the ClinGen ENIGMA BRCA1/2 Variant Curation Expert Panel (VCEP) specifications. Each criterion links the observation to its source, rule and review decision. This independent research tool does not issue clinical classifications on a laboratory’s behalf.

| 你想做什么 / Goal | 从这里开始 / Where to start |
|---|---|
| 先看示例 / Explore offline | 从 GitHub **Releases** 取产品 ZIP，完整解压后双击 `START_HERE.html`。Download the release asset, extract it, then open the HTML file. No install or account required. |
| 查一个位点 / Query a variant | 在解压目录运行 `python3 start.py`（Windows: `py -3 start.py`）。Choose the gene and original assembly. For example: **BRCA2, GRCh38, `chr13-32380106-G-T`**. Python 3.10+; no API key required for the supplied queries. |
| 看方法与来源 / Inspect the evidence | [逐项判读 / Criterion-by-criterion guide (中文)](../docs/Reviewer_Guide_ZH.md) · [APIs, GET/POST and curl (中文)](../METHODS_API_ZH.md) · [Clinical LR and literature (中文)](../docs/BRCA_Literature_and_Clinical_LR_ZH.md) |

**例 / Example:** `NM_000059.4:c.9217G>T` → BRCA2 p.Asp3073Tyr. 核对一个位点时，可以同时查看 VEP 的蛋白后果、gnomAD 的人群数据、功能实验原表和匹配的[专家个案 / expert curation](https://erepo.clinicalgenome.org/evrepo/ui/classification/237f2389-fa0c-4474-9f7a-50a8f9a3542d)。已发布的专家断言与本机勾选的证据分开显示。A missing response is not evidence that the allele is absent from a population.

**版本 / Version:** BRCA Review Desk 0.2.0rc8 uses the frozen ClinGen ENIGMA specifications [BRCA1 GN092 v1.2](https://cspec.genome.network/cspec/ui/svi/doc/GN092?version=1.2) and [BRCA2 GN097 v1.2](https://cspec.genome.network/cspec/ui/svi/doc/GN097?version=1.2). Live source versions may differ from the teaching snapshots. [版本说明 / Release notes](../RELEASE_NOTES.md) · [来源与使用条件 / Source terms](../SOURCE_TERMS.md)。

[API 与数据使用说明 / API and data use](../docs/API_SOURCE_TERMS.md)
