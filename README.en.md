# BRCA_clingen: BRCA1 and BRCA2 germline variant annotation and evidence review

[中文](README.zh-CN.md) · [Bilingual home](.github/README.md) · [Original rc8 guide (Chinese)](README.md)

BRCA Review Desk is an independent Python tool for molecular geneticists and clinical bioinformaticians reviewing germline variants in hereditary cancer predisposition. It uses Ensembl Variant Effect Predictor (VEP) to annotate transcript consequences and the gene-specific ClinGen ENIGMA BRCA1/2 Variant Curation Expert Panel (VCEP) specifications to review ACMG/AMP criteria. The report shows what was observed, which rule applies and what the reviewer actually selected. It does not sign out clinical classifications or assess somatic oncogenicity or treatment response.

## Walk through one public allele

1. **Start offline.** Once the release is public, find **Releases → v0.2.0-rc.8 → Assets** and download `brca-review-desk-0.2.0-rc.8.zip` and `SHA256SUMS.txt`. Extract the full product ZIP and open `START_HERE.html` in a browser. No installation is needed to read the examples. The application interface is primarily in Chinese.
2. **Run a local query.** With Python 3.10 or newer, run `python3 start.py` from the extracted directory (Windows: `py -3 start.py`). Choose **BRCA2 → GRCh38**, and enter `chr13-32380106-G-T`. The same allele on GRCh37 is `chr13-32954243-G-T`; the versioned transcript input is `NM_000059.4:c.9217G>T`. Expect CA10579818 and p.Asp3073Tyr. Enter one alternate allele at a time, using its original assembly. A protein change alone does not identify a nucleotide allele.
3. **Trace a criterion.** Open **PP4** in the example. The [expert curation](https://erepo.clinicalgenome.org/evrepo/ui/classification/237f2389-fa0c-4474-9f7a-50a8f9a3542d) cites a clinical likelihood ratio (LR) of 3.57 from [Li et al., PMID 31853058](https://pubmed.ncbi.nlm.nih.gov/31853058/). That value comes from family history, not segregation. The Parsons 2019 workbook linked on the same row is a separate lookup; it is not the source of this value.

## From annotation to a reviewable decision

| Question | What the software does | What to verify |
|---|---|---|
| Is this the same allele? | Accepts assembly-specific coordinates, a VCF row or versioned HGVS c. VEP annotates transcripts. The ClinGen Allele Registry (CAR) resolves the allele. The Linked Data Hub (LDH) can be queried by CAid independently of the VEP protein result. | Compare assembly, REF/ALT, transcript version and CAid. Cross-assembly results retain their source assembly. |
| Which criterion can the data support? | Shows gnomAD release, cohort, counts and coverage when supplied; ENIGMA Table 9 recommendations for available functional assays; and relevant ClinVar records or a matching eRepo expert curation. | Compare each observation with its rule. A missing response, an inapplicable rule and an observation below threshold are different results. A failed request does not imply AC=0. |
| What should I read first? | Reading cards point to a section, table row or supplement when available and link to PubMed or full text for manual review. PP4/BP5 also link to the [UCSC ENIGMA track information](https://hgdownload.soe.ucsc.edu/hubs/enigma/enigma.html) and the authors’ workbook. | Match the exact alternate allele; check whether studies share participants. Paper counts are not evidence strength. Do not count a combined LR and an included component twice. |
| What enters the calculation? | Expand all **19 applicable top-level criteria per gene** in the pinned specifications. Select a criterion and its allowed strength, then validate, save and calculate. | A published expert assertion is distinct from locally accepted evidence. Review both the subtotal and the allowed strength combinations; a sum alone does not establish a classification. |

## Versions and reproducibility

Release candidate 0.2.0rc8 uses [BRCA1 GN092 v1.2](https://cspec.genome.network/cspec/ui/svi/doc/GN092?version=1.2) and [BRCA2 GN097 v1.2](https://cspec.genome.network/cspec/ui/svi/doc/GN097?version=1.2). The offline examples are fixed snapshots. On-demand requests contact public services and send the queried nucleotide allele. Inspect the source version when supplied, assembly, transcript, retrieval time (UTC) and response SHA-256. Retrieval time is not a database version. A GRCh37 endpoint returning HTTP 500 or HTML is recorded as a service failure, not evidence of absence.

Local queries use the Python standard library; no Perl, Docker, VEP cache or API key is required for the included workflow. Technical references (in Chinese): [API methods and curl examples](METHODS_API_ZH.md), [reviewer guide](docs/Reviewer_Guide_ZH.md), [literature and clinical LR](docs/BRCA_Literature_and_Clinical_LR_ZH.md), [acceptance checks](docs/Acceptance_ZH.md) and [source terms](SOURCE_TERMS.md). This tool is not an official ClinGen or ENIGMA product.
