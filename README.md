# Matrix stiffness transcriptomic responses and osteoporosis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22726635.svg)](https://doi.org/10.5281/zenodo.22726635)

Reproducible code, dated protocols, systematic-search records, derived results and source figures for a study-clustered meta-analysis of matrix-stiffness transcriptomic responses and their relationship to human osteoblast-lineage programs and osteoporosis cohorts.

## Main findings

- Seven independent stiffness studies contributed 13 dependent contrasts.
- No gene among 12,241 eligible genes passed false-discovery-rate correction.
- Pooled-gene pathway signals did not survive inference with the independent study as the unit.
- Response magnitude localized to osteoblast-lineage programs, with strongest support for the spatial program.
- A human-only sensitivity analysis preserved the central conclusion.
- Three independent human osteoporosis cohorts did not support universal cross-cohort transfer of response magnitude, direction or a fixed signature.
- A model fixed on the earliest three studies failed to predict four later studies (median genome-wide Spearman rho = -0.050; median top-200 directional agreement = 0.493).
- In an independent 11-donor, four-tissue MSC atlas, the frozen response-magnitude set varied by tissue source but was not enriched in bone-marrow MSCs.

## Repository contents

- `work/upgrade/`: core meta-analysis and pathway scripts.
- `work/high_impact/`: localization, orthogonal validation, disease validation, sensitivity, figure and quality-control scripts.
- `outputs/研究升级/meta_analysis_v2/`: study-level effects and meta-analysis outputs.
- `outputs/研究升级/high_impact_v4/`: validation, sensitivity, search-audit and environment outputs.
- `outputs/研究升级/high_impact_v5/`: chronological transportability and independent multi-donor MSC atlas outputs.
- `outputs/final_analysis/`: processed matrices needed for the two additional osteoporosis cohorts.
- `figures/`: publication figures in PNG and PDF formats.
- `docs/`: PRISMA checklist and evidence-ceiling audit.

## Reproducibility

The analyses use Python 3.14 with package versions in `requirements.txt` and R 4.5.2 with versions recorded in `outputs/研究升级/high_impact_v4/R_package_versions.txt`. Large public source files are not redistributed. Their accessions and SHA-256 hashes are recorded in the input manifests.

The primary public sources include GSE193021, GSE181512, GSE226411, GSE288678, GSE55867, GSE255574 and GSE310514. Validation sources include GSE182158, GSE147390, GSE156508, GSE35958, GSE230665, GSE166824, GSE22011, GSE152708, GSE276529, GSE317531, GSE287556, GSE317069 and GSE299207, plus GWAS Catalog study GCST90726625.

Run the version-6 orchestration script after placing the documented public inputs, including the GSE182158 raw TAR, at the relative paths expected by the individual scripts:

```bash
python3 work/high_impact/run_full_pipeline_v60.py
```

Many downstream validation and audit results can be reproduced directly from the included derived matrices. Each protocol records whether it was frozen before its corresponding merge or added as a labelled sensitivity analysis.

## Data and code boundaries

This repository contains original analysis code and derived numerical results. It does not redistribute large raw public datasets. Results are for research use and do not constitute a diagnostic test, a measurement of tissue material stiffness or evidence of causality.

## Citation

For exact reproduction of the upgraded analyses, cite release v1.1.0 at [10.5281/zenodo.22839453](https://doi.org/10.5281/zenodo.22839453). The original v1.0.0 release remains available at [10.5281/zenodo.22726636](https://doi.org/10.5281/zenodo.22726636), and the concept DOI [10.5281/zenodo.22726635](https://doi.org/10.5281/zenodo.22726635) resolves to the latest archived release. Citation metadata are also provided in `CITATION.cff`.

## License

Original analysis code is released under the MIT License (`LICENSE`). Original derived tables, documentation and figures are released under the Creative Commons Attribution 4.0 International License (`LICENSE-DATA`). Source datasets retain the terms imposed by their originating repositories and data generators.

## Author

Zhengxiao Ouyang — ORCID: https://orcid.org/0000-0002-8997-0446  
Department of Orthopedics, The Second Xiangya Hospital, Central South University.
