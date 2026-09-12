# Matrix stiffness transcriptomic responses and osteoporosis

Reproducible code, dated protocols, systematic-search records, derived results and source figures for a study-clustered meta-analysis of matrix-stiffness transcriptomic responses and their relationship to human osteoblast-lineage programs and osteoporosis cohorts.

## Main findings

- Seven independent stiffness studies contributed 13 dependent contrasts.
- No gene among 12,241 eligible genes passed false-discovery-rate correction.
- Pooled-gene pathway signals did not survive inference with the independent study as the unit.
- Response magnitude localized to osteoblast-lineage programs, with strongest support for the spatial program.
- A human-only sensitivity analysis preserved the central conclusion.
- Three independent human osteoporosis cohorts did not support universal cross-cohort transfer of response magnitude, direction or a fixed signature.

## Repository contents

- `work/upgrade/`: core meta-analysis and pathway scripts.
- `work/high_impact/`: localization, orthogonal validation, disease validation, sensitivity, figure and quality-control scripts.
- `outputs/研究升级/meta_analysis_v2/`: study-level effects and meta-analysis outputs.
- `outputs/研究升级/high_impact_v4/`: validation, sensitivity, search-audit and environment outputs.
- `outputs/final_analysis/`: processed matrices needed for the two additional osteoporosis cohorts.
- `figures/`: publication figures in PNG and PDF formats.
- `docs/`: PRISMA checklist and evidence-ceiling audit.

## Reproducibility

The analyses use Python 3.14 with package versions in `requirements.txt` and R 4.5.2 with versions recorded in `outputs/研究升级/high_impact_v4/R_package_versions.txt`. Large public source files are not redistributed. Their accessions and SHA-256 hashes are recorded in the input manifests.

The primary public sources include GSE193021, GSE181512, GSE226411, GSE288678, GSE55867, GSE255574 and GSE310514. Validation sources include GSE147390, GSE156508, GSE35958, GSE230665, GSE166824, GSE22011, GSE152708, GSE276529, GSE317531, GSE287556, GSE317069 and GSE299207, plus GWAS Catalog study GCST90726625.

Run the version-5 orchestration script after placing the documented public inputs at the relative paths expected by the individual scripts:

```bash
python3 work/high_impact/run_full_pipeline_v50.py
```

Many downstream validation and audit results can be reproduced directly from the included derived matrices. Each protocol records whether it was frozen before its corresponding merge or added as a labelled sensitivity analysis.

## Data and code boundaries

This repository contains original analysis code and derived numerical results. It does not redistribute large raw public datasets. Results are for research use and do not constitute a diagnostic test, a measurement of tissue material stiffness or evidence of causality.

## Citation

Citation metadata will be finalized before the first public release and Zenodo archive. After release, cite the version DOI shown in `CITATION.cff`.
