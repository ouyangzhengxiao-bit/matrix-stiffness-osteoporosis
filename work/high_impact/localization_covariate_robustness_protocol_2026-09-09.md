# Post-review localization covariate robustness plan

Dated 2026-09-09, after the original localization results and therefore labelled post hoc.

The original Chai single-cell and spatial osteoblast-lineage unions will be retested after adding GENCODE v49 gene-length deciles to the original matching factors. The meta-level null jointly matches number of contributing studies, Hartung-Knapp standard-error decile and gene-length decile. The study-level null jointly matches study-effect variance decile and gene-length decile. Twenty thousand draws and seed 20260909 will be used.

Two union-level tests constitute one family and will receive Bonferroni correction. At the study level, the exact one-sided sign test is primary; one-sided t tests and t-based 95% confidence intervals are sensitivity summaries. Leave-one-study-out mean raw matched effects will be reported for every omission. These analyses assess robustness and do not replace the original prespecified tests.
