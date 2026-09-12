# Frozen protocol: multi-cohort osteoporosis transcriptome validation v1

Frozen: 2026-09-09, before computing stiffness–disease magnitude correlations for GSE35958 or GSE230665.

## Independent cohorts

1. GSE156508: primary femoral-head osteoblasts, 6 osteoporotic-fracture versus 6 severe-osteoarthritis donors.
2. GSE35958: cultured femoral-head marrow MSCs, 5 primary-osteoporosis versus 4 age-matched non-osteoporotic donors.
3. GSE230665: femoral tissue, 12 postmenopausal-osteoporosis versus 3 postmenopausal controls.

GSE35956 is excluded as an independent cohort because its five osteoporosis samples duplicate the five cases in GSE35958 and its controls are younger, confounding disease with age. GSE157322 is excluded from primary inference because low/high BMD is fully confounded with European/Polynesian ancestry.

## Primary endpoint

Within each cohort, compute Spearman correlation between absolute stiffness meta z score and absolute Welch disease t statistic among shared genes. Obtain an exact one-sided P value by enumerating every allocation with the observed case count.

For the combined test, standardize each observed cohort correlation against its complete within-cohort label-permutation distribution and average the three standardized values with equal cohort weight. Estimate a one-sided P value using 500,000 Monte Carlo draws from the Cartesian product of the three complete null distributions (seed 20260909). Equal weighting prevents the largest gene universe or most imbalanced cohort from dominating.

## Secondary endpoints

- Signed stiffness–disease correlation, two-sided exact within-cohort tests and a combined equal-weight standardized statistic.
- Donor-level fixed 200-gene signature separation, two-sided exact within-cohort tests and a combined equal-weight standardized statistic.
- Leave-one-donor-out direction audit for the primary correlation in each cohort.

No cohort is removed after results are inspected. Cohorts differ in cell source and comparator; synthesis is interpreted as triangulation of disease relevance, not a diagnostic accuracy analysis.
