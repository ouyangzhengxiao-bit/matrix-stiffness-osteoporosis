# Frozen protocol: transcriptomic prediction of femoral mechanical performance

Date frozen: 2026-09-07, before any gene-wise predictive model was fitted.

This analysis asks whether cortical-bone transcriptomic information improves prediction of measured femoral mechanical performance beyond animal characteristics and micro-CT geometry. It is an exploratory, prospectively specified reanalysis, not a registered study. The previously tested D3/D6 signatures and their null associations are not used for feature selection.

## Data split

- Eligible animals are the 192 GSE152708 mice with a unique expression identifier and matching phenotype record.
- Animals missing any required primary-model variable are excluded with reasons retained.
- Generation G29 is an untouched temporal/batch transport test set. Generations G23, G24, G27 and G28 form the development set.
- Hyperparameters are selected only in the development set by leave-one-generation-out cross-validation. G29 is evaluated once after the pipeline is fixed.

## Outcomes and model families

The primary outcome is four-point-bending maximum load (N). Structural stiffness (N/mm) is secondary.

Four prespecified ridge-regression models are compared:

1. Animal model: sex, body weight and age.
2. Geometry model: animal variables plus femoral length, cortical area, total area, cortical thickness, tissue mineral density, cortical porosity, Imin, Imax and pMOI.
3. Transcriptome model: animal variables plus log2(CPM+1) gene expression.
4. Combined model: animal, geometry and transcriptome variables.

Gene filtering and variance ranking occur independently inside every training fold. A gene must exceed 1 CPM in at least 20% of that fold's samples. At most the 2,000 most variable eligible genes are retained. Continuous predictors are standardized from the relevant training fold. Categorical sex is one-hot encoded. Ridge alpha is selected from 10^-4 to 10^6 on a logarithmic grid by mean cross-validation RMSE. No outcome-informed pathway or gene preselection is allowed.

## Evaluation

Primary evidence is the G29 difference in R-squared and RMSE between the combined and geometry models for maximum load. Transcriptome-versus-animal and geometry-versus-animal comparisons are descriptive. Secondary comparisons use structural stiffness.

Uncertainty for test-set RMSE differences is estimated by 10,000 paired bootstrap resamples of G29 animals with seed 20260907. The bootstrap interval is descriptive because there is only one held-out generation. Development cross-validation performance, test MAE, Pearson correlation and calibration slope/intercept are also reported. Negative test R-squared is retained.

## Interpretation rules

- Improvement is considered supported only if the combined model has lower G29 RMSE than the geometry model and the paired 95% bootstrap interval for `RMSE_combined - RMSE_geometry` excludes zero.
- A favorable development result without G29 improvement is classified as non-transportable.
- Gene coefficients are not interpreted mechanistically. Post hoc stability and pathway analyses are allowed only if the primary transport criterion is met; they are labeled exploratory.
- Generation is intentionally the transport axis and is not included as a predictor because G29 is unseen during training. Generation and micro-CT batch are perfectly confounded in the available records.
- The split reduces leakage but does not account for the complete genomic relatedness matrix. Results are predictive associations in mice, not causal effects or human clinical validation.
