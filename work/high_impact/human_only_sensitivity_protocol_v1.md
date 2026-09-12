# Frozen protocol: human-only stiffness sensitivity v1

Frozen: 2026-09-09, before inspecting human-only results.

## Scope
Exclude the two rat stiffness studies (GSE181512 and GSE310514). Retain five independent human MSC studies: GSE193021, GSE226411, GSE288678, GSE55867, and GSE255574. Reuse the already generated study-level effects without reprocessing or contrast reselection.

## Primary analysis
Fit the same REML random-effects model with modified Hartung-Knapp uncertainty as the seven-study analysis. Include genes observed in at least three of five human studies. Apply Benjamini-Hochberg FDR across all eligible genes. Report the number at q<0.05 and Spearman concordance with the seven-study pooled effects.

## Translational sensitivity
For GSE156508, GSE35958, and GSE230665, calculate (1) Spearman correlation between absolute human-only meta z scores and absolute case-control Welch t statistics and (2) a percentile-rank score based on the 100 largest and 100 smallest human-only pooled effects. Use every unique case-label allocation for cohort-level P values. Combine cohort standardized statistics with equal cohort weights and 500,000 Monte Carlo draws sampled from the exact cohort nulls (seed 20260909).

## Interpretation
This is a species sensitivity analysis, not a new discovery model. Lack of disease replication limits generalization; a positive result in one cohort cannot establish a diagnostic biomarker.
