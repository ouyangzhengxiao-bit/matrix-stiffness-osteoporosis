# Frozen multicohort stiffness-signature validation protocol (v1)

Frozen on 2026-09-07 after inspecting GEO study designs and before inspecting expression values from GSE193021, GSE181512, GSE226411, GSE288678, or GSE55867.

## Fixed discovery signature

The 100 positive and 100 negative human genes already frozen from GSE22011 are retained without alteration. No validation cohort may change gene membership, direction, or scoring.

For every sample, the score is the mean within-sample percentile rank of mapped positive genes minus the mean percentile rank of mapped negative genes. This rank score is used to limit platform-scale effects.

## Independent biological-context cohorts

- GSE193021: human bone-marrow MSCs, 0.5 versus 32 kPa, three samples per group at days 3 and 7. Analyze each time separately by exact two-sided allocation test; two contrasts.
- GSE181512: rat bone-marrow MSCs, 2 versus 18 kPa at 48 h, three samples per group. Map only NCBI ortholog-group-supported human-to-rat genes; exact two-sided allocation test.
- GSE226411: human bone-marrow MSCs in soft versus stiff hydrogels, three samples per group. Exact two-sided allocation test.
- GSE288678: human bone-marrow MSCs, 150 versus 900 kPa across three ECM formulations, four culture replicates per stiffness/formulation cell. Estimate the common high-minus-low score effect with ECM fixed effects. Obtain a two-sided randomization P value by permuting stiffness labels within each ECM stratum. Report formulation-specific effects as secondary estimates.
- GSE55867: human MSCs from three donors, paired 1.5 versus 50 kPa PDMS conditions. Exclude tissue-culture plastic. Use the mean paired high-minus-low difference and an exact two-sided sign-flip test over the three donor pairs.

GSE255574 remains the previously analyzed independent 3D human MSC context. Its poly(A) 24/48 h results and the separately frozen 48 h ribo-depleted technical replication are carried forward unchanged.

## Direction, multiplicity, and synthesis

Positive transfer means a higher frozen score in the stiffer condition. The six new biological contrast tests (two GSE193021 time points plus one test in each other cohort) receive a single Benjamini-Hochberg correction. Both raw and adjusted P values are reported; small designs are not declared validated from direction alone.

Across all biological contrasts, summarize the number of positive and negative effects and use an exact one-sided binomial sign test against 0.5 only as a direction-concordance summary. Do not pool effect magnitudes across platforms, species, dimensionality, or stiffness ranges into a single biological effect size.

## Integrity rules

- Use only official processed GEO matrices and official platform annotations; record SHA-256 hashes.
- No DEG-derived replacement signature, threshold changes, subgroup deletion, or outcome switching after values are seen.
- Technical/culture replicates are described according to source metadata and are not relabeled as independent donors.
- A null or reversed result is retained and interpreted as limited transferability/context dependence.

## Post-result estimation addition

After all primary exact tests above were completed, deterministic stratified bootstrap 95% intervals (seed 20260907; 10,000 resamples) were added for visualization and effect-size estimation. These intervals do not replace the prespecified exact tests and do not change multiplicity decisions.
