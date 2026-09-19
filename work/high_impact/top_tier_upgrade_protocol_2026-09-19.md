# Frozen protocol: transportability and independent multi-donor MSC validation

Protocol frozen: 19 September 2026, before inspecting either analysis endpoint.

## Objective

Test two weaknesses that cannot be resolved by adding more in-sample pathway analyses: whether stiffness-associated transcriptomic effects transport forward in publication time, and whether their magnitude localizes to bone-marrow mesenchymal stromal cells in an independent multi-donor human atlas.

## Analysis 1: chronological leave-future-out transportability

- Eligible studies are the seven studies already included in the study-unit meta-analysis. No study is added or removed after endpoint inspection.
- Studies are ordered by their GEO submission date. Sequential tests begin when three earlier studies are available. Ties are resolved by accession number.
- For each held-out study, training effects are inverse-variance pooled using only earlier studies and genes observed in at least three training studies.
- Primary diagnostic endpoint: Spearman correlation between the frozen training effect and the held-out study effect over all shared eligible genes.
- Secondary endpoint: directional agreement among the 200 genes with the largest absolute training effect. The number 200 is fixed before analysis and is not optimized.
- A fixed-signature endpoint uses the earliest three studies only and evaluates every later study without refitting.
- Results are interpreted as transportability diagnostics. Gene-wise nominal P values are not used to claim independence because genes are correlated.

## Analysis 2: independent multi-donor human MSC atlas

- Dataset: GSE182158, 11 donors from bone marrow (n=3), adipose (n=3), dermis (n=3), and umbilical cord (n=2).
- The donor is the inferential unit. Cells are summed to donor pseudobulk profiles; no cell is treated as an independent replicate.
- Counts are converted to counts per million, transformed as log2(CPM + 0.5), and ranked within donor.
- The frozen response-magnitude set is the 200 genes with the largest absolute meta-analysis z statistic in the existing seven-study analysis. It is not refit using GSE182158.
- Primary endpoint: mean within-donor percentile expression rank of the frozen response-magnitude set.
- Primary contrast: bone marrow versus the other three tissues combined. The P value is the exact two-sided permutation P value over all allocations of 3 of 11 donors to the bone-marrow group.
- Specificity endpoint: an exact four-group permutation test of the between-tissue sum of squares, enumerating all 11!/(3!3!3!2!) = 92,400 tissue-label assignments.
- Positive-direction and negative-direction frozen sets (top 100 each by signed meta z) are secondary and descriptive, with Holm correction across the two tests.
- Sensitivity analysis matches each signature gene to background genes by overall expression decile and reports the signature-minus-matched-background score. The random seed is 20260919.
- No cluster labels, cell-level P values, or post hoc gene-set changes are permitted.

## Decision rule

These analyses will be added to the manuscript regardless of direction if they complete successfully. Claims will be strengthened only when the prespecified effect direction and exact P value support them. Null or reversed findings will narrow the conclusion rather than be omitted.
