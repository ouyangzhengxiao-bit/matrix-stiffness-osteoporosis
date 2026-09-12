# Cross-study mechanotranscriptomic meta-analysis protocol v2 (frozen before analysis)

## Aim

Estimate which transcriptional responses to extracellular-matrix stiffness are reproducible across mesenchymal stromal-cell studies, quantify context dependence, and test whether any reproducible component relates to bone strength or osteoporosis.

## Eligibility and units of analysis

- Public transcriptomic studies with at least two experimentally imposed stiffness levels, an MSC/mesenchymal progenitor population, and at least two biological samples per stiffness condition.
- Main synthesis: human or rat bone-marrow MSC bulk or pseudobulk expression. Fibroblast GSE22011 is an external cell-type comparison, not part of the MSC primary meta-analysis.
- Each GEO series is one study cluster. Multiple times, ECM coatings, passages, or libraries within a series are dependent contrasts and will not be counted as independent studies.
- Ribosomal-RNA-depleted GSE255574 libraries are technical corroboration only.
- Studies are excluded before examining meta-analytic results when stiffness is confounded with another treatment or processed expression is unavailable and raw-read reprocessing is infeasible.

## Effect estimation

- Map features to official gene symbols and retain one-to-one human-rat orthologues for cross-species analysis.
- Within each contrast, transform expression to log2 scale when needed and standardize each gene across samples. Estimate stiff-minus-soft standardized mean differences (Hedges g) for two-level designs and standardized slopes per log2 stiffness for multi-level designs.
- For paired donors or blocked ECM designs, estimate effects within donor/stratum before aggregation.
- Combine dependent contrasts within each GEO series to one study-level effect per gene using equal weighting; sensitivity analyses use one prespecified representative contrast per series.

## Primary inference

- Random-effects meta-analysis by restricted maximum likelihood with Hartung-Knapp uncertainty where estimable; report pooled effect, 95% CI, tau-squared, I-squared, and prediction interval.
- Control FDR across genes with Benjamini-Hochberg. A conserved gene requires q < 0.05, at least four independent study clusters with the gene measured, no single study contributing more than 50% of inverse-variance weight, and the same pooled direction in every leave-one-study-out analysis.
- Primary pathway analysis uses rank-based enrichment of the complete gene-level meta-analysis, with BH correction across tested gene sets. Gene sets must contain 15-500 measured genes.

## Heterogeneity and validation

- Prespecified moderators: 2D versus 3D culture, human versus rat, exposure duration (<=72 h versus longer), and stiffness ratio on log2 scale. Moderator results are exploratory because study counts are small.
- Leave-one-study-out analyses assess dependence on individual series.
- External validation uses GSE166824 single-cell pseudobulk stiffness contrasts, GSE22011 fibroblast dose response, GSE152708 measured mouse-bone mechanics, and GSE276529 human glucocorticoid-induced osteoporosis tissue. Datasets used to estimate a signature are never reused as validation for that signature.
- Technical library replication is labeled separately from biological replication.

## Decision rule for journal tier

- Consider a broad high-impact computational/biological journal only if at least one coherent pathway passes FDR, remains directionally stable in leave-one-study-out analysis, and shows independent support in at least one bone-specific or human disease dataset.
- If the conserved signal is absent or context-dependent, frame the paper as a rigorous boundary-condition/resource study and select a transparent general-scope journal rather than claiming a universal stiffness program.
