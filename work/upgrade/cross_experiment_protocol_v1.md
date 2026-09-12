# Frozen protocol: cross-experiment matrix-response transport

Date frozen: 2026-09-07, before gene-level GSE22011 effects or GSE255574 expression values were examined.

## Question and study roles

This analysis asks whether a donor-consistent transcriptional response to increasing matrix stiffness in human lung fibroblasts transports to human marrow mesenchymal stromal cells (MSCs), and whether the resulting score associates with measured mouse-bone mechanics.

- GSE22011 is the discovery experiment: three independent human donors, each measured on 0.1, 0.4, 1.6, 6.4 and 25.6 kPa matrices after 48 hours.
- GSE255574 is the independent cell-context test: one selected human MSC donor, three reported culture replicates at 0.15, 0.5 and 2 kPa, at 24 and 48 hours. Only poly(A), DMSO samples are used. These are not independent donors.
- GSE152708 is the tissue-level test: 187 mice with complete covariates, cortical-bone expression and four-point-bending phenotypes.

The experiments differ in cell type, dimensional context, stiffness range and platform. Transport failure is therefore informative but does not establish absence of mechanotransduction.

## Discovery signature

GSE22011 processed RMA values are mapped using the official GPL6244 annotation. Missing or explicitly multigene symbols are excluded and repeated symbols are summarized by median expression. Within each donor, each gene is regressed on centered log2 matrix stiffness. The discovery effect is the mean donor-specific slope.

Eligible genes must have the same nonzero slope sign in all three donors. The 100 largest consistent positive slopes and 100 most negative consistent slopes form a fixed directional signature. Donor-level sign agreement and slopes for every mapped gene are retained. No P-value threshold or osteoporosis data are used for selection.

## Independent MSC test

GSE255574 transcript counts are summed to gene symbols using the deposited annotation, transformed to log2(CPM+1), and converted to within-sample percentile ranks. The signature score is mean rank of positive genes minus mean rank of negative genes. Separately at 24 and 48 hours, score is regressed on centered log2 stiffness.

For each time point, a two-sided exact permutation test enumerates all assignments of three samples to each of the three stiffness levels while retaining the observed values. The two time-point tests are adjusted by Benjamini-Hochberg. Replicates are culture preparations from one selected donor, so a significant result is cell-context transport, not donor replication.

### Amendment after the primary MSC result

The primary poly(A) analysis produced a significant response at 48 hours in the opposite direction from discovery. Before inspecting values from the ribosomal-RNA-depleted libraries, the nine deposited 48-hour DMSO-free/untreated samples are designated a technical/library-preparation replication. They use the same score and exact three-group permutation test. This single follow-up test is reported without combining it with the primary family. It can support robustness to library preparation, but it cannot add an independent donor or rescue the failed directional transport criterion.

## Mouse-bone test

One-to-one MGI homologs map the human signature to mouse genes. GSE152708 scores use within-sample expression ranks. Separate HC3 models for maximum load and structural stiffness adjust sex, age, body weight and generation. The two endpoints are one test family with Benjamini-Hochberg adjustment. The score is fixed before mouse outcomes are analyzed in this phase, but the mouse data were previously inspected for other scores; this is not a blind cohort.

## Interpretation

A generalizable response requires directional support at both MSC time points after adjustment. Bone relevance additionally requires an adjusted mouse association, but cannot convert cell-culture matrix modulus into whole-bone structural stiffness. Gene/pathway interpretation is undertaken only when the MSC transport criterion is met and remains exploratory. All null results and mapping losses are reported.
