# Frozen protocol: multi-donor primary osteoblast disease validation v1

Frozen: 2026-09-09, before inspection of gene-level disease effects or integration results.

## Dataset

GSE156508, primary osteoblasts grown from femoral-head explants: six women with osteoporotic hip fracture and six women with severe osteoarthritis. The comparison is interpreted as fracture-associated osteoporosis versus an osteoarthritis surgical control, not versus healthy bone.

## Processing

Use the deposited normalized series matrix. Map GPL16686 transcript-cluster identifiers to official symbols with Bioconductor `hugene20sttranscriptcluster.db` 8.8.0 and collapse multiple probes per symbol by the median. Retain genes shared with the stiffness meta-analysis.

## Primary test

Spearman correlation between absolute stiffness meta z score and absolute Welch disease t statistic. Use an exact one-sided sample-label permutation over all C(12,6)=924 balanced allocations. This asks whether genes most responsive to matrix stiffness are also the genes that most distinguish the two clinical osteoblast groups, without treating genes as independent replicates.

## Secondary tests

1. Two-sided exact label-permutation test for signed Spearman correlation between pooled stiffness effect and Welch disease t statistic.
2. Exact two-sided label-permutation test of a prespecified sample-level 200-gene stiffness signature score (100 positive minus 100 negative genes; within-sample ranks).
3. Exact one-sided label-permutation tests for excess mean absolute disease t statistic in the frozen GSE147390 osteoblast marker program and in the Chai spatial osteoblast-lineage union, each compared with its complement among shared genes.

The two program tests use Holm correction as one family. Secondary results do not override a null primary test. Missing genes are reported. No threshold is selected after looking at disease outcomes.
