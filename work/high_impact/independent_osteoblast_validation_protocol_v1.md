# Independent human osteoblast validation protocol v1

Frozen: 2026-09-09 Asia/Shanghai, before merging the marker table with stiffness meta-analysis outputs.

## Source and independence

Gong et al., *A systematic dissection of human primary osteoblasts in vivo at single-cell resolution* (Aging, 2021), GEO GSE147390, is independent of Chai et al. and of all seven stiffness studies. Supplementary Table 2 contains differential-expression results for author-identified osteoblast clusters C1/C2 versus contaminating clusters C3-C6. The experiment used freshly isolated femoral-head cells from one patient; donor-level replication is therefore unavailable and will be stated as a limitation.

## Primary gene program

The primary osteoblast program comprises genes in Supplementary Table 2 with source adjusted P < 0.05 and average fold change > 1.2, matching the source paper's threshold for biologically important cluster DEGs. Gene symbols are upper-cased and duplicates removed before merging.

## Primary endpoint and null

The endpoint is the absolute Hartung-Knapp meta z score from the seven-study stiffness synthesis. The observed mean is compared with 20,000 random non-program gene sets matched jointly on: (1) number of contributing stiffness studies, (2) meta-analysis standard-error decile, and (3) GENCODE v49 gene-length decile. The alternative is one-sided, greater enrichment. Report the observed-minus-null effect, standardized matched effect, 95% Monte Carlo null interval, and empirical P=(1+#null>=observed)/(N+1).

## Study-unit robustness

For each of seven GEO series, compare the program's mean absolute study effect with 20,000 background sets matched on study-effect variance decile and GENCODE gene-length decile. The primary cross-study inference is the exact one-sided sign test. Report the one-sided t test only as a sensitivity analysis, together with the mean raw matched effect and its t-based 95% confidence interval.

## Additional robustness

Repeat the study-unit synthesis leaving out each GEO series descriptively. Because the exact sign test with six remaining studies has limited resolution, report all leave-one-out mean effects and sign counts rather than treating each as an independent hypothesis.

## Decision rule

Support requires: (a) primary matched permutation P < 0.05; (b) positive study-unit effects in all seven studies with exact sign P < 0.05; and (c) positive leave-one-study-out mean effects for every omission. If any condition fails, describe the result as partial or unsupported. No threshold or program definition will be changed after outcome inspection.
