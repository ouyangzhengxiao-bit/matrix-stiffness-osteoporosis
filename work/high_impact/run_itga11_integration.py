#!/usr/bin/env python3
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

ROOT = Path(__file__).resolve().parents[2]
INP = ROOT / "work/high_impact/data"
META = ROOT / "outputs/研究升级/meta_analysis_v2"
OUT = ROOT / "outputs/研究升级/high_impact_v4"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(20260908)


def bh(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    q = np.empty(len(p), float)
    q[order] = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


meta = pd.read_csv(META / "gene_meta_analysis.csv")
meta = meta.dropna(subset=["gene", "pooled_effect", "se_hk"]).copy()
meta = meta[(meta.se_hk > 0) & np.isfinite(meta.se_hk)]
meta["z"] = meta.pooled_effect / meta.se_hk
meta = meta.drop_duplicates("gene").set_index("gene")

markers = pd.read_csv(INP / "MarkerGeneList.csv")
background = set(pd.read_csv(INP / "AllGene.csv").gene) & set(meta.index)
markers = markers[markers.Gene.isin(background)].copy()

# Primary competitive test: author-defined "Enriched" markers, five clusters.
rows = []
for cluster, g in markers[markers.Marker == "Enriched"].groupby("cluster"):
    genes = sorted(set(g.Gene))
    inside = meta.loc[genes, "z"].to_numpy()
    outside = meta.loc[list(background - set(genes)), "z"].to_numpy()
    for endpoint, a, b in [
        ("signed_stiffness_z", inside, outside),
        ("absolute_stiffness_z", np.abs(inside), np.abs(outside)),
    ]:
        stat, p = mannwhitneyu(a, b, alternative="two-sided")
        rows.append({"cluster": cluster, "endpoint": endpoint, "n_markers": len(inside),
                     "median_marker": float(np.median(a)), "median_background": float(np.median(b)),
                     "rank_biserial": float(2 * stat / (len(a) * len(b)) - 1), "p": float(p)})
cluster_tests = pd.DataFrame(rows)
cluster_tests["BH_q_within_endpoint"] = cluster_tests.groupby("endpoint")["p"].transform(bh)
cluster_tests.to_csv(OUT / "cluster_competitive_tests.csv", index=False)

# Orthogonal validation in paired ITGA11-high versus ITGA11-low donor cultures.
bulk = pd.read_csv(INP / "Bulk_Counts.txt", sep="\t", low_memory=False).drop_duplicates("Geneid").set_index("Geneid")
joined = meta[["pooled_effect", "z", "p_hk", "BH_q"]].join(
    bulk[["logFC_high_low", "padj_high_low", "logFC_Proli_high_low", "padj_Proli_high_low",
          "logFC_siITGA11_siCTR", "padj_siITGA11_siCTR"]], how="inner")

corr_rows = []
for col in ["logFC_high_low", "logFC_Proli_high_low", "logFC_siITGA11_siCTR"]:
    x = joined[["pooled_effect", col]].dropna()
    rho, p = spearmanr(x.pooled_effect, x[col])
    corr_rows.append({"comparison": col, "n_genes": len(x), "spearman_rho": rho, "p": p})
corr = pd.DataFrame(corr_rows)
corr["BH_q"] = bh(corr.p)
corr.to_csv(OUT / "genomewide_bulk_correlations.csv", index=False)

signature = json.loads((META / "meta_signature.json").read_text())
sig_rows = []
universe = np.array(sorted(set(joined.index)))
for col in ["logFC_high_low", "logFC_Proli_high_low", "logFC_siITGA11_siCTR"]:
    v = joined[col].dropna()
    pos = [g for g in signature["positive"] if g in v.index]
    neg = [g for g in signature["negative"] if g in v.index]
    obs = float(v.loc[pos].mean() - v.loc[neg].mean())
    vals = v.to_numpy()
    npos, nneg = len(pos), len(neg)
    null = np.empty(20000)
    for i in range(len(null)):
        ix = RNG.choice(len(vals), npos + nneg, replace=False)
        null[i] = vals[ix[:npos]].mean() - vals[ix[npos:]].mean()
    p = (1 + np.sum(np.abs(null) >= abs(obs))) / (len(null) + 1)
    sig_rows.append({"comparison": col, "n_positive": npos, "n_negative": nneg,
                     "directional_score": obs, "null_sd": float(null.std(ddof=1)), "permutation_p": p})
sig = pd.DataFrame(sig_rows)
sig["BH_q"] = bh(sig.permutation_p)
sig.to_csv(OUT / "stiffness_signature_in_itga11_bulk.csv", index=False)

# Context map: rank-normalized study effects averaged across each cluster's enriched markers.
study = pd.read_csv(META / "study_level_gene_effects.csv.gz")
study = study[study.gene.isin(background)].copy()
study["rank_z"] = study.groupby("study")["effect"].rank(pct=True).sub(0.5).mul(2)
context_rows = []
for st, ds in study.groupby("study"):
    series = ds.set_index("gene").rank_z
    for cluster, g in markers[markers.Marker == "Enriched"].groupby("cluster"):
        genes = list(set(g.Gene) & set(series.index))
        context_rows.append({"study": st, "cluster": cluster, "n_markers": len(genes),
                             "mean_rank_effect": float(series.loc[genes].mean())})
context = pd.DataFrame(context_rows)
context.to_csv(OUT / "study_by_cell_state_effects.csv", index=False)

summary = {
    "background_genes": len(background),
    "primary_hypothesis": "Higher-stiffness transcriptional effects are concentrated in the ITGA11-high osteogenic cluster 0 state.",
    "cluster0_tests": cluster_tests[cluster_tests.cluster == "Cl_0"].to_dict("records"),
    "bulk_signature_tests": sig.to_dict("records"),
    "genomewide_correlations": corr.to_dict("records"),
    "interpretation_guardrail": "All analyses use public data. Cells are not independent biological replicates; donor or study is the analysis unit where applicable."
}
(OUT / "itga11_integration_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
