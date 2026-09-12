#!/usr/bin/env python3
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom, norm, spearmanr

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "work/high_impact/data"
META_DIR = ROOT / "outputs/研究升级/meta_analysis_v2"
OUT = ROOT / "outputs/研究升级/high_impact_v4/genetics_validation"
ATAC = ROOT / "outputs/研究升级/high_impact_v4/atac_validation"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(20260908)
N_PERM = 20000


def bh(values):
    p = np.asarray(values, float)
    order = np.argsort(p)
    q = np.empty(len(p), float)
    q[order] = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def matched_null(df, genes, forbidden):
    members = df[df.gene.isin(genes)]
    candidates = df[~df.gene.isin(forbidden)]
    counts = members.stratum.value_counts().to_dict()
    out = np.zeros(N_PERM)
    replacement_strata = 0
    expanded_strata = 0
    for s, k in counts.items():
        pool = candidates.loc[candidates.stratum == s, "ZSTAT"].to_numpy()
        if len(pool) < k:
            lb, sb = map(int, s.split("_"))
            for radius in range(1, 19):
                near = candidates[(candidates.length_bin.sub(lb).abs() +
                                   candidates.snp_bin.sub(sb).abs()) <= radius]
                pool = near.ZSTAT.to_numpy()
                if len(pool) >= k:
                    expanded_strata += 1
                    break
        repl = len(pool) < k
        replacement_strata += int(repl)
        out += np.vstack([RNG.choice(pool, k, replace=repl) for _ in range(N_PERM)]).sum(axis=1)
    return out / len(members), len(members), replacement_strata, expanded_strata


def competitive(df, genes, forbidden, label):
    obs = float(df.loc[df.gene.isin(genes), "ZSTAT"].mean())
    null, n, repl, expanded = matched_null(df, genes, forbidden)
    p = float((1 + np.sum(null >= obs)) / (N_PERM + 1))
    return {"test": label, "n_genes": n, "mean_ZSTAT": obs, "matched_null_mean": float(null.mean()),
            "matched_null_sd": float(null.std(ddof=1)),
            "standardized_effect": float((obs - null.mean()) / null.std(ddof=1)),
            "one_sided_permutation_p": p, "strata_with_replacement": repl,
            "strata_expanded_to_neighbours": expanded}


def hc3_regression(df, member):
    d = df.copy()
    d["signature"] = d.gene.isin(member).astype(float)
    cov = pd.DataFrame({
        "intercept": 1.0,
        "signature": d.signature,
        "log_length": np.log1p(d.gene_length),
        "log_nsnps": np.log1p(d.NSNPS),
        "log_nparam": np.log1p(d.NPARAM),
    })
    chrom = pd.get_dummies(d.CHR.astype(str), prefix="chr", drop_first=True, dtype=float)
    X = pd.concat([cov, chrom], axis=1).to_numpy(float)
    y = d.ZSTAT.to_numpy(float)
    inv = np.linalg.pinv(X.T @ X)
    beta = inv @ X.T @ y
    residual = y - X @ beta
    leverage = np.einsum("ij,jk,ik->i", X, inv, X)
    adj = residual / np.maximum(1 - leverage, 1e-8)
    meat = X.T @ (X * (adj ** 2)[:, None])
    vcov = inv @ meat @ inv
    se = np.sqrt(np.diag(vcov))
    z = beta[1] / se[1]
    return {"signature_beta": float(beta[1]), "HC3_SE": float(se[1]), "z": float(z),
            "one_sided_p": float(norm.sf(z)), "n_genes": int(len(d))}


p = DATA / "NatureGenetics_2026_supp_tables_2_13.xlsx"
gwas = pd.read_excel(p, sheet_name="7e")
gwas = gwas.rename(columns={"GENE.SYMBOL": "gene"})
for c in ["ZSTAT", "NSNPS", "NPARAM", "START", "STOP", "CHR", "P_SNPWISE_MEAN"]:
    gwas[c] = pd.to_numeric(gwas[c], errors="coerce")
gwas = gwas.dropna(subset=["gene", "ZSTAT", "NSNPS", "NPARAM", "START", "STOP", "CHR"])
gwas["gene_length"] = gwas.STOP - gwas.START + 1
gwas = gwas.sort_values(["gene", "ZSTAT"], ascending=[True, False]).drop_duplicates("gene")

meta = pd.read_csv(META_DIR / "gene_meta_analysis.csv").drop_duplicates("gene")
df = gwas.merge(meta[["gene", "pooled_effect", "se_hk", "p_hk", "BH_q", "k"]], on="gene", how="inner")
df = df[(df.gene_length > 0) & (df.NSNPS > 0)].copy()
df["length_bin"] = pd.qcut(df.gene_length.rank(method="first"), 10, labels=False)
df["snp_bin"] = pd.qcut(df.NSNPS.rank(method="first"), 10, labels=False)
df["stratum"] = df.length_bin.astype(str) + "_" + df.snp_bin.astype(str)

sig = json.loads((META_DIR / "meta_signature.json").read_text())
positive = set(sig["positive"]) & set(df.gene)
negative = set(sig["negative"]) & set(df.gene)
union = positive | negative
primary = competitive(df, union, union, "200_gene_stiffness_signature")

secondary = pd.DataFrame([
    competitive(df, positive, union, "positive_stiffness_signature"),
    competitive(df, negative, union, "negative_stiffness_signature"),
])
secondary["BH_q"] = bh(secondary.one_sided_permutation_p)
secondary.to_csv(OUT / "signature_genetics_tests.csv", index=False)

markers = pd.read_csv(DATA / "MarkerGeneList.csv")
all_markers = set(markers.loc[markers.Marker == "Enriched", "Gene"])
cluster_rows = []
for cl, d in markers[markers.Marker == "Enriched"].groupby("cluster"):
    gs = set(d.Gene) & set(df.gene)
    cluster_rows.append(competitive(df, gs, all_markers, str(cl)))
clusters = pd.DataFrame(cluster_rows)
clusters["BH_q"] = bh(clusters.one_sided_permutation_p)
clusters.to_csv(OUT / "cell_state_genetics_tests.csv", index=False)

rho_abs = spearmanr(np.abs(df.pooled_effect), df.ZSTAT)
significant = set(df.loc[df.P_SNPWISE_MEAN < 2.5e-6, "gene"])
overlap = union & significant
M, K, n, x = len(df), len(significant), len(union), len(overlap)
overlap_p = float(hypergeom.sf(x - 1, M, K, n))
regression = hc3_regression(df, union)

atac = pd.read_csv(ATAC / "bidirectional_gene_promoter_results.csv.gz",
                   usecols=["gene", "atac_delta", "soft_delta"])
df = df.merge(atac, on="gene", how="left")
rho_atac_force = spearmanr(df.ZSTAT, np.abs(df.atac_delta), nan_policy="omit")
rho_atac_soft = spearmanr(df.ZSTAT, np.abs(df.soft_delta), nan_policy="omit")

df["in_signature"] = df.gene.isin(union)
df["signature_direction"] = np.where(df.gene.isin(positive), "positive",
                                      np.where(df.gene.isin(negative), "negative", "background"))
df["eBMD_significant"] = df.P_SNPWISE_MEAN < 2.5e-6
candidates = df[df.in_signature].copy().sort_values(["ZSTAT", "p_hk"], ascending=[False, True])
candidates.head(50).to_csv(OUT / "top_stiffness_eBMD_candidates.csv", index=False)
df.to_csv(OUT / "all_gene_genetics_integration.csv.gz", index=False, compression="gzip")

summary = {
    "source": "Nature Genetics 2026 Supplementary Table 7e; UK Biobank eBMD MAGMA gene tests, N=448010",
    "n_merged_genes": int(len(df)),
    "primary": primary,
    "adjusted_regression_sensitivity": regression,
    "secondary_signature_tests": secondary.to_dict("records"),
    "continuous_abs_meta_effect_vs_ZSTAT": {"rho": float(rho_abs.statistic), "p": float(rho_abs.pvalue)},
    "significant_eBMD_overlap": {"universe": M, "eBMD_significant": K, "signature": n,
                                 "overlap": x, "genes": sorted(overlap), "hypergeometric_p": overlap_p},
    "ATAC_absolute_delta_vs_ZSTAT_exploratory": {
        "force_plus": {"rho": float(rho_atac_force.statistic), "p": float(rho_atac_force.pvalue)},
        "soft": {"rho": float(rho_atac_soft.statistic), "p": float(rho_atac_soft.pvalue)}},
    "decision": "support" if primary["standardized_effect"] > 0 and primary["one_sided_permutation_p"] < 0.05
                else ("trend" if primary["standardized_effect"] > 0 else "not_support")
}
(OUT / "genetics_validation_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
