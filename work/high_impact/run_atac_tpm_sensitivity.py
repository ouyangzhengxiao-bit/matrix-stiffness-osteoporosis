#!/usr/bin/env python3
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent/'pydeps'))
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
META = ROOT / "outputs/研究升级/meta_analysis_v2"
ATAC = ROOT / "outputs/研究升级/high_impact_v4/atac_validation"
RNG = np.random.default_rng(20260908)
N = 20000


def analyze(df, baseline_cols, treated_cols, label):
    x = df.copy()
    x["baseline"] = np.log1p(x[baseline_cols]).mean(axis=1)
    x["treated"] = np.log1p(x[treated_cols]).mean(axis=1)
    x["delta_tpm_log"] = x.treated - x.baseline
    x["baseline_bin"] = pd.qcut(x.baseline.rank(method="first"), 10, labels=False)
    x["length_bin"] = pd.qcut(x.gene_length.rank(method="first"), 10, labels=False)
    x["stratum"] = x.baseline_bin.astype(str) + "_" + x.length_bin.astype(str)
    sig = json.loads((META / "meta_signature.json").read_text())
    pos, neg = set(sig["positive"]) & set(x.gene), set(sig["negative"]) & set(x.gene)
    forbidden = pos | neg

    def null(gs):
        m, c = x[x.gene.isin(gs)], x[~x.gene.isin(forbidden)]
        z = np.zeros(N)
        for s, k in m.stratum.value_counts().items():
            pool = c.loc[c.stratum == s, "delta_tpm_log"].to_numpy()
            z += np.vstack([RNG.choice(pool, k, replace=len(pool) < k) for _ in range(N)]).sum(axis=1)
        return z / len(m)

    po = float(x.loc[x.gene.isin(pos), "delta_tpm_log"].mean())
    no = float(x.loc[x.gene.isin(neg), "delta_tpm_log"].mean())
    pn, nn = null(pos), null(neg)
    d, dn = po - no, pn - nn
    alternative = "greater" if label == "force_plus" else "less"
    p = ((1 + np.sum(dn >= d)) if alternative == "greater" else (1 + np.sum(dn <= d))) / (N + 1)
    return x[["gene", "delta_tpm_log"]], {
        "comparison": label, "n_positive": len(pos), "n_negative": len(neg),
        "genomewide_mean_log1p_TPM_delta": float(x.delta_tpm_log.mean()),
        "positive_mean": po, "negative_mean": no, "directional_score": d,
        "matched_null_mean": float(dn.mean()), "matched_null_sd": float(dn.std(ddof=1)),
        "expected_alternative": alternative, "one_sided_permutation_p": float(p)}


old = pd.read_csv(ATAC / "gene_promoter_atac_results.csv.gz")
young = pd.read_csv(ATAC / "bidirectional_gene_promoter_results.csv.gz")
force_table, force = analyze(old, ["senescent_1", "senescent_2"],
                             ["force_plus_1", "force_plus_2"], "force_plus")
soft_table, soft = analyze(young, ["young_1", "young_2"],
                           ["force_minus_1", "force_minus_2"], "soft_force_minus")
merged = force_table.rename(columns={"delta_tpm_log": "force_plus_delta_tpm_log"}).merge(
    soft_table.rename(columns={"delta_tpm_log": "soft_delta_tpm_log"}), on="gene", how="outer")
merged.to_csv(ATAC / "posthoc_tpm_gene_results.csv.gz", index=False, compression="gzip")
summary = {"status": "post_hoc_method_sensitivity", "normalization": "author-provided TPM BigWig; log1p",
           "force_plus": force, "soft_force_minus": soft}
(ATAC / "posthoc_tpm_sensitivity_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
