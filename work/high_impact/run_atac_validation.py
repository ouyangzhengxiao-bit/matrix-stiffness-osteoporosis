#!/usr/bin/env python3
import gzip
import json
import re
import sys
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent/'pydeps'))

import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "work/high_impact/data"
META_DIR = ROOT / "outputs/研究升级/meta_analysis_v2"
OUT = ROOT / "outputs/研究升级/high_impact_v4/atac_validation"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(20260908)
N_PERM = 20000

SAMPLES = {
    "senescent_1": DATA / "GSM8747707_senescent-1.bw",
    "senescent_2": DATA / "GSM8747708_senescent-2.bw",
    "force_plus_1": DATA / "GSM8747709_Force_+_-1.bw",
    "force_plus_2": DATA / "GSM8747710_Force_+_-2.bw",
}


def attrs(s):
    return dict(re.findall(r'(\S+) "([^"]+)"', s))


def bh(values):
    p = np.asarray(values, float)
    order = np.argsort(p)
    q = np.empty(len(p), float)
    q[order] = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def read_genes(gtf):
    rows = []
    with gzip.open(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9 or f[2] != "gene":
                continue
            a = attrs(f[8])
            if a.get("gene_type") != "protein_coding":
                continue
            chrom = f[0].removeprefix("chr")
            if chrom not in {str(i) for i in range(1, 23)} | {"X"}:
                continue
            start, end = int(f[3]) - 1, int(f[4])
            tss = start if f[6] == "+" else end - 1
            rows.append({"gene": a["gene_name"], "gene_id": a["gene_id"].split(".")[0],
                         "chrom": chrom, "strand": f[6], "start": start, "end": end,
                         "tss": tss, "gene_length": end - start})
    g = pd.DataFrame(rows)
    # GENCODE may contain duplicated symbols; retain the longest annotated locus deterministically.
    return g.sort_values(["gene", "gene_length"], ascending=[True, False]).drop_duplicates("gene")


def promoter_signals(genes):
    frame = genes.copy()
    headers = {}
    for name, path in SAMPLES.items():
        bw = pyBigWig.open(str(path))
        chroms = bw.chroms()
        headers[name] = bw.header()
        values = []
        for r in frame.itertuples(index=False):
            a = max(0, r.tss - 2000)
            b = min(chroms[r.chrom], r.tss + 2001)
            v = bw.stats(r.chrom, a, b, type="mean", exact=True)[0]
            values.append(0.0 if v is None or not np.isfinite(v) else float(v))
        bw.close()
        frame[name] = values
        frame[name + "_rank"] = pd.Series(np.log1p(frame[name])).rank(method="average", pct=True).to_numpy()
        denom = max(float(headers[name]["sumData"]), 1.0)
        frame[name + "_scaled"] = np.log1p(frame[name] / denom * 1e9)
    frame["baseline_rank"] = frame[["senescent_1_rank", "senescent_2_rank"]].mean(axis=1)
    frame["force_plus_rank"] = frame[["force_plus_1_rank", "force_plus_2_rank"]].mean(axis=1)
    frame["atac_delta"] = frame.force_plus_rank - frame.baseline_rank
    frame["baseline_scaled"] = frame[["senescent_1_scaled", "senescent_2_scaled"]].mean(axis=1)
    frame["force_plus_scaled"] = frame[["force_plus_1_scaled", "force_plus_2_scaled"]].mean(axis=1)
    frame["atac_delta_scaled"] = frame.force_plus_scaled - frame.baseline_scaled
    return frame, headers


def add_strata(df):
    x = df.copy()
    x["baseline_bin"] = pd.qcut(x.baseline_rank.rank(method="first"), 10, labels=False)
    x["length_bin"] = pd.qcut(x.gene_length.rank(method="first"), 10, labels=False)
    x["stratum"] = x.baseline_bin.astype(str) + "_" + x.length_bin.astype(str)
    return x


def matched_null(df, group_genes, forbidden, value="atac_delta", n=N_PERM):
    members = df[df.gene.isin(group_genes)]
    candidates = df[~df.gene.isin(forbidden)]
    counts = members.stratum.value_counts().to_dict()
    pools = {s: candidates.loc[candidates.stratum == s, value].to_numpy() for s in counts}
    null = np.zeros(n, float)
    replacements = 0
    for s, k in counts.items():
        pool = pools[s]
        replace = len(pool) < k
        replacements += int(replace)
        draws = np.vstack([RNG.choice(pool, k, replace=replace) for _ in range(n)])
        null += draws.sum(axis=1)
    null /= max(len(members), 1)
    return null, len(members), replacements


def one_sided_greater(obs, null):
    return float((1 + np.sum(null >= obs)) / (len(null) + 1))


genes = read_genes(DATA / "gencode.v49.annotation.gtf.gz")
signal, headers = promoter_signals(genes)
meta = pd.read_csv(META_DIR / "gene_meta_analysis.csv").drop_duplicates("gene")
signal = signal.merge(meta[["gene", "pooled_effect", "se_hk", "p_hk", "BH_q"]], on="gene", how="inner")
signal = add_strata(signal)

signature = json.loads((META_DIR / "meta_signature.json").read_text())
positive = set(signature["positive"]) & set(signal.gene)
negative = set(signature["negative"]) & set(signal.gene)
forbidden = positive | negative

pos_obs = float(signal.loc[signal.gene.isin(positive), "atac_delta"].mean())
neg_obs = float(signal.loc[signal.gene.isin(negative), "atac_delta"].mean())
directional = pos_obs - neg_obs
pos_null, n_pos, repl_pos = matched_null(signal, positive, forbidden)
neg_null, n_neg, repl_neg = matched_null(signal, negative, forbidden)
directional_null = pos_null - neg_null
primary_p = one_sided_greater(directional, directional_null)

secondary = pd.DataFrame([
    {"test": "positive_vs_matched", "n": n_pos, "observed_mean_delta": pos_obs,
     "null_mean": float(pos_null.mean()), "effect": pos_obs - float(pos_null.mean()),
     "p_one_sided_expected": one_sided_greater(pos_obs, pos_null)},
    {"test": "negative_vs_matched", "n": n_neg, "observed_mean_delta": neg_obs,
     "null_mean": float(neg_null.mean()), "effect": neg_obs - float(neg_null.mean()),
     "p_one_sided_expected": float((1 + np.sum(neg_null <= neg_obs)) / (N_PERM + 1))},
])
secondary["BH_q"] = bh(secondary.p_one_sided_expected)
secondary.to_csv(OUT / "secondary_signature_tests.csv", index=False)

markers = pd.read_csv(DATA / "MarkerGeneList.csv")
cluster_rows = []
marker_forbidden = set(markers.loc[markers.Marker == "Enriched", "Gene"])
for cl, d in markers[markers.Marker == "Enriched"].groupby("cluster"):
    gs = set(d.Gene) & set(signal.gene)
    obs = float(signal.loc[signal.gene.isin(gs), "atac_delta"].mean())
    null, nn, replacements = matched_null(signal, gs, marker_forbidden)
    # Exploratory two-sided competitive test.
    centered = obs - float(null.mean())
    p = float((1 + np.sum(np.abs(null - null.mean()) >= abs(centered))) / (N_PERM + 1))
    cluster_rows.append({"cluster": cl, "n_genes": nn, "mean_delta": obs,
                         "matched_effect": centered, "permutation_p_two_sided": p,
                         "strata_with_replacement": replacements})
clusters = pd.DataFrame(cluster_rows)
clusters["BH_q"] = bh(clusters.permutation_p_two_sided)
clusters.to_csv(OUT / "cluster_atac_tests.csv", index=False)

sample_cols = list(SAMPLES)
qc_rows = []
for i, a in enumerate(sample_cols):
    for b in sample_cols[i + 1:]:
        rho, p = spearmanr(np.log1p(signal[a]), np.log1p(signal[b]))
        qc_rows.append({"sample_a": a, "sample_b": b, "spearman_rho": rho, "p": p})
qc = pd.DataFrame(qc_rows)
qc.to_csv(OUT / "sample_qc_correlations.csv", index=False)

rho, corr_p = spearmanr(signal.pooled_effect, signal.atac_delta)
scaled_directional = (signal.loc[signal.gene.isin(positive), "atac_delta_scaled"].mean()
                      - signal.loc[signal.gene.isin(negative), "atac_delta_scaled"].mean())
foxo = signal.loc[signal.gene == "FOXO1", ["gene", "atac_delta", "atac_delta_scaled"]].to_dict("records")

signal.to_csv(OUT / "gene_promoter_atac_results.csv.gz", index=False, compression="gzip")
summary = {
    "annotation": "GENCODE v49 comprehensive gene annotation, GRCh38, Ensembl 115",
    "n_protein_coding_genes_with_meta": int(len(signal)),
    "primary": {"n_positive": n_pos, "n_negative": n_neg, "positive_mean_delta": pos_obs,
                "negative_mean_delta": neg_obs, "directional_score": directional,
                "matched_null_mean": float(directional_null.mean()),
                "matched_null_sd": float(directional_null.std(ddof=1)),
                "one_sided_permutation_p": primary_p,
                "positive_strata_with_replacement": repl_pos,
                "negative_strata_with_replacement": repl_neg},
    "sensitivity_total_signal_scaled_directional_score": float(scaled_directional),
    "genomewide_spearman": {"rho": float(rho), "p": float(corr_p)},
    "within_group_qc": {
        "senescent_rho": float(qc.query("sample_a == 'senescent_1' and sample_b == 'senescent_2'").spearman_rho.iloc[0]),
        "force_plus_rho": float(qc.query("sample_a == 'force_plus_1' and sample_b == 'force_plus_2'").spearman_rho.iloc[0])},
    "nonzero_fraction": {s: float((signal[s] > 0).mean()) for s in SAMPLES},
    "FOXO1_descriptive": foxo,
    "bigwig_headers": headers,
    "decision": "support" if primary_p < 0.05 and directional > 0 and
                qc.query("sample_a == 'senescent_1' and sample_b == 'senescent_2'").spearman_rho.iloc[0] >= 0.70 and
                qc.query("sample_a == 'force_plus_1' and sample_b == 'force_plus_2'").spearman_rho.iloc[0] >= 0.70
                else ("trend" if directional > 0 else "not_support")
}
(OUT / "atac_validation_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
