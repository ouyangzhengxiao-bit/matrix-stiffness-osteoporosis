#!/usr/bin/env python3
"""Prespecified donor-unit validation in GSE182158."""
from pathlib import Path
import gzip, hashlib, itertools, json, tarfile
import numpy as np
import pandas as pd
from scipy.io import mmread

ROOT = Path(__file__).resolve().parents[2]
TAR = ROOT / "work/high_impact/data/GSE182158/GSE182158_RAW.tar"
META = ROOT / "outputs/研究升级/meta_analysis_v2/gene_meta_analysis.csv"
OUT = ROOT / "outputs/研究升级/high_impact_v5/multidonor_msc_atlas"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(20260919)

TISSUE = {
    "D01": "Dermis", "D02": "Dermis", "D03": "Dermis",
    "B01": "Bone marrow", "B02": "Bone marrow", "B03": "Bone marrow",
    "U01": "Umbilical cord", "U02": "Umbilical cord",
    "A01": "Adipose", "A02": "Adipose", "A03": "Adipose",
}


def donor_name(member):
    return Path(member).name.split("_")[1]


def load_pseudobulk():
    series = []
    with tarfile.open(TAR, "r") as tf:
        names = tf.getnames()
        matrices = sorted(x for x in names if x.endswith("_matrix.mtx.gz"))
        genes = {donor_name(x): x for x in names if x.endswith("_genes.tsv.gz")}
        for mname in matrices:
            donor = donor_name(mname)
            with gzip.GzipFile(fileobj=tf.extractfile(genes[donor])) as fh:
                g = pd.read_csv(fh, sep="\t", header=None, names=["ensembl", "gene"])
            with gzip.GzipFile(fileobj=tf.extractfile(mname)) as fh:
                mat = mmread(fh).tocsr()
            if mat.shape[0] != len(g):
                raise RuntimeError((donor, mat.shape, len(g)))
            counts = np.asarray(mat.sum(axis=1)).ravel()
            z = pd.Series(counts, index=g.gene.astype(str).str.upper(), name=donor)
            z = z.groupby(level=0).sum()
            series.append(z)
            print(donor, mat.shape, int(counts.sum()), flush=True)
    return pd.concat(series, axis=1).fillna(0)


def binary_exact(scores, bone_indices):
    scores = np.asarray(scores, float); n = len(scores)
    obs = scores[list(bone_indices)].mean() - np.delete(scores, bone_indices).mean()
    null = []
    for ix in itertools.combinations(range(n), len(bone_indices)):
        ix = np.array(ix); mask = np.ones(n, bool); mask[ix] = False
        null.append(scores[ix].mean() - scores[mask].mean())
    null = np.asarray(null)
    p = np.mean(np.abs(null) >= abs(obs) - 1e-15)
    return float(obs), float(p), null


def four_group_exact(scores):
    scores = np.asarray(scores, float); total = scores.mean(); labels = list(TISSUE.values())
    sizes = [3, 3, 3, 2]
    observed = sum(n * (scores[np.array(labels) == lab].mean() - total) ** 2
                   for lab, n in zip(["Bone marrow", "Adipose", "Dermis", "Umbilical cord"], sizes))
    vals = []
    allix = set(range(11))
    for b in itertools.combinations(range(11), 3):
        rem1 = sorted(allix - set(b))
        for a in itertools.combinations(rem1, 3):
            rem2 = sorted(set(rem1) - set(a))
            for d in itertools.combinations(rem2, 3):
                u = sorted(set(rem2) - set(d))
                val = (3*(scores[list(b)].mean()-total)**2 +
                       3*(scores[list(a)].mean()-total)**2 +
                       3*(scores[list(d)].mean()-total)**2 +
                       2*(scores[u].mean()-total)**2)
                vals.append(val)
    vals = np.asarray(vals)
    return float(observed), float(np.mean(vals >= observed - 1e-15)), vals


counts = load_pseudobulk()
lib = counts.sum(axis=0)
logcpm = np.log2(counts.div(lib, axis=1) * 1e6 + 0.5)
ranks = logcpm.rank(axis=0, pct=True, method="average")
meta = pd.read_csv(META)
meta["gene"] = meta.gene.astype(str).str.upper()
meta["meta_z"] = meta.pooled_effect / meta.se_hk
eligible = meta[meta.gene.isin(ranks.index)].drop_duplicates("gene")
sets = {
    "absolute_top200": eligible.loc[eligible.meta_z.abs().nlargest(200).index, "gene"].tolist(),
    "positive_top100": eligible.nlargest(100, "meta_z").gene.tolist(),
    "negative_top100": eligible.nsmallest(100, "meta_z").gene.tolist(),
}

rows = []
bone = [i for i, d in enumerate(ranks.columns) if TISSUE[d] == "Bone marrow"]
score_frame = pd.DataFrame({"donor": ranks.columns, "tissue": [TISSUE[x] for x in ranks.columns]})
for label, genes in sets.items():
    scores = ranks.loc[genes].mean(axis=0).to_numpy()
    diff, p, _ = binary_exact(scores, bone)
    four_stat, four_p, _ = four_group_exact(scores)
    rows.append({"gene_set": label, "mapped_genes": len(genes),
                 "bone_minus_other_mean_rank": diff, "exact_two_sided_p": p,
                 "four_tissue_between_ss": four_stat, "four_tissue_exact_p": four_p})
    score_frame[label] = scores

# Holm correction for the two prespecified signed secondary binary tests.
sec = [rows[1]["exact_two_sided_p"], rows[2]["exact_two_sided_p"]]
order = np.argsort(sec); adj = np.empty(2)
running = 0
for rank, ix in enumerate(order):
    running = max(running, sec[ix] * (2-rank)); adj[ix] = min(1, running)
rows[1]["holm_p_signed_tests"] = float(adj[0]); rows[2]["holm_p_signed_tests"] = float(adj[1])

# Prespecified expression-decile-matched sensitivity for the absolute set.
sig = sets["absolute_top200"]
overall = logcpm.mean(axis=1)
bins = pd.qcut(overall.rank(method="first"), 10, labels=False)
background = [g for g in ranks.index if g not in set(sig)]
bin_pools = {b: np.array([g for g in background if bins[g] == b], dtype=object) for b in range(10)}
matched = np.zeros((1000, ranks.shape[1]))
for rep in range(1000):
    chosen = [RNG.choice(bin_pools[int(bins[g])]) for g in sig]
    matched[rep] = ranks.loc[chosen].mean(axis=0).to_numpy()
sens_scores = ranks.loc[sig].mean(axis=0).to_numpy() - matched.mean(axis=0)
sens_diff, sens_p, _ = binary_exact(sens_scores, bone)
sens_four, sens_four_p, _ = four_group_exact(sens_scores)
score_frame["absolute_top200_minus_matched_background"] = sens_scores

pd.DataFrame(rows).to_csv(OUT / "multidonor_gene_set_tests.csv", index=False)
score_frame.to_csv(OUT / "donor_scores.csv", index=False)
pd.DataFrame({"gene_set": [k for k, v in sets.items() for _ in v],
              "gene": [g for v in sets.values() for g in v]}).to_csv(OUT / "frozen_gene_sets.csv", index=False)

primary = rows[0]
digest = hashlib.sha256()
with TAR.open("rb") as fh:
    for chunk in iter(lambda: fh.read(8 * 1024 * 1024), b""):
        digest.update(chunk)
summary = {
    "dataset": "GSE182158", "donors": 11, "cells_reported_by_source": ">130000",
    "primary": primary,
    "matched_background_sensitivity": {
        "bone_minus_other": sens_diff, "exact_two_sided_p": sens_p,
        "four_tissue_between_ss": sens_four, "four_tissue_exact_p": sens_four_p,
        "matched_sets": 1000, "seed": 20260919,
    },
    "interpretation": "bone_marrow_enriched" if primary["bone_minus_other_mean_rank"] > 0 and primary["exact_two_sided_p"] < .05 else "not_bone_marrow_specific",
    "input": {"bytes": TAR.stat().st_size, "sha256": digest.hexdigest()},
}
(OUT / "multidonor_msc_atlas_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(pd.DataFrame(rows).to_string(index=False))
print(json.dumps(summary, indent=2))
