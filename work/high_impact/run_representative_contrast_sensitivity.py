#!/usr/bin/env python3
"""Post-review sensitivity using one prespecified representative contrast per GEO series."""
import sys, json, gzip, re
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/upgrade'))
import run_meta_analysis_v2 as m
OUT=ROOT/'outputs/\u7814\u7a76\u5347\u7ea7/high_impact_v4/representative_contrast_sensitivity'
OUT.mkdir(parents=True,exist_ok=True)
original_aggregate=m.aggregate_study
m.aggregate_study=lambda parts: parts[0]
studies={}
for name,loader in m.LOADERS.items(): studies[name]=loader()
studies['GSE310514'],_=m.load_gse310514()
meta=m.meta_table(studies)
meta.to_csv(OUT/'representative_contrast_gene_meta.csv',index=False)
long=[]
for name,z in studies.items():
    q=z.copy(); q.insert(0,'gene',q.index); q.insert(0,'study',name); long.append(q.reset_index(drop=True))
pd.concat(long).to_csv(OUT/'representative_contrast_study_effects.csv.gz',index=False)
primary=pd.read_csv(ROOT/'outputs/\u7814\u7a76\u5347\u7ea7/meta_analysis_v2/gene_meta_analysis.csv')
common=primary[['gene','pooled_effect']].merge(meta[['gene','pooled_effect']],on='gene',suffixes=('_all','_representative'))
rho,p=spearmanr(common.pooled_effect_all,common.pooled_effect_representative)
summary={
 'design':'one representative contrast per GEO series; first ordered time/coating/passage contrast for multi-contrast series',
 'representative_choices':{
  'GSE193021':'first deposited comparison (samples 1-3 vs 4-6)',
  'GSE181512':'only eligible contrast','GSE226411':'only eligible contrast',
  'GSE288678':'first ECM stratum','GSE55867':'only paired contrast',
  'GSE255574':'24-hour 2000 Pa vs 150 Pa','GSE310514':'passage 1, 100 kPa vs 2 kPa'},
 'independent_studies':len(studies),'genes_meta_analyzed':len(meta),
 'genes_fdr':int((meta.BH_q<.05).sum()),'minimum_p':float(meta.p_hk.min()),
 'minimum_q':float(meta.BH_q.min()),'pooled_effect_spearman_vs_all_contrasts':float(rho),
 'pooled_effect_spearman_p':float(p)}
(OUT/'representative_contrast_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
