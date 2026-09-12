#!/usr/bin/env python3
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]; DATA=ROOT/'work/high_impact/data'; META=ROOT/'outputs/研究升级/meta_analysis_v2'; OUT=ROOT/'outputs/研究升级/high_impact_v4/human_bone_localization'
R=np.random.default_rng(20260908); N=20000
p=DATA/'NatureGenetics_2026_supp_tables_2_13.xlsx'
sc=pd.read_excel(p,sheet_name='2c').rename(columns={'GENE':'gene','CLUSTER NUMBER AND CELL TYPE':'cluster'})
sp=pd.read_excel(p,sheet_name='2d').rename(columns={'GENE':'gene','CLUSTER NUMBER AND CELL TYPE':'cluster'})
sc.cluster=sc.cluster.astype(str).str.upper(); sp.cluster=sp.cluster.astype(str).str.upper()
a=set(sc.loc[sc.cluster.isin({'0. MSCS (CXCL12+GALNT17+)','1. MSCS (CXCL12+STMN2+)','3. PRE-OSTEOBLASTS','12. MATURE OSTEOBLASTS'}),'gene'].str.upper())
b=set(sp.loc[sp.cluster.isin({'0. MSCS','9. OSTEO-MSCS','1. OSTEOBLASTS'}),'gene'].str.upper())
unique=b-a
m=pd.read_csv(META/'gene_meta_analysis.csv'); m=m[(m.se_hk>0)&np.isfinite(m.se_hk)].drop_duplicates('gene'); m.gene=m.gene.str.upper(); m['abs_z']=(m.pooled_effect/m.se_hk).abs(); m['sebin']=pd.qcut(m.se_hk.rank(method='first'),10,labels=False); m['stratum']=m.k.astype(int).astype(str)+'_'+m.sebin.astype(str)
present=m.gene.isin(unique); obs=float(m.loc[present,'abs_z'].mean()); candidates=m[~present]; null=np.zeros(N)
for s,k in m.loc[present,'stratum'].value_counts().items():
 pool=candidates.loc[candidates.stratum==s,'abs_z'].to_numpy(); null += np.vstack([R.choice(pool,k,replace=len(pool)<k) for _ in range(N)]).sum(1)
null/=int(present.sum()); result={'status':'post_hoc_overlap_removed_sensitivity','spatial_total_genes':len(b),'overlap_removed':len(a&b),'spatial_unique_source_genes':len(unique),'spatial_unique_meta_genes':int(present.sum()),'observed_mean_abs_meta_z':obs,'matched_null_mean':float(null.mean()),'effect':float(obs-null.mean()),'standardized_effect':float((obs-null.mean())/null.std(ddof=1)),'one_sided_permutation_p':float((1+np.sum(null>=obs))/(N+1))}
(OUT/'overlap_removed_spatial_sensitivity.json').write_text(json.dumps(result,indent=2),encoding='utf-8'); print(json.dumps(result,indent=2))
