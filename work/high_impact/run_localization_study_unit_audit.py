#!/usr/bin/env python3
import itertools,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import binomtest,ttest_1samp,t
ROOT=Path(__file__).resolve().parents[2];DATA=ROOT/'work/high_impact/data';META=ROOT/'outputs/研究升级/meta_analysis_v2';OUT=ROOT/'outputs/研究升级/high_impact_v4/human_bone_localization';RNG=np.random.default_rng(20260908);N=10000

def programs():
 p=DATA/'NatureGenetics_2026_supp_tables_2_13.xlsx';sc=pd.read_excel(p,sheet_name='2c').rename(columns={'GENE':'gene','CLUSTER NUMBER AND CELL TYPE':'cluster'});sp=pd.read_excel(p,sheet_name='2d').rename(columns={'GENE':'gene','CLUSTER NUMBER AND CELL TYPE':'cluster'});sc.cluster=sc.cluster.astype(str).str.upper();sp.cluster=sp.cluster.astype(str).str.upper();sl={'0. MSCS (CXCL12+GALNT17+)','1. MSCS (CXCL12+STMN2+)','3. PRE-OSTEOBLASTS','12. MATURE OSTEOBLASTS'};pl={'0. MSCS','9. OSTEO-MSCS','1. OSTEOBLASTS'};return set(sc.loc[sc.cluster.isin(sl),'gene'].str.upper()),set(sp.loc[sp.cluster.isin(pl),'gene'].str.upper())

def audit_one(d,genes):
 d=d.copy();d.gene=d.gene.str.upper();d=d.drop_duplicates('gene');d=d[np.isfinite(d.effect)&np.isfinite(d.variance)&(d.variance>0)];present=d.gene.isin(genes);d['vbin']=pd.qcut(d.variance.rank(method='first'),10,labels=False);obs=float(np.abs(d.loc[present,'effect']).mean());null=np.zeros(N)
 for b,k in d.loc[present,'vbin'].value_counts().items():
  pool=np.abs(d.loc[(~present)&(d.vbin==b),'effect'].to_numpy());null += RNG.choice(pool,size=(N,int(k)),replace=True).sum(1)
 null/=int(present.sum());raw=obs-float(null.mean());return int(present.sum()),obs,float(null.mean()),float(null.std(ddof=1)),float(raw/null.std(ddof=1)),raw

def signflip_p(x):
 x=np.asarray(x,float);obs=x.mean();null=np.array([(x*np.array(s)).mean() for s in itertools.product([-1,1],repeat=len(x))]);return float(np.mean(null>=obs-1e-15))

def summarize(table):
 x=table.raw_matched_effect.to_numpy();z=table.standardized_effect.to_numpy();n=len(x);tt=ttest_1samp(x,0,alternative='greater');se=x.std(ddof=1)/np.sqrt(n);crit=t.ppf(.975,n-1)
 return {'n_studies':n,'mean_raw_matched_effect':float(x.mean()),'raw_effect_ci_low':float(x.mean()-crit*se),'raw_effect_ci_high':float(x.mean()+crit*se),'mean_standardized_effect':float(z.mean()),'positive_studies':int((x>0).sum()),'one_sided_exact_sign_p':float(binomtest(int((x>0).sum()),n,.5,alternative='greater').pvalue),'one_sided_exact_signflip_p':signflip_p(x),'sensitivity_one_sided_t_p':float(tt.pvalue)}
sc_genes,sp_genes=programs();effects=pd.read_csv(META/'study_level_gene_effects.csv.gz');rows=[]
for study,d in effects.groupby('study'):
 for source,gs in [('scRNA',sc_genes),('spatial',sp_genes)]:
  n,obs,nm,ns,z,raw=audit_one(d,gs);rows.append({'study':study,'program_source':source,'n_program_genes':n,'observed_mean_abs_effect':obs,'matched_null_mean':nm,'matched_null_sd':ns,'standardized_effect':z,'raw_matched_effect':raw})
result=pd.DataFrame(rows);result.to_csv(OUT/'study_unit_localization_effects.csv',index=False);summary={s:summarize(d) for s,d in result.groupby('program_source')};loo=[]
for s,g in result.groupby('program_source'):
 for omit in g.study:
  y=g.loc[g.study!=omit,'raw_matched_effect'];loo.append({'program_source':s,'omitted_study':omit,'mean_remaining':float(y.mean()),'positive_remaining':int((y>0).sum())})
pd.DataFrame(loo).to_csv(OUT/'study_unit_localization_leave_one_out.csv',index=False);summary['all_leave_one_out_means_positive']={s:bool((pd.DataFrame(loo).query('program_source==@s').mean_remaining>0).all()) for s in ['scRNA','spatial']};summary['post_review_inference']='exact sign-flip test of raw matched study effects; exact sign count and t test retained as robustness summaries';(OUT/'study_unit_localization_audit_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary,indent=2))
