#!/usr/bin/env python3
import gzip,json,re,itertools
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import binomtest,ttest_1samp,t
ROOT=Path(__file__).resolve().parents[2]; DATA=ROOT/'work/high_impact/data'; META=ROOT/'outputs/研究升级/meta_analysis_v2'; OUT=ROOT/'outputs/研究升级/high_impact_v4/human_bone_localization'; N=20000; SEED=20260909

def lengths():
 o={}
 with gzip.open(DATA/'gencode.v49.annotation.gtf.gz','rt') as f:
  for l in f:
   if l.startswith('#'): continue
   q=l.rstrip().split('\t')
   if len(q)>=9 and q[2]=='gene':
    m=re.search(r'gene_name "([^"]+)"',q[8])
    if m:o[m.group(1).upper()]=int(q[4])-int(q[3])+1
 return pd.Series(o)

def programs():
 p=DATA/'NatureGenetics_2026_supp_tables_2_13.xlsx'
 def rd(sheet):
  d=pd.read_excel(p,sheet_name=sheet).rename(columns={'GENE':'gene','CLUSTER NUMBER AND CELL TYPE':'cluster'});d.gene=d.gene.astype(str).str.upper();d.cluster=d.cluster.astype(str).str.upper();return d
 sc=rd('2c');sp=rd('2d')
 sl={'0. MSCS (CXCL12+GALNT17+)','1. MSCS (CXCL12+STMN2+)','3. PRE-OSTEOBLASTS','12. MATURE OSTEOBLASTS'}
 pl={'0. MSCS','9. OSTEO-MSCS','1. OSTEOBLASTS'}
 return {'scRNA':set(sc.loc[sc.cluster.isin(sl),'gene']),'spatial':set(sp.loc[sp.cluster.isin(pl),'gene'])}

def bins(d,spec):
 d=d.copy()
 for new,old in spec:d[new]=pd.qcut(d[old].rank(method='first'),10,labels=False,duplicates='drop').astype(int)
 return d

def null(d,gs,value,strata,rng):
 mem=d[d.gene.isin(gs)];bg=d[~d.gene.isin(gs)];z=np.zeros(N)
 for key,k in mem.groupby(strata).size().items():
  if not isinstance(key,tuple):key=(key,)
  dist=np.zeros(len(bg))
  exact=np.ones(len(bg),bool)
  for c,v in zip(strata,key):
   exact &= bg[c].eq(v).to_numpy();dist += (2 if c=='k' else 1)*np.abs(bg[c].to_numpy(float)-float(v))
  a=bg.loc[exact,value].to_numpy()
  if len(a)<max(k,20): a=bg.iloc[np.argsort(dist)[:min(max(int(k*4),50),len(bg))]][value].to_numpy()
  z += rng.choice(a,size=(N,int(k)),replace=True).sum(1)
 return mem,z/len(mem)
L=lengths();P=programs();M=pd.read_csv(META/'gene_meta_analysis.csv').drop_duplicates('gene');M.gene=M.gene.str.upper();M['length']=M.gene.map(L);M=M[M.length.notna()&(M.se_hk>0)].copy();M['abs_z']=(M.pooled_effect/M.se_hk).abs();M=bins(M,[('sebin','se_hk'),('lenbin','length')])
meta_rows=[]
for i,(name,gs) in enumerate(P.items()):
 m,n=null(M,gs,'abs_z',['k','sebin','lenbin'],np.random.default_rng(SEED+i));obs=m.abs_z.mean();p=(1+(n>=obs).sum())/(N+1)
 meta_rows.append({'program':name,'n_genes':len(m),'observed':obs,'null_mean':n.mean(),'raw_matched_effect':obs-n.mean(),'standardized_effect':(obs-n.mean())/n.std(ddof=1),'permutation_p':p})
meta_res=pd.DataFrame(meta_rows);meta_res['bonferroni_p']=np.minimum(1,2*meta_res.permutation_p)
D=pd.read_csv(META/'study_level_gene_effects.csv.gz');D.gene=D.gene.str.upper();D['length']=D.gene.map(L);study_rows=[]
for pi,(name,gs) in enumerate(P.items()):
 for si,(study,d) in enumerate(sorted(D.groupby('study'))):
  d=d[d.length.notna()&(d.variance>0)].drop_duplicates('gene').copy();d['abs_effect']=d.effect.abs();d=bins(d,[('vbin','variance'),('lenbin','length')]);m,n=null(d,gs,'abs_effect',['vbin','lenbin'],np.random.default_rng(SEED+100+pi*20+si));obs=m.abs_effect.mean();study_rows.append({'program':name,'study':study,'n_genes':len(m),'observed':obs,'null_mean':n.mean(),'raw_matched_effect':obs-n.mean(),'standardized_effect':(obs-n.mean())/n.std(ddof=1)})
study=pd.DataFrame(study_rows);summ=[];loos=[]
for name,g in study.groupby('program'):
 x=g.raw_matched_effect.to_numpy();n=len(x);tt=ttest_1samp(x,0,alternative='greater');crit=t.ppf(.975,n-1);se=x.std(ddof=1)/np.sqrt(n)
 summ.append({'program':name,'n_studies':n,'positive_studies':int((x>0).sum()),'mean_raw_effect':x.mean(),'ci_low':x.mean()-crit*se,'ci_high':x.mean()+crit*se,'exact_sign_p':binomtest(int((x>0).sum()),n,.5,alternative='greater').pvalue,'exact_signflip_p':float(np.mean(np.array([(x*np.array(z)).mean() for z in itertools.product([-1,1],repeat=n)])>=x.mean()-1e-15)),'sensitivity_t_p':tt.pvalue})
 for omit in g.study:
  y=g.loc[g.study!=omit,'raw_matched_effect'];loos.append({'program':name,'omitted_study':omit,'mean_remaining':y.mean(),'positive_remaining':int((y>0).sum())})
meta_res.to_csv(OUT/'covariate_matched_union_tests.csv',index=False);study.to_csv(OUT/'covariate_matched_study_effects.csv',index=False);pd.DataFrame(summ).to_csv(OUT/'covariate_matched_study_summary.csv',index=False);pd.DataFrame(loos).to_csv(OUT/'covariate_matched_leave_one_out.csv',index=False)
summary={'meta_level':meta_res.to_dict('records'),'study_unit':summ,'all_leave_one_out_means_positive':{n:bool((pd.DataFrame(loos).query('program==@n').mean_remaining>0).all()) for n in P}}
(OUT/'covariate_matched_robustness_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary,indent=2))
