#!/usr/bin/env python3
import itertools,json,sys
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata,spearmanr,false_discovery_control
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/upgrade'))
from run_meta_analysis_v2 import reml
OUT=ROOT/'outputs/研究升级/high_impact_v4/human_only_sensitivity';OUT.mkdir(parents=True,exist_ok=True)
LONG=ROOT/'outputs/研究升级/meta_analysis_v2/study_level_gene_effects.csv.gz'
HUMAN=['GSE193021','GSE226411','GSE288678','GSE55867','GSE255574']
long=pd.read_csv(LONG);long=long[long.study.isin(HUMAN)]
rows=[]
for gene,d in long.groupby('gene',sort=True):
 if len(d)<3: continue
 mu,se,p,tau,i2,pl,ph,w=reml(d.effect.to_numpy(),d.variance.to_numpy())
 rows.append({'gene':gene,'k':len(d),'pooled_effect':mu,'se_hk':se,'p_hk':p,'tau2':tau,'I2_percent':i2,'prediction_low':pl,'prediction_high':ph,'max_weight':w.max(),'studies':';'.join(d.study)})
hm=pd.DataFrame(rows);hm['BH_q']=false_discovery_control(hm.p_hk);hm=hm.sort_values('p_hk');hm.to_csv(OUT/'human_only_gene_meta.csv',index=False)
full=pd.read_csv(ROOT/'outputs/研究升级/meta_analysis_v2/gene_meta_analysis.csv')[['gene','pooled_effect','se_hk']].rename(columns={'pooled_effect':'full_effect','se_hk':'full_se'})
comp=hm.merge(full,on='gene');rho,p=spearmanr(comp.pooled_effect,comp.full_effect)
sig={'positive':hm.nlargest(100,'pooled_effect').gene.str.upper().tolist(),'negative':hm.nsmallest(100,'pooled_effect').gene.str.upper().tolist()}
(OUT/'human_only_signature.json').write_text(json.dumps(sig,indent=2),encoding='utf-8')
cohorts=[('GSE156508',ROOT/'outputs/研究升级/high_impact_v4/osteoporosis_osteoblast_validation/GSE156508_gene_expression.csv',6,6),('GSE35958',ROOT/'outputs/final_analysis/GSE35958_gene_expression.csv.gz',4,5),('GSE230665',ROOT/'outputs/final_analysis/GSE230665_gene_expression.csv.gz',3,12)]
def load(path):
 d=pd.read_csv(path);g='gene' if 'gene' in d.columns else d.columns[0];d=d.set_index(g);d.index=d.index.astype(str).str.upper();d=d[~d.index.duplicated(keep='first')];return d.apply(pd.to_numeric,errors='coerce')
def wt(x,hi):
 m=np.zeros(x.shape[1],bool);m[list(hi)]=1;a=x[:,m];b=x[:,~m];den=np.sqrt(a.var(1,ddof=1)/a.shape[1]+b.var(1,ddof=1)/b.shape[1]);return (a.mean(1)-b.mean(1))/den
def score(e):
 r=e.rank(axis=0,pct=True);u=r.index.intersection(sig['positive']);d=r.index.intersection(sig['negative']);return (r.loc[u].mean()-r.loc[d].mean()).to_numpy(),len(u),len(d)
mi=hm.set_index(hm.gene.str.upper());nulls={};cres=[]
for name,path,nc,ncase in cohorts:
 e=load(path);common=e.index.intersection(mi.index);e=e.loc[common];z=(mi.loc[common].pooled_effect/mi.loc[common].se_hk).to_numpy();x=e.to_numpy(float);ok=np.isfinite(z)&np.isfinite(x).all(1);e=e.iloc[np.where(ok)[0]];z=z[ok];x=x[ok];zr=rankdata(z);zar=rankdata(abs(z));alloc=list(itertools.combinations(range(nc+ncase),ncase));obs=alloc.index(tuple(range(nc,nc+ncase)));av=[];sv=[]
 scores,nu,nd=score(e);dv=[]
 for hi in alloc:
  t=wt(x,hi);o=np.isfinite(t);av.append(np.corrcoef(zar[o],rankdata(abs(t[o])))[0,1]);sv.append(np.corrcoef(zr[o],rankdata(t[o]))[0,1]);m=np.zeros(len(scores),bool);m[list(hi)]=1;dv.append(scores[m].mean()-scores[~m].mean())
 av=np.array(av);sv=np.array(sv);dv=np.array(dv);nulls[name]={'abs':av,'signed':sv,'signature':dv}
 cres.append({'cohort':name,'shared_genes':len(z),'allocations':len(alloc),'abs_rho':av[obs],'abs_exact_one_sided_p':np.mean(av>=av[obs]-1e-15),'signed_rho':sv[obs],'signed_exact_two_sided_p':np.mean(abs(sv)>=abs(sv[obs])-1e-15),'signature_difference':dv[obs],'signature_exact_two_sided_p':np.mean(abs(dv)>=abs(dv[obs])-1e-15),'signature_up':nu,'signature_down':nd})
res=pd.DataFrame(cres);res.to_csv(OUT/'human_only_disease_results.csv',index=False)
rng=np.random.default_rng(20260909);B=500000;combined={}
for endpoint,alt,col in [('abs','greater','abs_rho'),('signed','two-sided','signed_rho'),('signature','two-sided','signature_difference')]:
 zs=[];draw=np.zeros(B)
 for name,_,_,_ in cohorts:
  v=nulls[name][endpoint];mu=v.mean();sd=v.std(ddof=1);o=res.loc[res.cohort==name,col].iloc[0];zs.append((o-mu)/sd);draw+=(rng.choice(v,B)-mu)/sd
 draw/=len(cohorts);obs=float(np.mean(zs));pcomb=float((1+(draw>=obs).sum())/(B+1)) if alt=='greater' else float((1+(abs(draw)>=abs(obs)).sum())/(B+1));combined[endpoint]={'statistic':obs,'p':pcomb,'draws':B,'alternative':alt,'cohort_z':zs}
summary={'human_studies':HUMAN,'genes_k_ge_3':len(hm),'genes_fdr':int((hm.BH_q<.05).sum()),'full_vs_human_effect_spearman':rho,'full_vs_human_p_descriptive':p,'cohorts':res.to_dict('records'),'combined':combined}
(OUT/'human_only_sensitivity_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
