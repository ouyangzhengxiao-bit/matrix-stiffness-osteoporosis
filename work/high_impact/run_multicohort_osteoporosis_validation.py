#!/usr/bin/env python3
import itertools,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'outputs/研究升级/high_impact_v4/multicohort_osteoporosis_validation';OUT.mkdir(parents=True,exist_ok=True)
meta=pd.read_csv(ROOT/'outputs/研究升级/meta_analysis_v2/gene_meta_analysis.csv').set_index('gene');meta.index=meta.index.astype(str).str.upper();meta=meta[~meta.index.duplicated(keep='first')];meta['z']=meta.pooled_effect/meta.se_hk
sig=json.loads((ROOT/'outputs/研究升级/meta_analysis_v2/meta_signature.json').read_text())
cohorts=[('GSE156508',ROOT/'outputs/研究升级/high_impact_v4/osteoporosis_osteoblast_validation/GSE156508_gene_expression.csv',6,6),('GSE35958',ROOT/'outputs/final_analysis/GSE35958_gene_expression.csv.gz',4,5),('GSE230665',ROOT/'outputs/final_analysis/GSE230665_gene_expression.csv.gz',3,12)]

def load(path):
 d=pd.read_csv(path);g='gene' if 'gene' in d.columns else d.columns[0];d=d.set_index(g);d.index=d.index.astype(str).str.upper();d=d[~d.index.duplicated(keep='first')];return d.apply(pd.to_numeric,errors='coerce')
def disease_t(x,hi):
 mask=np.zeros(x.shape[1],bool);mask[list(hi)]=1;a=x[:,mask];b=x[:,~mask];den=np.sqrt(a.var(1,ddof=1)/a.shape[1]+b.var(1,ddof=1)/b.shape[1]);return (a.mean(1)-b.mean(1))/den
def score_diff(e,nctrl,ncase):
 r=e.rank(axis=0,pct=True);up=r.index.intersection([x.upper() for x in sig['positive']]);dn=r.index.intersection([x.upper() for x in sig['negative']]);s=(r.loc[up].mean()-r.loc[dn].mean()).to_numpy();obs=tuple(range(nctrl,nctrl+ncase));alloc=list(itertools.combinations(range(len(s)),ncase));vals=[]
 for hi in alloc:
  m=np.zeros(len(s),bool);m[list(hi)]=1;vals.append(s[m].mean()-s[~m].mean())
 return np.array(vals),alloc.index(obs),len(up),len(dn)
rows=[];nulls={};lodo=[]
for name,path,nctrl,ncase in cohorts:
 e=load(path);common=e.index.intersection(meta.index);e=e.loc[common];x=e.to_numpy(float);z=meta.loc[common,'z'].to_numpy(float);base_ok=np.isfinite(z)&np.isfinite(x).all(axis=1);e=e.iloc[np.where(base_ok)[0]];x=x[base_ok];z=z[base_ok];assert x.shape[0]==z.shape[0]==e.shape[0];obs=tuple(range(nctrl,nctrl+ncase));alloc=list(itertools.combinations(range(nctrl+ncase),ncase));oi=alloc.index(obs);a=[];ss=[]
 zr=rankdata(z);zar=rankdata(np.abs(z))
 for hi in alloc:
  t=disease_t(x,hi);ok=np.isfinite(t);a.append(np.corrcoef(zar[ok],rankdata(np.abs(t[ok])))[0,1]);ss.append(np.corrcoef(zr[ok],rankdata(t[ok]))[0,1])
 a=np.array(a);ss=np.array(ss);sg,soi,nup,ndn=score_diff(e,nctrl,ncase)
 row={'cohort':name,'n_control':nctrl,'n_case':ncase,'shared_genes':len(common),'allocations':len(alloc),'abs_rho':a[oi],'abs_exact_one_sided_p':np.mean(a>=a[oi]-1e-15),'signed_rho':ss[oi],'signed_exact_two_sided_p':np.mean(np.abs(ss)>=abs(ss[oi])-1e-15),'signature_difference':sg[soi],'signature_exact_two_sided_p':np.mean(np.abs(sg)>=abs(sg[soi])-1e-15),'signature_up':nup,'signature_down':ndn};rows.append(row);nulls[name]={'abs':a,'signed':ss,'signature':sg}
 # Leave one donor out: direction only, plus exact p.
 groups=np.array([0]*nctrl+[1]*ncase)
 for omit in range(len(groups)):
  keep=np.arange(len(groups))!=omit;xx=x[:,keep];gg=groups[keep];nc=int(gg.sum());obs2=tuple(np.where(gg==1)[0]);al2=list(itertools.combinations(range(len(gg)),nc));oj=al2.index(obs2);vals=[]
  for hi in al2:
   t=disease_t(xx,hi);ok=np.isfinite(t);vals.append(np.corrcoef(zar[ok],rankdata(np.abs(t[ok])))[0,1])
  vals=np.array(vals);lodo.append({'cohort':name,'omitted_sample':e.columns[omit],'omitted_group':'case' if groups[omit] else 'control','rho':vals[oj],'exact_p':np.mean(vals>=vals[oj]-1e-15),'allocations':len(vals)})
res=pd.DataFrame(rows);res.to_csv(OUT/'cohort_results.csv',index=False);pd.DataFrame(lodo).to_csv(OUT/'leave_one_donor_out.csv',index=False)
rng=np.random.default_rng(20260909);B=500000;combined={}
for endpoint,side in [('abs','greater'),('signed','two-sided'),('signature','two-sided')]:
 obs=[];draw=np.zeros(B)
 for name,_,_,_ in cohorts:
  v=nulls[name][endpoint];mu=v.mean();sd=v.std(ddof=1);o=res.loc[res.cohort==name,{'abs':'abs_rho','signed':'signed_rho','signature':'signature_difference'}[endpoint]].iloc[0];obs.append((o-mu)/sd);draw+=(rng.choice(v,B,replace=True)-mu)/sd
 observed=float(np.mean(obs));draw/=len(cohorts);p=float((1+(draw>=observed).sum())/(B+1)) if side=='greater' else float((1+(np.abs(draw)>=abs(observed)).sum())/(B+1));combined[endpoint]={'equal_weight_mean_standardized_statistic':observed,'monte_carlo_p':p,'draws':B,'alternative':side,'cohort_standardized_statistics':obs}
summary={'cohorts':res.to_dict('records'),'combined':combined,'leave_one_donor_out':{name:{'all_positive':bool((pd.DataFrame(lodo).query('cohort==@name').rho>0).all()),'rho_min':float(pd.DataFrame(lodo).query('cohort==@name').rho.min()),'rho_max':float(pd.DataFrame(lodo).query('cohort==@name').rho.max())} for name,_,_,_ in cohorts},'guardrails':['GSE156508 comparator is severe osteoarthritis, not healthy bone.','GSE230665 has only three controls.','Cell source differs across cohorts.','Magnitude association is not directional transfer or diagnostic accuracy.']}
(OUT/'multicohort_osteoporosis_validation_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
for name,d in nulls.items():pd.DataFrame(d).to_csv(OUT/f'{name}_exact_nulls.csv.gz',index=False)
print(json.dumps(summary,indent=2))
