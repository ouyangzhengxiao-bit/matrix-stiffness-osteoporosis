#!/usr/bin/env python3
import itertools,json,hashlib
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata,spearmanr
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/研究升级/high_impact_v4/osteoporosis_osteoblast_validation';OUT.mkdir(parents=True,exist_ok=True)
META=ROOT/'outputs/研究升级/meta_analysis_v2'
expr=pd.read_csv(OUT/'GSE156508_gene_expression.csv').drop_duplicates('gene').set_index('gene')
expr.index=expr.index.astype(str).str.upper();expr=expr.apply(pd.to_numeric,errors='coerce')
meta=pd.read_csv(META/'gene_meta_analysis.csv').drop_duplicates('gene').set_index('gene');meta.index=meta.index.str.upper();meta['meta_z']=meta.pooled_effect/meta.se_hk
common=expr.index.intersection(meta.index);x=expr.loc[common].to_numpy(float);m=meta.loc[common]
alloc=list(itertools.combinations(range(12),6));obs_hi=tuple(range(6,12));obs_i=alloc.index(obs_hi)

def welch(a,hi):
 mask=np.zeros(12,bool);mask[list(hi)]=True
 z1=a[:,mask];z0=a[:,~mask]
 den=np.sqrt(z1.var(1,ddof=1)/6+z0.var(1,ddof=1)/6)
 return (z1.mean(1)-z0.mean(1))/den
stiff=m.meta_z.to_numpy(float);stiff_abs_rank=rankdata(np.abs(stiff));stiff_rank=rankdata(stiff)
# Fixed programs before outcomes.
gong=pd.read_excel(ROOT/'work/high_impact/data/independent_bone/Gong2021_SuppTable2.xlsx',header=1)
gong.columns=['p_val','avg_FC','p_val_adj','gene'];gong['gene']=gong.gene.astype(str).str.upper()
gong_set=set(gong.loc[(pd.to_numeric(gong.p_val_adj,errors='coerce')<.05)&(pd.to_numeric(gong.avg_FC,errors='coerce')>1.2),'gene'])
sp=pd.read_excel(ROOT/'work/high_impact/data/NatureGenetics_2026_supp_tables_2_13.xlsx',sheet_name='2d').rename(columns={'GENE':'gene','CLUSTER NUMBER AND CELL TYPE':'cluster'})
sp.gene=sp.gene.astype(str).str.upper();sp.cluster=sp.cluster.astype(str).str.upper();sp_set=set(sp.loc[sp.cluster.isin({'0. MSCS','9. OSTEO-MSCS','1. OSTEOBLASTS'}),'gene'])
gong_mask=np.array([g in gong_set for g in common]);sp_mask=np.array([g in sp_set for g in common])
rows=[];null_abs=[];null_signed=[];null_g=[];null_sp=[]
for i,hi in enumerate(alloc):
 t=welch(x,hi);ok=np.isfinite(t)&np.isfinite(stiff)
 ar=rankdata(np.abs(t[ok]));sr=rankdata(t[ok])
 null_abs.append(np.corrcoef(stiff_abs_rank[ok],ar)[0,1])
 null_signed.append(np.corrcoef(stiff_rank[ok],sr)[0,1])
 for mask,dest in [(gong_mask,null_g),(sp_mask,null_sp)]:
  q=ok&mask;bg=ok&~mask;dest.append(np.mean(np.abs(t[q]))-np.mean(np.abs(t[bg])))
 if i==obs_i: obs_t=t
null_abs=np.array(null_abs);null_signed=np.array(null_signed);null_g=np.array(null_g);null_sp=np.array(null_sp)
obs_abs=null_abs[obs_i];obs_signed=null_signed[obs_i];obs_g=null_g[obs_i];obs_sp=null_sp[obs_i]
p_abs=float(np.mean(null_abs>=obs_abs-1e-15));p_signed=float(np.mean(np.abs(null_signed)>=abs(obs_signed)-1e-15))
p_g=float(np.mean(null_g>=obs_g-1e-15));p_sp=float(np.mean(null_sp>=obs_sp-1e-15))
# Holm for two program tests.
ps=[p_g,p_sp];order=np.argsort(ps);adj=[0,0];running=0
for rank,idx in enumerate(order):running=max(running,(2-rank)*ps[idx]);adj[idx]=min(1,running)
# Fixed 200-gene signature, sample as inference unit.
sig=json.loads((META/'meta_signature.json').read_text());up=expr.index.intersection([g.upper() for g in sig['positive']]);dn=expr.index.intersection([g.upper() for g in sig['negative']])
ranks=expr.rank(axis=0,pct=True);score=ranks.loc[up].mean()-ranks.loc[dn].mean();sv=score.to_numpy()
def diff(hi):
 mask=np.zeros(12,bool);mask[list(hi)]=True;return sv[mask].mean()-sv[~mask].mean()
null_sig=np.array([diff(a) for a in alloc]);obs_sig=null_sig[obs_i];p_sig=float(np.mean(np.abs(null_sig)>=abs(obs_sig)-1e-15))
res=pd.DataFrame({'gene':common,'stiffness_meta_z':stiff,'disease_welch_t':obs_t,'abs_stiffness_meta_z':np.abs(stiff),'abs_disease_t':np.abs(obs_t)})
res.to_csv(OUT/'GSE156508_stiffness_disease_gene_results.csv.gz',index=False)
pd.DataFrame({'sample':expr.columns,'group':['osteoarthritis']*6+['osteoporotic_fracture']*6,'signature_score':sv}).to_csv(OUT/'GSE156508_signature_scores.csv',index=False)
pd.DataFrame({'permutation':range(len(alloc)),'abs_correlation':null_abs,'signed_correlation':null_signed,'gong_program_effect':null_g,'spatial_program_effect':null_sp,'signature_group_difference':null_sig}).to_csv(OUT/'exact_label_permutation_distribution.csv.gz',index=False)
summary={'source':'GSE156508','design':'six osteoporotic hip-fracture versus six severe-osteoarthritis primary femoral-head osteoblast donors','shared_genes':len(common),'exact_balanced_allocations':len(alloc),'primary':{'test':'Spearman absolute stiffness meta z versus absolute disease Welch t','rho':float(obs_abs),'one_sided_exact_label_permutation_p':p_abs},'secondary_signed':{'rho':float(obs_signed),'two_sided_exact_label_permutation_p':p_signed},'fixed_signature':{'positive_genes':len(up),'negative_genes':len(dn),'fracture_minus_OA_score':float(obs_sig),'two_sided_exact_label_permutation_p':p_sig},'program_disease_magnitude':[{'program':'GSE147390_osteoblast_markers','mapped_genes':int(gong_mask.sum()),'effect_vs_complement':float(obs_g),'one_sided_exact_p':p_g,'Holm_p':adj[0]},{'program':'Chai_spatial_osteoblast_union','mapped_genes':int(sp_mask.sum()),'effect_vs_complement':float(obs_sp),'one_sided_exact_p':p_sp,'Holm_p':adj[1]}],'interpretation_guardrail':'Osteoporotic fracture versus severe osteoarthritis; not healthy controls and not causal.'}

# Leave-one-donor-out audit of the primary absolute correlation.
lodo=[]
groups=np.array([0]*6+[1]*6)
for omit in range(12):
 keep=np.arange(12)!=omit;xx=x[:,keep];gg=groups[keep];nhi=int(gg.sum());obs_hi2=tuple(np.where(gg==1)[0]);alloc2=list(itertools.combinations(range(11),nhi));oi=alloc2.index(obs_hi2)
 vals=[]
 for hi in alloc2:
  mask=np.zeros(11,bool);mask[list(hi)]=True
  z1=xx[:,mask];z0=xx[:,~mask];den=np.sqrt(z1.var(1,ddof=1)/z1.shape[1]+z0.var(1,ddof=1)/z0.shape[1]);tt=(z1.mean(1)-z0.mean(1))/den;ok=np.isfinite(tt)&np.isfinite(stiff)
  vals.append(np.corrcoef(rankdata(np.abs(stiff[ok])),rankdata(np.abs(tt[ok])))[0,1])
 vals=np.array(vals);lodo.append({'omitted_sample':expr.columns[omit],'omitted_group':'fracture' if groups[omit] else 'osteoarthritis','rho':float(vals[oi]),'exact_p':float(np.mean(vals>=vals[oi]-1e-15)),'allocations':len(vals)})
pd.DataFrame(lodo).to_csv(OUT/'leave_one_donor_out_primary.csv',index=False)
summary['leave_one_donor_out']={'all_rho_positive':bool(all(z['rho']>0 for z in lodo)),'rho_min':min(z['rho'] for z in lodo),'rho_max':max(z['rho'] for z in lodo),'significant_at_0.05':sum(z['exact_p']<.05 for z in lodo)}

(OUT/'osteoporosis_osteoblast_validation_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
manifest={}
for p in [ROOT/'work/high_impact/data/GSE156508_series_matrix.txt.gz',ROOT/'work/high_impact/osteoporosis_osteoblast_validation_protocol_v1.md']:
 manifest[str(p.relative_to(ROOT))]={'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
(OUT/'input_protocol_manifest.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(summary,indent=2))
