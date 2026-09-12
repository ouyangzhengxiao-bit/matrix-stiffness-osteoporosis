#!/usr/bin/env python3
import gzip, json, re, itertools
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import binomtest, ttest_1samp, t

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'work/high_impact/data'
META=ROOT/'outputs/研究升级/meta_analysis_v2'
OUT=ROOT/'outputs/研究升级/high_impact_v4/independent_osteoblast_validation'
OUT.mkdir(parents=True,exist_ok=True)
N=20000
SEED=20260909

def read_program():
    p=DATA/'independent_bone/Gong2021_SuppTable2.xlsx'
    d=pd.read_excel(p,header=1)
    d.columns=['p_val','avg_FC','p_val_adj','gene']
    for c in ['p_val','avg_FC','p_val_adj']: d[c]=pd.to_numeric(d[c],errors='coerce')
    d['gene']=d.gene.astype(str).str.upper()
    return set(d.loc[(d.p_val_adj<.05)&(d.avg_FC>1.2),'gene'].dropna())

def gene_lengths():
    out={}
    p=DATA/'gencode.v49.annotation.gtf.gz'
    with gzip.open(p,'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            q=line.rstrip().split('\t')
            if len(q)<9 or q[2]!='gene': continue
            m=re.search(r'gene_name "([^"]+)"',q[8])
            if m: out[m.group(1).upper()]=int(q[4])-int(q[3])+1
    return pd.Series(out,name='gene_length')

def add_bins(d, cols):
    x=d.copy()
    for name,source in cols:
        x[name]=pd.qcut(x[source].rank(method='first'),10,labels=False,duplicates='drop').astype(int)
    return x

def matched_null(d, member, value, strata, rng):
    mem=d[d.gene.isin(member)].copy(); bg=d[~d.gene.isin(member)].copy()
    null=np.zeros(N)
    for key,k in mem.groupby(strata,dropna=False).size().items():
        if not isinstance(key,tuple): key=(key,)
        pool=bg.copy()
        exact=np.ones(len(pool),dtype=bool)
        for col,val in zip(strata,key): exact &= pool[col].eq(val).to_numpy()
        arr=pool.loc[exact,value].to_numpy(float)
        if len(arr)<k:
            dist=np.zeros(len(pool),float)
            for col,val in zip(strata,key):
                if col=='k': dist += 2*np.abs(pool[col].to_numpy(float)-float(val))
                else: dist += np.abs(pool[col].to_numpy(float)-float(val))
            take=max(int(k*4),50)
            arr=pool.iloc[np.argsort(dist)[:min(take,len(pool))]][value].to_numpy(float)
        # Sampling with replacement avoids exhausting small exact strata.
        null += rng.choice(arr,size=(N,int(k)),replace=True).sum(axis=1)
    return mem, null/len(mem)

program=read_program()
length=gene_lengths()
meta=pd.read_csv(META/'gene_meta_analysis.csv').drop_duplicates('gene')
meta['gene']=meta.gene.str.upper(); meta['gene_length']=meta.gene.map(length)
meta=meta[np.isfinite(meta.se_hk)&(meta.se_hk>0)&meta.gene_length.notna()].copy()
meta['meta_z']=meta.pooled_effect/meta.se_hk; meta['abs_meta_z']=meta.meta_z.abs()
meta=add_bins(meta,[('se_bin','se_hk'),('length_bin','gene_length')])
rng=np.random.default_rng(SEED)
mem,null=matched_null(meta,program,'abs_meta_z',['k','se_bin','length_bin'],rng)
obs=float(mem.abs_meta_z.mean()); effect=obs-float(null.mean())
primary={'source':'GSE147390_Gong2021','source_program_genes':len(program),'mapped_genes':len(mem),
         'observed_mean_abs_meta_z':obs,'null_mean':float(null.mean()),
         'raw_matched_effect':effect,'standardized_matched_effect':float(effect/null.std(ddof=1)),
         'null_95_low':float(np.quantile(null,.025)),'null_95_high':float(np.quantile(null,.975)),
         'permutation_p':float((1+(null>=obs).sum())/(N+1))}

long=pd.read_csv(META/'study_level_gene_effects.csv.gz')
long['gene']=long.gene.str.upper(); long['gene_length']=long.gene.map(length)
rows=[]
for i,(study,d) in enumerate(sorted(long.groupby('study'))):
    d=d[np.isfinite(d.effect)&np.isfinite(d.variance)&(d.variance>0)&d.gene_length.notna()].drop_duplicates('gene').copy()
    d=add_bins(d,[('vbin','variance'),('length_bin','gene_length')])
    d['abs_effect']=d.effect.abs()
    m,nul=matched_null(d,program,'abs_effect',['vbin','length_bin'],np.random.default_rng(SEED+100+i))
    ob=float(m.abs_effect.mean()); raw=ob-float(nul.mean())
    rows.append({'study':study,'n_program_genes':len(m),'observed_mean_abs_effect':ob,
                 'matched_null_mean':float(nul.mean()),'matched_null_sd':float(nul.std(ddof=1)),
                 'raw_matched_effect':raw,'standardized_effect':float(raw/nul.std(ddof=1))})
res=pd.DataFrame(rows)
x=res.raw_matched_effect.to_numpy(); n=len(x); tm=ttest_1samp(x,0,alternative='greater')
crit=t.ppf(.975,n-1); ci=(float(x.mean()-crit*x.std(ddof=1)/np.sqrt(n)),float(x.mean()+crit*x.std(ddof=1)/np.sqrt(n)))
sign=binomtest(int((x>0).sum()),n,.5,alternative='greater')
null_signflip=np.array([(x*np.array(z)).mean() for z in itertools.product([-1,1],repeat=n)])
signflip_p=float(np.mean(null_signflip>=x.mean()-1e-15))
loo=[]
for omit in res.study:
    y=res.loc[res.study!=omit,'raw_matched_effect']
    loo.append({'omitted_study':omit,'n_remaining':len(y),'mean_raw_matched_effect':float(y.mean()),
                'positive_remaining':int((y>0).sum())})
summary={'primary_gene_set_test':primary,
         'study_unit':{'n_studies':n,'positive_studies':int((x>0).sum()),'mean_raw_matched_effect':float(x.mean()),
                       't_95_ci_low':ci[0],'t_95_ci_high':ci[1],
                       'one_sided_exact_sign_p':float(sign.pvalue),'one_sided_exact_signflip_p':signflip_p,
                       'sensitivity_one_sided_t_p':float(tm.pvalue)},
         'leave_one_study_out_all_means_positive':bool(all(z['mean_raw_matched_effect']>0 for z in loo)),
         'decision':'support' if primary['permutation_p']<.05 and sign.pvalue<.05 and all(z['mean_raw_matched_effect']>0 for z in loo) else 'partial_or_not_support'}
res.to_csv(OUT/'independent_osteoblast_study_effects.csv',index=False)
pd.DataFrame(loo).to_csv(OUT/'independent_osteoblast_leave_one_out.csv',index=False)
pd.DataFrame([primary]).to_csv(OUT/'independent_osteoblast_primary_test.csv',index=False)
(OUT/'independent_osteoblast_validation_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
