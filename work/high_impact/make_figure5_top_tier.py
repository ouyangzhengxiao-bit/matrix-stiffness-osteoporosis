#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
V5=ROOT/'outputs/研究升级/high_impact_v5'
OUT=V5/'figures';OUT.mkdir(parents=True,exist_ok=True)
t=pd.read_csv(V5/'temporal_transportability/chronological_transportability_results.csv')
d=pd.read_csv(V5/'multidonor_msc_atlas/donor_scores.csv')

plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.spines.top':False,
                     'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42})
fig,ax=plt.subplots(1,3,figsize=(7.4,2.55),gridspec_kw={'wspace':.52})

# a: chronological transportability
z=t[t['mode']=='fixed_earliest_three'].copy(); x=np.arange(len(z))
colors=['#B24A4A' if v<0 else '#3A78A8' for v in z.spearman_rho]
ax[0].bar(x,z.spearman_rho,color=colors,width=.68)
ax[0].axhline(0,color='black',lw=.7)
ax[0].set_xticks(x,z.held_out_study.str.replace('GSE',''),rotation=40,ha='right')
ax[0].set_ylabel(r'Held-out Spearman $\rho$')
ax[0].set_xlabel('Later GEO study')
ax[0].set_title('Fixed earliest-three model',fontsize=9)
ax[0].text(.02,.96,'Median = −0.050',transform=ax[0].transAxes,va='top',fontsize=7)

# b: directional agreement
ax[1].bar(x,z.top200_directional_agreement,color='#6B8E6B',width=.68)
ax[1].axhline(.5,color='black',lw=.8,ls='--')
ax[1].set_ylim(0,1)
ax[1].set_xticks(x,z.held_out_study.str.replace('GSE',''),rotation=40,ha='right')
ax[1].set_ylabel('Top-200 direction agreement')
ax[1].set_xlabel('Later GEO study')
ax[1].set_title('Chronological transfer',fontsize=9)
ax[1].text(.02,.96,'Median = 0.493',transform=ax[1].transAxes,va='top',fontsize=7)

# c: donor-unit independent atlas
order=['Bone marrow','Adipose','Dermis','Umbilical cord']
palette={'Bone marrow':'#C75B39','Adipose':'#D8A329','Dermis':'#4C8F69','Umbilical cord':'#5678A6'}
rng=np.random.default_rng(20260919)
for i,g in enumerate(order):
    vals=d.loc[d.tissue==g,'absolute_top200_minus_matched_background'].to_numpy()
    xx=i+np.linspace(-.08,.08,len(vals))
    ax[2].scatter(xx,vals,s=25,color=palette[g],edgecolor='white',lw=.5,zorder=3)
    ax[2].plot([i-.18,i+.18],[vals.mean(),vals.mean()],color='black',lw=1.3)
ax[2].axhline(0,color='#777777',lw=.7)
ax[2].set_xticks(range(4),['BM','Adipose','Dermis','UC'])
ax[2].set_ylabel('Matched expression-rank difference')
ax[2].set_title('Independent 11-donor atlas',fontsize=9)
ax[2].text(.02,.96,'Bone vs others: exact $P$=0.218\nFour tissues: exact $P$=0.00299',
           transform=ax[2].transAxes,va='top',fontsize=6.6)

for i,a in enumerate(ax):
    a.text(-.2,1.08,chr(97+i),transform=a.transAxes,fontweight='bold',fontsize=11)
fig.savefig(OUT/'Figure_5.png',dpi=600,bbox_inches='tight')
fig.savefig(OUT/'Figure_5.tif',dpi=600,bbox_inches='tight',pil_kwargs={'compression':'tiff_lzw'})
fig.savefig(OUT/'Figure_5.pdf',bbox_inches='tight')
plt.close(fig)
t.to_csv(OUT/'Figure_5a_b_source_data.csv',index=False)
d.to_csv(OUT/'Figure_5c_source_data.csv',index=False)
print(OUT)
