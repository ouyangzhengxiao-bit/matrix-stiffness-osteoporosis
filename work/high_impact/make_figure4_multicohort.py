#!/usr/bin/env python3
from pathlib import Path
import numpy as np,pandas as pd,matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];D=ROOT/'outputs/研究升级/high_impact_v4/multicohort_osteoporosis_validation';O=ROOT/'outputs/研究升级/Communications_Biology_v4'
res=pd.read_csv(D/'cohort_results.csv'); names=res.cohort.tolist(); labels=['Primary osteoblasts\nGSE156508','BM-MSCs\nGSE35958','Femur tissue\nGSE230665']
null={n:pd.read_csv(D/f'{n}_exact_nulls.csv.gz') for n in names}
plt.rcParams.update({'font.size':9,'font.family':'DejaVu Sans','axes.linewidth':0.8})
fig,axs=plt.subplots(1,3,figsize=(10.2,3.25),constrained_layout=True)
# a cohort absolute correlations against exact-null intervals
ax=axs[0];y=np.arange(3)[::-1]
for i,(n,lab) in enumerate(zip(names,labels)):
 v=null[n]['abs'].to_numpy();lo,hi=np.quantile(v,[.025,.975]);obs=res.loc[res.cohort==n,'abs_rho'].iloc[0];yy=y[i]
 ax.plot([lo,hi],[yy,yy],color='#6b7280',lw=3,solid_capstyle='round');ax.scatter(obs,yy,s=45,color='#d95f02',zorder=3);ax.text(max(hi,obs)+.002,yy,f"P={res.loc[res.cohort==n,'abs_exact_one_sided_p'].iloc[0]:.3f}",va='center',fontsize=8)
ax.axvline(0,color='black',lw=.7,ls=':');ax.set_yticks(y,labels);ax.set_xlabel('Absolute-effect Spearman $\\rho$');ax.set_title('a  Cohort-level exact tests',loc='left',fontweight='bold');ax.set_xlim(-.045,.065)
# b combined exact-null Monte Carlo
rng=np.random.default_rng(20260909);B=500000;draw=np.zeros(B);zs=[]
for n in names:
 v=null[n]['abs'].to_numpy();mu=v.mean();sd=v.std(ddof=1);obs=res.loc[res.cohort==n,'abs_rho'].iloc[0];zs.append((obs-mu)/sd);draw+=(rng.choice(v,B)-mu)/sd
draw/=3;obs=np.mean(zs);ax=axs[1];ax.hist(draw,bins=80,color='#9ecae1',edgecolor='none',density=True);ax.axvline(obs,color='#d95f02',lw=2);ax.text(obs+.05,ax.get_ylim()[1]*.88,f'Observed={obs:.2f}\nP=0.259',color='#a33f00',fontsize=8);ax.set_xlabel('Equal-cohort standardized statistic');ax.set_ylabel('Density');ax.set_title('b  Cross-cohort synthesis',loc='left',fontweight='bold')
# c signed and signature exact P
ax=axs[2];x=np.arange(3);w=.34;ps=res.signed_exact_two_sided_p.to_numpy();pg=res.signature_exact_two_sided_p.to_numpy();ax.bar(x-w/2,-np.log10(ps),w,label='Signed effects',color='#4c78a8');ax.bar(x+w/2,-np.log10(pg),w,label='Fixed signature',color='#f2a65a');ax.axhline(-np.log10(.05),color='black',lw=.8,ls='--');ax.set_xticks(x,['GSE156508','GSE35958','GSE230665'],rotation=30,ha='right');ax.set_ylabel('$-\\log_{10}$(exact P)');ax.set_title('c  Directional transfer',loc='left',fontweight='bold');ax.legend(frameon=False,fontsize=8)
for ext in ['png','pdf','tif']:
 fig.savefig(O/f'Figure_4.{ext}',dpi=400 if ext!='pdf' else None,bbox_inches='tight')
plt.close(fig)
