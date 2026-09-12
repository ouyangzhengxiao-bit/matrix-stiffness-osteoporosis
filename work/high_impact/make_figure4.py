#!/usr/bin/env python3
from pathlib import Path
import json
import pandas as pd,numpy as np
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];O=ROOT/'outputs/研究升级/Communications_Biology_v4';D=ROOT/'outputs/研究升级/high_impact_v4/osteoporosis_osteoblast_validation'
s=json.loads((D/'osteoporosis_osteoblast_validation_summary.json').read_text());g=pd.read_csv(D/'GSE156508_stiffness_disease_gene_results.csv.gz');n=pd.read_csv(D/'exact_label_permutation_distribution.csv.gz');sc=pd.read_csv(D/'GSE156508_signature_scores.csv')
plt.rcParams.update({'font.family':'Arial','font.size':8.5,'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(1,3,figsize=(7.2,2.65),gridspec_kw={'width_ratios':[1,1.15,1]})
ax=axs[0];ax.hist(n.abs_correlation,bins=28,color='#9AA1A6',edgecolor='white');ax.axvline(s['primary']['rho'],color='#D87A2C',lw=2);ax.set(xlabel='Permuted Spearman rho',ylabel='Balanced allocations',title='a  Exact donor-label test');ax.text(.04,.94,f"Observed rho={s['primary']['rho']:.4f}\nP={s['primary']['one_sided_exact_label_permutation_p']:.4f}",transform=ax.transAxes,va='top')
ax=axs[1];hb=ax.hexbin(g.abs_stiffness_meta_z,g.abs_disease_t,gridsize=36,mincnt=1,cmap='Blues',bins='log');ax.set(xlabel='Absolute stiffness meta z',ylabel='Absolute fracture–OA Welch t',title='b  Distributed magnitude association');fig.colorbar(hb,ax=ax,label='Gene density',fraction=.05,pad=.03)
ax=axs[2];groups=['osteoarthritis','osteoporotic_fracture'];labels=['Severe OA','Osteoporotic\nfracture'];colors=['#2C6EAA','#D87A2C']
for i,(gr,lab,c) in enumerate(zip(groups,labels,colors)):
 y=sc.loc[sc.group==gr,'signature_score'].to_numpy();j=np.linspace(-.08,.08,len(y));ax.scatter(i+j,y,color=c,s=24,zorder=3);ax.plot([i-.2,i+.2],[y.mean(),y.mean()],color='black',lw=1.4)
ax.set_xticks([0,1],labels);ax.set_ylabel('Fixed signature score');ax.set_title('c  No fixed-signature separation');ax.text(.5,.96,f"Exact two-sided P={s['fixed_signature']['two_sided_exact_label_permutation_p']:.3f}",ha='center',va='top',transform=ax.transAxes)
for i,ax in enumerate(axs):ax.title.set_fontweight('bold');ax.title.set_fontsize(9)
fig.tight_layout()
for ext in ['png','pdf','tif']:
 kw={'dpi':320,'bbox_inches':'tight'} if ext!='pdf' else {'bbox_inches':'tight'}
 if ext=='tif':kw['pil_kwargs']={'compression':'tiff_lzw'}
 fig.savefig(O/f'Figure_4.{ext}',**kw)
plt.close(fig)
