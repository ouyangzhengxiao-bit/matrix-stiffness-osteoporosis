from pathlib import Path
import numpy as np, pandas as pd, json
from scipy.stats import ttest_1samp, false_discovery_control

DATA=Path("work/upgrade/data"); OUT=Path("outputs/研究升级/meta_analysis_v2")
long=pd.read_csv(OUT/"study_level_gene_effects.csv.gz")
studies=sorted(long.study.unique())
effects={s:g.set_index("gene").effect for s,g in long.groupby("study")}
sets=[]
for p in [DATA/"h.all.v2025.1.Hs.symbols.gmt",DATA/"c2.cp.reactome.v2025.1.Hs.symbols.gmt"]:
    with open(p) as f:
        for line in f:
            q=line.rstrip().split("\t")
            sets.append((p.name.split(".")[0],q[0],set(x.upper() for x in q[2:])))
rows=[]; detail=[]
for coll,name,gs in sets:
    zs=[]
    for s in studies:
        e=effects[s]; mask=e.index.isin(gs); m=int(mask.sum()); n=len(e)
        if m<15 or m>500 or n-m<15: continue
        ranks=e.rank().to_numpy(); rank_sum=ranks[mask].sum()
        expected=m*(n+1)/2; sd=np.sqrt(m*(n-m)*(n+1)/12)
        z=(rank_sum-expected)/sd
        zs.append(z); detail.append({"collection":coll,"pathway":name,"study":s,"n_genes":m,"enrichment_z":z})
    if len(zs)<4: continue
    stat,p=ttest_1samp(zs,0)
    signs=np.sign(zs)
    rows.append({"collection":coll,"pathway":name,"k_studies":len(zs),"mean_study_z":np.mean(zs),
                 "sd_study_z":np.std(zs,ddof=1),"t":stat,"p_study_unit":p,
                 "same_direction_fraction":max(np.mean(signs>0),np.mean(signs<0))})
res=pd.DataFrame(rows);res["BH_q_study_unit"]=false_discovery_control(res.p_study_unit)
res["loo_direction_stable"]=False
for i in res.index[res.BH_q_study_unit<.05]:
    zs=pd.DataFrame(detail)
    z=zs[(zs.collection==res.at[i,"collection"])&(zs.pathway==res.at[i,"pathway"])].enrichment_z.to_numpy()
    full=np.sign(z.mean()); res.at[i,"loo_direction_stable"]=all(np.sign(np.delete(z,j).mean())==full for j in range(len(z)))
pooled=pd.read_csv(OUT/"pathway_meta_analysis.csv")[["collection","pathway","BH_q"]].rename(columns={"BH_q":"BH_q_pooled_gene_rank"})
res=res.merge(pooled,on=["collection","pathway"],how="left")
res["robust_pathway"]=(res.BH_q_study_unit<.05)&(res.BH_q_pooled_gene_rank<.05)&res.loo_direction_stable&(res.k_studies>=5)
res=res.sort_values("p_study_unit")
res.to_csv(OUT/"pathway_study_unit_meta.csv",index=False)
pd.DataFrame(detail).to_csv(OUT/"pathway_study_effects.csv.gz",index=False)
summary={"pathways_tested":len(res),"study_unit_fdr":int((res.BH_q_study_unit<.05).sum()),
         "robust_pathways":int(res.robust_pathway.sum()),
         "robust_names":res.loc[res.robust_pathway,"pathway"].tolist()}
(OUT/"pathway_robustness_summary.json").write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2));print(res.head(25).to_string(index=False))
