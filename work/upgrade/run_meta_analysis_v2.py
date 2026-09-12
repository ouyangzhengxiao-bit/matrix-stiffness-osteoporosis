from pathlib import Path
import gzip, io, itertools, json, hashlib, re, sys
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import t as tdist, mannwhitneyu, false_discovery_control, spearmanr
import statsmodels.api as sm

sys.path.insert(0, "work/upgrade")
from run_cross_experiment_v1 import gencode_transcript_map, read_gse22011, score_samples
from run_prediction_v1 import prepare_data

def read_series_matrix(path):
    with gzip.open(path, "rt") as handle:
        text = handle.read()
    table = text.split("!series_matrix_table_begin\n")[1].split("!series_matrix_table_end")[0]
    frame = pd.read_csv(io.StringIO(table), sep="\t", index_col=0)
    frame.index = frame.index.astype(str)
    return frame

def collapse_by_symbol(expression, mapping):
    mapping = mapping.dropna().drop_duplicates()
    unique = mapping.groupby("PROBEID").SYMBOL.agg(lambda x: sorted(set(x)))
    unique = unique[unique.map(len) == 1].map(lambda x: x[0])
    expression = expression.copy()
    expression["gene"] = unique.reindex(expression.index)
    return expression.dropna(subset=["gene"]).groupby("gene").median()

DATA = Path("work/upgrade/data")
OUT = Path("outputs/研究升级/meta_analysis_v2")
OUT.mkdir(parents=True, exist_ok=True)

def full_human_rat_map():
    hi = pd.read_csv(DATA/"Homo_sapiens.gene_info.gz", sep="\t", dtype=str, low_memory=False)[["GeneID","Symbol"]]
    ri = pd.read_csv(DATA/"Rattus_norvegicus.gene_info.gz", sep="\t", dtype=str, low_memory=False)[["GeneID","Symbol"]]
    pairs=[]
    with gzip.open(DATA/"gene_orthologs.complete.gz","rt") as f:
        next(f)
        for line in f:
            tax,g,rel,otax,og=line.rstrip().split("\t")
            if tax=="9606" and otax=="10116": pairs.append((g,og))
    p=pd.DataFrame(pairs,columns=["human_id","rat_id"]).drop_duplicates()
    hc=p.groupby("human_id").rat_id.nunique(); rc=p.groupby("rat_id").human_id.nunique()
    p=p[p.human_id.map(hc).eq(1)&p.rat_id.map(rc).eq(1)]
    p["human"]=p.human_id.map(hi.set_index("GeneID").Symbol)
    p["rat"]=p.rat_id.map(ri.set_index("GeneID").Symbol)
    p=p.dropna()[["human","rat"]].drop_duplicates()
    p=p[~p.human.str.upper().duplicated(False)&~p.rat.str.upper().duplicated(False)]
    p["human"]=p.human.str.upper(); p["rat"]=p.rat.str.upper()
    p.to_csv(OUT/"human_rat_one_to_one_full.csv",index=False)
    return p.set_index("rat").human

HR = full_human_rat_map()

def prep_expr(x, kind="continuous", species="human"):
    x=x.apply(pd.to_numeric,errors="coerce").dropna(how="all")
    if kind=="counts":
        keep=(x.div(x.sum(axis=0),axis=1)*1e6 >= 1).sum(axis=1)>=2
        x=np.log2(x.loc[keep].div(x.loc[keep].sum(axis=0),axis=1)*1e6+0.5)
    elif kind=="linear":
        x=x.loc[(x>1).sum(axis=1)>=2]
        x=np.log2(x+0.5)
    if species=="rat":
        idx=pd.Series(x.index.astype(str).str.upper(),index=x.index).map(HR)
        x=x.loc[idx.notna()].copy(); x.index=idx[idx.notna()].values
    else:
        x.index=x.index.astype(str).str.upper()
    return x.groupby(level=0).median()

def unpaired_effect(x, high, low):
    a=x.iloc[:,list(high)].to_numpy(float); b=x.iloc[:,list(low)].to_numpy(float)
    n1,n0=a.shape[1],b.shape[1]; df=n1+n0-2
    sp=np.sqrt(((n1-1)*a.var(1,ddof=1)+(n0-1)*b.var(1,ddof=1))/df)
    J=1-3/(4*df-1)
    g=J*(a.mean(1)-b.mean(1))/sp
    v=(n1+n0)/(n1*n0)+g*g/(2*df)
    ok=np.isfinite(g)&np.isfinite(v)&(v>0)
    return pd.DataFrame({"effect":g[ok],"variance":v[ok]},index=x.index[ok])

def paired_effect(x, high, low):
    d=x.iloc[:,list(high)].to_numpy(float)-x.iloc[:,list(low)].to_numpy(float)
    n=d.shape[1]; J=1-3/(4*(n-1)-1)
    g=J*d.mean(1)/d.std(1,ddof=1)
    v=1/n+g*g/(2*n)
    ok=np.isfinite(g)&np.isfinite(v)&(v>0)
    return pd.DataFrame({"effect":g[ok],"variance":v[ok]},index=x.index[ok])

def aggregate_study(parts):
    genes=sorted(set().union(*[set(z.index) for z in parts]))
    ys=pd.concat([z.effect for z in parts],axis=1)
    vs=pd.concat([z.variance for z in parts],axis=1)
    m=ys.notna().sum(1)
    effect=ys.mean(1)
    within=vs.mean(1)
    between=ys.var(1,ddof=1).fillna(0)
    variance=within+between
    out=pd.DataFrame({"effect":effect,"variance":variance,"n_contrasts":m}).dropna()
    return out

def load_gse193021():
    rows=[]
    with gzip.open(DATA/"GSE193021_family.soft.gz","rt",errors="replace") as f:
        for line in f:
            if line.startswith("!platform_table_begin"):
                h=next(f).rstrip().split("\t"); ii=h.index("ID"); ai=h.index("gene_assignment")
                for row in f:
                    if row.startswith("!platform_table_end"): break
                    q=row.rstrip().split("\t")
                    if len(q)<=ai: continue
                    sy=set()
                    for a in q[ai].split(" /// "):
                        p=a.split(" // ")
                        if len(p)>1 and re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]*",p[1]): sy.add(p[1])
                    if len(sy)==1: rows.append((q[ii],sy.pop()))
                break
    mp=pd.DataFrame(rows,columns=["PROBEID","SYMBOL"])
    x=prep_expr(collapse_by_symbol(read_series_matrix(DATA/"GSE193021_series_matrix.txt.gz"),mp))
    return aggregate_study([unpaired_effect(x,range(0,3),range(3,6)),unpaired_effect(x,range(6,9),range(9,12))])

def load_gse181512():
    z=pd.read_csv(DATA/"GSE181512_gene_fpkm.txt.gz",sep="\t")
    x=z.set_index("gene_name")[[f"Soft_{i}" for i in range(1,4)]+[f"Stiff_{i}" for i in range(1,4)]]
    x=prep_expr(x,"linear","rat")
    return aggregate_study([unpaired_effect(x,range(3,6),range(0,3))])

def load_gse226411():
    z=pd.read_csv(DATA/"GSE226411_rpkm.txt.gz",sep="\t")
    x=prep_expr(z.set_index("symbol").iloc[:,:6],"linear")
    return aggregate_study([unpaired_effect(x,range(3,6),range(0,3))])

def load_gse288678():
    z=pd.read_excel(DATA/"GSE288678_gene_fpkm.xlsx")
    cols=[c for c in z if c.startswith("A")]
    x=prep_expr(z.set_index("gene_name")[cols],"linear")
    parts=[]
    for start in [0,4,8]: parts.append(unpaired_effect(x,range(12+start,16+start),range(start,start+4)))
    return aggregate_study(parts)

def load_gse55867():
    mp=pd.read_csv(DATA/"hugene10stprobeset_symbol.tsv",sep="\t",dtype=str)
    x=prep_expr(collapse_by_symbol(read_series_matrix(DATA/"GSE55867_series_matrix.txt.gz"),mp))
    return aggregate_study([paired_effect(x,[1,4,7],[2,5,8])])

def load_gse255574():
    tm=gencode_transcript_map(); arr=[]
    meta=[]
    for p in sorted((DATA/"GSE255574_DMSO").glob("*.txt.gz")):
        m=re.search(r"_(150|500|t20)-(24|48)-d_(\d)-3",p.name)
        stiff={"150":150,"500":500,"t20":2000}[m.group(1)]; tmh=int(m.group(2))
        z=pd.read_csv(p,sep="\t",usecols=["Name","NumReads"]); z["gene"]=tm.reindex(z.Name).to_numpy()
        s=z.dropna(subset=["gene"]).groupby("gene").NumReads.sum(); s.name=p.name.split("_")[0]
        arr.append(s); meta.append((s.name,stiff,tmh))
    x=prep_expr(pd.concat(arr,axis=1).fillna(0),"counts")
    md=pd.DataFrame(meta,columns=["sample","stiffness","time"]).set_index("sample").loc[x.columns]
    parts=[]
    for tmh in [24,48]:
        lo=np.where((md.time==tmh)&(md.stiffness==150))[0]
        hi=np.where((md.time==tmh)&(md.stiffness==2000))[0]
        parts.append(unpaired_effect(x,hi,lo))
    return aggregate_study(parts)

def load_gse310514():
    z=pd.read_csv(DATA/"GSE310514_RawCountMatrix.csv.gz")
    mapz=pd.read_csv(DATA/"GSE181512_gene_fpkm.txt.gz",sep="\t",usecols=["gene_id","gene_name"]).drop_duplicates("gene_id").set_index("gene_id").gene_name
    z["symbol"]=mapz.reindex(z.Geneid).to_numpy()
    x=prep_expr(z.dropna(subset=["symbol"]).set_index("symbol").drop(columns=["Geneid"]),"counts","rat")
    parts=[]
    for prefix in ["P1_","P3_","P5_"]:
        lo=[i for i,c in enumerate(x.columns) if c.startswith(prefix+"2-")]
        hi=[i for i,c in enumerate(x.columns) if c.startswith(prefix+"100-")]
        parts.append(unpaired_effect(x,hi,lo))
    return aggregate_study(parts),x

LOADERS={"GSE193021":load_gse193021,"GSE181512":load_gse181512,"GSE226411":load_gse226411,
         "GSE288678":load_gse288678,"GSE55867":load_gse55867,"GSE255574":load_gse255574}

def reml(y,v):
    k=len(y)
    def obj(tau):
        w=1/(v+tau); mu=np.sum(w*y)/np.sum(w)
        return .5*(np.log(v+tau).sum()+np.log(w.sum())+np.sum(w*(y-mu)**2))
    upper=max(10,np.var(y)*10)
    opt=minimize_scalar(obj,bounds=(0,upper),method="bounded",options={"xatol":1e-8})
    tau=max(0,opt.x if opt.fun<obj(0)-1e-8 else 0)
    w=1/(v+tau); mu=np.sum(w*y)/np.sum(w); Q=np.sum(w*(y-mu)**2)
    scale=max(1,Q/(k-1)); se=np.sqrt(scale/w.sum()); crit=tdist.ppf(.975,k-1)
    p=2*tdist.sf(abs(mu/se),k-1)
    q0=np.sum((1/v)*(y-np.sum(y/v)/np.sum(1/v))**2)
    i2=max(0,(q0-(k-1))/q0)*100 if q0>0 else 0
    pi=crit*np.sqrt(tau+se*se)
    return mu,se,p,tau,i2,mu-pi,mu+pi,w/w.sum()

def meta_table(studies):
    genes=sorted(set().union(*[set(x.index) for x in studies.values()]))
    rows=[]
    for g in genes:
        obs=[(name,z.at[g,"effect"],z.at[g,"variance"]) for name,z in studies.items() if g in z.index]
        if len(obs)<4: continue
        names=[o[0] for o in obs]; y=np.array([o[1] for o in obs]); v=np.array([o[2] for o in obs])
        mu,se,p,tau,i2,pl,ph,w=reml(y,v)
        rows.append({"gene":g,"k":len(obs),"pooled_effect":mu,"se_hk":se,"ci_low":mu-tdist.ppf(.975,len(obs)-1)*se,
                     "ci_high":mu+tdist.ppf(.975,len(obs)-1)*se,"p_hk":p,"tau2":tau,"I2_percent":i2,
                     "prediction_low":pl,"prediction_high":ph,"max_weight":w.max(),"studies":";".join(names)})
    out=pd.DataFrame(rows)
    out["BH_q"]=false_discovery_control(out.p_hk)
    out=out.sort_values("p_hk")
    return out

def pathway_test(meta,path):
    score=meta.set_index("gene").pooled_effect
    allgenes=set(score.index); rank=score.rank().to_numpy(); genes=score.index.to_numpy()
    rows=[]
    with open(path) as f:
        for line in f:
            q=line.rstrip().split("\t"); name=q[0]; gs=allgenes.intersection(x.upper() for x in q[2:])
            if not 15<=len(gs)<=500: continue
            mask=np.isin(genes,list(gs))
            u,p=mannwhitneyu(rank[mask],rank[~mask],alternative="two-sided")
            expected=mask.sum()*(len(rank)+1)/2
            sd=np.sqrt(mask.sum()*(~mask).sum()*(len(rank)+1)/12)
            z=(rank[mask].sum()-expected)/sd
            rows.append({"collection":Path(path).name.split(".")[0],"pathway":name,"n_genes":mask.sum(),"rank_enrichment_z":z,"p":p})
    return pd.DataFrame(rows)

def external_validation(meta,signature,memx):
    out=[]
    # GSE166824: one-donor pseudobulk, descriptive only
    with gzip.open("work/data/GSE166824_sc.txt.gz","rt",encoding="latin1") as f:
        next(f); labels=np.array(next(f).rstrip("\n").split("\t")[2:]); df=pd.read_csv(f,sep="\t",header=None)
    x=df.iloc[:,2:].to_numpy(float); valid=labels!=""; x=x[:,valid]; labels=labels[valid]
    genes=df.iloc[:,1].astype(str).str.upper().to_numpy(); norm=np.log1p(x/x.sum(0)*10000)
    for day in ["D3","D6"]:
        hi=np.array([s.startswith("25kPa") and day in s for s in labels]); lo=np.array([s.startswith("2kPa") and day in s for s in labels])
        delta=pd.Series(norm[:,hi].mean(1)-norm[:,lo].mean(1),index=genes).groupby(level=0).mean()
        common=meta.gene[meta.gene.isin(delta.index)]
        rho,p=spearmanr(meta.set_index("gene").loc[common].pooled_effect,delta.loc[common])
        out.append({"dataset":"GSE166824","contrast":day+"_stiff_minus_soft","metric":"gene_effect_spearman","estimate":rho,"p_descriptive":p,"n":len(common)})
    # Fibroblast dose response
    fib=read_gse22011(); st=np.log2(np.array([100,400,1600,6400,25600],float)); st-=st.mean()
    slopes=[]
    for d in range(3):
        a=fib.iloc[:,d*5:(d+1)*5].to_numpy()@st/(st@st)
        slopes.append(pd.Series(a,index=fib.index.str.upper()))
    fs=pd.concat(slopes,axis=1).mean(1)
    common=meta.gene[meta.gene.isin(fs.index)]; rho,p=spearmanr(meta.set_index("gene").loc[common].pooled_effect,fs.loc[common])
    out.append({"dataset":"GSE22011","contrast":"fibroblast_dose_response","metric":"gene_effect_spearman","estimate":rho,"p_descriptive":p,"n":len(common)})
    # Mechanical memory scores
    ranks=memx.rank(axis=0,pct=True); up=ranks.index.intersection(signature["positive"]); dn=ranks.index.intersection(signature["negative"])
    scores=ranks.loc[up].mean()-ranks.loc[dn].mean()
    for c,val in scores.items(): out.append({"dataset":"GSE310514","contrast":c,"metric":"meta_signature_score","estimate":val,"p_descriptive":np.nan,"n":len(up)+len(dn)})
    return pd.DataFrame(out)

def bone_validation(signature):
    pheno,logcpm,_=prepare_data()
    hom=pd.read_csv("outputs/研究升级/one_to_one_homology.csv").set_index("human").mouse
    mapped={d:hom.reindex(gs).dropna().tolist() for d,gs in signature.items()}
    score,nup,ndn=score_samples(logcpm,mapped); pheno=pheno.copy(); pheno["score"]=score.reindex(pheno.index)
    rows=[]
    for outcome in ["maximum_load","structural_stiffness"]:
        cols=["score",outcome,"sex_binary","age","weight","generation"]; d=pheno[cols].dropna()
        X=pd.DataFrame(index=d.index)
        for c in ["score","age","weight"]: X[c]=(d[c]-d[c].mean())/d[c].std(ddof=1)
        X["sex_binary"]=d.sex_binary
        X=pd.concat([X,pd.get_dummies(d.generation,prefix="generation",drop_first=True,dtype=float)],axis=1)
        X=sm.add_constant(X,has_constant="add"); y=(d[outcome]-d[outcome].mean())/d[outcome].std(ddof=1)
        fit=sm.OLS(y,X).fit(cov_type="HC3",use_t=True); ci=fit.conf_int().loc["score"]
        rows.append({"dataset":"GSE152708","outcome":outcome,"n":len(d),"beta":fit.params["score"],"ci_low":ci.iloc[0],"ci_high":ci.iloc[1],"p":fit.pvalues["score"],"up_mapped":nup,"down_mapped":ndn})
    # Human GIOP exact permutation
    y=pd.read_csv(DATA/"GSE276529_rsem_count.tsv.gz",sep="\t",index_col=0); y.index=y.index.str.split(".").str[0]
    mp=pd.read_csv("outputs/final_analysis/source_gene_map.csv").drop_duplicates("ensembl").set_index("ensembl").gene
    y["symbol"]=mp.reindex(y.index); y=y.dropna().groupby("symbol").sum(); ranks=y.rank(axis=0,pct=True)
    up=ranks.index.intersection(signature["positive"]); dn=ranks.index.intersection(signature["negative"]); z=ranks.loc[up].mean()-ranks.loc[dn].mean()
    case=np.array([False]*4+[True]*5); obs=z[case].mean()-z[~case].mean(); null=[]
    for ix in itertools.combinations(range(9),5):
        m=np.zeros(9,bool);m[list(ix)]=True;null.append(z.to_numpy()[m].mean()-z.to_numpy()[~m].mean())
    p=np.mean(np.abs(null)>=abs(obs)-1e-15)
    rows.append({"dataset":"GSE276529","outcome":"GIOP_case_control","n":9,"beta":obs,"ci_low":np.nan,"ci_high":np.nan,"p":p,"up_mapped":len(up),"down_mapped":len(dn)})
    out=pd.DataFrame(rows); out["BH_q"]=false_discovery_control(out.p)
    return out

def main():
    studies={}
    for name,loader in LOADERS.items():
        studies[name]=loader(); print(name,studies[name].shape,flush=True)
    studies["GSE310514"],memx=load_gse310514(); print("GSE310514",studies["GSE310514"].shape,flush=True)
    long=[]
    for name,z in studies.items():
        q=z.copy();q.insert(0,"gene",q.index);q.insert(0,"study",name);long.append(q.reset_index(drop=True))
    pd.concat(long).to_csv(OUT/"study_level_gene_effects.csv.gz",index=False)
    meta=meta_table(studies)
    # leave-one-study-out stability for nominal and FDR candidates
    cand=meta.loc[meta.BH_q<.05,"gene"].tolist()
    meta["loo_direction_stable"]=False
    loo=[]
    for gene in cand:
        full=np.sign(meta.set_index("gene").at[gene,"pooled_effect"])
        stable=True
        for omit in studies:
            obs=[(z.at[gene,"effect"],z.at[gene,"variance"]) for n,z in studies.items() if n!=omit and gene in z.index]
            if len(obs)<3: continue
            y=np.array([a for a,b in obs]);v=np.array([b for a,b in obs]);mu,*_=reml(y,v)
            loo.append({"gene":gene,"omitted":omit,"pooled_effect":mu,"same_direction":np.sign(mu)==full})
            stable &= np.sign(mu)==full
        meta.loc[meta.gene==gene,"loo_direction_stable"]=stable
    meta["conserved"]=((meta.BH_q<.05)&(meta.k>=4)&(meta.max_weight<=.5)&meta.loo_direction_stable)
    meta.to_csv(OUT/"gene_meta_analysis.csv",index=False); pd.DataFrame(loo).to_csv(OUT/"gene_leave_one_out.csv",index=False)
    paths=pd.concat([pathway_test(meta,DATA/"h.all.v2025.1.Hs.symbols.gmt"),pathway_test(meta,DATA/"c2.cp.reactome.v2025.1.Hs.symbols.gmt")],ignore_index=True)
    paths["BH_q"]=false_discovery_control(paths.p); paths=paths.sort_values("p"); paths.to_csv(OUT/"pathway_meta_analysis.csv",index=False)
    pos=meta.nlargest(100,"pooled_effect").gene.tolist(); neg=meta.nsmallest(100,"pooled_effect").gene.tolist()
    sig={"positive":pos,"negative":neg}; (OUT/"meta_signature.json").write_text(json.dumps(sig,indent=2))
    external_validation(meta,sig,memx).to_csv(OUT/"external_transcriptomic_validation.csv",index=False)
    bone=bone_validation(sig);bone.to_csv(OUT/"bone_validation.csv",index=False)
    # Study correlation matrix
    wide=pd.concat({n:z.effect for n,z in studies.items()},axis=1)
    wide.corr(method="spearman",min_periods=1000).to_csv(OUT/"study_effect_spearman.csv")
    manifest={}
    for p in list(DATA.glob("GSE*"))+list(DATA.glob("*.gmt")):
        if p.is_file(): manifest[str(p)]={"bytes":p.stat().st_size,"sha256":hashlib.sha256(p.read_bytes()).hexdigest()}
    (OUT/"input_manifest.json").write_text(json.dumps(manifest,indent=2))
    summary={"independent_studies":len(studies),"genes_meta_analyzed":len(meta),"genes_fdr":int((meta.BH_q<.05).sum()),
             "conserved_genes":int(meta.conserved.sum()),"pathways_fdr":int((paths.BH_q<.05).sum()),
             "pathways_fdr_by_collection":paths[paths.BH_q<.05].groupby("collection").size().to_dict(),
             "bone_any_fdr":bool((bone.BH_q<.05).any())}
    (OUT/"analysis_summary.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2));print(meta.head(15).to_string(index=False));print(paths.head(15).to_string(index=False));print(bone.to_string(index=False))

if __name__=="__main__": main()
