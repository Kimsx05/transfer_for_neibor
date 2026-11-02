"""Integrate spot programs; retain per-section associations without pooling spots."""
from figures import *
def load_scores():
    coverage=[]
    for s in secs:
        path=out/"spots"/(s+"_programs.tsv.gz")
        u=pd.read_csv(path,sep="\t")
        d=D[s].drop(columns=[c for c in D[s] if c.startswith(("HALLMARK_","REACTOME_","GSE15659_","AMS_"))])
        d=d.merge(u,on=key,validate="one_to_one")
        assert len(d)==len(D[s]) and all(k in d for k in focus)
        D[s]=d
        d.to_parquet(out/"spots"/(s+".parquet"),index=False)
        coverage.append(pd.read_csv(tab/(s+"_match.tsv"),sep="\t"))
    save(pd.concat(coverage),"gene_set_coverage.tsv")
def tables():
    effects=[];checks=[];consistent=[];state=[]
    for s in secs:
        d=D[s]
        for name in [c for c in d if c.startswith(("HALLMARK_","REACTOME_","GSE15659_"))]:
            v=d[name]
            for t in targets:
                dist=d[t+"_P80_all_distance_um"];rr,rho=corr(dist,v)
                effects.append(dict(section_id=s,target=t,program=name,n=(dist.notna()&v.notna()).sum(),pearson_r=rr,spearman_rho=rho))
            if "CHOLESTEROL" in name:
                rr,rho=corr(d[raw].where(d["query"]),v);consistent.append(dict(section_id=s,program=name,n=d["query"].sum(),pearson_r=rr,spearman_rho=rho))
            if name.startswith("GSE15659") or "CTLA4" in name:
                rr,rho=corr(d.Treg_fraction,v)
                state.append(dict(section_id=s,program=name,n=len(d),pearson_r=rr,spearman_rho=rho))
            if "AMS_"+name in d:
                other=d["AMS_"+name];rr,rho=corr(v,other)
                for t in targets:
                    x=d[t+"_P80_all_distance_um"]
                    checks.append(dict(section_id=s,program=name,target=t,all_tissue_pearson=rr,all_tissue_spearman=rho,ucell_distance_r=corr(x,v)[0],AMS_distance_r=corr(x,other)[0]))
    save(effects,"program_distance.tsv");save(consistent,"cholesterol_MP10.tsv");save(state,"Treg_program_composition.tsv");save(checks,"method_comparison.tsv")
def maps(s):
    d=D[s]
    for page,cols in enumerate([focus[:6],focus[6:12],["Treg_fraction","Treg_abundance",focus[-2],focus[-1],focus[11],focus[6]]]):
        title=["Growth, stress and inflammation","Immune, lipid and matrix programs","Treg context across the full tissue"][page]
        f,a,l=layout(2,3,s+"  |  "+title,"Whole-spot RNA programs (UCell)   •   display scales fixed per feature across pilot sections")
        for j,c in enumerate(cols):
            unit="UCell score";title=short(c)
            if c=="Treg_abundance":unit="Posterior mean estimate";title="Treg estimated abundance"
            if c=="Treg_fraction":unit="Fraction of all 30";title="Treg composition"
            field(a[j],l[j],d,c,title,unit)
        export(f,a,l,s,f"0{5+page}_program_maps")
    if s in secs[:2]:
        f,a,l=layout(3,2,s+"  |  UCell and AddModuleScore","Each row uses one complete standard set   •   method-specific color scales",w=5)
        for i,name in enumerate([focus[9],focus[4],focus[2]]):
            field(a[2*i],l[2*i],d,name,short(name),"UCell")
            c="AMS_"+name
            combined=pd.concat([D[x][c] for x in secs[:2]])
            limits=tuple(combined.quantile([.01,.99]).to_numpy())
            field(a[2*i+1],l[2*i+1],d,c,short(name),"AddModuleScore",limits,cmap="coolwarm")
        export(f,a,l,s,"10_method_comparison")
def trends(s):
    d=D[s]; modules=[focus[i] for i in [4,2,0,1,9,10]]
    f,a,l=layout(2,3,s+"  |  Programs by target distance","Epi queries   •   medians at observed distance bins   •   lines stop where a bin has fewer than 5 spots",h=3,w=4.7)
    for j,c in enumerate(modules):
        ax=a[j]
        for i,t in enumerate(targets):
            x=d[t+"_P80_all_distance_um"].round(6)
            b=pd.cut(x,[-1,0,200,500,1000,np.inf],labels=False)
            g=d[c].groupby(b).agg(["median","count"]).reindex(range(5))
            ax.plot(range(5),g["median"].where(g["count"]>=5),color=colors_target[i],marker="o",ms=4,label=t)
        ax.set_title(short(c));ax.set_xticks(range(5),["0","≤200","≤500","≤1000",">1000"],rotation=15)
        ax.set_xlabel("Distance-bin upper bound (μm)");ax.set_ylabel("Median spot UCell")
        l[j].legend(handles=ax.get_lines(),loc="center",ncol=3,frameon=False,handlelength=1)
    export(f,a,l,s,"08_program_trends")
    z=pd.read_csv(tab/"program_distance.tsv",sep="\t")
    z=z[z.section_id==s].pivot(index="program",columns="target",values="pearson_r")
    f,a,l=layout(1,2,s+"  |  Program–distance associations","All 55 UCell programs   •   one section   •   including zero distances   •   display cap ±0.6",w=8,h=12,left=.27,wspace=1.85)
    for j,mat in enumerate([z.iloc[:28],z.iloc[28:]]):
        mat=mat.copy();mat.index=[short(x) for x in mat.index]
        heat(a[j],l[j],mat,"Whole-spot programs in Epi queries",-.6,.6,"RdBu_r","Pearson r")
    export(f,a,l,s,"09_program_effects")
if __name__=="__main__":
    load_scores();tables()
    for s in secs:print(s,flush=True);maps(s);trends(s)

