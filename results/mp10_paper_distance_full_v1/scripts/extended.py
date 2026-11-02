"""Same-pair sensitivities, spatial program summaries, and coverage."""
from distance import *
read=lambda name:pd.read_csv(tab/name,sep="\t")
def effect(x,y):
    keep=np.isfinite(x)&np.isfinite(y)
    u,v=x[keep],y[keep]
    why="insufficient_pairs" if len(u)<3 else "constant_distance_or_predictor" if u.nunique()<2 else "constant_score" if v.nunique()<2 else "ok"
    rr,rho=corr(u,v) if why=="ok" else (np.nan,np.nan)
    return dict(n_used=len(u),pearson_r=rr,spearman_rho=rho,reason=why)
def adjusted(d,t):
    dist=d[t+"_P80_all_distance_um"]
    c=pd.DataFrame({"MP10":d[raw],"distance":dist,"epithelial_fraction":d.Epithelial_fraction,
                    "log1p_nCount_RNA":np.log1p(d.nCount_RNA.where(d.nCount_RNA>=0)),
                    "lesion_location":d.lesion_location}).replace([np.inf,-np.inf],np.nan).dropna()
    e=effect(c.distance,c.MP10)
    cols=pd.concat([c[["epithelial_fraction","log1p_nCount_RNA"]],pd.get_dummies(c.lesion_location,prefix="location",drop_first=True,dtype=float)],axis=1)
    used=[k for k in cols if cols[k].nunique()>1]
    design=np.column_stack([np.ones(len(c)),cols[used].to_numpy(dtype=float)])
    reason=e["reason"];partial=np.nan
    if reason=="ok":
        values=c[["distance","MP10"]].to_numpy()
        coefficients,_,rank,_=np.linalg.lstsq(design,values,rcond=None)
        resid=values-design@coefficients
        if len(c)-rank<2:reason="insufficient_residual_df"
        elif np.any(np.linalg.norm(resid,axis=0)<1e-10*np.linalg.norm(values-values.mean(axis=0),axis=0)):reason="constant_residual"
        else:partial=pearsonr(*resid.T).statistic
    return dict(n_used=len(c),unadjusted_paired_r=e["pearson_r"],adjusted_r=partial,actual_covariates=";".join(used) or "intercept_only",location_levels=";".join(sorted(c.lesion_location.unique())),reason=reason)
def main():
    manifest=read("input_manifest.tsv")
    modules=sorted(read("gene_sets.tsv.gz").gs_name.unique())
    adj=[];model=[];fx=[];regional=[];concord=[];tregs=[];covers=[];rep=[]
    for sec in manifest.section_id:
        path=out/"spots"/(sec+".parquet")
        d=pd.read_parquet(path)
        for t in targets:
            adj.append(dict(section_id=sec,target=t,**adjusted(d,t)))
            x=d[t+"_P80_all_distance_um"].where(d.epi_score_valid.eq(True)&d.spot_MP10_model_universe.notna()&d.epi_MP10_model_universe.notna())
            for layer in ("spot_MP10_model_universe","epi_MP10_model_universe"):
                model.append(dict(section_id=sec,target=t,layer=layer,**effect(x,d[layer])))
        scorefile=out/"scores"/(sec+"_programs.tsv.gz")
        old=[k for k in d if k in modules or k.startswith("AMS_")]
        if scorefile.exists():
            score=pd.read_csv(scorefile,sep="\t")
            assert len(score)==len(d) and not score.duplicated(key).any()
            d=d.drop(columns=old).merge(score,on=key,how="left",validate="1:1",indicator=True)
            assert d.pop("_merge").eq("both").all()
            coverage=pd.read_csv(out/"scores"/(sec+"_match.tsv"),sep="\t")
        else:
            for name in modules:d[name]=np.nan
            coverage=read("gene_sets.tsv.gz").groupby("gs_name").gene_symbol.nunique().rename("n_original").reset_index()
            coverage=coverage.assign(section_id=sec,n_matched=np.nan,n_universe=np.nan,db_version="2026.1.Hs",status="missing_score_input")
        coverage["match_fraction"]=coverage.n_matched/coverage.n_original
        covers.append(coverage)
        for name in modules:
            v=d[name]
            if scorefile.exists() and v.notna().any() and v.nunique()==1:
                coverage.loc[coverage.gs_name==name,"status"]="constant"
                d[name]=np.nan;v=d[name]
            for t in targets:
                e=effect(d[t+"_P80_all_distance_um"],v)
                if not v.notna().any(): e["reason"]=str(coverage.loc[coverage.gs_name==name,"status"].iloc[0])
                fx.append(dict(section_id=sec,target=t,program=name,**e))
            if "CHOLESTEROL" in name:
                concord.append(dict(section_id=sec,program=name,**effect(d[raw].where(d["query"]),v)))
            if name.startswith("GSE15659") or "CO_INHIBITION_BY_CTLA4" in name:
                tregs.append(dict(section_id=sec,program=name,**effect(d.Treg_fraction,v)))
        for axis in ["lesion_location","distance_band_200um"]:
            for level,g in d.groupby(axis):
                for population,idx in [("tissue",g.center_qc & g.lesion_location.ne("quality_excluded")),("epi_query",g["query"])]:
                    v=g.loc[idx,modules]
                    q=v.quantile([.25,.5,.75]).T
                    q.columns=["q25","median","q75"]
                    q["n_valid"]=v.count()
                    regional.append(q.rename_axis("program").reset_index().assign(section_id=sec,axis=axis,region=level,population=population,n_population=len(v)))
        d.to_parquet(path,index=False)
        rep.append(dict(section_id=sec,n_sets_scored=coverage.status.eq("scored").sum(),n_sets_constant=coverage.status.eq("constant").sum(),
                        n_depth=d.nCount_RNA.notna().sum(),n_MP10_query_zero=(d["query"]&d[raw].eq(0)).sum()))
    for rows,name in [(adj,"adjusted"),(model,"matched_layers"),(fx,"program_distance"),(pd.concat(regional),"region_programs"),(concord,"cholesterol_MP10"),(tregs,"Treg_program_composition"),(pd.concat(covers),"gene_set_coverage"),(rep,"extended_coverage")]:
        save(rows,name+".tsv")
    z=read("region_niche.tsv")
    z["n_niche_across_regions"]=z.groupby(["section_id","version","axis","niche"]).n.transform("sum")
    z["fraction_of_niche_in_region"]=z.n/z.n_niche_across_regions.where(z.n_niche_across_regions.gt(0))
    z.loc[z.n_labelled.eq(0),"fraction_of_niche_in_region"]=np.nan
    save(z,"region_niche.tsv")
if __name__=="__main__":main()

