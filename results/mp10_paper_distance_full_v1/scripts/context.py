from distance import *
p=pd.read_csv(tab/"input_manifest.tsv",sep="\t")
rows=[];variants=[];support=[]
for s in p.section_id:
    d=pd.read_parquet(out/"spots"/(s+".parquet"))
    d=d[d["query"]].copy()
    for t in targets:
        x=d[t+"_P80_all_distance_um"]
        for field in ["lesion_location","distance_band_200um","accepted_label"]:
            for k,g in d.groupby(field):
                v=x[g.index];rr,rho=corr(v,g[raw]) if v.notna().sum()>=20 else (np.nan,np.nan)
                rows.append(dict(section_id=s,target=t,axis=field,level=k,n_query=len(g),n_used=v.notna().sum(),distance_median=v.median(),n_zero=v.eq(0).sum(),pearson_r=rr,spearman_rho=rho))
        for k in ["common_anchor","E0"]:
            g=d[d[k]];n_used=x[g.index].notna().sum();rr,rho=corr(x[g.index],g[raw]) if n_used>=20 else (np.nan,np.nan);variants.append(dict(section_id=s,target=t,query_subset=k,n_used=n_used,reason=("fewer_than_20_pairs" if n_used<20 else "constant_or_unavailable" if not np.isfinite(rr) else "ok"),pearson_r=rr,spearman_rho=rho))
        support.append(dict(section_id=s,target=t,n_query_no_boundary_with_target_distance=(d.signed_distance_um.isna()&x.notna()).sum(),n_without_component_target=d[t+"_P80_all_reason"].eq("no_target_in_component").sum()))
save(rows,"context_distance.tsv");save(variants,"anchor_subset.tsv");save(support,"distance_coverage.tsv")
gc=pd.read_csv(root/"spatial_ecology_top20_v2/00_inputs/gene_coverage.tsv",sep="\t")
save(gc[(gc.MP=="EpiMP10")&gc.section_id.isin(p.section_id)],"raw_MP10_coverage.tsv")

