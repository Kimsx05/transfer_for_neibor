"""Run from the project root, using Python 3.12. Only manifest sections are selected."""
import os, sys, shutil
from pathlib import Path
for k in ["OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"]: os.environ[k]="1"
root=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(root/"spatial_ecology_top20_v2/runtime"))
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import spearmanr,pearsonr
out=root/"mp10_paper_distance_pilot_v1"
tab=out/"tables"
key=["section_id","barcode"]
raw="EpiMP10_Top20_AUCell"
targets={"CAF":["Matrix_CAF","iCAF","rCAF","vCAF"],"TAM":["C1QC TAM","SPP1 TAM"],"Treg":["Treg"]}
def save(rows,file):
    pd.DataFrame(rows).to_csv(tab/file,sep="\t",index=False,na_rep="NA")
def subset(path,sections):
    return pd.concat([g[g.section_id.isin(sections)] for g in pd.read_csv(path,sep="\t",chunksize=20000)],ignore_index=True)
def corr(x,y):
    ok=x.notna()&y.notna()
    if x[ok].nunique()<2 or y[ok].nunique()<2 or ok.sum()<3:return np.nan,np.nan
    return pearsonr(x[ok],y[ok]).statistic,spearmanr(x[ok],y[ok]).statistic
def search(d,query,target):
    xy=d[["x_um","y_um"]].to_numpy()
    distance=pd.Series(np.nan,index=d.index)
    reason=pd.Series("not_epi_query",index=d.index)
    reason[query]="missing_component_or_coordinates"
    for _,idx in d[d.tissue_component_id.notna()&np.isfinite(xy).all(1)].groupby("tissue_component_id").groups.items():
        ti=idx[target.loc[idx]]; qi=idx[query.loc[idx]]
        if not len(qi):continue
        if not len(ti):
            reason.loc[qi]="no_target_in_component"
            continue
        distance.loc[qi]=cKDTree(xy[ti]).query(xy[qi],workers=1)[0]
        for i in qi[:2]:assert np.isclose(distance[i],np.linalg.norm(xy[ti]-xy[i],axis=1).min())
        reason.loc[qi]="ok"
    return distance,reason
def main():
    pilots=pd.read_csv(out/"tables/input_manifest.tsv",sep="\t")
    assert len(pilots)==7 and pilots.dataset_id.is_unique
    secs=pilots.section_id
    labels=subset(root/"results/mp10_pathology_regions_v3_4_final/final_spot_labels.tsv.gz",secs)
    feat=subset(root/"mp10_cnv_regions_interface_v1/input_review/position_features.tsv.gz",secs)
    feat["barcode"]=feat.position_id
    maps=pd.read_csv(root/"spatial_ecology_top20_v2/00_inputs/major_map.tsv",sep="\t")
    src=root/"neighborhood_mp10_source_v4/04_niche_comparison/anchor_MP10_all_layers.parquet"
    flags=["spot_MP10_frozen_status","frozen_definition_id","signature_coverage_status","missing_reason","n_signature_available_raw","assignment_status"]
    r12=root/"neighborhood_ring_harmony_v3/03_main/SAMPLE_R12_LEIDEN"
    shutil.copyfile(r12/"profiles_R1_R2.tsv",tab/"source_R12_profiles.tsv")
    cov=[]; effects=[]; regions=[]; niches=[]; counts=[]; bins=[]; localization=[]
    for r in pilots.itertuples():
        sec=r.section_id
        d=labels[labels.section_id==sec].merge(feat[feat.section_id==sec],on=["dataset_id"]+key,validate="1:1")
        assert len(d)==r.n_spots
        nodes=pd.read_parquet(root/"neighborhood_spottedpy_v8/01_geometry"/sec/"nodes.parquet",columns=key+["spot_uid","common_anchor","center_qc","accepted_label","E0","score_eligible"])
        d=d.merge(nodes,on=key,validate="1:1",how="left")
        d=d.merge(pd.read_parquet(src,columns=key+flags,filters=[("section_id","==",sec)]),on=key,validate="1:1",how="left")
        assert d.spot_uid.notna().all()
        d["accepted_label"]=d.accepted_label.where(d.accepted_label.ge(0)).astype("Int64")
        d["niche_version"]="SAMPLE_R12_LEIDEN"
        d["main_score_layer"]=raw
        d["epi_anchor"]=d.analysis_eligible_legacy
        mean=pd.read_csv(r.c2l_mean_path,index_col=0)
        assert mean.index.is_unique and set(mean.index)==set(d.barcode) and set(mean)==set(maps.subtype)
        mean=mean.loc[d.barcode].reset_index(drop=True)
        total=mean.sum(axis=1,min_count=len(mean.columns))
        valid=(mean.ge(0)&np.isfinite(mean)).all(axis=1)&total.gt(0)&d.lesion_location.ne("quality_excluded")&d.center_qc
        prop=mean.div(total.where(valid),axis=0)
        d=pd.concat([d,mean.add_prefix("mean__"),prop.add_prefix("fraction__")],axis=1)
        for group,g in maps.groupby("major"):d[group+"_fraction"]=prop[g.subtype].sum(axis=1,min_count=len(g))
        d["Treg_over_T_NK"]=prop.Treg/d.T_NK_fraction.where(d.T_NK_fraction.gt(0))
        query=d.epi_anchor & d[raw].notna() & d.lesion_location.ne("quality_excluded")
        d["query"]=query
        for t,factors in targets.items():
            a=mean.reindex(columns=factors).sum(axis=1,min_count=len(factors))
            v=a/total.where(valid)
            d[t+"_fraction"]=v; d[t+"_abundance"]=a
            for q,scope in [(80,"all"),(90,"all"),(80,"margin")]:
                tag=f"{t}_P{q}_{scope}"
                threshold=v.quantile(q/100)
                defined=v.nunique()>1 and np.isfinite(threshold)
                target=(v.ge(threshold)&v.gt(0)) if defined else pd.Series(False,index=d.index)
                if scope=="margin":target &= d.distance_band_200um.isin(["inner_margin","outer_margin"])
                d[tag+"_target"]=target.astype("boolean").where(valid & defined)
                dist,reason=search(d,query,target)
                if not defined:reason[query]="undefined_relative_target"
                d[tag+"_distance_um"]=dist;d[tag+"_reason"]=reason
                d[tag+"_query_is_target"]=d[tag+"_target"].where(query)
                common=dict(section_id=sec,target=t,quantile=q,scope=scope)
                cov.append(dict(**common,factors=" + ".join(factors),n_eligible=valid.sum(),threshold=threshold,defined=defined,n_target=target.sum(),target_fraction=target.sum()/valid.sum(),eligible_min=v.min(),eligible_max=v.max(),abundance_min=a[valid].min(),abundance_max=a[valid].max(),target_composition_min=v[target].min(),target_composition_median=v[target].median(),target_composition_max=v[target].max(),target_abundance_min=a[target].min(),target_abundance_median=a[target].median(),target_abundance_max=a[target].max(),low_abundance_exploratory=a[target].median()<1))
                for positive in ([False,True] if (q,scope)==(80,"all") else [False]):
                    x=dist.where(dist.gt(0)) if positive else dist
                    rr,rho=corr(x,d[raw])
                    effects.append(dict(**common,positive_only=positive,layer=raw,n_query=query.sum(),n_target=target.sum(),n_used=x.notna().sum(),n_zero=x.eq(0).sum(),zero_fraction=x.eq(0).sum()/max(1,x.notna().sum()),distance_min=x.min(),distance_max=x.max(),pearson_r=rr,spearman_rho=rho))
                if (q,scope)==(80,"all"):
                    groups=pd.cut(dist.round(6),[-1,0,200,500,1000,np.inf],labels=["0","(0,200]","(200,500]","(500,1000]",">1000"])
                    for b in groups.cat.categories:
                        vals=d.loc[groups.eq(b),raw];bins.append(dict(section_id=sec,target=t,bin_um=b,n=len(vals),median=vals.median(),q25=vals.quantile(.25),q75=vals.quantile(.75)))
                    for field in ["lesion_location","distance_band_200um","accepted_label"]:
                        for level,z in d.groupby(field,dropna=False):
                            ix=z.index
                            localization.append(dict(section_id=sec,target=t,axis=field,level=level,n=len(ix),n_target=target[ix].sum(),target_fraction_of_level=target[ix].sum()/max(1,valid[ix].sum()),share_of_all_targets=target[ix].sum()/max(1,target.sum()),n_queries=query[ix].sum(),median_distance_um=dist[ix].median()))
        for field in ["lesion_location","distance_band_200um"]:
            for level,z in d.groupby(field):
                common=dict(section_id=sec,axis=field,region=level)
                cols=[g+"_fraction" for g in maps.major.unique()]+[t+"_fraction" for t in targets]+["Treg_over_T_NK"]
                regions.append(dict(**common,n=len(z),n_valid=valid[z.index].sum(),**z[cols].mean().to_dict()))
                for version,col,ok in [("SAMPLE_R12_LEIDEN","accepted_label",z["epi_anchor"]),("dataset_local","local_niche_id",valid[z.index])]:
                    u=z[ok & z[col].notna()]; den=len(u)
                    levels=list(range(9)) if col=="accepted_label" else d.loc[valid,col].dropna().unique()
                    for k in levels:
                        g=u[u[col].eq(k)];n=len(g)
                        niches.append(dict(**common,version=version,niche=k,n_eligible=ok.sum(),n_labelled=den,coverage=den/ok.sum() if ok.sum() else np.nan,n=n,fraction=n/den if den else np.nan,MP10_median=g[raw].median() if n>=5 else np.nan))
        counts.append(dict(section_id=sec,n_spots=len(d),n_mean=valid.sum(),n_raw=d[raw].notna().sum(),n_queries=query.sum(),n_R12=d.accepted_label.notna().sum(),n_local=d.local_niche_id.notna().sum(),n_epi_valid=d.epi_score_valid.eq(True).sum(),n_signed=d.signed_distance_um.notna().sum()))
        d.to_parquet(out/"spots"/(sec+".parquet"),index=False)
    for a,b in [(cov,"targets"),(effects,"associations"),(regions,"region_composition"),(niches,"region_niche"),(counts,"coverage"),(bins,"distance_bins"),(localization,"target_localization")]:save(a,b+".tsv")
if __name__=="__main__":main()

