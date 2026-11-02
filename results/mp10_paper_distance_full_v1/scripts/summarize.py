"""Three descriptive scales. Equal section weights within units; datasets remain separate by source."""
from distance import *
def read(file):return pd.read_csv(tab/(file+".tsv"),sep="\t")
ID=["source_context","dataset_id","family","condition","target","feature","statistic"]
def effects():
    parts=[]
    def add(d,family,condition,feature,target,field,stat):
        parts.append(pd.DataFrame(dict(section_id=d.section_id,family=family,condition=condition,target=target,feature=feature,statistic=stat,effect=d[field],n_used=d.n_used)))
    a=read("associations")
    condition="P"+a["quantile"].astype(str)+"_"+a.scope+np.where(a.positive_only,"_positive","")
    for metric in ["pearson_r","spearman_rho"]:add(a,"MP10",condition,raw,a.target,metric,metric)
    a=read("adjusted")
    for c in ["adjusted_r","unadjusted_paired_r"]:add(a,"adjustment",c,raw,a.target,c,"pearson_r")
    a=read("matched_layers")
    for metric in ["pearson_r","spearman_rho"]:add(a,"model_layers","matched_P80_all",a.layer,a.target,metric,metric)
    for nm,family,target in [("program_distance","program_distance",None),("cholesterol_MP10","MP10_concordance","MP10"),("Treg_program_composition","Treg_composition","Treg")]:
        a=read(nm)
        for metric in ["pearson_r","spearman_rho"]:add(a,family,"P80_all" if target is None else "within_section",a.program,a.target if target is None else target,metric,metric)
    a=read("anchor_subset").rename(columns={"query_subset":"subset"})
    for metric in ["pearson_r","spearman_rho"]:add(a,"query_subset",a.subset,raw,a.target,metric,metric)
    return pd.concat(parts,ignore_index=True)
def summary():
    m=read("input_manifest")
    identity=["section_id","source_context","dataset_id","analysis_unit_uid","patient_uid","patient_id_status","eligible_for_patient_balanced_analysis","eligible_under_revised_analysis_unit_policy"]
    sec=effects().merge(m[identity],on="section_id",validate="m:1")
    save(sec,"section_effects.tsv.gz")
    d=sec.dropna(subset=["effect"]).copy()
    d["verified"]=d.patient_uid.where(d.eligible_for_patient_balanced_analysis.eq(True))
    d["unit"]=d.analysis_unit_uid.where(d.eligible_under_revised_analysis_unit_policy.eq(True))
    d["z"]=np.arctanh(d.effect.clip(-1+1e-12,1-1e-12))
    eligible=d[d.unit.notna()]
    unit=eligible.groupby(ID+["unit"],as_index=False).agg(n_sections=("section_id","nunique"),sections=("section_id","|".join),patient_uid=("verified","first"),mean_z=("z","mean"))
    unit["effect"]=np.tanh(unit.mean_z)
    unit["patient_policy_eligible"]=unit.patient_uid.notna()
    save(unit,"unit_effects.tsv.gz")
    ds=d.groupby(ID,as_index=False).agg(n_sections=("section_id","nunique"),n_units=("unit","nunique"),n_verified_patients=("verified","nunique"),median=("effect","median"),q25=("effect",lambda a:a.quantile(.25)),q75=("effect",lambda a:a.quantile(.75)),n_positive=("effect",lambda a:(a>0).sum()),n_negative=("effect",lambda a:(a<0).sum()))
    ds["population"]="all_sections"
    ds["dataset_effect"]=ds["median"]
    unitparts=[]
    for name,g in [("analysis_units",unit),("verified_patients",unit[unit.patient_policy_eligible])]:
        # Verified patient units follow the crosswalk; do not invent an ID for missing individuals.
        z=g.groupby(ID,as_index=False).agg(n_sections=("n_sections","sum"),n_units=("unit","nunique"),n_verified_patients=("patient_uid","nunique"),mean_z=("mean_z","mean"),median=("effect","median"),q25=("effect",lambda a:a.quantile(.25)),q75=("effect",lambda a:a.quantile(.75)),n_positive=("effect",lambda a:(a>0).sum()),n_negative=("effect",lambda a:(a<0).sum()))
        z["dataset_effect"]=np.tanh(z.pop("mean_z"));z["population"]=name;unitparts.append(z)
    ds=pd.concat([ds,*unitparts],ignore_index=True)
    save(ds,"dataset_effects.tsv.gz")
    ds["z"]=np.arctanh(ds.dataset_effect.clip(-1+1e-12,1-1e-12))
    equal=ds.groupby([k for k in ID if k!="dataset_id"]+["population"],as_index=False).agg(n_sections=("n_sections","sum"),n_units=("n_units","sum"),n_verified_patients=("n_verified_patients","sum"),n_datasets=("dataset_id","nunique"),datasets=("dataset_id","|".join),mean_z=("z","mean"),median_dataset_effect=("dataset_effect","median"))
    equal["dataset_equal_r"]=np.tanh(equal.pop("mean_z"))
    save(equal,"dataset_equal.tsv")
    # Direct equal-unit summary in addition to the explicitly dataset-balanced description.
    direct=[]
    for name,g in [("analysis_units",unit),("verified_patients",unit[unit.patient_policy_eligible])]:
        z=g.groupby([k for k in ID if k!="dataset_id"],as_index=False).agg(n_sections=("n_sections","sum"),n_units=("unit","nunique"),n_verified_patients=("patient_uid","nunique"),n_datasets=("dataset_id","nunique"),mean_z=("mean_z","mean"))
        z["equal_unit_r"]=np.tanh(z.pop("mean_z"));z["population"]=name;direct.append(z)
    save(pd.concat(direct),"unit_equal.tsv")
    # Module-level denominators include all registered sections, even when no effect is estimable.
    records=[]
    for field,file,col in [("distance","associations","pearson_r"),("adjusted","adjusted","adjusted_r"),("model","matched_layers","pearson_r"),("program_distance","program_distance","pearson_r"),("query_subset","anchor_subset","pearson_r"),("cholesterol","cholesterol_MP10","pearson_r"),("Treg_composition","Treg_program_composition","pearson_r")]:
        a=read(file)
        for s,g in a.groupby("section_id"):
            key=dict(section_id=s,module=field)
            for _,v in g[g[col].isna()].iterrows():
                records.append(dict(**key,target=v.get("target",field),condition=(str(v.get("quantile",""))+" "+str(v.get("scope",""))+" "+str(v.get("positive_only",v.get("query_subset","")))).strip(),feature=v.get("program",v.get("layer",raw)),reason=v.reason))
    save(records,"unavailable_effects.tsv")
    reg=m[identity+["pilot_representative","counts_exists","image_exists","genes_exists","barcodes_exists","c2l_mean_exists"]].merge(read("coverage"),on="section_id",validate="1:1").merge(read("extended_coverage"),on="section_id",validate="1:1")
    for name,g in [("P80_all",read("associations").query('quantile==80 and scope=="all" and not positive_only')),("P80_positive",read("associations").query('positive_only')),("P80_margin",read("associations").query('scope=="margin"')),("P90_all",read("associations").query('quantile==90')),("matched_layers",read("matched_layers")),("program_distance",read("program_distance"))]:
        reg=reg.merge(g.groupby("section_id").pearson_r.count().rename("n_valid_"+name),on="section_id",how="left",validate="1:1")
    reg=reg.merge(read("adjusted").groupby("section_id").adjusted_r.count().rename("n_valid_adjusted"),on="section_id",validate="1:1")
    reg.to_csv(out/"full_sections.tsv",sep="\t",index=False,na_rep="NA")
if __name__=="__main__":summary()

