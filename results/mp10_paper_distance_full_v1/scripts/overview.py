"""Queue coverage, across-section effects, spatial summaries and seven prespecified examples."""
from figures import *
import subprocess
read=lambda name:pd.read_csv(tab/name,sep="\t")
tumor="tumor_study_specimen_not_region_diagnosis"
normal="adjacent_normal_source"
meta=manifest.set_index("section_id")
ds_names={x:x.replace("__UNVERIFIED","*") for x in manifest.dataset_id.unique()}
def with_meta(df):
    return df.drop(columns=["dataset_id","source_context"],errors="ignore").merge(manifest[["section_id","source_context","dataset_id"]],on="section_id",validate="m:1")
def sorted_regions(index):
    return [c for c in region_order if c in index]+[c for c in index if c not in region_order]
def describe_coverage():
    reg=pd.read_csv(out/"full_sections.tsv",sep="\t")
    rows=[]
    for (source,ds),g in reg.groupby(["source_context","dataset_id"]):
        rows.append(dict(source_context=source,dataset_id=ds,sections=len(g),spots=g.n_spots.sum(),queries=g.n_queries.sum(),units=g.analysis_unit_uid.nunique(),verified=g.loc[g.eligible_for_patient_balanced_analysis,"patient_uid"].nunique(),P80=g.n_valid_P80_all.eq(3).sum(),margin=g.n_valid_P80_margin.eq(3).sum(),adjusted=g.n_valid_adjusted.eq(3).sum(),matched=g.n_valid_matched_layers.eq(6).sum(),UCell55=g.n_sets_scored.eq(55).sum(),R12_sections=g.n_R12.gt(0).sum()))
    z=pd.DataFrame(rows);save(z,"cohort_coverage.tsv")
    f,ax=plt.subplots(figsize=(15,6));ax.axis("off")
    f.text(.04,.90,"Epithelial MP10   |   Cohort and estimability",fontsize=20,fontweight="bold")
    f.text(.04,.80,"A section counts for a distance module when all three routes are estimable. Source groups are kept separate.")
    for_print=z.copy();for_print["dataset_id"]=[ds_names[v]+(" (normal)" if s==normal else "") for v,s in zip(z.dataset_id,z.source_context)]
    for_print=for_print.drop(columns=["source_context","R12_sections"]).rename(columns={"dataset_id":"Dataset / source"})
    widths=[.20]+[.075]*(len(for_print.columns)-1)
    tb=ax.table(cellText=for_print.to_numpy(),colLabels=for_print.columns,loc="center",cellLoc="center",colWidths=widths);tb.auto_set_font_size(False);tb.set_fontsize(12);tb.scale(1,1.7)
    for (i,j),cell in tb.get_celld().items():
        cell.set_edgecolor("white");cell.set_facecolor("#D8E6E8" if not i else "#F4F6F6" if i%2==1 else "white")
    f.text(.04,.08,"* Original dataset alias unverified. Units include eligible anonymous sample units; verified counts are policy-eligible patient IDs.",fontsize=11)
    dest=out/"figures/summary";dest.mkdir(parents=True,exist_ok=True)
    for ext in ["pdf","png"]:f.savefig(dest/("00_coverage."+ext),dpi=300,facecolor="white")
    plt.close(f)
def section_pages():
    adj=read("adjusted.tsv");mod=read("matched_layers.tsv")
    blocks=[]
    for q,sc,pos,label in [(80,"all",False,"P80"),(80,"all",True,"Positive"),(90,"all",False,"P90"),(80,"margin",False,"Margin")]:
        z=A[(A["quantile"]==q)&(A.scope==sc)&(A.positive_only==pos)]
        blocks.append((label,z.pivot(index="section_id",columns="target",values="pearson_r")))
    contrasts=[(k,adj.pivot(index="section_id",columns="target",values=v)) for k,v in [("Paired raw","unadjusted_paired_r"),("Adjusted","adjusted_r")]]
    contrasts+= [(k,mod[mod.layer.eq(v)].pivot(index="section_id",columns="target",values="pearson_r")) for k,v in [("Model spot","spot_MP10_model_universe"),("Model Epi","epi_MP10_model_universe")]]
    for source in [tumor,normal]:
        ids=manifest[manifest.source_context.eq(source)].sort_values(["dataset_id","section_id"]).section_id.tolist()
        for start in range(0,len(ids),20):
            ss=ids[start:start+20];tag=("normal" if source==normal else "tumor")+f"_{start//20+1}"
            for group,name in [(blocks,"conditions"),(contrasts,"sensitivity")]:
                f,a,l=layout(1,4,tag.replace("_"," ")+"   |   MP10–distance effects","Within-section Pearson r   •   positive = higher score farther from target   •   grey = unavailable",w=4.0,h=max(3.2,.28*len(ss)),left=.14,wspace=.48)
                for j,(label,mat) in enumerate(group):
                    mat=mat.reindex(ss)
                    if j:mat.index=["" for _ in ss]
                    heat(a[j],l[j],mat,label,-.6,.6,"RdBu_r","Pearson r")
                export(f,a,l,"summary",tag+"_"+name)
def program_and_region():
    frac=[c for c in R if c.endswith("_fraction") and c.split("_fraction")[0] not in targets]
    reg=with_meta(R)
    z=reg.melt(id_vars=["source_context","dataset_id","section_id","axis","region"],value_vars=frac,var_name="feature",value_name="fraction").dropna()
    z["feature"]=z.feature.str.removesuffix("_fraction")
    by=["source_context","axis","region","feature"]
    ds=z.groupby(by+["dataset_id"]).fraction.mean().reset_index()
    summary=ds.groupby(by).fraction.agg(["mean","count"]).reset_index().rename(columns={"mean":"fraction","count":"n_datasets"})
    save(summary,"overview_region_lineages.tsv")
    niche=N[N.version=="SAMPLE_R12_LEIDEN"].merge(manifest[["section_id","source_context"]],on="section_id",validate="m:1")
    nr=niche.groupby(["source_context","axis","region","niche","dataset_id"])[["fraction","MP10_median"]].mean().reset_index()
    nsum=nr.groupby(["source_context","axis","region","niche"])[["fraction","MP10_median"]].mean().reset_index()
    save(nsum,"overview_R12.tsv")
    eq=read("dataset_equal.tsv")
    pr=read("region_programs.tsv");pr=with_meta(pr[pr.program.isin(focus)&pr.axis.eq("distance_band_200um")])
    grp=["source_context","region","program","population"]
    pr=pr.groupby(grp+["dataset_id"])["median"].mean().reset_index().groupby(grp)["median"].mean().reset_index()
    save(pr,"overview_region_programs.tsv")
    for src in [tumor,normal]:
        lab="tumor" if src==tumor else "normal"
        f,a,l=layout(1,2,lab.title()+" study sections   |   R12 by v3.4 band","Mean across sections within dataset, then equal datasets   •   regional MP10 medians require ≥5 Epi queries   •   grey = NA",w=7,h=4.4,left=.17,wspace=.60)
        for j,col in enumerate(["fraction","MP10_median"]):
            z=nsum[(nsum.source_context==src)&nsum.axis.eq("distance_band_200um")].pivot(index="region",columns="niche",values=col)
            z=z.reindex(sorted_regions(z.index))
            heat(a[j],l[j],z,["Fraction among labelled Epi anchors","Mixed-spot MP10 in Epi queries"][j],0,1 if j==0 else .15,"viridis",["Mean fraction","Mean regional median AUCell"][j])
        export(f,a,l,"summary",lab+"_R12")
        f,a,l=layout(1,3,lab.title()+" study sections   |   Prespecified programs","Distance effects: dataset-equal summary of section medians   •   region maps: equal dataset means of section medians",w=6.0,h=7,left=.16,wspace=.80)
        z=eq[(eq.source_context==src)&eq.family.eq("program_distance")&eq.statistic.eq("pearson_r")&eq.population.eq("all_sections")].pivot(index="feature",columns="target",values="dataset_equal_r").reindex(focus)
        z.index=[short(k) for k in z.index];heat(a[0],l[0],z,"Epi-query distance effects",-.6,.6,"RdBu_r","Pearson r")
        for j,pop in [(1,"tissue"),(2,"epi_query")]:
            z=pr[(pr.source_context==src)&pr.population.eq(pop)].pivot(index="program",columns="region",values="median").reindex(focus)
            z=z[sorted_regions(z.columns)].rename(columns={"candidate_inside_beyond_w":"Inside >200","candidate_outside_beyond_w":"Outside >200","inside_outer_unobserved":"Inside, no edge","outside_outer_unobserved":"Outside, no edge","outside_no_external_route":"Outside, no route","band_unresolved":"Unresolved","unlocalized":"Unlocalized"})
            z.index=["" for _ in focus];heat(a[j],l[j],z,["","Tissue program medians","Epi-query program medians"][j],0,.2,"viridis","Spot UCell (display cap 0.2)")
        export(f,a,l,"summary",lab+"_programs")
    subprocess.run(["Rscript",str(Path(__file__).with_name("composition.R")),str(out)],check=True)
def dataset_page():
    z=read("dataset_effects.tsv.gz")
    z=z[(z.source_context==tumor)&z.family.eq("MP10")&z.statistic.eq("pearson_r")&z.condition.eq("P80_all")]
    f,a,l=layout(1,3,"Tumor study sections   |   Main effects by dataset","Left: median section effect. Middle and right: equal-unit Fisher-z means within dataset. Grey indicates no eligible units.",w=4.9,h=4,left=.15,wspace=.70)
    for j,(pop,title) in enumerate([("all_sections","All sections"),("analysis_units","Analysis units"),("verified_patients","Verified patient subset")]):
        mat=z[z.population==pop].pivot(index="dataset_id",columns="target",values="dataset_effect").reindex(sorted(manifest.dataset_id.unique()))
        mat=mat.rename(index=ds_names)
        heat(a[j],l[j],mat,title,-.5,.5,"RdBu_r","Dataset effect")
    export(f,a,l,"summary","dataset_weights")
if __name__=="__main__":
    describe_coverage();section_pages();dataset_page();program_and_region()

