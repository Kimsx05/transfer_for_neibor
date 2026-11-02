from figures import *
import subprocess,zipfile
def summary():
    f,a,l=layout(2,2,"MP10 distance pilot  |  Seven sections","Section-level Pearson correlations. Positive r means higher mixed-spot MP10 farther from target.",w=6.4,h=3.9,left=.19,wspace=.80)
    for j,(q,scope,pos,title) in enumerate([(80,"all",False,"P80   /   All targets"),(80,"all",True,"P80   /   Positive distances"),(90,"all",False,"P90   /   All targets"),(80,"margin",False,"P80   /   Margin targets")]):
        z=A[(A["quantile"]==q)&(A.scope==scope)&(A.positive_only==pos)].pivot(index="section_id",columns="target",values="pearson_r").reindex(secs)
        heat(a[j],l[j],z,title,-.5,.5,"RdBu_r","Pearson r")
        for row in range(len(z)):
            for col in range(len(z.columns)):
                v=z.iloc[row,col]
                a[j].text(col+.5,row+.5,f"{v:+.2f}",ha="center",va="center",color="white" if abs(v)>.31 else "#111")
    export(f,a,l,"summary","00_section_effects")
def supporting():
    profile=pd.read_csv(tab/"source_R12_profiles.tsv",sep="\t")
    f,a,l=layout(1,2,"Accepted R12 environments  |  Reference ring profiles","Original SAMPLE_R12_LEIDEN neighborhood profile   •   mean relative q05 composition in R1 and R2",w=5.6,h=5,left=.10,wspace=.35)
    mats=[]
    for ring in ["R1","R2"]:
        vals=profile[profile.ring==ring].pivot(index="cluster",columns="cell_type",values="mean_relative_composition")
        mats.append(pd.DataFrame({t:vals[fac].sum(axis=1) for t,fac in targets.items()}))
    cap=max(x.max().max() for x in mats)
    for j,m in enumerate(mats):
        m.index=["R12 "+str(k) for k in m.index];heat(a[j],l[j],m,["First ring","Second ring"][j],0,cap,"viridis","Mean ring composition")
    export(f,a,l,"summary","01_R12_ring_profiles")
    z=pd.read_csv(tab/"program_distance.tsv",sep="\t")
    f,a,l=layout(2,2,"Spot transcriptional programs  |  Distance gradients","Pearson r within Epi queries   •   full-tissue P80 targets   •   color cap ±0.6",w=6.4,h=3.9,left=.19,wspace=.80)
    for j,name in enumerate([focus[4],focus[10],focus[2],focus[1]]):
        mat=z[z.program==name].pivot(index="section_id",columns="target",values="pearson_r").reindex(secs)
        heat(a[j],l[j],mat,short(name),-.6,.6,"RdBu_r","Pearson r")
        for row in range(len(mat)):
            for col in range(len(mat.columns)):
                v=mat.iloc[row,col];a[j].text(col+.5,row+.5,f"{v:+.2f}",ha="center",va="center",color="white" if abs(v)>.36 else "#111")
    export(f,a,l,"summary","02_program_gradients")
def package():
    folders=[out/"figures"/s for s in ["summary"]+secs]
    index=[]
    for folder in folders:
        for p in sorted(folder.glob("*.png")):
            base=p.with_suffix("")
            record={"section_id":folder.name,"figure":p.stem}
            for col,ext in [("full_pdf",".pdf"),("body_pdf","_body.pdf"),("legend_pdf","_legend.pdf"),("png",".png")]:
                file=Path(str(base)+ext)
                assert file.exists(),file
                record[col]=str(file.relative_to(out))
            index.append(record)
    assert len(index)==68 and len({(r["section_id"],r["figure"]) for r in index})==68
    pd.DataFrame(index).to_csv(out/"figure_index.tsv",sep="\t",index=False)
    pdfs=[out/"figures/summary"/p for p in ["00_section_effects.pdf","02_program_gradients.pdf","01_R12_ring_profiles.pdf"]]+[out/"figures/NPJ_sample6"/name for name in ["01_context.pdf","02_targets.pdf","03_distance.pdf","04_regions.pdf","07_program_maps.pdf","08_program_trends.pdf"]]+[out/"figures/GSE285715_T2/03_distance.pdf",out/"figures/GSE171351_B1/03_distance.pdf"]
    subprocess.run(["Rscript","-e","args=commandArgs(TRUE);qpdf::pdf_combine(args[-1],args[1])",str(out/"pilot_overview.pdf"),*map(str,pdfs)],check=True)
    with zipfile.ZipFile(out/"review.zip","w",compression=zipfile.ZIP_DEFLATED,compresslevel=6) as pack:
        files=[out/n for n in ["pilot_overview.pdf","METHODS.md","RESULTS.md","figure_index.tsv","pilot_sections.tsv"]]
        files+=list(tab.glob("*.tsv"))
        files+=[out/r["full_pdf"] for r in index]
        for p in files:
            if p.name.endswith("_match.tsv"):continue
            pack.write(p,p.relative_to(out))
if __name__=="__main__":summary();supporting();package()

