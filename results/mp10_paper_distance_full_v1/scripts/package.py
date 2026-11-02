"""Index full figures. Vector-preserving PDF assembly uses the installed qpdf R binding."""
from distance import *
import subprocess,zipfile
def combine(inputs,output):
    assert inputs and all(p.is_file() for p in inputs)
    subprocess.run(["Rscript","-e","a=commandArgs(TRUE); qpdf::pdf_combine(a[-1],a[1]); stopifnot(qpdf::pdf_length(a[1])==length(a)-1)",str(output),*map(str,inputs)],check=True)
def main():
    m=pd.read_csv(tab/"input_manifest.tsv",sep="\t")
    rows=[];atlas_rows=[];summ=out/"figures/summary"
    for ds,group in m.groupby("dataset_id"):
        pdf=[]
        for sec in sorted(group.section_id):
            folder=out/"figures"/ds/sec
            for name in ["01_context","02_targets","03_distance","04_regions","05_program_maps","06_program_maps","07_program_maps","08_program_trends","09_program_effects"]:
                base=folder/name
                record=dict(dataset_id=ds,section_id=sec,figure=name,source_context=group.loc[group.section_id==sec,"source_context"].item(),atlas_page=len(pdf)+1,atlas=f"atlases/{ds}.pdf")
                for col,ext in [("pdf",".pdf"),("png",".png"),("body_pdf","_panel_only.pdf"),("legend_pdf","_legend_only.pdf")]:
                    p=Path(str(base)+ext);assert p.stat().st_size>0,p;record[col]=str(p.relative_to(out))
                pdf.append(Path(str(base)+".pdf"));rows.append(record)
        combine(pdf,out/"atlases"/(ds+".pdf"))
        atlas_rows.append(dict(dataset_id=ds,sections=len(group),pages=len(pdf),path=f"atlases/{ds}.pdf"))
    figs=sorted(summ.glob("*.png"))
    for p in figs:
        record=dict(dataset_id="overview",section_id="overview",figure=p.stem,pdf=str(p.with_suffix(".pdf").relative_to(out)),png=str(p.relative_to(out)))
        for col,end in [("body_pdf","_panel_only.pdf"),("legend_pdf","_legend_only.pdf")]:
            aux=p.with_name(p.stem+end)
            if p.stem!="00_coverage":assert aux.is_file(),aux
            record[col]=str(aux.relative_to(out)) if aux.exists() else None
        rows.append(record)
    pd.DataFrame(rows).to_csv(out/"figure_index.tsv",sep="\t",index=False,na_rep="NA")
    save(atlas_rows,"atlas_index.tsv")
    pages=[summ/"00_coverage.pdf",summ/"dataset_weights.pdf"]+[p.with_suffix(".pdf") for p in figs if "conditions" in p.name or "sensitivity" in p.name]+[summ/(s+".pdf") for s in ["lineage_composition","tumor_R12","tumor_programs","normal_R12","normal_programs"]]
    # Fixed pilot representatives supply the tissue maps in this cohort overview.
    reps=m[m.pilot_representative]
    pages.extend(out/"figures"/r.dataset_id/r.section_id/"01_context.pdf" for r in reps.itertuples())
    combine(pages,out/"full_overview.pdf")
    with zipfile.ZipFile(out/"review.zip","w",zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        files=[out/f for f in ["README.md","METHODS.md","RESULTS.md","parameters.tsv","full_sections.tsv","figure_index.tsv","full_overview.pdf"]]
        files+=list(tab.glob("*.tsv"))+list(tab.glob("*.tsv.gz"))
        files += [out/"figures"/r.dataset_id/r.section_id/"03_distance.pdf" for r in reps.itertuples()]
        for p in files:archive.write(p,p.relative_to(out))
    with zipfile.ZipFile(out/"review.zip") as archive:assert archive.testzip() is None
if __name__=="__main__":main()

