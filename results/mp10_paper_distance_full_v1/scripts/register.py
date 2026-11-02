import os, shutil
from pathlib import Path
import pandas as pd
root=Path(os.environ.get("MP10_PROJECT",Path(__file__).resolve().parents[2]))
out=Path(os.environ.get("MP10_OUTPUT",Path(__file__).resolve().parents[1]))
pilot=root/"mp10_paper_distance_pilot_v1"
for folder in ["tables","spots","scores","logs","figures","atlases"]: (out/folder).mkdir(parents=True,exist_ok=True)
read=lambda path:pd.read_csv(path,sep="\t")
m=read(os.environ.get("MP10_MANIFEST",root/"mp10_pathology_regions_v2/tables/section_manifest.tsv"))
assert m.section_id.is_unique
ids=read(root/"mp10_cnv_regions_interface_v1/input_review/section_identity_crosswalk.tsv").drop(columns=["n_spots"])
m=m.merge(ids,on="section_id",how="left",validate="1:1")
pilot_sections=read(pilot/"pilot_sections.tsv").section_id
m["pilot_representative"]=m.section_id.isin(pilot_sections)
assert set(pilot_sections)<=set(m.section_id)
for c in ["counts","genes","barcodes","c2l_mean","image"]:m[c+"_exists"]=m[c+"_path"].map(lambda v:Path(v).is_file())
m.to_csv(out/"tables/input_manifest.tsv",sep="\t",index=False,na_rep="NA")
for f in ["gene_sets.tsv.gz","method_comparison.tsv"]:shutil.copyfile(pilot/"tables"/f,out/"tables"/f)
for s in pilot_sections:
 for src in [pilot/"spots"/(s+"_programs.tsv.gz"),pilot/"tables"/(s+"_match.tsv")]:
  shutil.copyfile(src,out/"scores"/src.name)
