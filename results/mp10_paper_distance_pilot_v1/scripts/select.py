"""Freeze the pre-outcome v8 pilot IDs and verify physical input availability."""
from pathlib import Path
import pandas as pd
root=Path(__file__).resolve().parents[2]
out=root/"mp10_paper_distance_pilot_v1"
read=lambda p:pd.read_csv(p,sep="\t")
v8=read(root/"neighborhood_spottedpy_v8/00_audit/pilot_sections.tsv")
manifest=v8.merge(read(root/"mp10_pathology_regions_v2/tables/section_manifest.tsv"),on=["section_id","dataset_id"],validate="one_to_one")
assert len(manifest)==7 and manifest.dataset_id.is_unique
for k in ["counts","barcodes","genes","image","c2l_mean"]:
    manifest[k+"_exists"]=manifest[k+"_path"].map(lambda x:Path(x).is_file())
manifest.to_csv(out/"tables/input_manifest.tsv",sep="\t",index=False)
p=manifest[["section_id","dataset_id","source_context","n_spots","n_E0"]].copy()
p=p.merge(read(root/"mp10_cnv_regions_interface_v1/input_review/section_identity_crosswalk.tsv")[["section_id","patient_uid","patient_id_status","analysis_unit_uid"]],on="section_id",validate="one_to_one")
p["selection_reason"]="Frozen v8 representative: closest to dataset median eligible center count; lexical tie break"
p["input_coverage"]=manifest.filter(regex="_exists$").apply(lambda row:", ".join(k.replace("_exists","") for k in row.index if row[k]),axis=1)
p.to_csv(out/"pilot_sections.tsv",sep="\t",index=False,na_rep="NA")

