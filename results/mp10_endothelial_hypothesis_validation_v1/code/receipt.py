"""Close the run with resource provenance and reproducible file checksums."""
import csv, gzip, hashlib, sys
from pathlib import Path
from xml.etree import ElementTree
from zipfile import ZipFile

out=Path(__file__).resolve().parents[1]
project=out.parents[1]
cache=project/"cache/mp10_public_resources_v1"

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1048576), b""): h.update(block)
    return h.hexdigest()

def save(name, columns, rows):
    with (out/"tables"/name).open("w") as fp:
        w=csv.writer(fp,delimiter="\t",lineterminator="\n")
        w.writerow(columns);w.writerows(rows)

with (out/"logs/resource_attempts.tsv").open() as f:
    attempts=list(csv.DictReader(f,delimiter="\t"))
rows=[]
for acc,name in sorted({(x["resource"],x["file"]) for x in attempts}):
    path=cache/acc/name; history=[x for x in attempts if (x["resource"],x["file"])==(acc,name)]
    urls=[x["url"] for x in history if x["status"]=="downloaded"]
    if not path.is_file():
        rows.append((acc,name,history[-1]["url"],"unavailable",0,"","not used"));continue
    checksum=sha(path)
    assert checksum in {x["sha256"] for x in history if x["status"] in ("downloaded","reused")}
    if acc.startswith("GSE"):
        role="parsed SOFT library mapping" if "soft" in name else "expression read into external contrast models"
    elif acc.startswith("GENCODE"):
        role="transcript-to-gene mapping read into ATF4 contrast models"
    elif name.endswith(".xml"):
        root=ElementTree.parse(path);role="article text and supplementary index; "+str(len(root.findall(".//p")))+" paragraphs"
    elif name.endswith(".xlsx"):
        with ZipFile(path) as z:
            sheets=[i.attrib["name"] for i in ElementTree.fromstring(z.read("xl/workbook.xml")).iter() if i.tag.endswith("}sheet")]
        role=("original signature input; " if acc in ("PMC8432967","eLife.63326") else "literature supplementary reference; ")+" | ".join(sheets)
    elif name.endswith(".pdf"):role="literature supplementary PDF; see supplement_pdf_readout.tsv"
    else:role="reference"
    rows.append((acc,name,urls[-1] if urls else history[0]["url"],"verified",path.stat().st_size,checksum,role))
save("public_resource_manifest.tsv",["accession","file","source_url","status","bytes","sha256","read_role"],rows)

v1=project/"results/mp10_immune_endothelial_mechanisms_v1";v2=project/"results/mp10_hypothesis_followup_v2"
paths=[v1/"objects/pseudobulk.qs",v1/"plot_data/receiver_Hallmark.tsv.gz",v1/"plot_data/Hallmark_sensitivity_Capillary_EC.tsv.gz",v1/"config/run_config.json"]
paths += [v1/"tables"/k for k in ("de/M0073.tsv.gz","gsea_results.tsv.gz","gsea_coverage.tsv.gz","gsea_sensitivity_results.tsv.gz","specimen_exposure_eligibility.tsv","full_celltype_counts.tsv","patient_relationship_audit.tsv","activity_associations.tsv.gz","scenic_Endothelial_regulon_coverage.tsv","scenic_second_seed_coverage.tsv","prior_molecules_per_profile.tsv.gz","nichenet_lr_prior.tsv.gz","nichenet_ligand_activities_annotated.tsv.gz","pseudobulk_profiles.tsv","final_label_dictionary.tsv","input_manifest.tsv")]
paths += [v2/"tables"/k for k in ("EC_M20_Hallmark_complete.tsv.gz","model_results.tsv.gz","sample_level_scores.tsv.gz","mTORC1_M0073_membership.tsv","mTORC1_spatial_source_evidence.tsv","input_manifest.tsv")]
paths += [project/"resources/mp10_regulatory_priors"/k for k in ("hallmark.tsv.gz","ligand_target_matrix_nsga2r_final.rds","hgnc_complete_set_20261118.txt")]
paths += [project/"mp10_spatial_communication_v1/epi_state_type_contrasts.tsv",project/"neighborhood_mp10_calibration_v5/01_score_fidelity/top20.txt"]
paths += list((project/"mp10_cellchat_high_definition_sensitivity_v1/tables/lr").glob("C10__ALL__*__raw_sparse.tsv.gz"))
save("local_input_manifest.tsv",["path","bytes","sha256"],[(str(f.relative_to(project)),f.stat().st_size,sha(f)) for f in paths])
files=[p for p in sorted(out.rglob("*")) if p.is_file() and "__pycache__" not in p.parts and p.name!="delivery_manifest.tsv"]
save("delivery_manifest.tsv",["path","bytes","sha256"],[(str(p.relative_to(out)),p.stat().st_size,sha(p)) for p in files])
