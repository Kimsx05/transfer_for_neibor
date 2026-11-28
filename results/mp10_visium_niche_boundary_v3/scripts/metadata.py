"""Create the outcome-free metadata table before importing the analysis module."""
import csv
from pathlib import Path
root=Path(__file__).resolve().parents[1]
columns=['section_id','dataset_id','analysis_unit_uid','patient_uid','patient_id_status','source_context','image_path','evidence_method']
source=root.parent/'mp10_visium_niche_boundary_v2/audit/input_manifest.tsv'
target=root/'cache/immune_metadata.tsv'
with source.open() as inp:
    rows=list(csv.DictReader(inp,delimiter='\t'))
target.parent.mkdir(exist_ok=True)
with target.open('w') as out:
    w=csv.DictWriter(out,fieldnames=columns,delimiter='\t',lineterminator='\n',extrasaction='ignore')
    w.writeheader();w.writerows(rows)

