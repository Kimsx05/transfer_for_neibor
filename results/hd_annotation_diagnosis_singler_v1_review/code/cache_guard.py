from pathlib import Path
import json,hashlib
root=Path(__file__).resolve().parents[1];old=root.parent/'hd_major_annotation_v1'
paths=[old/'reference/reference.rds',old/'reference/marker_panels.rds']
for s in ['Cervilla','D1_Johnson']:
 paths += [old/s/'seurat.rds']+[old/s/'tables'/f for f in ['analysis_metadata.tsv.gz','marker_evidence.tsv.gz','annotations.tsv.gz','confidence_support_bands.tsv','gene_identity.tsv']]
sig={'input_files':{str(p.relative_to(root.parent)):{'bytes':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns} for p in paths},'core_scripts':{f:hashlib.sha256((root/'code'/f).read_bytes()).hexdigest() for f in ['diagnose_rules.py','train_blca.R','train_hpca.R','predict_section.R','compare.py','marker_summary.R','plots.R','validate.R','focused_evidence.R','review_evidence.py']},'parameters':{'seed':123,'threads':4,'chunk_size':5000,'BLCA_wilcox_de_n':50,'BLCA_aggr_ref':True,'HPCA_classic':True,'prune_after_whole_section_merge':True,'prune':'nmads=3,min.diff.med=-Inf,min.diff.next=0'}}
p=root/'cache_manifest.json'
if p.exists():
 if json.loads(p.read_text())!=sig:raise SystemExit('Input or core parameter/code changed: do not reuse cached predictions; inspect and rebuild affected outputs deliberately.')
else:p.write_text(json.dumps(sig,indent=2))
print('CACHE_INPUTS_PARAMETERS_MATCH')
