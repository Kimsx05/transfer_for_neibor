from pathlib import Path
import csv,json,hashlib,datetime,shutil
R=Path('/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_k32_validation_v4');V=R.parent/'mp10_resolution_scan_v3';O=R.parent/'mp_ucell_state_v1'
def read(p):
 with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows):
 with open(p,'w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
paths=[V/x for x in ['README.md','CONFIG.yaml','COMPLETED.json','00_audit/INPUT_MANIFEST.tsv','00_audit/cell_metadata.qs','00_audit/global_scaling.tsv','00_audit/frozen_high_labels.tsv.gz','00_audit/high_thresholds.tsv','00_audit/high_thresholds.qs','00_audit/gene_annotation.tsv.gz','00_audit/RNA_counts_check.tsv','00_audit/RNA_data_check.tsv','01_partitions/FULL_Z_r1.00.qs','03_candidates/candidates_preDE.tsv','03_candidates/comparisons_preDE.tsv','04_pseudobulk_markers/comparison_model_index.tsv','04_pseudobulk_markers/aggregate_FULL_Z_r1.00.qs','04_pseudobulk_markers/sample_totals.qs','05_function/frozen_GO_BP_Reactome.qs','05_function/resource_provenance.txt','07_figures/FULL_Z_shared_embedding.tsv.gz','07_figures/RNA_shared_embedding.tsv.gz','code/04_models.R','code/03_aggregate.R']]
paths += [O/x for x in ['02_ucell_scores/M_RAW.qs','02_ucell_scores/M_Z.qs','02_ucell_scores/original_usage_raw.qs','01_signatures/top50_signatures.tsv']]
paths.append(Path('/home/data/t070721/Metabolistic_scRNA/Multi2/3.4 malignant cell process/PCA_ALLGENE1/epi_cycling_allgene_harmony_clustered_final_after_small_removal.qs'))
for r in read(V/'04_pseudobulk_markers/comparison_model_index.tsv'):
 if r['candidate_id']=='K32':
  paths += [f for f in (V/r['model_dir']).glob('*') if f.is_file()]
  paths += [f for f in (V/'05_function'/r['model_key']).glob('*') if f.is_file()]
rows=[]
for p in sorted(set(paths)):
 st=p.stat();rows.append(dict(path=str(p),bytes=st.st_size,mtime_ns=st.st_mtime_ns,mtime=datetime.datetime.fromtimestamp(st.st_mtime,datetime.timezone(datetime.timedelta(hours=8))).isoformat(),sha256=hashlib.file_digest(p.open('rb'),'sha256').hexdigest() if st.st_size<150*1024**2 else 'not_rehashed_large_file',version_evidence='completed v3 / original frozen v1 and NMF inputs; exact K32 membership and sparse aggregate checks in 00_audit'))
f=R/'00_audit/INPUT_MANIFEST.tsv'
if f.exists():
 old=read(f);assert [(x['path'],x['bytes'],x['mtime_ns']) for x in old]==[(x['path'],str(x['bytes']),str(x['mtime_ns'])) for x in rows]
else:write(f,rows)
for x in ['resource_provenance.txt','pathway_manifest.tsv']:shutil.copy2(V/'05_function'/x,R/'00_audit'/x)
for x in ['RNA_counts_check.tsv','RNA_data_check.tsv','software_versions.tsv']:shutil.copy2(V/'00_audit'/x,R/'00_audit'/('v3_'+x))
shutil.copy2(O/'01_signatures/top50_signatures.tsv',R/'00_audit/top50_signatures.tsv')
# Distinguish filtered pathways from fgsea numerical estimation NA.
f=R/'06_supplement/focus_pathways_all_models.tsv';rs=read(f)
import gzip
cache={}
for r in rs:
 key=(r['branch'],r['control'])
 if key not in cache:
  with gzip.open(R/'02_rna_DEG'/key[0]/key[1]/'function/enrichment_all.tsv.gz','rt') as h:cache[key]={x['pathway']:x for x in csv.DictReader(h,delimiter='\t')}
 r['status']='NOT_RETURNED_BY_FIXED_SIZE_MAPPING_FILTER' if r['pathway'] not in cache[key] else ('ESTIMATION_NA' if r['NES']=='NA' else 'TESTED')
write(f,rs)
print('Input manifest frozen:',len(rows),'files')
