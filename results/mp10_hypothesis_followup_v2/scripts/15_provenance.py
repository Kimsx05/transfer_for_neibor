import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import pathlib,pandas as pd,hashlib,json,re,datetime,importlib.metadata,shutil
P=pathlib.Path.cwd();B=P/'results/mp10_hypothesis_followup_v2';T=B/'tables';OLD=P/'results/mp10_immune_endothelial_mechanisms_v1';V1=P/'results/mp10_hypothesis_followup_v1';REF=P/'resources/mp10_regulatory_priors'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for c in iter(lambda:f.read(8*1024**2),b''):h.update(c)
 return h.hexdigest()
base=pd.read_csv(T/'input_manifest.tsv',sep='\t');paths=set(pathlib.Path(p) for p in base.path)
paths|={P/'subtype/v2/final_annotation_pre_reference/final_annotation_v2_pre_reference.qs',P/'neighborhood_mp10_calibration_v5/01_score_fidelity/top20.txt'}
paths|={pathlib.Path(p) for p in json.loads((B/'config/scenic_mapping_correction.json').read_text())['sources']}
paths|={OLD/'objects/model_parts/M0073.qs'}
for n in ['coverage.tsv.gz','treg_cycle_composition.tsv','original_source_input_manifest.tsv']:paths.add(V1/'tables'/n)
for n in ['pseudobulk_profiles.tsv','model_status.tsv','epi_existing_cluster_composition.tsv','activity_associations.tsv.gz','scenic_membership_Epi.tsv.gz','scenic_genes_Epi.tsv','scenic_Epi_aucell_parts.tsv','scenic_Epi_M20-masked_aucell_parts.tsv','scenic_second_seed_coverage.tsv','specimen_exposure_eligibility.tsv']:paths.add(OLD/'tables'/n)
for n in ['collectri.tsv','dorothea_AB.tsv']:paths.add(REF/n)
for id in ['M0001','M0016','M0031','M0073']:
 paths.add(OLD/'tables/de'/f'{id}.tsv.gz')
 for mask in ['Full','M20-masked']:
  p=(V1/'tables/gsea_repair' if id=='M0031' else OLD/'tables/gsea')/f'{id}_Hallmark_{mask}.tsv.gz'
  if p.exists():paths.add(p)
for fn in ['scenic_Epi_aucell_parts.tsv','scenic_Epi_M20-masked_aucell_parts.tsv']:
 paths.add(pathlib.Path(pd.read_csv(OLD/'tables'/fn,sep='\t').file.iloc[0]))
# Exact spatial source audit already records used paths.
sp=pd.read_csv(T/'spatial_input_audit.tsv',sep='\t')
for col in sp.columns:
 if 'path' in col:
  for x in sp[col].dropna():
   if pathlib.Path(str(x)).is_file():paths.add(pathlib.Path(x))
for folder,files in [('mp10_dataset_distance_heterogeneity_v1',['tables/discussion_sensitivity.tsv.gz','tables/unit_endpoints.tsv.gz','tables/unit_metadata.tsv']),('mp10_p95_front_programs_v1',['tables/context85_effects.parquet','scripts/analyze.py'])]:
 for f in files:
  p=P/'results'/folder/f
  if p.exists():paths.add(p)
for f in ['mp10_immune_endothelial_mechanisms_v1/scripts/13_scenic_ctx.py','mp10_immune_endothelial_mechanisms_v1/scripts/35_reference_coverage.py','mp10_immune_endothelial_mechanisms_v1/scripts/02_models_worker.R','mp10_hypothesis_followup_v1/scripts/04_scores.R']:
 p=P/f
 if p.exists():paths.add(p)
a=[]
for p in sorted(paths):
 if not p.exists():raise FileNotFoundError(p)
 st=p.stat();a.append(dict(path=str(p),bytes=st.st_size,mtime_utc=datetime.datetime.fromtimestamp(st.st_mtime,datetime.timezone.utc).isoformat(),sha256=sha(p),reuse='read_only_original_dependency'))
pd.DataFrame(a).to_csv(T/'input_reuse_audit.tsv',sep='\t',index=False)
v=[]
for n in ['pyscenic','ctxcore','arboreto','numpy','pandas','scipy','numba','pyarrow','matplotlib','Pillow','PyMuPDF']:
 try:v.append(dict(package=n,version=importlib.metadata.version(n)))
 except importlib.metadata.PackageNotFoundError:pass
pd.DataFrame(v).to_csv(B/'logs/environment_versions.tsv',sep='\t',index=False)
resources=[]
for p in (B/'logs').glob('*.log'):
 s=p.read_text(errors='replace');m=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',s)
 if m:resources.append(dict(stage=p.name,peak_RSS_KiB=int(m[1]),peak_RSS_GiB=int(m[1])/1024**2))
resources.append(dict(stage='delivery_disk_increment',bytes=sum(p.stat().st_size for p in B.rglob('*') if p.is_file()),description='Own result directory only; original objects and databases not copied. Code directory additional small scripts.'))
pd.DataFrame(resources).to_csv(T/'resource_usage.tsv',sep='\t',index=False)
(B/'scripts').mkdir(exist_ok=True)
for p in (P/'mp10_hypothesis_followup_v2/scripts').iterdir():
 if p.is_file():shutil.copy2(p,B/'scripts'/p.name)
print('Provenance',len(a),'inputs, peak RSS',max(x.get('peak_RSS_GiB',0) for x in resources))
