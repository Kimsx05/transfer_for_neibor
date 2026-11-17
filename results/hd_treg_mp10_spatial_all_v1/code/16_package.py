from common import *
import shutil,zipfile
checks=pd.read_csv(R/'tables/validation_checks.tsv',sep='\t');assert checks.passed.all()
assert len(pd.read_csv(R/'tables/association_primary.tsv',sep='\t'))==10
assert (R/'logs/R_loader_validation.log').read_text().count('PASS')==5
for s in SECTIONS:
 for k in ['T32','T16','Epi64']:
  d=pd.read_parquet(R/s/f'{k}_regions.parquet');assert d.region_id.is_unique and 'measurement_eligible'in d
# Delete only regenerable intermediates produced by this task. Retain80-column native raw check.
for p in (R/'temporary').glob('*_ucell_*'):p.unlink()
review=R/'review'
for sub in ['code','tables','figures','regional_tables']:(review/sub).mkdir(parents=True,exist_ok=True)
for p in list(R.glob('*.md'))+list(R.glob('*.yaml'))+list(R.glob('*.json'))+list(R.glob('*.sh')):shutil.copy2(p,review/p.name)
for p in (R/'code').iterdir():
 if p.suffix in ['.py','.R']:shutil.copy2(p,review/'code'/p.name)
for p in (R/'tables').glob('*.tsv'):
 if '_gene_mapping'in p.name or 'UCell_native_raw'in p.name:continue
 if p.stat().st_size<3_000_000:shutil.copy2(p,review/'tables'/p.name)
for name in ['spottedpy_hotspot_distances.tsv.gz']:
 shutil.copy2(R/'tables'/name,review/'tables'/name)
shutil.copy2(R/'sources/literature_signature.tsv',review/'tables/literature_original309.tsv')
for p in (R/'figures').glob('*.png'):shutil.copy2(p,review/'figures'/p.name)
for s in SECTIONS:
 for k in ['T32','T16','Epi64']:
  d=pd.read_parquet(R/s/f'{k}_regions.parquet');cols=[c for c in d if not c.startswith(('raw_','rank_'))];d[cols].to_parquet(review/'regional_tables'/f'{s}_{k}.parquet',index=False)
status={'all_five_full_supplied_inputs_processed':True,'main_results_rows':10,'main_testable_rows':7,'inference_unavailable_main':['Cervilla50','MPBC1_50','MPBC2_50'],'formal_communications_run':False,'all_inputs_read_only':True,'R_loader_sections_passed':5,'validation_passed':int(checks.passed.sum()),'regions_reloaded_after_final_fields':True,'scope_finished':True}
(R/'completion_status.json').write_text(json.dumps(status,indent=2));shutil.copy2(R/'completion_status.json',review/'completion_status.json')

def archive():
 with zipfile.ZipFile(R/'review.zip','w',zipfile.ZIP_DEFLATED,compresslevel=6)as z:
  for p in sorted(review.rglob('*')):
   if p.is_file():z.write(p,Path('review')/p.relative_to(review))
 with zipfile.ZipFile(R/'review.zip')as z:assert z.testzip()is None
archive()
files=[p for p in R.rglob('*')if p.is_file()];logical=sum(p.stat().st_size for p in files);allocated=sum(p.stat().st_blocks*512 for p in files)
rows=pd.read_csv(R/'tables/resource_usage.tsv',sep='\t');rss=float(rows.peak_RSS_KiB.max())/1024**2
resource=f'''# Resource and delivery size

- Total CPU affinity cap: 4 logical CPUs; BLAS/OpenMP/Numba threads1, R sequential. No shared environment upgrades.
- Largest recorded single-process RSS: {rss:.2f} GiB. Not a simultaneous aggregate RSS measurement.
- Native all-region UCell:14:08.87 elapsed. Per-step wall/RSS: tables/resource_usage.tsv; original logs retained.
- User explicitly authorized daytime continuation: runtime_authorization.json. No queued job claimed complete.
- Result logical bytes including review directory+zip: {logical:,} ({logical/1024**3:.3f} GiB).
- Allocated bytes at packaging: {allocated:,} ({allocated/1024**3:.3f} GiB).
- review.zip: {(R/'review.zip').stat().st_size:,} bytes ({(R/'review.zip').stat().st_size/1024**2:.2f} MiB).
- Approximately167 MiB regenerable UCell MatrixMarket intermediates removed after validation. Original large expression objects are referenced, not copied.
- Size is a packaging-time snapshot; this small resource note and final log updates can change it slightly.
'''
(R/'RESOURCE_USAGE.md').write_text(resource);shutil.copy2(R/'RESOURCE_USAGE.md',review/'RESOURCE_USAGE.md');archive()
print(resource)
