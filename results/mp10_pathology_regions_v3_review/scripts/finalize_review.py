from run_v3 import *
from zipfile import ZipFile,ZIP_DEFLATED
# Check actual deliverables, not just completion markers.
s=pd.read_csv(O/'tables/section_summary.tsv',sep='\t');cohort=pd.read_parquet(O/'tables/all_spots_v3.parquet');fm=pd.read_csv(R/'spatial_ecology_top20_v2/00_inputs/section_manifest.tsv',sep='\t').set_index('section_id');checks=[]
for z in s.itertuples():
 p=O/'by_dataset'/z.dataset_id/z.section_id;d=pd.read_parquet(p/'spot_regions_v3.parquet');old=pd.read_parquet(Path(M.loc[z.section_id,'output_dir'])/'spot_region_labels.parquet').set_index('barcode').reindex(d.barcode);mp=pd.read_csv(fm.loc[z.section_id,'score_path'],index_col=0).reindex(d.barcode)
 assert np.array_equal(d.MP10_original.to_numpy(),mp.EpiMP10_Top20_AUCell.to_numpy())
 assert np.array_equal(d.malignant_fraction.to_numpy(),old.malignant_fraction.to_numpy())
 assert len(d)==z.n_spots and set(d.barcode)==set(cohort.loc[cohort.section_id==z.section_id,'barcode'])
 for stem in ['01_v3_core','02_v2_SpaCET_v3_comparison','03_bandwidth_sensitivity','04_domains','05_main_regions']:
  with Image.open(p/(stem+'.png')) as im:assert min(im.info.get('dpi',(0,0)))>299
  assert (p/(stem+'.pdf')).stat().st_size>1000
 for fn in ['lesion_rois_um.geojson','lesion_rois_image_px.geojson','lesion_outer_boundaries_um.geojson','tumor_nest_contact_boundaries_um.geojson']:
  js=json.loads((p/fn).read_text());assert js['type']=='FeatureCollection'
 checks.append({'section_id':z.section_id,'n_spots':len(d),'figures':'5 PNG at 300 dpi + 5 PDF','MP10_original_numeric_equality':True,'SpaCET_original_numeric_equality':True,'geometry':'present'})
pd.DataFrame(checks).to_csv(O/'tables/final_delivery_checks.tsv',sep='\t',index=False)
# Portable review package: all H&E/ROI and main maps, pilot comparisons, small tables and code; no large matrices/full objects.
files=[O/'METHODS.md',O/'RESULTS.md',O/'config.json'];files+=list((O/'scripts').glob('*.py'))+list((O/'scripts').glob('*.R'));files+=list((O/'tables').glob('*.tsv'))+list((O/'tables').glob('*.tsv.gz'));files+=list((O/'figures').glob('*'));files+=list((O/'comparators').glob('*/comparison.*'))
rows=[]
for z in s.itertuples():
 p=O/'by_dataset'/z.dataset_id/z.section_id;base=str(p.relative_to(O));files += [p/'01_v3_core.png',p/'05_main_regions.png',p/'05_main_regions.pdf',p/'coordinate_scale.json',p/'lesion_rois_um.geojson',p/'lesion_boundary_edges.tsv.gz']
 extra=''
 if z.section_id in C['pilot_sections']:
  files += [p/'02_v2_SpaCET_v3_comparison.png',p/'02_v2_SpaCET_v3_comparison.pdf',p/'03_bandwidth_sensitivity.png'];extra=f' · <a href="{base}/02_v2_SpaCET_v3_comparison.png">pilot comparison</a>'
 rows.append(f'<tr><td>{z.section_id}</td><td><a href="{base}/01_v3_core.png">H&E/ROI/distance</a> · <a href="{base}/05_main_regions.png">main spatial regions</a>{extra}</td></tr>')
review='<!doctype html><meta charset="utf-8"><title>v3 review</title><style>body{font:15px sans-serif;margin:30px}td{padding:8px;border-bottom:1px solid #ddd}</style><h1>BLCA lesion regions v3 — review</h1><p>All 63 sections. Computational candidates, no documented pathologist confirmation. This compact archive excludes full per-section evidence tables and large objects; use the original project directory for full reproducibility inputs.</p><p><a href="RESULTS.md">Results</a> · <a href="METHODS.md">Methods</a> · <a href="tables/all_spot_labels_v3.tsv.gz">265,002 spot core labels</a></p><table>'+''.join(rows)+'</table>'
with ZipFile(O/'review.zip','w',ZIP_DEFLATED,compresslevel=5) as archive:
 for p in sorted(set(files)):
  if p.is_file():archive.write(p,str(p.relative_to(O)))
 archive.writestr('index.html',review)
size=sum(p.stat().st_size for p in O.rglob('*') if p.is_file());(O/'logs/delivery_size.json').write_text(json.dumps({'result_total_GiB':size/1024**3,'review_zip_MiB':(O/'review.zip').stat().st_size/1024**2,'n_sections':len(s),'n_spots':len(cohort)},indent=2));print('PASS final delivery',len(s),len(cohort),'total GiB',size/1024**3,'zip MiB',(O/'review.zip').stat().st_size/1024**2)
