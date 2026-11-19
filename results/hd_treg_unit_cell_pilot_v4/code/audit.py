from base import *
from datetime import datetime
manifest=[]; checks=[]
roi=rd(Z/'tables/roi_manifest.tsv');roi=roi[roi.section.isin(sections)].copy();roi['evaluation']='256um half-open core; nucleus or fixed-grid geometric center';roi['read_scope']='complete source members of selected IDs; full section Epi64 query';put(roi,'roi_manifest.tsv')
ref=Path(js(V/'models/reference_prepared.json')['reference'])
for s in sections:
 (R/s).mkdir(exist_ok=True)
 audit=js(N/s/'metadata/input_audit.json');trans=js(N/s/'metadata/transform.json')
 paths=[ref,V/'tables/reference_membership.parquet',Path(audit['input_matrix']),Path(audit['source_image']),N/s/'metadata/transform.json',N/s/'inputs/model_image.npy',S/s/'inputs/nuclei.npz',S/s/'bin_to_cell.parquet',S/s/'object_qc.parquet',S/s/'cells.h5ad',N/s/'cells.h5ad',V/s/'bin_evidence.parquet',Z/'tables/MP10_high_thresholds.tsv',V/'frozen_config.yaml',F/'tables/paired_effects.tsv',O/'tables/bin8_rctd_results.tsv.gz',O/'tables/conservative4um_ROI_cells.tsv.gz']
 for kind in ['T16','T32','Epi64']:paths.extend([V/s/(kind+'_regions.parquet'),V/s/(kind+'_members.parquet')])
 for p in paths:manifest.append(dict(section=s,path=str(p),readable=os.access(p,os.R_OK),size=p.stat().st_size if p.exists()else None,native_mpp=trans['native_mpp'],crop_xyxy=str(trans['fullres_crop_xyxy']),coords='fullres*metadata mpp; x right y down',counts='supplied raw integers; zero retained'))
 fields=['region_id','component','x_um','y_um','measurement_eligible','MP10','H05','H10','n_bins']
 e=pd.read_parquet(V/s/'Epi64_regions.parquet',columns=fields).set_index('region_id');z=pd.read_parquet(Z/s/'Epi_regions.parquet').set_index('region_id');m=pd.read_parquet(V/s/'Epi64_members.parquet');assert e.index.is_unique and m.bin_id.is_unique and set(e.index)==set(m.region_id)
 c=pd.DataFrame(index=e.index);c['section']=s;c['eligible_matches_v3_presence']=e.measurement_eligible==e.index.isin(z.index);c['member_count_exact']=e.n_bins==m.region_id.value_counts().reindex(e.index);c['member_source']=str(V/s/'Epi64_members.parquet')
 for k in ['component','x_um','y_um','MP10','H05','H10']:
  c[k+'_v3_exact']=pd.Series(pd.NA,index=e.index,dtype='boolean');c.loc[z.index,k+'_v3_exact']=e.loc[z.index,k].eq(z[k]).values
 assert c.select_dtypes('bool').fillna(True).all().all()
 th=rd(Z/'tables/MP10_high_thresholds.tsv').query('section==@s')
 for h in ['H05','H10']:
  cutoff=th[th.definition.eq(h)].threshold.item();assert np.array_equal(z[h],z.MP10>cutoff);c[h+'_frozen_cutoff']=cutoff
 checks.append(c.reset_index())
put(pd.concat(checks),'mp10_frozen_check.tsv');put(pd.DataFrame(manifest),'input_manifest.tsv')
config=dict(frozen=datetime.now().astimezone().isoformat(),roi_source=str(Z/'tables/roi_manifest.tsv'),reference=str(ref),membership=str(V/'tables/reference_membership.parquet'),groups='sample grouped, 3 fold within train; historical sample_holdout evaluation only',features='top 24 positive mean log-expression contrasts per class in each training fold; common symbols only; sum duplicates before normalization',excluded=exclude,model='multinomial L2 LBFGS, C=.1 or 1 chosen by grouped out-of-fold macro recall',threshold='minimum grid cutoff 0.3..0.95 step0.05 meeting Treg OOF precision >=0.8, FPR<=0.01 and >=10 true positives at 100 UMI thinning; require raw T support and predicted Treg. Otherwise no supported tier',depth_check=[25,100,500,'native'],seed=123,threads=1,roi_count=len(roi),MP10='direct v1 full frozen Epi64 read; v3 eligible values checked; no rescoring')
if not (R/'method_config.yaml').exists():(R/'method_config.yaml').write_text(yaml.safe_dump(config,sort_keys=False))
print(pd.DataFrame(manifest).groupby('section').readable.agg(['sum','count']),len(roi))
