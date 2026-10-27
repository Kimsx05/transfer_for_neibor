from common import *
import shutil
s=pd.read_csv(V/'tables/samples.tsv',sep='\t');m=pd.read_csv(V/'tables/cells.tsv.gz',sep='\t',low_memory=False)
# Existing frozen annotation and contamination decisions; never select old malignant labels.
m=m.drop(columns=['communication_type','malignant','MP10_cell'])
s=s.drop(columns=[c for c in s if c.startswith('MP10_') or c in ['n_malignant','n_MP10_scored','malignant_source']])
m=m.merge(s[['sample_uid','tissue_status','patient_uid_recorded']],on='sample_uid',validate='many_to_one')
f=P/'neighborhood_mp10_calibration_v5/01_score_fidelity'
a=pd.read_parquet(f/'cell_level_fidelity.parquet',columns=['cell_id','Top20_score']).rename(columns={'cell_id':'barcode','Top20_score':'MP10_cell'})
m=m.merge(a,on='barcode',how='left',validate='one_to_one');m.loc[m.major.ne('Epithelial'),'MP10_cell']=np.nan
m['score_source']=np.where(m.major.eq('Epithelial')&m.MP10_cell.notna(),'reused',np.where(m.major.eq('Epithelial'),'unavailable','not_applicable'))
save(m,'cells_pre_score.tsv.gz');save(s,'samples_base.tsv')
e=m[m.major.eq('Epithelial')].copy();save(e[e.MP10_cell.isna()][['barcode','sample_uid','dataset']],'score_fill_cells.tsv')
# Spread validation across datasets and original score range; compare counts against old frozen cache later.
v=e[e.MP10_cell.notna()].sort_values('MP10_cell').groupby('dataset',group_keys=False).apply(lambda d:d.iloc[np.unique(np.linspace(0,len(d)-1,12).astype(int))],include_groups=False).reset_index()
# pandas apply removed dataset: reconstruct from original m
v=v.drop(columns=['index'],errors='ignore').merge(m[['barcode','dataset']],on='barcode') if 'dataset' not in v else v
save(v[['barcode','sample_uid','dataset','MP10_cell']],'score_validation_cells.tsv')
for name in ['subtype_parent.tsv','cell_colors.tsv','CellChat_DB_interactions.tsv','CellChat_DB_complexes.tsv','CellChat_gene_coverage.tsv','package_versions.tsv','triMean_acceleration_validation.tsv']:
 shutil.copy2(V/'tables'/name,T/name)
shutil.copy2(V/'tables/proportions.tsv.gz',T/'proportions_v2.tsv.gz')
for name in ['fast_trimean.cpp','speckle_getTransformedProps_official.R']:
 shutil.copy2(V/'scripts'/name,O/'scripts'/name)
(O/'scripts/fast_trimean.R').write_text((V/'scripts/fast_trimean.R').read_text().replace(str(V),str(O)))
save(pd.DataFrame({'gene':TOP}),'frozen_MP10_Top20.tsv')
print('All epi',len(e),'reuse',e.MP10_cell.notna().sum(),'fill',e.MP10_cell.isna().sum(),flush=True)
