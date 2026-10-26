from pathlib import Path
import sys,json
P=Path('/home/data/t070721/codex_workspace/Bladder Metabolism')
sys.path=[q for q in sys.path if 'dist-packages' not in q]
sys.path[:0]=[str(P/'spatial_ecology_top20_v2/runtime'),str(P/'cell2location_setup_v1/cell2location_env/lib/python3.12/site-packages')]
import pandas as pd
import numpy as np
O=P/'mp10_subtype_screen_round1'; T=O/'tables'
base=Path('/home/data/t070721/Metabolistic_scRNA/Multi2')
def save(x,name): x.to_csv(T/name,sep='\t',index=False,na_rep='NA')
m=pd.read_csv(P/'subtype/v2/final_annotation_pre_reference/frozen_annotation_barcode_mapping_v2.csv',low_memory=False)
assert m.barcode.is_unique
save(m.groupby(['dataset','annotation_status_v2']).size().reset_index(name='n_cells'),'annotation_inclusion.tsv')
# Existing contamination removals apply; reference-training exclusions do NOT remove cells.
m=m[m.annotation_status_v2.ne('excluded_before_round2')].copy()
m['sample_uid']=m.dataset+'::'+m['orig.ident'];m['major']=m.major_celltype.replace({'Cycling':'Epithelial'})
m['subtype']=m.subtype_v2.fillna(m.major+'_unclassified')
m.loc[m.major.eq('Epithelial'),'subtype']='Epithelial_unclassified'
cl=pd.read_csv(P/'cellchat/sample_by_clinical_group.csv')
assert cl.specimen_id.is_unique
cl=cl.rename(columns={'patient_uid':'patient_uid_recorded','fine_stage':'T_raw'})
s=m[['dataset','orig.ident','sample_uid']].drop_duplicates().merge(cl.drop(columns=['dataset','orig.ident']),left_on='orig.ident',right_on='specimen_id',how='left',validate='one_to_one')
assert s.patient_uid_recorded.notna().all()
s['patient_uid_verified']=np.nan
s['tissue_status']=s.clinical_group.map({'Normal':'Normal','MIBC':'Tumor','NMIBC':'Tumor'})
s['T_category']=s.T_raw.str.extract(r'^(?:p|c)?(Ta|Tis|T[1-4])(?:[abc])?$')[0]
s['T_basis']=np.where(s.T_raw.fillna('').str.startswith('pT'),'pT',np.where(s.T_raw.fillna('').str.startswith('cT'),'cT',np.where(s.T_category.notna(),'T_unspecified_basis','unknown')))
s['normal_origin']=np.nan;s['treatment_status']=np.nan;s['enrichment_status']='not_documented'
m=m.merge(s[['sample_uid','tissue_status']],on='sample_uid',validate='many_to_one')
# Reuse the old per-cell rough-cluster assignment and frozen malignant-like-main map.
cm=pd.read_csv(base/'3.3 recognize malignant/rough_cluster_infercnv_score_cell_metadata.csv',usecols=['cell','rough_cluster'])
status=pd.read_csv(base/'3.4 malignant cell process/NMF1/rough_cluster_malignant_like_status_summary.csv')
cm=cm.merge(status[['rough_cluster','malignant_like_status']],on='rough_cluster',how='left',validate='many_to_one')
m=m.merge(cm.rename(columns={'cell':'barcode'}),on='barcode',how='left',validate='one_to_one')
m['malignant']=m.major.eq('Epithelial') & m.tissue_status.eq('Tumor') & m.malignant_like_status.eq('malignant_like_main')
m['communication_type']=m.subtype
m.loc[m.major.eq('Epithelial'),'communication_type']='Other_epithelial_unresolved'
m.loc[m.major.eq('Epithelial') & m.tissue_status.eq('Normal'),'communication_type']='Normal_epithelial'
m.loc[m.malignant,'communication_type']='Malignant_like_epithelial'
fidelity=P/'neighborhood_mp10_calibration_v5/01_score_fidelity'
expected='HMGCS1 MSMO1 INSIG1 FDFT1 SQLE FDPS IDI1 ACSS2 ELOVL5 LDLR MVD SCD HMGCR ACAT2 SLC20A1 ACSL3 ACSL1 ADIPOR2 MXD1 NPC1'.split()
actual=(fidelity/'top20.txt').read_text().split();assert set(expected)==set(actual) and len(actual)==20
score=pd.read_parquet(fidelity/'cell_level_fidelity.parquet',columns=['cell_id','Top20_score'])
assert score.cell_id.is_unique
m=m.merge(score.rename(columns={'cell_id':'barcode','Top20_score':'MP10_cell'}),on='barcode',how='left',validate='one_to_one')
counts=m.groupby('sample_uid').agg(n_total=('barcode','size'),n_malignant=('malignant','sum'))
counts['n_TME']=m[m.major.ne('Epithelial')].groupby('sample_uid').size()
sc=m[m.malignant].groupby('sample_uid').MP10_cell.agg(n_MP10_scored='count',MP10_sample='mean',MP10_median='median')
s=s.merge(counts,on='sample_uid').merge(sc,on='sample_uid',how='left')
s['n_MP10_scored']=s.n_MP10_scored.fillna(0).astype(int)
print(s[['dataset','orig.ident','n_malignant','n_MP10_scored']].to_string(index=False))
print('Coverage quantiles',s.loc[s.tissue_status.eq('Tumor'),'n_MP10_scored'].quantile([0,.1,.25,.5,.75,.9,1]).to_dict())
# Fixed BEFORE inspecting any association: retain project precedent of >=30 target cells.
s['MP10_eligible']=s.tissue_status.eq('Tumor') & s.n_MP10_scored.ge(30)
s['MP10_threshold']=s.dataset.map(s[s.MP10_eligible].groupby('dataset').MP10_sample.median())
s['MP10_group']=pd.Series(index=s.index,dtype='object')
s.loc[s.MP10_eligible,'MP10_group']=np.where(s.loc[s.MP10_eligible,'MP10_sample']>s.loc[s.MP10_eligible,'MP10_threshold'],'high','low')
s['MP10_source']='cell_level_fidelity.parquet:Top20_score; frozen average-tie AUCell-like; mean of existing malignant_like_main cells'
s['MP10_coverage_rule']='>=30 malignant-like epithelial cells with finite stored score'
s['malignant_source']='frozen rough_cluster_malignant_like_status_summary.csv + per-cell rough_cluster; tumor only'
parents=m.groupby(['sample_uid','major']).size().unstack(fill_value=0).add_prefix('n_parent_')
s=s.merge(parents,on='sample_uid')
map_=m[m.major.ne('Epithelial')][['major','subtype']].drop_duplicates().sort_values(['major','subtype'])
assert map_.subtype.is_unique
save(s,'samples.tsv');save(map_,'subtype_parent.tsv')
save(m[['barcode','sample_uid','dataset','major','subtype','communication_type','malignant','MP10_cell','annotation_status_v2','spatial_reference_include']],'cells.tsv.gz')
save(pd.DataFrame({'gene':expected,'rank':range(1,21)}),'frozen_MP10_Top20.tsv')
rows=[]
for comp,col,groups in [('Normal_Tumor','tissue_status',['Normal','Tumor']),('T1_T2','T_category',['T1','T2']),('MP10_low_high','MP10_group',['low','high'])]:
 z=s[s[col].isin(groups)]
 for (ds,g),a in z.groupby(['dataset',col]):rows.append(dict(comparison=comp,dataset=ds,group=g,n_specimens=len(a),n_patient_ids_recorded=a.patient_uid_recorded.nunique(),n_patients_verified=0,n_missing_T=a.T_category.isna().sum(),n_missing_MP10=(~a.MP10_eligible).sum(),n_missing_treatment=a.treatment_status.isna().sum()))
save(pd.DataFrame(rows),'group_coverage.tsv')
# Full zero-completed compositions; absent denominator is NA.
rows=[]
for _,a in s.iterrows():
 z=m[m.sample_uid.eq(a.sample_uid)]
 for denom,types in [('all_cells',m.major.unique()),('TME',m.loc[m.major.ne('Epithelial'),'major'].unique())]:
  den=a.n_total if denom=='all_cells' else a.n_TME
  for ct in types:rows.append(dict(sample_uid=a.sample_uid,dataset=a.dataset,patient_uid=a.patient_uid_recorded,level='major',major=ct,subtype=ct,denominator_type=denom,numerator=int(z.major.eq(ct).sum()),denominator=den))
 for _,b in map_.iterrows():
  for denom,den in [('parent',int(z.major.eq(b.major).sum())),('TME',a.n_TME)]:rows.append(dict(sample_uid=a.sample_uid,dataset=a.dataset,patient_uid=a.patient_uid_recorded,level='subtype',major=b.major,subtype=b.subtype,denominator_type=denom,numerator=int(z.subtype.eq(b.subtype).sum()),denominator=den))
p=pd.DataFrame(rows);p['proportion']=p.numerator/p.denominator.replace(0,np.nan)
save(p,'proportions.tsv.gz')
print('samples',len(s),'recorded patients',s.patient_uid_recorded.nunique(),'cells',len(m),'MP10 eligible',s.MP10_eligible.sum())
print(s.groupby('dataset').agg(specimens=('sample_uid','size'),patients_recorded=('patient_uid_recorded','nunique'),MP10_eligible=('MP10_eligible','sum')))
