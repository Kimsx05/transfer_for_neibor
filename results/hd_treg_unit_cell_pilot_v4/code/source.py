from base import *
from scipy.spatial import cKDTree
kinds=['T8_from2','T16_from2','cell_nucleus','cell_conservative','existing_cell8']
def contextual(s):
 root=R/s;src=pd.read_parquet(root/'sources.parquet');native=pd.read_parquet(V/s/'bin_evidence.parquet',columns=['bin_id','x_um','y_um','component','array_row','array_col','type_blca','type_hpca','eligible_pool'])
 e=pd.read_parquet(V/s/'Epi64_regions.parquet',columns=['region_id','component','x_um','y_um','measurement_eligible','MP10','H05','H10']).query('measurement_eligible')
 em=pd.read_parquet(V/s/'Epi64_members.parquet').set_index('bin_id').region_id
 e_trees={c:(v.reset_index(drop=True),cKDTree(v[['x_um','y_um']])) for c,v in e.groupby('component')}
 cov_trees={c:cKDTree(v[['x_um','y_um']])for c,v in native.groupby('component')}
 p=native[native.bin_id.isin(src.parent8.dropna())].copy();p['epi64']=p.bin_id.map(em)
 for size in ['T16','T32']:
  mem=pd.read_parquet(V/s/(size+'_members.parquet')).set_index('bin_id').region_id
  p[size+'_id']=p.bin_id.map(mem)
  z=pd.read_parquet(V/s/(size+'_regions.parquet'),columns=['region_id','measurement_eligible','Treg_main','evidence']).set_index('region_id')
  for col in z.columns:p[size+'_'+col]=p[size+'_id'].map(z[col])
 v3=pd.read_parquet(Z/s/'T_regions.parquet').set_index('region_id')
 for col in ['strict','FOXP3_supported']:p['T32_'+col]=p.T32_id.map(v3[col])
 if s=='Cervilla':
  prior=rd(O/'tables/bin8_rctd_results.tsv.gz').set_index('bin_id')
  for col in ['spot_class','first_type','second_type']:p['RCTD_'+col]=p.bin_id.map(prior[col])
 p.to_parquet(root/'historical_native8_context.parquet',index=False)
 old=pd.read_csv(A/'hd_annotation_diagnosis_singler_v1'/s/'diagnosis.tsv.gz',sep='\t',usecols=['cell_id','singler_blca_labels','singler_hpca_labels','final_major'])
 if s=='Cervilla':
  fine=pd.read_parquet(O/'tables/Cervilla_cell_annotation.parquet',columns=['cell_id','fine_prediction','fine_pruned']);old=old.merge(fine,on='cell_id',how='left')
 nuc=[];res=[]
 for kind in kinds:
  obj=pd.read_parquet(root/(kind+'.parquet'));mem=pd.read_parquet(root/(kind+'_members.parquet'));mem['epi_id']=mem.parent8.map(em)
  assert mem.barcode.is_unique
  if kind.startswith('cell') or kind=='existing_cell8':
   if not nuc:ids=obj[['id','cell_id','roi','x_um','y_um']].merge(old,on='cell_id',how='left')
   nuc.append(obj.set_index('id')[['nCount','nFeature','T_UMI','FOXP3','IL2RA','CTLA4','IKZF2','Epi_UMI','Fib_UMI','Myeloid_UMI']].add_suffix('_'+kind))
  group=mem.groupby('id');comps=group.component.agg(lambda v: v.dropna().iloc[0] if v.nunique()==1 else np.nan);nmiss=group.parent8.apply(lambda v:v.isna().sum())
  shared=mem.dropna(subset=['epi_id']).groupby('id').epi_id.agg(set).to_dict()
  for c,units in obj.groupby(obj.id.map(comps).fillna(-1)):
   xy=units[['x_um','y_um']].to_numpy();v,tree=e_trees.get(c,(None,None));hits=tree.query_ball_point(xy,200)if tree is not None else [[] for _ in xy]
   coverage={r:cov_trees[c].query_ball_point(xy,r,return_length=True)if c in cov_trees else np.zeros(len(xy)) for r in [50,200]}
   if v is not None:coords=v[['x_um','y_um']].to_numpy();ei=v.region_id.to_numpy();h05=v.H05.to_numpy();h10=v.H10.to_numpy()
   for j,(q,h) in enumerate(zip(units.itertuples(),hits)):
    h=np.array(h,int);dist=np.linalg.norm(coords[h]-xy[j],axis=1)if len(h)else np.empty(0);bad=np.array([i in shared.get(q.id,set())for i in ei[h]],bool)if len(h)else np.empty(0,bool)
    for r in [50,200]:
     take=dist<=r;ok=take&~bad;unknown=c==-1 or nmiss.get(q.id,0)>0
     res.append(dict(section=s,branch=kind,id=q.id,radius_um=r,component=c if c!=-1 else np.nan,source_uncertain=unknown,eligible_epi_before=int(take.sum())if c!=-1 else np.nan,shared_excluded=int((bad&take).sum())if c!=-1 else np.nan,independent_epi=int(ok.sum())if not unknown else np.nan,H05_neighbors=int(h05[h[ok]].sum())if len(h)and not unknown else (np.nan if unknown else 0),H10_neighbors=int(h10[h[ok]].sum())if len(h)and not unknown else (np.nan if unknown else 0),nearest_independent_um=dist[ok].min()if ok.any()and not unknown else np.nan,supplied_same_component_centers=coverage[r][j]if c!=-1 else np.nan,independent_region_ids='|'.join(ei[h[ok]])if len(h)and not unknown else '',unmapped_parent_bins=nmiss.get(q.id,0)))
 audit=ids.set_index('id').join(pd.concat(nuc,axis=1)).reset_index();audit['section']=s
 for b,k in [('nucleus','cell_nucleus'),('4um','cell_conservative'),('8um','existing_cell8')]:
  audit['FOXP3_and_T_'+b]=(audit['FOXP3_'+k]>0)&(audit['T_UMI_'+k]>0)
 for col in ['FOXP3','T_UMI','nCount','Epi_UMI','Fib_UMI']:
  assert (audit[col+'_cell_nucleus']<=audit[col+'_cell_conservative']).all() and (audit[col+'_cell_conservative']<=audit[col+'_existing_cell8']).all()
 return audit,pd.DataFrame(res)
if __name__=='__main__':
 audit=[];cov=[]
 for s in sections:
  a,c=contextual(s);audit.append(a);cov.append(c)
 put(pd.concat(audit),'cell_assignment_audit.tsv');out=pd.concat(cov);out.to_parquet(R/'source_overlap_and_coverage.parquet',index=False)
 put(out.groupby(['section','branch','radius_um']).agg(units=('id','size'),source_uncertain=('source_uncertain','sum'),shared_pairs=('shared_excluded','sum'),independent_pairs=('independent_epi','sum'),units_with_opportunity=('independent_epi',lambda x:(x>0).sum())).reset_index(),'source_overlap_and_coverage.tsv')
