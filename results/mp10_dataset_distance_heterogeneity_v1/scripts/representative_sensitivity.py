from heterogeneity import *
def main():
 reg=read('endpoints.tsv.gz');u=read('unit_endpoints.tsv.gz');u=u[u.identity_level.eq('analysis_unit')&u.value.notna()];mats,meta=matrix_inputs();out=[]
 name='C1_CAF_full';base,kind=mats[name];datasets=sorted(meta.loc[base.index,'dataset_id'].unique())
 variants=[('UCell_RNA','common_support_0_1000','r',10),('AddModuleScore_SCT_data','common_support','r',10),('AddModuleScore_SCT_data','full_range','r',10),('AddModuleScore_SCT_data','common_support_0_1000','partial_r',10),('AddModuleScore_SCT_data','common_support_0_1000','r',30)]
 for method,ran,st,nmin in variants:
  r=reg[reg.module.eq('C1')&reg.score_method.eq(method)&reg.range_type.eq(ran)&reg.statistic.eq(st)&reg.signature.eq('full')&reg.core&reg.candidate_block.eq('CAF')];v=u[u.endpoint_id.isin(r.endpoint_id)&u.n_min.ge(nmin)].merge(r[['endpoint_id','canonical_program','candidate_type']],on='endpoint_id');v['feature']=v.canonical_program+' / '+v.candidate_type;alt=v.pivot(index='record_id',columns='feature',values='value');idx=base.index.intersection(alt.index);features=base.columns.intersection(alt.columns);pa=processed(base.loc[idx,features],'r',datasets,meta);pb=processed(alt.loc[idx,features],'r',datasets,meta)
  if pa is None or pb is None:continue
  ii=pa['index'].intersection(pb['index']);ff=pa['features'].intersection(pb['features']);pa=processed(base.loc[ii,ff],'r',datasets,meta);pb=processed(alt.loc[ii,ff],'r',datasets,meta)
  if pa is None or pb is None:continue
  k=max(cluster(pa['X']),key=lambda x:x[1])[0];la=cut_tree(linkage(pa['X'],method='ward'),n_clusters=k).ravel();lb=cut_tree(linkage(pb['X'],method='ward'),n_clusters=k).ravel()
  if nmin==30:
   original=processed(base,'r',datasets,meta);k=max(cluster(original['X']),key=lambda x:x[1])[0];original_lab=cut_tree(linkage(original['X'],method='ward'),n_clusters=k).ravel();la=original_lab[original['index'].get_indexer(pb['index'])];lb=cut_tree(linkage(pb['X'],method='ward'),n_clusters=k).ravel();pa['coords']=original['coords'][original['index'].get_indexer(pb['index'])]
  out.append(dict(input_family=name,variant=method+'|'+ran+'|'+st+'|n'+str(nmin),n_paired_units=len(pa['index']),n_common_features=len(pa['features']),K=k,adjusted_rand=ari(la,lb),absolute_PC1_correlation=abs(stats.pearsonr(pa['coords'][:,0],pb['coords'][:,0]).statistic)))
 # Top100 replaces Top30, not an extra duplicate theme in a full matrix.
 for name in ['B1_R0_partial_r_full','B3_Epi_distance_partial_r_full']:
  base,kind=mats[name];mo=name.split('_partial')[0];r=reg[reg.module.eq(mo)&reg.statistic.eq('partial_r')&reg.signature.eq('Top100')];v=u[u.endpoint_id.isin(r.endpoint_id)];alt=base.copy();key='PURAM_PEMT_TOP30 / none'
  if len(v) and key in alt:
   alt[key]=v.set_index('record_id').value.reindex(alt.index);ds=sorted(meta.loc[alt.index,'dataset_id'].unique());pa=processed(base,kind,ds,meta);pb=processed(alt,kind,ds,meta)
   if pa and pb:
    k=max(cluster(pa['X']),key=lambda x:x[1])[0];la=cut_tree(linkage(pa['X'],method='ward'),n_clusters=k).ravel();lb=cut_tree(linkage(pb['X'],method='ward'),n_clusters=k).ravel();out.append(dict(input_family=name,variant='Puram_Top100_replaces_Top30',n_paired_units=len(pa['index']),n_common_features=len(pa['features']),K=k,adjusted_rand=ari(la,lb),absolute_PC1_correlation=abs(stats.pearsonr(pa['coords'][:,0],pb['coords'][:,0]).statistic)))
 save(out,'representative_method_range_sensitivity.tsv');print('representative sensitivity complete',flush=True)
if __name__=='__main__':main()
