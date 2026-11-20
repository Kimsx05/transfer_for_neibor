from base6 import *
parts=[];masks=[]
for s in sections:
 o=R/s;d=pd.read_parquet(o/'cells_scored.parquet')
 ov=pd.read_parquet(o/'nuclear_overlap.parquet').query('substantial');dominant=ov.sort_values(['old_id','intersection','new_id']).drop_duplicates('old_id',keep='last')
 pp=pd.read_parquet(o/'priority_old.parquet').rename(columns={'id':'old_id'});z=pp.merge(dominant,on='old_id',how='left').merge(d[['id','outcome','mixed_v6','A','B','extended','p_T','p_Treg_given_T','nCount','FOXP3','direct_width','other_core_width','nuclear_action','review_action']],left_on='new_id',right_on='id',how='left',suffixes=('_v5','_v6'));z['outcome']=z.outcome.fillna('no_substantial_successor');put(z,o/'priority_dominant_successor.tsv')
 for tag,sel in [('former_FOXP3_mixed',z.mixed&z.FOXP3_v5.gt(0)),('former_FOXP3_nonmixed',~z.mixed&z.FOXP3_v5.gt(0)),('former_auxiliary',z.group.eq('RNA_candidate'))]:
  v=z[sel].groupby('outcome').size().reset_index(name='old_objects');v['section']=s;v['old_group']=tag;parts.append(v)
 old=pd.read_parquet(o/'nuclear_audit_old.parquet');v=old.action.value_counts().rename_axis('action').reset_index(name='objects');v['denominator']='old_masks';v['section']=s;masks.append(v)
 v=d.nuclear_action.value_counts().rename_axis('action').reset_index(name='objects');v['denominator']='final_masks';v['section']=s;masks.append(v)
 L=pd.read_parquet(o/'assignment.parquet',columns=['source_bin_index','id','usable']);org=pd.read_parquet(S/s/'bin_to_cell.parquet',columns=['source_bin_index','cell_label','raw_counts']).sort_values('source_bin_index')
 assert np.array_equal(org.source_bin_index,L.source_bin_index)
 L['old_id']=org.cell_label.to_numpy();L['UMI']=org.raw_counts.to_numpy();L.groupby(['old_id','id','usable']).agg(bins=('UMI','size'),UMI=('UMI','sum')).reset_index().to_parquet(o/'RNA_old_new_by_source.parquet',index=False)
tab(pd.concat(parts,ignore_index=True),'old_candidate_disposition.tsv');tab(pd.concat(masks,ignore_index=True),'mask_actions.tsv')
