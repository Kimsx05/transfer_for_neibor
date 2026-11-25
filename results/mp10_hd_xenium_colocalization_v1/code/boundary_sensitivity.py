from base import *
rows=[]
for s in XE:
 d=O/'xenium'/s;z=pd.read_parquet(d/'labels.parquet');qc=pd.read_parquet(d/'boundary_geometry.parquet').reindex(z.index)
 assert qc.review_boundary.notna().all();ok=~qc.review_boundary.to_numpy()
 # Maintain original segmentation; exclude flagged polygon endpoints on original graph for sensitivity analysis.
 e=z.epi_eligible.to_numpy()&ok;score=z.Chol_biosynthesis_8.to_numpy();q=np.quantile(score[e],.95);high=(score>q);old=z.chol_group.eq('Chol-high Epi').fillna(False).to_numpy()
 node=pd.read_parquet(d/'graph_nodes.parquet');assert np.array_equal(node.id,z.index)
 for r in [50,200]:
  g=sparse.load_npz(d/f'graph_{r}.npz')[e];degree=np.asarray(g@ok.astype('i4')).ravel();eligible=(degree>0)&(node[f'coverage_{r}'].to_numpy()[e]>=.5)
  orig=pd.read_csv(d/f'effects_{r}.tsv',sep='\t').set_index('neighbor')
  for label,mask in [('Treg',z.supported=='Treg')]+[(k+' (all)',z.major==k) for k in ['Endothelial','Fibroblast','Myeloid']]:
   cnt=np.asarray(g@(ok&mask).astype('i4')).ravel();frac=np.divide(cnt,degree,out=np.full(len(degree),np.nan),where=eligible);delta=np.nanmean(frac[high[e]])-np.nanmean(frac[~high[e]])
   rows.append((s,r,label,int(e.sum()),q,int(high[e].sum()),int(np.sum(high[e]!=old[e])),orig.loc[label,'delta'],delta,delta-orig.loc[label,'delta']))
put(pd.DataFrame(rows,columns=['sample','radius','neighbor','epi_after_exclusion','q95','n_high_after_exclusion','rank_group_changes_among_retained','primary_delta','exclude_boundary_delta','difference']),O/'reports/boundary_sensitivity.tsv')

