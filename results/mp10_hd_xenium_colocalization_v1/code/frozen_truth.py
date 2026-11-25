from base import *
MEM=P/'results/mp10_immune_endothelial_mechanisms_v1/tables/pseudobulk_membership.tsv.gz'
def frozen_truth():
 m=pd.read_csv(MEM,sep='\t',usecols=['barcode','sample_uid','profile_type'])
 e=m[m.profile_type.eq('Epi_all')].set_index('barcode').copy()
 h=m[m.profile_type.eq('H05')].set_index('barcode')
 assert e.index.is_unique and h.index.is_unique and h.index.isin(e.index).all()
 e['H05_truth']=np.where(e.index.isin(h.index),'H05','nonH05');return e[['sample_uid','H05_truth']]
if __name__=='__main__':
 e=frozen_truth();base=pd.read_csv(P/'mp10_sc_c2l_final_annotation_v1/tables/analysis_cells.tsv.gz',sep='\t').set_index('barcode');q=pd.read_csv(P/'mp10_cellchat_high_definition_sensitivity_v1/tables/specimen_quantiles_coverage.tsv',sep='\t').set_index('sample_uid')
 v=pd.read_csv(O/'shared/reference/membership.tsv.gz',sep='\t').set_index('barcode');expected=e.H05_truth.reindex(v.index)
 assert (expected.fillna('n/a')==v.H05_truth.fillna('n/a')).all()
 e['MP10_cell']=base.loc[e.index,'MP10_cell'];rows=[]
 for k,t in e.groupby('sample_uid'):
  n=int(q.loc[k,'n_H05']);cut=t.MP10_cell.nlargest(n).min()
  disagreement=t.H05_truth.eq('H05')!=(t.MP10_cell>=cut)
  assert all(abs(t.loc[disagreement,'MP10_cell']-q.loc[k,'Q95'])<1e-12)
  rows.append((k,len(t),n,int(t.H05_truth.eq('H05').sum()),q.loc[k,'Q95'],cut,int(disagreement.sum())))
 e.to_parquet(O/'shared/frozen_epithelial_truth.parquet')
 put(pd.DataFrame(rows,columns=['sample','eligible_epi','original_n_H05','member_n_H05','frozen_Q95','minimum_rank_reconstructed_high','threshold_serialization_disagreements']),O/'shared/H05_membership_check.tsv')
 print('Exact full-set and model member agreement',len(e),(e.H05_truth=='H05').sum(),len(v),flush=True)

