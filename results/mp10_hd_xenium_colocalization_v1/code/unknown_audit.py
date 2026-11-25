from base import *
rules=runpy.run_path(str(A/'cell_annotation_v2_1/code/programs.py'));tx=pd.read_csv(O/'shared/taxonomy.tsv',sep='\t')
rows=[]
for s in XE:
 d=O/'xenium'/s;z=pd.read_parquet(d/'labels.parquet');x,ids,features=tenx(next((P/'highres_ST/highres_spatial_audit_v2_external/00_local_input_inventory/GSE326226_h5').glob('*_'+s+'_*h5')))
 assert np.array_equal(ids,z.index);keep=features.feature_type.eq('Gene Expression').to_numpy();x=x[:,keep];genes=pd.Index(features.name[keep])
 for major in ['Epithelial','T/NK','Endothelial','Fibroblast','Myeloid']:
  markers=set()
  if major=='Epithelial':
   for k in ['Epithelial_basal','Epithelial_luminal','Epithelial_squamous']:markers.update(rules['STATE'][k])
  else:
   for v in tx.loc[tx.major==major,'anchors']:markers.update(v.split())
  markers=sorted((markers-set(rules['MAJOR'][major]))&set(genes));n=x[:,genes.get_indexer(markers)].getnnz(axis=1)
  unknown=(z.major=='Unknown')&z.major_candidate.eq(major)
  confident=unknown&(z.major_support>=.8)&((z.major_support-z.major_runner)>=.3)
  rows.append((s,major,int(unknown.sum()),int(confident.sum()),int((confident&(n>=1)).sum()),int((confident&(n>=2)).sum()),len(markers),','.join(markers)))
put(pd.DataFrame(rows,columns=['sample','candidate_major','unknown_with_candidate','candidate_ge_0.8_margin_ge_0.3','one_extra_marker','two_extra_markers','measured_extra_markers','markers']),O/'reports/unknown_marker_review.tsv')

