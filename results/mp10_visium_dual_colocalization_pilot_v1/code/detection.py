from shared import *
from scipy.sparse import csr_matrix
rows=[]
g=read(O/'gene_sets.tsv')
for s in S:
 p=O/'.work'/s
 genes=pd.Index((p/'genes').read_text().splitlines());ids=pd.Index((p/'ids').read_text().splitlines())
 x=csr_matrix((np.memmap(p/'x',dtype='<i4',mode='r'),np.memmap(p/'i',dtype='<i4',mode='r'),np.fromfile(p/'p',dtype='<i4')),shape=(len(ids),len(genes)))
 f=join(read(O/f'spots/{s}.base.tsv.gz')[KEY+['query']],read(O/f'spots/{s}.function.tsv.gz')).set_index('barcode').loc[ids]
 ok=f['query'].to_numpy(bool);x=x[ok].tocsc()
 for (mod,v),members in g.groupby(['module','variant']):
  used=genes.intersection(members.gene_symbol);mat=x[:,genes.get_indexer(used)]
  ng=mat.getnnz(axis=1);scores=f.loc[ok,mod+'__'+v]
  rows.append(dict(section_id=s,module=mod,variant=v,n=len(ng),zero_UMI_spots=(ng==0).sum(),median_detected_members=np.median(ng),zero_score_spots=(scores==0).sum(),zero_score_fraction=(scores==0).mean(),median_UMI=np.median(np.asarray(mat.sum(axis=1)).ravel())))
save(rows,'tables/functional_detection.tsv')

