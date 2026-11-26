from shared import *
from scipy.stats import rankdata
sets=read(O/'gene_sets.tsv')
for sec in S:
 p=O/'.work'/sec
 genes=pd.Index((p/'genes').read_text().splitlines());ids=(p/'ids').read_text().splitlines()
 ptr=np.fromfile(p/'p','<i4');index=np.memmap(p/'i',dtype='<i4',mode='r');v=np.memmap(p/'x',dtype='<i4',mode='r')
 K=int(np.ceil(.05*len(genes)))
 group={f'{mod}__{variant}':genes.get_indexer(genes.intersection(g.gene_symbol)) for (mod,variant),g in sets.groupby(['module','variant'])}
 denominator={k:sum(K-r for r in range(1,min(len(ix),K)+1)) for k,ix in group.items()}
 w=np.zeros(len(genes));rows=[]
 for j,id in enumerate(ids):
  sl=slice(ptr[j],ptr[j+1]);indices=index[sl];values=v[sl]
  w[indices]=np.maximum(0,K-rankdata(-values,method='average'))
  # A tied zero has mean rank above K, including when no gene is detected.
  result={k:w[ix].sum()/denominator[k] if len(ix)>=5 else np.nan for k,ix in group.items()}
  if j in [0,len(ids)//2,len(ids)-1]:
   dense=np.zeros(len(genes));dense[indices]=values
   ranks=rankdata(-dense,method='average')
   for k,ix in group.items():
    if len(ix)>=5:assert abs(result[k]-np.maximum(0,K-ranks[ix]).sum()/denominator[k])<1e-12
  rows.append(dict(section_id=sec,barcode=id,**result));w[indices]=0
 save(rows,f'spots/{sec}.function.tsv.gz')

