from base import *
from scipy.stats import rankdata
sets=json.loads((O/'shared/programs.json').read_text());rows=[]
for s in XE:
 d=O/'xenium'/s;v=pd.read_csv(d/'ucell.tsv.gz',sep='\t').sample(n=64,random_state=123).set_index('id')
 x,ids,g=tenx(next((P/'highres_ST/highres_spatial_audit_v2_external/00_local_input_inventory/GSE326226_h5').glob('*_'+s+'_*h5')))
 jj=g.feature_type.eq('Gene Expression').to_numpy();bg=pd.Index(g.name[jj]);a=x[pd.Index(ids).get_indexer(v.index)][:,jj].toarray().astype(float)
 ranks=rankdata(-a,axis=1,method='average'); cap=min(1500,len(bg))
 for k,genes in sets.items():
  ix=bg.get_indexer(genes);assert ix.min()>=0;n=len(genes);r=np.minimum(ranks[:,ix],cap);minimum=n*(n+1)/2
  expected=np.maximum(0,1-(r.sum(1)-minimum)/(n*cap-minimum));error=np.max(np.abs(v[k].to_numpy()-expected));rows.append((s,k,len(v),len(bg),n,cap,error))
res=pd.DataFrame(rows,columns=['sample','program','cells_checked','background','members','maxRank','maximum_abs_error'])
put(res,O/'reports/UCell_native_check.tsv')
assert res.maximum_abs_error.max()<1e-10,res

