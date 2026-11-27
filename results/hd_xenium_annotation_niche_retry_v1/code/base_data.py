from pathlib import Path
import sys,json,os
import numpy as np,pandas as pd
from scipy import sparse
import anndata as ad
O=Path(__file__).resolve().parents[1];P=O.parents[1];V=P/'results/mp10_hd_xenium_colocalization_v1';A=P/'highres_VisiumHD/analysis';H=A/'cell_annotation_v2'
sys.path.append(str(V/'code'))
from base import tenx,read_export
HD=['Cervilla','D1_Johnson'];XE=['A1R2','A1R4','A1R5','A1R6','A2R2']
M=json.loads((V/'shared/major_markers.json').read_text());tax=pd.read_csv(V/'shared/taxonomy.tsv',sep='\t').drop_duplicates('label').set_index('label')
def write(df,name):df.to_csv(O/'tables'/name,sep='\t',index=False)
def src(s):
 return list((P/'highres_VisiumHD/data').glob(('Cervilla' if s=='Cervilla' else 'GSE300692')+'/**/square_002um/filtered_feature_bc_matrix.h5'))[0]
def paths(s):return O/('hd' if s in HD else 'xenium')/s
def dataset(s,image_qc=True):
 if s in HD:
  a=ad.read_h5ad(paths(s)/'cells.h5ad')
  if image_qc and (paths(s)/'secondary_image_support.parquet').exists():
   v=pd.read_parquet(paths(s)/'secondary_image_support.parquet').reindex(a.obs_names)
   a.obs['secondary_image_suspect']=v.image_suspect.eq(True)
   a.obs['geometry_flag']=a.obs.geometry_flag|a.obs.secondary_image_suspect
  return a
 z=pd.read_parquet(V/'xenium'/s/'labels.parquet')
 f=pd.read_csv(V/'shared/xenium_inputs.tsv',sep='\t').query('sample==@s and role=="expression"').path.item()
 x,ids,g=tenx(f);ok=g.feature_type.eq('Gene Expression').values
 assert np.array_equal(ids,z.index) and g.name[ok].is_unique
 return ad.AnnData(x[:,ok].tocsr(),obs=z,var=pd.DataFrame(index=g.name[ok].values))
def norm(x):
 x=x.astype('f4');x=x.multiply((1e4/np.maximum(1,x.sum(1).A1))[:,None]).tocsr();x.data=np.log1p(x.data);return x
def eligible(a,s):
 if s in HD:return (a.obs.nFeature>=10)&(a.obs.n_counts>=10)&~a.obs.geometry_flag
 return (a.obs.nFeature>=10)&~a.obs.boundary_degenerate&~a.obs.area_outlier

