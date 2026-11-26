from common import *
from scipy.stats import rankdata
from itertools import combinations
import json
P=R/'results/mp10_visium_h05_mapping_pilot_v1'
T=O/'tables'
KEY=['section_id','barcode']
MP='EpiMP10_Top20_AUCell'
H='H05_separate_transfer_score'
N='nonH05_separate_transfer_score'
GENES=R/'mp10_paper_distance_full_v1/tables/gene_sets.tsv.gz'
PATHWAYS={'Synthesis':'REACTOME_CHOLESTEROL_BIOSYNTHESIS','Homeostasis':'HALLMARK_CHOLESTEROL_HOMEOSTASIS','Fatty_acid':'HALLMARK_FATTY_ACID_METABOLISM'}
def join(a,b):
 assert not a.duplicated(KEY).any() and not b.duplicated(KEY).any()
 assert set(map(tuple,a[KEY].values))==set(map(tuple,b[KEY].values))
 return a.merge(b,on=KEY,validate='one_to_one',sort=False)
def corr(x,y):
 if len(x)<3 or np.ptp(x)==0 or np.ptp(y)==0:return np.nan
 return np.corrcoef(x,y)[0,1]
def stat(x,y,e,z):
 v=np.column_stack((x,y,e,z)).astype(float);v=v[np.isfinite(v).all(1)]
 n=len(v);missing=dict(n=n,pearson=np.nan,spearman=np.nan,resid_pearson=np.nan,resid_spearman=np.nan)
 if n<4:return missing
 a=v[:,:2];design=np.column_stack((np.ones(n),v[:,2:]))
 r=a-design@np.linalg.lstsq(design,a,rcond=None)[0]
 ok=np.all(np.ptp(r,axis=0)>1e-12*np.maximum(1,np.ptp(a,axis=0)))
 return dict(n=n,pearson=corr(*a.T),spearman=corr(*rankdata(a,axis=0).T),resid_pearson=corr(*r.T) if ok else np.nan,resid_spearman=corr(*rankdata(r,axis=0).T) if ok else np.nan)
def relate(d,x,y):
 return stat(d[x],d[y],d.Epithelial_fraction,np.log1p(d.nCount_RNA))
def square(d,s,method,cols):
 return [dict(section_id=s,method=method,left=x,right=y,**relate(d,x,y)) for x,y in combinations(cols,2)]
def ccol(t):return 'c2l_mean__'+t
def lcol(t):return 'log1p_c2l_mean__'+t
def scol(t):return 'TME_score__'+t

