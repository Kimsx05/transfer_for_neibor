from common import *
from niche import data
import squidpy as sq
def run(s):
 v,d,z,i,ct,xy,a,p=data(s);a.obsp['spatial_connectivities']=sparse.load_npz(v/'graph_50.npz')
 # HD component attribute has its own graph statistics. Mutually exclusive subtype identity remains v1.
 lab=z.supported.copy()
 if s in XE:
  lab[lab=='Treg']='T';lab[i.grade.eq('supported')]='Treg'
 a.obs['subtype']=pd.Categorical(lab);res=sq.gr.nhood_enrichment(a,'subtype',n_perms=1000,n_jobs=1,rng=CFG['seed'],copy=True,handle_nan='keep',show_progress_bar=False)
 k=a.obs.subtype.cat.categories;ix=pd.MultiIndex.from_product([k,k],names=['from','to']);pd.DataFrame({'observed':res.counts.ravel(),'Z':res.zscore.ravel()},index=ix).to_parquet(d/'squidpy_50.parquet')
 # Native coordinate co-occurrence in bounded complete spatial field, preselected at tissue median.
 mid=xy[np.argmin(((xy-np.median(xy,axis=0))**2).sum(1))];sel=(abs(xy-mid)<=320).all(1);a=a[sel].copy()
 c=sq.gr.co_occurrence(a,'subtype',spatial_key='spatial',interval=np.array([0,25,50,100,150,200],dtype=float),copy=True)
 np.savez_compressed(d/'squidpy_native_distance.npz',values=c[0],interval=c[1],types=a.obs.subtype.cat.categories.astype(str).to_numpy(dtype=str),center=mid)
if __name__=='__main__':
 for s in sys.argv[1:]:run(s)
