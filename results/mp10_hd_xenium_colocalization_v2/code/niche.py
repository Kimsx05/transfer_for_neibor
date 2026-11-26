from common import *
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.decomposition import TruncatedSVD
from scipy.spatial import cKDTree
from geometry import edges,pixel
import anndata as ad
def data(s):
 v,d=folders(s);z=labels(s);i=pd.read_parquet(d/'identity.parquet');ct=z.supported.copy()
 tnames=set(TAX.index[(TAX.major=='T/NK') & (TAX.index!='NK')])
 ct[ct.isin(tnames)]='T';ct[(ct=='T/NK')&i.reliable_T]='T';ct[(ct=='T/NK')&~i.reliable_T]='T/NK unresolved'
 if s in XE:ct[i.grade.eq('supported')]='T'
 xy=z[['x_um','y_um']].values
 a=ad.AnnData(sparse.csr_matrix((len(z),1),dtype='f4'),obs=pd.DataFrame({'ct':pd.Categorical(ct),'sample':pd.Categorical([s]*len(z))},index=z.index))
 a.obsm['spatial']=xy
 p=OUT/'niche'/s;p.mkdir(exist_ok=True);return v,d,z,i,ct,xy,a,p
def mender(s):
 from MENDER.MENDER import MENDER_single
 v,d,z,i,ct,xy,a,p=data(s);w=np.load(v/'window.npz');rc=pixel(xy,w['affine']);comp=pd.read_parquet(v/'graph_nodes.parquet').component.to_numpy('i4')
 def obj(a):
  m=MENDER_single(a);m.set_MENDER_para(nn_mode='radius',nn_para=25,n_scales=3,include_self=False,count_rep='s');return m
 if s in ['Cervilla','A1R2']:
  mid=xy[np.argmin(((xy-np.median(xy,axis=0))**2).sum(1))];sel=(abs(xy-mid)<180).all(1);o=obj(a[sel].copy());o.run_representation(group_norm=True);assert o.adata_MENDER.n_obs==sel.sum();(p/'native_mender_pilot.txt').write_text(str((mid,sel.sum(),o.adata_MENDER.shape)))
 model=obj(a);cats=ct.astype('category').cat.categories;cm=np.eye(len(cats),dtype='f4')[ct.astype('category').cat.codes];last=np.zeros_like(cm)
 for l,r in enumerate([25,50,75]):
  G=sparse.load_npz(v/'graph_50.npz') if r==50 else edges(xy,xy,rc,rc,w['mask'],comp,comp,r,True)
  accum=G@cm;ring=accum-last;assert ring.min()>=0;last=accum
  a.obsm['scale'+str(l)]=ring/np.maximum(1,ring.sum(1,keepdims=True))
 model.generate_ct_representation();X=model.adata_MENDER.X;metric=[];best=None
 for k in [4,6]:
  ys=[MiniBatchKMeans(n_clusters=k,n_init=3,batch_size=2048,random_state=sd).fit_predict(X) for sd in CFG['seeds']]
  ar=np.median([adjusted_rand_score(u,vv) for u,vv in itertools.combinations(ys,2)]);metric.append((k,ar))
  if best is None or ar>best[0]:best=(ar,ys[0])
 pd.Series(best[1],index=z.index,name='niche').to_frame().to_parquet(p/'MENDER.parquet');put(pd.DataFrame(metric,columns=['K','median_seed_ARI']),p/'MENDER_stability.tsv')
 tbl=pd.DataFrame({'id':z.index,'x':xy[:,0],'y':xy[:,1],'ct':ct});put(tbl,p/'niche_input.tsv.gz')
 kd=cKDTree(xy).query(xy,k=20)[0][:,-1];put(pd.DataFrame({'quantile':[.05,.5,.95],'k20_um':np.quantile(kd,[.05,.5,.95])}),p/'k20.tsv')
def ccfit(s):
 import torch,cellcharter as cc
 torch.set_num_threads(1)
 if torch.get_num_interop_threads()!=1:torch.set_num_interop_threads(1)
 v,d,z,i,ct,xy,a,p=data(s);x,genes,ids=expression(s);assert np.array_equal(ids,a.obs_names)
 a.obsp['spatial_connectivities']=sparse.load_npz(v/'graph_50.npz')
 drop=set(TOP+CORE+PAN+list(itertools.chain.from_iterable(json.loads((V/'shared/programs.json').read_text()).values())))
 rng=np.random.default_rng(CFG['seed']);train=np.sort(rng.choice(len(z),min(5000,len(z)),replace=False));rows=[]
 for name in ['all','masked']:
  lib=x.sum(1).A1
  ok=~np.isin(genes,list(drop)) if name=='masked' else np.ones(len(genes),bool);pop=x.getnnz(0);jj=np.flatnonzero(ok)[np.argsort(-pop[ok],kind='stable')[:1000]]
  xs=x[:,jj].astype('f4');xs=xs.multiply((1e4/np.maximum(1,lib))[:,None]).tocsr();xs.data=np.log1p(xs.data);svd=TruncatedSVD(12,random_state=CFG['seed']).fit(xs[train]);a.obsm['latent']=svd.transform(xs)
  cc.gr.aggregate_neighbors(a,n_layers=[0,1,2],sample_key='sample',connectivity_key='spatial_connectivities',use_rep='latent');X=a.obsm['X_cellcharter'].astype('f4')
  ys=[]
  for seed in CFG['seeds']:
   f=cc.tl.GaussianMixture(n_clusters=4,random_state=seed,batch_size=1024,trainer_params=dict(accelerator='cpu',devices=1,max_epochs=40,enable_progress_bar=False))
   f.fit(X[train]);ys.append(np.asarray(f.predict(X[train])));rows.append((name,seed,int(f.num_iter_),bool(f.converged_)))
   if seed==CFG['seed']:
    yp=np.concatenate([np.asarray(f.predict(X[h:h+4096])) for h in range(0,len(z),4096)]);pd.Series(yp,index=z.index,name='niche').to_frame().to_parquet(p/('CellCharter_'+name+'.parquet'))
  (p/('CellCharter_'+name+'_ARI.txt')).write_text(str(np.median([adjusted_rand_score(a1,b) for a1,b in itertools.combinations(ys,2)])))
 put(pd.DataFrame(rows,columns=['genes','seed','iterations','converged']),p/'CellCharter_stability.tsv')
if __name__=='__main__':
 for s in sys.argv[2:]:({'mender':mender,'cellcharter':ccfit}[sys.argv[1]])(s)
