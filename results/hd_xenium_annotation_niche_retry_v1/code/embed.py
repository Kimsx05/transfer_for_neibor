from base_data import *
import torch,scvi
torch.set_num_threads(3);torch.set_num_interop_threads(1)
platform=sys.argv[1];samples=XE if platform=='xenium' else HD
scvi.settings.seed=1729;scvi.settings.dl_num_workers=0
xs=[];names=[];ss=[];lookup=[];common=None
for s in samples:
 a=dataset(s,image_qc=False);ok=eligible(a,s).to_numpy()
 g=pd.Index(a.var_names)
 common=g if common is None else common.intersection(g,sort=False)
 idx=np.flatnonzero(ok); xs.append((a.X[idx].tocsr(),g));names.extend(s+'::'+a.obs_names[idx]);ss.extend([s]*len(idx));lookup.append((s,idx,a.n_obs));del a
xs=[x[:,g.get_indexer(common)] for x,g in xs];x=sparse.vstack(xs,format='csr');del xs
# Frequent, variable genes selected using sparse first/second moments after depth normalization.
y=norm(x);mu=y.mean(0).A1;v=y.power(2).mean(0).A1-mu**2
j=np.sort(np.argsort(v/(.1+mu))[-1200:]);del y
a=ad.AnnData(x[:,j],obs=pd.DataFrame({'sample':pd.Categorical(ss)},index=names),var=pd.DataFrame(index=common[j]));del x
train_object=a
if platform=='hd':
 rng=np.random.default_rng(1729);positions=a.obs.groupby('sample',observed=True).indices
 idx=np.sort(np.concatenate([rng.choice(i,min(35000,len(i)),replace=False) for i in positions.values()]))
 train_object=a[idx].copy()
 pd.Series(train_object.obs_names).to_csv(O/platform/'scvi_balanced_fit_ids.tsv.gz',header=False,index=False)
scvi.model.SCVI.setup_anndata(train_object,batch_key='sample');m=scvi.model.SCVI(train_object,n_latent=12,n_hidden=64,n_layers=1,gene_likelihood='nb')
m.train(max_epochs=24,batch_size=512,accelerator='cpu',devices=1,train_size=.9,early_stopping=True,early_stopping_patience=5,check_val_every_n_epoch=1,enable_progress_bar=False,enable_model_summary=False)
out=O/platform;m.save(str(out/'scvi_model'),overwrite=True,save_anndata=False)
for k,v in m.history.items():v.to_csv(out/(k+'.tsv'),sep='\t')
z=m.get_latent_representation(adata=a,batch_size=512);start=0
for s,i,n in lookup:
 d=paths(s);d.mkdir(exist_ok=True);arr=np.full((n,z.shape[1]),np.nan,dtype='f4');arr[i]=z[start:start+len(i)];np.save(d/'latent.npy',arr);pd.Series(a.obs_names[start:start+len(i)]).to_csv(d/'latent_eligible_ids.txt.gz',index=False,header=False);start+=len(i)
pd.Series(common[j]).to_csv(out/'latent_genes.tsv',index=False,header=False)

