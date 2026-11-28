from data import *
import torch,scvi
torch.set_num_threads(2);torch.set_num_interop_threads(1);scvi.settings.seed=1729;scvi.settings.dl_num_workers=0
p=sys.argv[1];d=O/p;x,g,q=combine(p,save=True)
y=x.astype('f4').multiply((1e4/q.n_counts.to_numpy())[:,None]).tocsr().astype('f4');y.data=np.log1p(y.data);mu=y.mean(0).A1;disp=(y.power(2).mean(0).A1-mu**2)/(.1+mu);j=np.sort(np.argsort(disp)[-1200:]);del y
gg=g.assign(dispersion=disp,selected=np.isin(np.arange(len(g)),j));write(gg,d/'model_features.tsv')
a=ad.AnnData(x[:,j],obs=pd.DataFrame({'sample':pd.Categorical(q['sample'].values)},index=q.cell_key),var=g.iloc[j].set_index('feature_id'));del x
rng=np.random.default_rng(1729);tr=np.sort(np.concatenate([rng.choice(v,min(len(v),35000),False) for v in a.obs.groupby('sample',observed=True).indices.values()])) if p=='hd' else np.arange(len(q))
b=a[tr].copy() if p=='hd' else a;scvi.model.SCVI.setup_anndata(b,batch_key='sample')
model=scvi.model.SCVI(b,n_latent=12,n_hidden=64,n_layers=1,gene_likelihood='nb',dispersion='gene')
model.train(max_epochs=24,batch_size=512,accelerator='cpu',devices=1,train_size=.9,early_stopping=True,early_stopping_patience=5,check_val_every_n_epoch=1,enable_progress_bar=False,enable_model_summary=False)
model.save(str(d/'scvi_model'),overwrite=False,save_anndata=False)
for k,df in model.history.items():df.to_csv(d/(k+'.tsv'),sep='\t')
splits=pd.DataFrame({'cell_key':b.obs_names,'role':'fit_frame'});splits.iloc[model.train_indices,1]='training';splits.iloc[model.validation_indices,1]='validation';write(splits,d/'fit_members.tsv.gz')
z=model.get_latent_representation(adata=a,batch_size=512);assert np.isfinite(z).all();pd.DataFrame(z,index=q.cell_key,columns=['z'+str(i) for i in range(12)]).to_parquet(d/'latent.parquet')

