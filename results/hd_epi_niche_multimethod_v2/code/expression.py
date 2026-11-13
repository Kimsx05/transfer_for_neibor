from common import *
import anndata as ad,scanpy as sc,scvi,torch,gc
from scipy import sparse
scvi.settings.logging_dir=str(R/'logs/scvi');scvi.settings.seed=123;torch.set_num_threads(4);torch.set_num_interop_threads(1)

def data(x,ids,gn,s):return ad.AnnData(x.astype(np.float32),obs=pd.DataFrame({'section':s},index=pd.Index(ids.astype(str))),var=pd.DataFrame(index=pd.Index(gn.astype(str))))
def prep():
 out=R/'models/expression';out.mkdir(exist_ok=True)
 if (out/'T_expr.h5ad').exists():return
 rows=[];xs=[];obs=[];genes=None;symbols=None;rng=np.random.default_rng(123)
 for s in SECTIONS[:4]:
  n=pd.read_parquet(V/s/'nodes.parquet');n=n[n.analysis_qc];label=n.type_blca.fillna('Unknown');chosen=[]
  for _,q in n.groupby(label):chosen.extend(rng.choice(q.index,min(1000,len(q)),replace=False))
  rest=np.setdiff1d(n.index,np.asarray(chosen));chosen.extend(rng.choice(rest,50000-len(chosen),replace=False));q=n.loc[chosen].sort_values('bin_id');x,ids,gn,sy=counts(s,ids=q.bin_id.to_numpy())
  if genes is None:genes=gn;symbols=sy
  assert np.array_equal(gn,genes);xs.append(x);obs.append(pd.DataFrame({'section':s},index=s+'::'+ids));rows.append(q[['section','bin_id','type_blca','type_hpca','nCount_RNA','nFeature_RNA']]);del n;gc.collect()
 a=ad.AnnData(sparse.vstack(xs).tocsr().astype(np.float32),obs=pd.concat(obs),var=pd.DataFrame({'symbol':symbols},index=genes));del xs
 pd.concat(rows).to_parquet(R/'tables/T_expr.parquet',index=False)
 # HVG selection only from balanced reference data; temporary log view, raw counts retained.
 raw=a.X.copy();sc.pp.normalize_total(a,target_sum=10000);sc.pp.log1p(a);sc.pp.highly_variable_genes(a,flavor='seurat',n_top_genes=3000,batch_key='section');var=a.var.copy();a.X=raw;a=a[:,var.highly_variable].copy();a.write_h5ad(out/'T_expr.h5ad',compression='gzip');var.to_csv(R/'tables/expression_features.tsv',sep='\t')
 sig=R.parents[2]/'neighborhood_v1/00_inputs_and_audit/mp_completion/frozen_MP_Top20_gene_sets.csv';d=pd.read_csv(sig);d=d[d.MP=='EpiMP10'].copy();assert len(d)==20;d['in_model_features']=d.gene.isin(a.var.symbol);d['source_path']=str(sig);d.to_csv(R/'tables/Top20_feature_exclusion.tsv',sep='\t',index=False)
 del a,raw;gc.collect();print('EXPRESSION_TRAINING_INPUT_COMPLETE',flush=True)

def main():
 gate();prep();out=R/'models/expression';ref=ad.read_h5ad(out/'T_expr.h5ad');top=pd.read_csv(R/'tables/Top20_feature_exclusion.tsv',sep='\t');device='gpu' if torch.cuda.is_available() else 'cpu';timings=[]
 for variant in ['scVI','scVI_without_Top20']:
  dest=out/variant;dest.mkdir(exist_ok=True)
  if (dest/'complete.json').exists():continue
  if variant!='scVI' and not top.in_model_features.any():
   (dest/'complete.json').write_text(json.dumps({'identical_to':'scVI','reason':'no Top20 model features'}));
   for f in (out/'scVI').glob('*_latent.npz'):(dest/f.name).symlink_to(Path('../scVI')/f.name)
   continue
  keep=~ref.var.symbol.isin(top.gene) if variant!='scVI' else np.ones(ref.n_vars,bool);a=ref[:,keep].copy();gn=a.var_names.to_numpy();pd.DataFrame({'feature_id':gn,'symbol':a.var.symbol.to_numpy()}).to_csv(dest/'features.tsv',sep='\t',index=False)
  scvi.settings.seed=123
  if (dest/'reference/model.pt').exists():model=scvi.model.SCVI.load(dest/'reference',adata=a)
  else:
   scvi.model.SCVI.setup_anndata(a,batch_key='section');model=scvi.model.SCVI(a,n_latent=10,n_layers=2,use_layer_norm='both',use_batch_norm='none',encode_covariates=True,dropout_rate=.2)
   start=time.time();model.train(max_epochs=100,accelerator=device,devices=1,batch_size=1024,early_stopping=True,early_stopping_patience=10,train_size=.9,plan_kwargs={"n_epochs_kl_warmup":20},enable_progress_bar=False);model.save(dest/'reference',overwrite=True);timings.append(dict(variant=variant,step='reference',seconds=time.time()-start,epochs=int(model.history['elbo_train'].shape[0])));model.history['elbo_train'].to_csv(dest/'reference_elbo.tsv',sep='\t')
  model.history['elbo_validation'].to_csv(dest/'reference_validation_elbo.tsv',sep='\t')
  assert float(model.history['kl_weight'].iloc[-1,0])==1.0
  original_latent=model.get_latent_representation(a[:1000].copy())
  for s in SECTIONS[:4]:
   f=dest/f'{s}_latent.npz'
   if f.exists():continue
   x,ids,_,_=counts(s,genes=gn);allad=data(x,ids,gn,s);latent=model.get_latent_representation(allad,batch_size=2048);np.savez_compressed(f,bin_id=ids.astype(str),latent=latent);del x,allad,latent;gc.collect();print('ENCODED',variant,s,flush=True)
  s='MPBC2';x,ids,_,_=counts(s,genes=gn);q=data(x,ids,gn,s);rng=np.random.default_rng(123);qi=np.sort(rng.choice(q.n_obs,min(50000,q.n_obs),replace=False));qtrain=q[qi].copy();pd.DataFrame({'section':s,'bin_id':ids[qi]}).to_csv(dest/'query_adaptation_ids.tsv',sep='\t',index=False)
  if (dest/'query/model.pt').exists():qm=scvi.model.SCVI.load(dest/'query',adata=qtrain)
  else:
   qm=scvi.model.SCVI.load_query_data(qtrain,model);start=time.time();qm.train(max_epochs=50,accelerator=device,devices=1,batch_size=1024,plan_kwargs={'weight_decay':0.0,'n_epochs_kl_warmup':10},enable_progress_bar=False);qm.save(dest/'query',overwrite=True);timings.append(dict(variant=variant,step='MPBC2_query',seconds=time.time()-start,epochs=50));qm.history['elbo_train'].to_csv(dest/'query_elbo.tsv',sep='\t')
  # Original reference model remains separate and exactly unchanged; adaptation mapping
  # of reference points is also checked, not silently assumed invariant.
  assert float(qm.history['kl_weight'].iloc[-1,0])==1.0
  after=model.get_latent_representation(a[:1000].copy());query_ref=qm.get_latent_representation(a[:1000].copy());assert np.allclose(original_latent,after)
  np.savez_compressed(dest/f'{s}_latent.npz',bin_id=ids.astype(str),latent=qm.get_latent_representation(q,batch_size=2048))
  (dest/'complete.json').write_text(json.dumps(dict(variant=variant,device=device,reference_train=a.n_obs,features=a.n_vars,reference_kl_warmup_epochs=20,query_kl_warmup_epochs=10,reference_final_KL_weight=float(model.history['kl_weight'].iloc[-1,0]),query_final_KL_weight=float(qm.history['kl_weight'].iloc[-1,0]),query_adaptation_n=len(qi),query_encoded_n=len(q),reference_model_unchanged_max_abs=float(np.max(abs(after-original_latent))),adapted_model_reference_max_abs=float(np.max(abs(query_ref-original_latent))),scope='four-section reference T_expr; 50000 MPBC2 query adaptation; all QC neighbor bins encoded'),indent=2));del q,qtrain,x,a,model,qm;gc.collect();torch.cuda.empty_cache();print('SCVI_COMPLETE',variant,flush=True)
 pd.DataFrame(timings).to_csv(R/'tables/expression_runtime.tsv',sep='\t',index=False)
if __name__=='__main__':main()
