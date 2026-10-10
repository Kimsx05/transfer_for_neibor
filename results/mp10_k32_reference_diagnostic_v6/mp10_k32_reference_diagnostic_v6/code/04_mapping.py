from common import *
import sys,time
import anndata as ad,torch,scvi
from cell2location.models import Cell2location
from lightning.pytorch.callbacks import Callback
from scipy.stats import pearsonr
labelset,collection=sys.argv[1:3];assert labelset in ['S30','S31'];assert collection in ['FULL_DEV','TARGET_UMI','FIXED_CAPTURE'];assert not (labelset=='S31' and collection=='FULL_DEV'),'Reuse existing final S31 full DEV; no refit'
O=R/'05_deconvolution'/(labelset+'_'+collection);O.mkdir(exist_ok=True);cfg=json.loads((R/'CONFIG.yaml').read_text());c=cfg['mapping'];labels=json.loads((R/f'00_manifest/labels_{labelset}.json').read_text());I=S/'04_pseudospots/DEV' if collection=='FULL_DEV' else R/'04_paired_depth_control'/collection;ref=S/'03_reference/MAIN' if labelset=='S31' else R/'03_merged_reference';st=json.loads((ref/'REFERENCE_DONE.json').read_text())
if st['converged'] is not True:js(O/'STATUS.json',{'status':'UNAVAILABLE','reason':'required reference NOT_CONVERGED'});raise SystemExit
if collection!='FULL_DEV':assert (R/'04_paired_depth_control/READY.json').exists()
fitseed=seed('mapping_DEV' if collection=='FULL_DEV' else 'v6_control_mapping');postseed=seed('mapping_posterior_DEV' if collection=='FULL_DEV' else 'v6_control_posterior');torch.set_num_threads(1);torch.set_num_interop_threads(1);scvi.settings.seed=fitseed
a=ad.read_h5ad(I/'counts.h5ad');assert list(a.obs.columns)==['sample'];a.obs['sample']=pd.Categorical(a.obs['sample']);sig=pd.read_csv(ref/'signatures.tsv.gz',sep='\t',index_col=0);assert list(sig.index)==list(a.var_names);sig=sig[labels];Cell2location.setup_anndata(a,batch_key='sample')
finger={'labelset':labelset,'collection':collection,'counts_path':str(I/'counts.h5ad'),'counts_sha256':sha(I/'counts.h5ad'),'signature_sha256':sha(ref/'signatures.tsv.gz'),'gene_order_hash':key(a.var_names),'spot_order_hash':key(a.obs_names),'labels':labels,'mapping_config':c,'fit_seed':fitseed,'posterior_seed':postseed,'control_hash':sha(R/'CONTROL_FROZEN.yaml') if collection!='FULL_DEV' else None}
if (O/'run_fingerprint.json').exists():assert json.loads((O/'run_fingerprint.json').read_text())==finger
else:js(O/'run_fingerprint.json',finger)
if (O/'MAPPING_DONE.json').exists():done=json.loads((O/'MAPPING_DONE.json').read_text());assert sha(O/'w_sf_draws.npy')==done['draws_sha256'];say('VALIDATED_DONE',labelset,collection);raise SystemExit
cp=O/'model';history=[];epochs=0
if (cp/'model.pt').exists():m=Cell2location.load(cp,adata=a,accelerator='gpu',device=1);history=pd.read_csv(O/'training_history.tsv',sep='\t').loss.tolist();epochs=len(history)
else:m=Cell2location(a,cell_state_df=sig,N_cells_per_location=c['N_cells_per_location'],detection_alpha=c['detection_alpha'])
assert list(m.factor_names_)==labels
class Progress(Callback):
 def on_train_epoch_end(self,trainer,pl_module):
  e=epochs+trainer.current_epoch+1
  if e%500==0:js(O/'progress.json',{'stage':'TRAINING','epochs':e,'time':datetime.datetime.now().isoformat()});say('EPOCH',labelset,collection,e)
bpath=O/'abundance_checkpoint_mean.tsv.gz';baseline=pd.read_csv(bpath,sep='\t',index_col=0) if bpath.exists() else None;sh=O/'abundance_convergence_history.tsv';stability=pd.read_csv(sh,sep='\t').to_dict('records') if sh.exists() else [];converged=False
if epochs>0 and (O/'convergence.json').exists():converged=json.loads((O/'convergence.json').read_text())['converged']
while epochs<c['max_epochs'] and not converged:
 n=min(c['max_epochs']-epochs,c['initial_epochs'] if epochs==0 else c['extension_epochs']);scvi.settings.seed=fitseed;say('MAP_START',labelset,collection,n,a.shape);m.train(max_epochs=n,batch_size=None,train_size=1,lr=c['lr'],accelerator='gpu',device=1,enable_progress_bar=False,enable_checkpointing=False,callbacks=[Progress()]);h=np.asarray(m.history_['elbo_train']).ravel().astype(float)[-n:];assert np.isfinite(h).all();history.extend(h.tolist());epochs=len(history);m.save(cp,overwrite=True,save_anndata=False);wt(pd.DataFrame({'epoch':np.arange(1,epochs+1),'loss':history}),O/'training_history.tsv')
 w=min(500,epochs//4);losschange=abs(np.mean(history[-w:])-np.mean(history[-2*w:-w]))/abs(np.mean(history[-w:]));scvi.settings.seed=postseed;ps=m.sample_posterior(num_samples=c['abundance_stability_samples'],batch_size=min(1000,a.n_obs),accelerator='gpu',device=1,return_sites=['w_sf'],return_samples=False,summary_frequency=1)['post_sample_means']['w_sf'];current=pd.DataFrame(ps,index=a.obs_names,columns=labels);current.index.name='pseudospot_id';stable=False
 if baseline is not None:
  passed=[]
  for lab in (['Epi_K32','Epi_other','total_epithelial'] if labelset=='S31' else ['Epithelial']):
   old=baseline.Epi_K32+baseline.Epi_other if lab=='total_epithelial' else baseline[lab];now=current.Epi_K32+current.Epi_other if lab=='total_epithelial' else current[lab];delta=abs(now-old)/(old+.1);med=float(delta.median());q95=float(delta.quantile(.95));ok=med<c['abundance_median_relative_change_max'] and q95<c['abundance_q95_relative_change_max'];passed.append(ok);stability.append({'epoch':epochs,'abundance':lab,'median_relative_change':med,'q95_relative_change':q95,'stable':ok})
  stable=all(passed);wt(pd.DataFrame(stability),sh)
 current.to_csv(bpath,sep='\t');baseline=current;converged=bool(epochs>=c['minimum_epochs'] and stable and losschange<c['convergence_relative_window_change']);js(O/'convergence.json',{'epochs':epochs,'relative_window_change':losschange,'abundance_stable':stable,'converged':converged,'S30_stability_target':'Epithelial replaces S31 K32/other/sum; same numerical rule'});say('CONVERGENCE',labelset,collection,epochs,converged)
scvi.settings.seed=postseed;draw=m.sample_posterior(num_samples=cfg['posterior_samples'],batch_size=min(1000,a.n_obs),accelerator='gpu',device=1,return_sites=['w_sf'],return_samples=True,summary_frequency=1)['posterior_samples']['w_sf'];assert draw.shape==(1000,a.n_obs,len(labels)) and np.isfinite(draw).all();np.save(O/'w_sf_draws.npy',draw.astype(np.float32))
for name,fun in [('mean',np.mean),('median',np.median),('q05',lambda x,axis:np.quantile(x,.05,axis=axis)),('q95',lambda x,axis:np.quantile(x,.95,axis=axis))]:
 d=pd.DataFrame(fun(draw,axis=0),index=a.obs_names,columns=labels);d.index.name='pseudospot_id';d.to_csv(O/f'abundance_{name}.tsv.gz',sep='\t')
epi=draw[:,:,labels.index('Epithelial')] if labelset=='S30' else draw[:,:,labels.index('Epi_K32')]+draw[:,:,labels.index('Epi_other')];total=draw.sum(2);ratio=pd.DataFrame({'pseudospot_id':a.obs_names,'epithelial_denominator_unstable':(epi<.1).mean(0)>.05,'fraction_draws_epi_below_min':(epi<.1).mean(0)});arrays={'epithelial':epi,'total_cells':total,'epithelial_fraction_all':epi/total}
if labelset=='S31':
 k=draw[:,:,labels.index('Epi_K32')];rk=np.divide(k,epi,out=np.full_like(k,np.nan),where=epi>=.1);rk[:,ratio.epithelial_denominator_unstable]=np.nan;arrays.update(K32_fraction_epithelial=rk,K32_fraction_all=k/total)
for lab,x in arrays.items():
 for name,fn in [('mean',np.nanmean),('median',np.nanmedian),('q05',lambda z,axis:np.nanquantile(z,.05,axis=axis)),('q95',lambda z,axis:np.nanquantile(z,.95,axis=axis))]:ratio[lab+'_'+name]=fn(x,axis=0)
wt(ratio,O/'posterior_derived.tsv.gz')
# Low-dimensional technical terms and posterior means only; no spot x gene x draw export.
scvi.settings.seed=seed('v6_technical_posterior');tech=m.sample_posterior(num_samples=100,batch_size=min(1000,a.n_obs),accelerator='gpu',device=1,return_sites=['m_g','s_g_gene_add','detection_y_s','detection_mean_y_e','alpha_g_inverse'],return_samples=False,summary_frequency=1);rp=tech['post_sample_means'];np.savez_compressed(O/'technical_parameter_means.npz',**rp);rp['w_sf']=draw.mean(0);batch=a.obs['_scvi_batch'].to_numpy();det=rp['detection_y_s'].reshape(-1);bg=rp['s_g_gene_add'].sum(1)[batch]*det;celltotal=(rp['w_sf']@(sig.to_numpy()*rp['m_g'].reshape(-1,1)).sum(0))*det;td=pd.DataFrame({'pseudospot_id':a.obs_names,'detection_y_s_mean':det,'detection_batch_mean':rp['detection_mean_y_e'].reshape(-1)[batch],'background_expected_UMI_from_means':bg,'cellular_expected_UMI_from_means':celltotal,'expected_UMI_from_means':bg+celltotal,'observed_G0_UMI':np.asarray(a.X.sum(1)).ravel()});wt(td,O/'technical_terms.tsv.gz')
ix=np.linspace(0,a.n_obs-1,min(256,a.n_obs),dtype=int);mu=m.module.model.compute_expected(rp,m.adata_manager,ind_x=ix)['mu'];obs=a.X[ix].toarray();assert np.allclose(mu.sum(1),td.expected_UMI_from_means.iloc[ix],rtol=2e-5);wt(pd.DataFrame({'gene':a.var_names,'observed_mean':obs.mean(0),'expected_mean':mu.mean(0)}),O/'reconstruction_genes.tsv');js(O/'reconstruction_QC.json',{'log1p_gene_mean_pearson':float(pearsonr(np.log1p(obs.mean(0)),np.log1p(mu.mean(0))).statistic),'expected_observed_UMI_ratio':float(mu.sum()/obs.sum()),'all_finite':bool(np.isfinite(mu).all())});wt(a.obs.reset_index(),O/'registered_obs.tsv.gz');m.save(cp,overwrite=True,save_anndata=False);js(O/'MAPPING_DONE.json',{'status':'COMPLETED' if converged else 'NOT_CONVERGED','converged':converged,'epochs':epochs,'n_spots':a.n_obs,'n_labels':len(labels),'counts_sha256':finger['counts_sha256'],'draws_sha256':sha(O/'w_sf_draws.npy'),'input_path':str(I/'counts.h5ad'),'K32_prediction_available':labelset=='S31'});say('MAPPING_DONE',labelset,collection,epochs,converged)
