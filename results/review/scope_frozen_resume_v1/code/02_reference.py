from common import *
import sys,time,gc,inspect
import anndata as ad,torch,scvi
from cell2location.models import RegressionModel
from lightning.pytorch.callbacks import Callback
from scipy.stats import pearsonr,spearmanr
branch=sys.argv[1];O=R/'03_reference'/branch;cfg=json.loads((R/'CONFIG.yaml').read_text());freeze=json.loads((R/'02_inputs/INPUTS_READY.json').read_text());labels=json.loads((R/'02_inputs/labels.json').read_text())
torch.set_num_threads(1);torch.set_num_interop_threads(1);scvi.settings.seed=seed('reference_fit');np.random.seed(seed('reference_fit'));torch.manual_seed(seed('reference_fit'))
a=ad.read_h5ad(R/'02_inputs/TRAIN_counts.h5ad');a=a[a.obs['sample']!=P07].copy() if branch=='NO_P07T' else a
assert a.obs.new_reference_label.nunique()==31;a.obs['new_reference_label']=pd.Categorical(a.obs.new_reference_label,categories=labels);a.obs['sample']=pd.Categorical(a.obs['sample']);RegressionModel.setup_anndata(a,labels_key='new_reference_label',batch_key='sample')
fingerprint={**freeze,'branch':branch,'member_order':key(a.obs_names),'reference_settings':cfg['NB'],'seed':seed('reference_fit')};cp=O/'model';history=[];epochs=0
if (O/'run_fingerprint.json').exists():assert json.loads((O/'run_fingerprint.json').read_text())==fingerprint
else:js(O/'run_fingerprint.json',fingerprint)
if (O/'REFERENCE_DONE.json').exists():st=json.loads((O/'REFERENCE_DONE.json').read_text());assert sha(O/'signatures.tsv.gz')==st['signature_sha256'] and (cp/'model.pt').exists();say('VALIDATED_DONE',branch);raise SystemExit
if (cp/'model.pt').exists():
 m=RegressionModel.load(cp,adata=a,accelerator='gpu',device=1);history=pd.read_csv(O/'training_history.tsv',sep='\t').loss.tolist();epochs=len(history)
else:m=RegressionModel(a)
assert list(m.factor_names_)==labels
class Progress(Callback):
 def on_train_epoch_end(self,trainer,pl_module):
  e=epochs+trainer.current_epoch+1
  if e%25==0:js(O/'progress.json',{'stage':'TRAINING','epochs':e,'time':datetime.datetime.now().isoformat()});say(branch,'epoch',e)
converged=False;t0=time.time()
while epochs<cfg['NB']['max_epochs']:
 n=cfg['NB']['initial_epochs'] if epochs==0 else cfg['NB']['extension_epochs'];say('TRAIN',branch,'shape',a.shape,'epochs',n)
 m.train(max_epochs=n,batch_size=2500,train_size=1,lr=.002,accelerator='gpu',device=1,enable_progress_bar=False,enable_checkpointing=False,callbacks=[Progress()])
 h=np.asarray(m.history_['elbo_train']).ravel().astype(float);h=h[-n:];assert np.isfinite(h).all();history.extend(h.tolist());epochs=len(history);wt(pd.DataFrame({'epoch':np.arange(1,epochs+1),'loss':history}),O/'training_history.tsv');m.save(cp,overwrite=True,save_anndata=False)
 w=min(50,len(history)//4);rel=abs(np.mean(history[-w:])-np.mean(history[-2*w:-w]))/abs(np.mean(history[-w:]));converged=epochs>=250 and rel<cfg['NB']['convergence_relative_window_change'];js(O/'convergence.json',{'epochs':epochs,'window':w,'relative_window_change':rel,'plateau_pass':converged,'elapsed_seconds_this_process':time.time()-t0});say('CONVERGENCE',branch,epochs,rel)
 if converged:break
# Sample only parameters needed for signatures and reconstruction; no cell-level draws retained.
scvi.settings.seed=seed('reference_posterior');samples=m.sample_posterior(num_samples=1000,batch_size=2500,accelerator='gpu',device=1,return_sites=['per_cluster_mu_fg','detection_mean_y_e','s_g_gene_add','alpha_g_inverse'],return_samples=False)
means=samples['post_sample_means'];np.savez_compressed(O/'posterior_parameter_means.npz',**means)
batch_codes=a.obs['_scvi_batch'].to_numpy();scale=float(np.mean(means['detection_mean_y_e'][batch_codes]));sig=pd.DataFrame(means['per_cluster_mu_fg'].T*scale,index=a.var_names,columns=labels);sig.index.name='gene';assert np.isfinite(sig.to_numpy()).all() and (sig.to_numpy()>=0).all();sig.to_csv(O/'signatures.tsv.gz',sep='\t')
# Numerical reconstruction QC on 1000 fixed TRAIN cells; no held-out truth used.
ix=np.sort(np.random.default_rng(seed('reference_reconstruction_QC')).choice(a.n_obs,min(1000,a.n_obs),replace=False));mu=m.module.model.compute_expected(means,m.adata_manager,ind_x=ix)['mu'];obs=a.X[ix].toarray();rec=pd.DataFrame({'cell_id':a.obs_names[ix],'observed_UMI':obs.sum(1),'expected_UMI':mu.sum(1),'mean_absolute_gene_error':np.mean(abs(obs-mu),axis=1)});wt(rec,O/'reconstruction_cells.tsv')
gene_rec=pd.DataFrame({'gene':a.var_names,'observed_mean':obs.mean(0),'expected_mean':mu.mean(0)});wt(gene_rec,O/'reconstruction_genes.tsv')
qc={'cell_total_pearson':float(pearsonr(rec.observed_UMI,rec.expected_UMI).statistic),'log1p_gene_mean_pearson':float(pearsonr(np.log1p(obs.mean(0)),np.log1p(mu.mean(0))).statistic),'observed_total':float(obs.sum()),'expected_total':float(mu.sum()),'all_finite':bool(np.isfinite(mu).all()),'n_QC_cells':len(ix),'signature_detection_scale':scale}
js(O/'reconstruction_QC.json',qc);a.write_h5ad(O/'training_registered.h5ad',compression='lzf');m.save(cp,overwrite=True,save_anndata=False)
js(O/'REFERENCE_DONE.json',{'training_completed':True,'converged':converged,'epochs':epochs,'n_cells':a.n_obs,'n_genes':a.n_vars,'n_labels':len(labels),'signature_sha256':sha(O/'signatures.tsv.gz'),'QC':qc,'biological_identifiability':'not_determined_by_reference_training'})
say('REFERENCE_EXPORTED',branch,converged)
