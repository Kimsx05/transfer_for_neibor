from common import *
import torch,scvi,anndata as ad,time
from cell2location.models import Cell2location
from lightning.pytorch.callbacks import Callback
O=R/'05_deconvolution/ARCHIVE_MAIN_DEV_3000' if (R/'05_deconvolution/ARCHIVE_MAIN_DEV_3000').exists() else R/'05_deconvolution/MAIN_DEV';D=R/'05_deconvolution/DEV_convergence_probe_5000';D.mkdir(exist_ok=True)
assert len(pd.read_csv(O/'training_history.tsv',sep='\t'))==3000,'Probe requires the actual3000epoch checkpoint'
if (D/'PROBE_DONE.json').exists():
 assert json.loads((D/'PROBE_CONFIG.json').read_text())['baseline_model_sha256']==sha(O/'model/model.pt');say('VALIDATED_EXISTING_PROBE');raise SystemExit
torch.set_num_threads(1);torch.set_num_interop_threads(1);scvi.settings.seed=seed('mapping_DEV');a=ad.read_h5ad(O/'mapping_registered.h5ad');m=Cell2location.load(O/'model',adata=a,accelerator='gpu',device=1)
js(D/'PROBE_CONFIG.json',{'purpose':'DEV-only abundance stability check; loss plateau alone may not establish abundance convergence','baseline_epochs':3000,'additional_epochs':2000,'counts_sha256':sha(R/'04_pseudospots/DEV/counts.h5ad'),'baseline_model_sha256':sha(O/'model/model.pt'),'TEST_predictions_seen':False,'fit_seed':seed('mapping_DEV'),'posterior_seed':seed('mapping_posterior_DEV'),'same_priors':True})
class Progress(Callback):
 def on_train_epoch_end(self,trainer,pl_module):
  e=3000+trainer.current_epoch+1
  if e%250==0:say('DEV_PROBE_EPOCH',e)
say('DEV_PROBE_START')
m.train(max_epochs=2000,batch_size=None,train_size=1,lr=.002,accelerator='gpu',device=1,enable_progress_bar=False,enable_checkpointing=False,callbacks=[Progress()]);m.save(D/'model',overwrite=True,save_anndata=False);h=np.asarray(m.history_['elbo_train']).ravel().astype(float)[-2000:];wt(pd.DataFrame({'epoch':np.arange(3001,5001),'loss':h}),D/'training_history.tsv')
scvi.settings.seed=seed('mapping_posterior_DEV');means=m.sample_posterior(num_samples=256,batch_size=1000,accelerator='gpu',device=1,return_sites=['w_sf'],return_samples=False,summary_frequency=1)['post_sample_means']['w_sf'];labels=list(m.factor_names_);p=pd.DataFrame(means,index=a.obs_names,columns=labels);p.index.name='pseudospot_id';p.to_csv(D/'abundance_mean.tsv.gz',sep='\t');old=pd.read_csv(O/'abundance_mean.tsv.gz',sep='\t',index_col=0);delta=p-old;rel=abs(delta)/(old+.1);wt(pd.DataFrame({'subtype':labels,'mean_abs_change':abs(delta).mean().values,'median_relative_change':rel.median().values,'q95_relative_change':rel.quantile(.95).values,'mean_3000':old.mean().values,'mean_5000':p.mean().values}),D/'abundance_stability_by_subtype.tsv');js(D/'PROBE_DONE.json',{'epochs':5000,'posterior_samples':256,'loss_window_relative_change':float(abs(h[-500:].mean()-h[-1000:-500].mean())/abs(h[-500:].mean())),'K32_mean_absolute_change':float(abs(delta.Epi_K32).mean()),'K32_median_relative_change':float(rel.Epi_K32.median()),'K32_q95_relative_change':float(rel.Epi_K32.quantile(.95))});say('DEV_PROBE_DONE')
