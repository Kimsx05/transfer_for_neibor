from common import *
import anndata as ad,torch,scvi
from cell2location.models import Cell2location
from lightning.pytorch.callbacks import Callback
class Progress(Callback):
 def on_train_batch_end(self,trainer,pl_module,outputs,batch,batch_idx):
  if trainer.global_step and trainer.global_step%500==0:log(f"synthetic optimizer step {trainer.global_step}; epoch {trainer.current_epoch}")

def main(out):
 torch.set_num_threads(4);torch.set_num_interop_threads(1);scvi.settings.seed=123;P=out/'03_attribution_calibration';dest=P/'synthetic_refits/model1';dest.mkdir(exist_ok=True);a=ad.read_h5ad(P/'synthetic_refits/input.h5ad');axis=json.loads((P/'panel_axes.json').read_text());ref=pd.read_csv(P/'frozen_reference.tsv',sep='\t',index_col=0);assert list(ref.index)==list(a.var_names) and list(ref.columns)==axis['types'];assert 'state' not in a.obs and 'true_score' not in a.obs
 if (dest/'posterior.h5ad').exists():log('synthetic posterior exists; no extra model fit');return
 attempt_file=dest/'attempts.json'
 attempts=json.loads(attempt_file.read_text()) if attempt_file.exists() else [{'attempt':1,'status':'FAILED_AFTER_5000_STEPS_LOGGER_HISTORY_MISSING_NO_SAVED_MODEL'}]
 assert len(attempts)<2,'Two fresh fitting attempts already consumed; use saved checkpoint, no third fit'
 attempts.append({'attempt':2,'status':'RUNNING_FIXED_LOGGER_AND_CHECKPOINT_ORDER'});write_json(out,attempt_file,attempts)
 Cell2location.setup_anndata(a,batch_key='technical_batch');m=Cell2location(a,cell_state_df=ref,N_cells_per_location=20,detection_alpha=20);history=[];total=0;converged=False
 for steps in [5000,5000,10000]:
  m.train(max_epochs=steps,batch_size=128,train_size=1,accelerator='cpu',device='auto',max_steps=steps,enable_checkpointing=False,enable_progress_bar=False,default_root_dir=str(dest),callbacks=[Progress()]);total+=steps;m.save(str(safe(out,dest/'model')),overwrite=True)
  h=np.array([float(x) for x in np.asarray(m.history_['elbo_train']).ravel()]);assert np.isfinite(h).all();change=abs(np.median(h[-500:])-np.median(h[-1000:-500]))/max(abs(np.median(h[-1000:-500])),1) if len(h)>=1000 else np.inf;converged=change<.01;history.append(dict(total_optimizer_step_cap=total,n_history_epochs=len(h),relative_last500_change=change,converged=converged));write_tsv(out,dest/'convergence_history.tsv',history);m.save(str(safe(out,dest/'model')),overwrite=True);np.save(safe(out,dest/'elbo.npy'),h);log(f'synthetic checkpoint {total}, objective change {change}, converged {converged}')
  if converged:break
 a=m.export_posterior(a,sample_kwargs={'num_samples':100,'batch_size':128,'accelerator':'cpu','device':'auto'});a.write_h5ad(safe(out,dest/'posterior.h5ad'));write_json(out,dest/'fit_status.json',dict(status='CONVERGED_OBJECTIVE_RULE' if converged else 'TECHNICALLY_UNRESOLVED',total_optimizer_step_cap=total,n_bags=a.n_obs,n_genes=a.n_vars,device='cpu',threads=4,posterior_samples=100,new_models=1,truth_not_model_input=True));log('synthetic refit posterior saved')
if __name__=='__main__':main(OUT)
