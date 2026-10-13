from common import *
import anndata as ad,torch,scvi
from cell2location.models import Cell2location
from scvi.module.base import PyroBaseModuleClass

def main(out):
 P=out/'03_attribution_calibration/synthetic_refits/model1'
 if (P/'posterior.h5ad').exists() and (P/'fit_status.json').exists():log('posterior already complete');return
 assert (P/'model/model.pt').exists(),'Saved synthetic checkpoint required; no fresh fit permitted'
 torch.set_num_threads(4);torch.set_num_interop_threads(1);scvi.settings.seed=123
 a=ad.read_h5ad(out/'03_attribution_calibration/synthetic_refits/input.h5ad')
 original=PyroBaseModuleClass.on_load
 def cpu_load(self,model,**kwargs):
  self.on_load_kwargs.update(accelerator='cpu',device='auto',enable_progress_bar=False,default_root_dir=str(P))
  return original(self,model,**kwargs)
 PyroBaseModuleClass.on_load=cpu_load
 m=Cell2location.load(str(P/'model'),adata=a,accelerator='cpu',device='auto')
 h=np.array([float(x) for x in np.asarray(m.history_['elbo_train']).ravel()]);assert np.isfinite(h).all();change=abs(np.median(h[-500:])-np.median(h[-1000:-500]))/max(abs(np.median(h[-1000:-500])),1);converged=change<.01
 write_tsv(out,P/'convergence_history.tsv',[dict(total_optimizer_step_cap=5000,n_history_epochs=len(h),relative_last500_change=change,converged=converged,history_conversion='float(scalar tensor), no new fresh fit')])
 assert converged,'Checkpoint not converged; do not silently declare success'
 np.save(safe(out,P/'elbo.npy'),h)
 a=m.export_posterior(a,sample_kwargs={'num_samples':100,'batch_size':128,'accelerator':'cpu','device':'auto'})
 # New anndata forbids slash-bearing DataFrame column keys in obsm; all full posterior arrays remain in uns.mod.
 for key in list(a.obsm.keys()):del a.obsm[key]
 a.write_h5ad(safe(out,P/'posterior.partial.h5ad'));(P/'posterior.partial.h5ad').replace(P/'posterior.h5ad')
 write_json(out,P/'fit_status.json',dict(status='CONVERGED_OBJECTIVE_RULE',total_optimizer_steps=5000,n_history_epochs=len(h),n_bags=a.n_obs,n_genes=a.n_vars,relative_objective_change=change,device='cpu',threads=4,posterior_samples=100,fresh_fit_attempts=2,usable_models=1,truth_not_model_input=True,loader_initialization='scvi one-step parameter shape initialization before restoring checkpoint; original checkpoint unchanged'))
 at=json.loads((P/'attempts.json').read_text());at[-1]['status']='COMPLETED_FROM_SAVED_5000_STEP_CHECKPOINT';write_json(out,P/'attempts.json',at);log('Saved synthetic checkpoint posterior exported')
if __name__=='__main__':main(OUT)
