from common import *
import sys,time
import anndata as ad,torch,scvi
from cell2location.models import Cell2location
from lightning.pytorch.callbacks import Callback
from scipy.stats import pearsonr
branch,part=sys.argv[1:3];smoke='--smoke' in sys.argv
if branch=='NO_P07T' and part=='DEV':
 while (R/'05_deconvolution/DEV_convergence_probe_5000/PROBE_CONFIG.json').exists() and not (R/'06_evaluation/DEV_CONVERGENCE_REVIEW.json').exists():time.sleep(10)
O=R/'05_deconvolution'/('SMOKE_DEV' if smoke else branch+'_'+part);O.mkdir(exist_ok=True);cfg=json.loads((R/'CONFIG.yaml').read_text());labels=json.loads((R/'02_inputs/labels.json').read_text());refstat=json.loads((R/f'03_reference/{branch}/REFERENCE_DONE.json').read_text());assert refstat['converged'] in (True,'True'),'Reference not converged; mapping cannot support valid evaluation'
if not smoke and (R/'06_evaluation/DEV_CONVERGENCE_REVIEW.json').exists():cfg['mapping'].update(json.loads((R/'06_evaluation/DEV_CONVERGENCE_REVIEW.json').read_text()).get('mapping_overrides',{}))
if not smoke:
 for k,v in {'minimum_epochs':5000,'abundance_stability_samples':256,'abundance_median_relative_change_max':.05,'abundance_q95_relative_change_max':.20}.items():cfg['mapping'].setdefault(k,v)
if part=='TEST':assert (R/'06_evaluation/DEV_FROZEN.json').exists(),'Freeze DEV before TEST mapping'
torch.set_num_threads(1);torch.set_num_interop_threads(1);fitseed=seed('mapping_'+part);scvi.settings.seed=fitseed;np.random.seed(fitseed);torch.manual_seed(fitseed)
a=ad.read_h5ad(R/f'04_pseudospots/{part}/counts.h5ad');sig=pd.read_csv(R/f'03_reference/{branch}/signatures.tsv.gz',sep='\t',index_col=0).loc[a.var_names,labels]
if smoke:
 ix=np.sort(np.random.default_rng(seed('smoke_DEV')).choice(a.n_obs,64,replace=False));a=a[ix].copy()
assert set(a.obs.columns)=={'sample'};a.obs['sample']=pd.Categorical(a.obs['sample']);Cell2location.setup_anndata(a,batch_key='sample')
finger={'reference':branch,'split':part,'smoke':smoke,'counts_sha256':sha(R/f'04_pseudospots/{part}/counts.h5ad'),'signatures_sha256':sha(R/f'03_reference/{branch}/signatures.tsv.gz'),'gene_order':key(a.var_names),'labels':labels,'config':cfg['mapping'],'fit_seed':fitseed,'spot_order':key(a.obs_names)}
if (O/'run_fingerprint.json').exists():assert json.loads((O/'run_fingerprint.json').read_text())==finger
else:js(O/'run_fingerprint.json',finger)
if (O/'MAPPING_DONE.json').exists():st=json.loads((O/'MAPPING_DONE.json').read_text());assert sha(O/'w_sf_draws.npy')==st['w_sf_draws_sha256'] and (O/'model/model.pt').exists();say('VALIDATED_DONE',branch,part);raise SystemExit
cp=O/'model';history=[];epochs=0
if (cp/'model.pt').exists():
 m=Cell2location.load(cp,adata=a,accelerator='gpu',device=1);history=pd.read_csv(O/'training_history.tsv',sep='\t').loss.tolist();epochs=len(history)
else:m=Cell2location(a,cell_state_df=sig,N_cells_per_location=20,detection_alpha=20)
assert list(m.factor_names_)==labels
class Progress(Callback):
 def on_train_epoch_end(self,trainer,pl_module):
  e=epochs+trainer.current_epoch+1
  if e%250==0:js(O/'progress.json',{'stage':'TRAINING','epochs':e,'time':datetime.datetime.now().isoformat()});say('MAPPING_EPOCH',branch,part,e)
converged=False;t0=time.time();maxep=50 if smoke else cfg['mapping']['max_epochs']
baseline_path=O/'abundance_checkpoint_mean.tsv.gz';baseline=pd.read_csv(baseline_path,sep='\t',index_col=0) if baseline_path.exists() else None
stab_path=O/'abundance_convergence_history.tsv';stability_history=pd.read_csv(stab_path,sep='\t').to_dict('records') if stab_path.exists() else []
if '--export-only' in sys.argv:
 assert epochs>=cfg['mapping'].get('minimum_epochs',3000) and not smoke
 st=json.loads((O/'convergence.json').read_text());converged=bool(st['plateau_pass'] and st.get('abundance_stable',False));assert converged;maxep=epochs
while epochs<maxep:
 n=50 if smoke else min(maxep-epochs,cfg['mapping']['initial_epochs'] if epochs==0 else cfg['mapping']['extension_epochs']);say('MAP_START',branch,part,a.shape,n)
 if epochs>0:scvi.settings.seed=fitseed
 m.train(max_epochs=n,batch_size=None,train_size=1,lr=.002,accelerator='gpu',device=1,enable_progress_bar=False,enable_checkpointing=False,callbacks=[Progress()]);h=np.asarray(m.history_['elbo_train']).ravel().astype(float)[-n:];assert np.isfinite(h).all();history.extend(h.tolist());epochs=len(history);wt(pd.DataFrame({'epoch':np.arange(1,epochs+1),'loss':history}),O/'training_history.tsv');m.save(cp,overwrite=True,save_anndata=False)
 w=min(500,len(history)//4);rel=abs(np.mean(history[-w:])-np.mean(history[-2*w:-w]))/abs(np.mean(history[-w:]));loss_stable=bool(rel<cfg['mapping']['convergence_relative_window_change']);abundance_stable=False
 if not smoke:
  scvi.settings.seed=seed('mapping_posterior_'+part);st=m.sample_posterior(num_samples=cfg['mapping']['abundance_stability_samples'],batch_size=1000,accelerator='gpu',device=1,return_sites=['w_sf'],return_samples=False,summary_frequency=1)['post_sample_means']['w_sf'];current=pd.DataFrame(st,index=a.obs_names,columns=labels);current.index.name='pseudospot_id'
  if baseline is not None:
   tests=[]
   for label in ['Epi_K32','Epi_other','total_epithelial']:
    old=baseline.Epi_K32+baseline.Epi_other if label=='total_epithelial' else baseline[label];now=current.Epi_K32+current.Epi_other if label=='total_epithelial' else current[label];delta=abs(now-old);change=delta/(old+.1);med=float(change.median());q95=float(change.quantile(.95));passed=bool(med<cfg['mapping']['abundance_median_relative_change_max'] and q95<cfg['mapping']['abundance_q95_relative_change_max']);tests.append(passed);stability_history.append({'epoch':epochs,'abundance':label,'mean_absolute_change':float(delta.mean()),'median_relative_change':med,'q95_relative_change':q95,'stable':passed})
   abundance_stable=all(tests);wt(pd.DataFrame(stability_history),stab_path)
  current.to_csv(baseline_path,sep='\t');baseline=current
 converged=bool(not smoke and epochs>=cfg['mapping'].get('minimum_epochs',5000) and loss_stable and abundance_stable);js(O/'convergence.json',{'epochs':epochs,'relative_window_change':rel,'window':w,'plateau_pass':loss_stable,'abundance_stable':abundance_stable,'converged':converged,'technical_smoke':smoke,'elapsed_this_process':time.time()-t0});say('MAP_CONVERGENCE',branch,part,epochs,rel,'abundance_stable',abundance_stable)
 if converged or smoke:break
scvi.settings.seed=seed('mapping_posterior_'+part);ns=16 if smoke else cfg['posterior_samples'];p=m.sample_posterior(num_samples=ns,batch_size=min(1000,a.n_obs),accelerator='gpu',device=1,return_sites=['w_sf'],return_samples=True,summary_frequency=1);draw=p['posterior_samples']['w_sf'];assert draw.shape==(ns,a.n_obs,len(labels)) and np.isfinite(draw).all();np.save(O/'w_sf_draws.npy',draw.astype(np.float32))
for name,fun in [('mean',lambda x:np.mean(x,axis=0)),('median',lambda x:np.median(x,axis=0)),('q05',lambda x:np.quantile(x,.05,axis=0)),('q95',lambda x:np.quantile(x,.95,axis=0))]:
 d=pd.DataFrame(fun(draw),index=a.obs_names,columns=labels);d.index.name='pseudospot_id';d.to_csv(O/f'abundance_{name}.tsv.gz',sep='\t')
k=draw[:,:,labels.index('Epi_K32')];other=draw[:,:,labels.index('Epi_other')];epi=k+other;total=draw.sum(2);unstable=(epi<cfg['ratio_min_epithelial_count']).mean(0)>.05;pr=np.divide(k,epi,out=np.full_like(k,np.nan),where=epi>=cfg['ratio_min_epithelial_count']);pr[:,unstable]=np.nan
ratio=pd.DataFrame({'pseudospot_id':a.obs_names,'epithelial_denominator_unstable':unstable,'fraction_draws_epi_below_min':(epi<cfg['ratio_min_epithelial_count']).mean(0)})
for label,x in [('K32_fraction_all',k/total),('K32_fraction_epithelial',pr),('total_epithelial',epi)]:
 for name,fun in [('mean',lambda z:np.nanmean(z,axis=0)),('median',lambda z:np.nanmedian(z,axis=0)),('q05',lambda z:np.nanquantile(z,.05,axis=0)),('q95',lambda z:np.nanquantile(z,.95,axis=0))]:ratio[label+'_'+name]=fun(x)
wt(ratio,O/'posterior_ratios_and_epithelial.tsv.gz')
# Independent reconstruction summaries, no abundance truth/scenario covariates.
scvi.settings.seed=seed('mapping_reconstruction_'+part);rp=m.sample_posterior(num_samples=16 if smoke else 100,batch_size=min(1000,a.n_obs),accelerator='gpu',device=1,return_sites=['m_g','s_g_gene_add','detection_y_s','alpha_g_inverse'],return_samples=False,summary_frequency=1)['post_sample_means'];rp['w_sf']=draw.mean(0);np.savez_compressed(O/'reconstruction_parameter_means.npz',**rp);ix=np.linspace(0,a.n_obs-1,min(256,a.n_obs),dtype=int);mu=m.module.model.compute_expected(rp,m.adata_manager,ind_x=ix)['mu'];obs=a.X[ix].toarray();rec=pd.DataFrame({'pseudospot_id':a.obs_names[ix],'observed_UMI':obs.sum(1),'expected_UMI':mu.sum(1),'gene_mean_absolute_error':abs(obs-mu).mean(1)});wt(rec,O/'reconstruction_spots.tsv');wt(pd.DataFrame({'gene':a.var_names,'observed_mean':obs.mean(0),'expected_mean':mu.mean(0)}),O/'reconstruction_genes.tsv');js(O/'reconstruction_QC.json',{'log1p_gene_mean_pearson':float(pearsonr(np.log1p(obs.mean(0)),np.log1p(mu.mean(0))).statistic),'total_UMI_ratio':float(mu.sum()/obs.sum()),'all_finite':bool(np.isfinite(mu).all())})
a.write_h5ad(O/'mapping_registered.h5ad',compression='lzf');m.save(cp,overwrite=True,save_anndata=False);js(O/'MAPPING_DONE.json',{'training_completed':True,'converged':converged,'technical_smoke':smoke,'epochs':epochs,'posterior_samples':ns,'n_spots':a.n_obs,'n_genes':a.n_vars,'counts_hash':finger['counts_sha256'],'w_sf_draws_sha256':sha(O/'w_sf_draws.npy'),'biological_identifiability':'requires_truth_evaluation'});say('MAPPING_EXPORTED',branch,part,epochs)
