from common import *
import sys,anndata as ad,torch,scvi
from cell2location.models import Cell2location
from lightning.pytorch.callbacks import Callback
from scipy.stats import pearsonr
name=sys.argv[1];mode=sys.argv[2] if len(sys.argv)>2 else 'block';O=outmap(name);O.mkdir(parents=True,exist_ok=True);cfg=config();c=cfg['mapping'];ref=S/'03_reference/MAIN' if name=='DEV_OLD' else R/'01_reference';st=json.loads((ref/'REFERENCE_DONE.json').read_text());inputpath=R/'02_pseudospots/counts.h5ad' if name.startswith('DEV') else O/'counts.h5ad';a=ad.read_h5ad(inputpath);assert list(a.obs.columns)==['sample'];sig=pd.read_csv(ref/'signatures.tsv.gz',sep='\t',index_col=0);labels=list(sig.columns);kl='Epi_K32' if name=='DEV_OLD' else K;ll='Epi_other' if name=='DEV_OLD' else L;assert list(sig.index)==list(a.var_names) and len(labels)==31 and np.isfinite(sig.to_numpy()).all();prior={'N_cells_per_location':20,'detection_alpha':20} if name.startswith('DEV') else json.loads((O/'INPUT.json').read_text());torch.set_num_threads(1);torch.set_num_interop_threads(1);fitseed=seed('v7_mapping_'+('DEV' if name.startswith('DEV') else name));postseed=seed('v7_posterior_'+('DEV' if name.startswith('DEV') else name));scvi.settings.seed=fitseed;a.obs['sample']=pd.Categorical(a.obs['sample']);Cell2location.setup_anndata(a,batch_key='sample');finger={'counts_sha256':sha(inputpath),'signature_sha256':sha(ref/'signatures.tsv.gz'),'gene_hash':key(a.var_names),'spot_hash':key(a.obs_names),'labels':labels,'fit_seed':fitseed,'post_seed':postseed,'N_cells_per_location':prior['N_cells_per_location'],'detection_alpha':prior['detection_alpha'],'lr':.002,'mapping_settings':dict(c),'optimizer_resume':'Pyro optim get_state/set_state if actual training plan exposes it; saved per checkpoint'}
if (O/'run_fingerprint.json').exists():assert json.loads((O/'run_fingerprint.json').read_text())==finger
else:js(O/'run_fingerprint.json',finger)
if '--resume-seconds' in sys.argv:
 seconds=float(sys.argv[sys.argv.index('--resume-seconds')+1]);assert 0<seconds<=3600;cfg['training_deadline']=time.time()+seconds;js(O/('explicit_manual_resume_'+str(int(time.time()))+'.json'),{'authorized_by':'explicit --resume-seconds command','seconds':seconds,'start_epoch_time':time.time()})
cp=O/'model';history=pd.read_csv(O/'training_history.tsv',sep='\t').loss.tolist() if (O/'training_history.tsv').exists() else [];epochs=len(history);m=Cell2location.load(cp,adata=a,accelerator='gpu',device=[0]) if (cp/'model.pt').exists() else Cell2location(a,cell_state_df=sig,N_cells_per_location=prior['N_cells_per_location'],detection_alpha=prior['detection_alpha']);assert list(m.factor_names_)==labels;conv=json.loads((O/'convergence.json').read_text()) if (O/'convergence.json').exists() else {'epochs':epochs,'converged':False};start=time.time()
def save(trainer=None,plan=None,newhist=None):
 m.save(cp,overwrite=True,save_anndata=False)
 if newhist is not None:wt(pd.DataFrame({'epoch':np.arange(1,len(newhist)+1),'loss':newhist}),O/'training_history.tsv')
 if plan is not None and hasattr(plan,'optim') and hasattr(plan.optim,'get_state'):
  state=plan.optim.get_state();torch.save(state,O/'pyro_optimizer.pt');steps=[float(v['step']) for opt in state.values() for v in opt.get('state',{}).values() if 'step' in v];js(O/'OPTIMIZER_AUDIT.json',{'n_parameter_states':len(state),'min_step':min(steps) if steps else None,'max_step':max(steps) if steps else None,'saved_epoch':len(newhist) if newhist is not None else None,'API':'PyroOptim.get_state/set_state'})
class Progress(Callback):
 def __init__(self):self.done=0;self.lastsave=time.time();self.loss=[];self.batch_losses=[]
 def on_train_batch_end(self,trainer,pl_module,outputs,batch,batch_idx):
  loss=outputs["loss"] if isinstance(outputs,dict) else outputs;self.batch_losses.append(float(loss.detach().cpu()))
 def on_train_start(self,trainer,pl_module):
  if (O/'pyro_optimizer.pt').exists():pl_module.optim.set_state(torch.load(O/'pyro_optimizer.pt',map_location='cuda:0',weights_only=False));say('OPTIMIZER_RESTORED',name)
 def on_train_epoch_end(self,trainer,pl_module):
  self.done=trainer.current_epoch+1;e=epochs+self.done;assert self.batch_losses;self.loss.append(float(np.mean(self.batch_losses)));self.batch_losses=[];assert np.isfinite(self.loss[-1]);now=time.time()
  if e%100==0:js(O/'progress.json',{'stage':'TRAINING','actual_epochs':e,'elapsed_block_seconds':now-start,'epoch_seconds':(now-start)/self.done,'time_epoch':now,'pid':os.getpid()})
  if now-self.lastsave>=c['checkpoint_seconds']:save(trainer,pl_module,history+self.loss);self.lastsave=now;say('CHECKPOINT',name,e)
  if now-start>=c['block_seconds'] or now>=cfg['training_deadline']:trainer.should_stop=True
if mode=='block' and not conv.get('converged',False) and epochs<c['max_epochs'] and time.time()<cfg['training_deadline']:
 target=3000 if epochs<3000 else 3000+((epochs-3000)//2000+1)*2000;n=min(target-epochs,c['max_epochs']-epochs);cb=Progress();say('TRAIN_START',name,epochs,n,a.shape)
 m.train(max_epochs=n,batch_size=None,train_size=1,lr=.002,accelerator='gpu',device=[0],enable_progress_bar=False,enable_checkpointing=False,callbacks=[cb]);assert cb.done>0;history.extend(cb.loss);epochs=len(history);save(m.trainer,m.trainer.lightning_module,history);say('BLOCK_FINISHED',name,epochs,'seconds',time.time()-start)
 # Formal stability only at target milestones, never shortened comparison intervals.
 if epochs>=3000 and ((epochs-3000)%2000==0 or epochs==c['max_epochs']):
  scvi.settings.seed=postseed;ps=m.sample_posterior(num_samples=256,batch_size=min(1000,a.n_obs),accelerator='gpu',device=[0],return_sites=['w_sf'],return_samples=True,summary_frequency=1)['posterior_samples']['w_sf'];assert np.isfinite(ps).all();current=pd.DataFrame(ps.mean(0),index=a.obs_names,columns=labels);current.to_csv(O/f'checkpoint_mean_{epochs}.tsv.gz',sep='\t');np.save(O/'latest_256_draws.npy',ps.astype(np.float32));previous=epochs-2000
  if not (O/f'checkpoint_mean_{previous}.tsv.gz').exists():
   eligible=[int(p.name.split('_')[-1].split('.')[0]) for p in O.glob('checkpoint_mean_*.tsv.gz') if int(p.name.split('_')[-1].split('.')[0])<=epochs-2000];previous=max(eligible) if eligible else previous
  bp=O/f'checkpoint_mean_{previous}.tsv.gz';passed=False;rec=[]
  if bp.exists():
   old=pd.read_csv(bp,sep='\t',index_col=0);assert old.index.equals(current.index)
   for lab in [kl,ll,'epithelial']:
    before=old[kl]+old[ll] if lab=='epithelial' else old[lab];after=current[kl]+current[ll] if lab=='epithelial' else current[lab];d=abs(after-before)/(before+.1);rec.append({'epoch':epochs,'previous_epoch':previous,'label':lab,'median_relative_change':float(d.median()),'q95_relative_change':float(d.quantile(.95)),'passed':bool(d.median()<.05 and d.quantile(.95)<.2)})
   passed=all(x['passed'] for x in rec);hp=O/'abundance_convergence.tsv';allrec=pd.concat([pd.read_csv(hp,sep='\t'),pd.DataFrame(rec)],ignore_index=True) if hp.exists() else pd.DataFrame(rec);wt(allrec,hp)
  losschange=abs(np.mean(history[-500:])-np.mean(history[-1000:-500]))/abs(np.mean(history[-500:]));conv={'epochs':epochs,'converged':bool(epochs>=5000 and passed and losschange<.005),'loss_relative_change':losschange,'abundance_stable':passed,'comparison_epoch_gap':epochs-previous if bp.exists() else None,'reference_numerical_status':st.get('status'),'wall_seconds_last_block':time.time()-start};js(O/'convergence.json',conv)
 else:conv.update(epochs=epochs,converged=False);js(O/'convergence.json',conv)
 block={'start_epoch':epochs-cb.done,'end_epoch':epochs,'requested_epochs':n,'actual_epochs_this_block':cb.done,'elapsed_seconds':time.time()-start,'time_epoch':time.time(),'pid':os.getpid(),'converged':conv['converged'],'fit_seed':fitseed,'wall_limited':cb.done<n,'optimizer_state_preserved':(O/'pyro_optimizer.pt').exists()};js(O/'BLOCK_DONE.json',block)
 bp=O/'TRAINING_BLOCKS.tsv';blocks=pd.concat([pd.read_csv(bp,sep='\t'),pd.DataFrame([block])],ignore_index=True) if bp.exists() else pd.DataFrame([block]);wt(blocks,bp)
# Export a 256-draw provisional overview at every first stage; final 1000 only when stable or final request.
final=mode=='final' or conv.get('converged',False);n_draws=1000 if final else 256
if final or epochs>=3000:
 scvi.settings.seed=postseed;draw=m.sample_posterior(num_samples=n_draws,batch_size=min(1000,a.n_obs),accelerator='gpu',device=[0],return_sites=['w_sf'],return_samples=True,summary_frequency=1)['posterior_samples']['w_sf'];assert draw.shape==(n_draws,a.n_obs,len(labels)) and np.isfinite(draw).all() and (draw>=0).all();np.save(O/'w_sf_draws.npy',draw.astype(np.float32));summary={}
 for nm,fn in [('mean',np.mean),('median',np.median),('q05',lambda x,axis:np.quantile(x,.05,axis=axis)),('q95',lambda x,axis:np.quantile(x,.95,axis=axis))]:
  d=pd.DataFrame(fn(draw,axis=0),index=a.obs_names,columns=labels);d.index.name='spot_id';d.to_csv(O/f'abundance_{nm}.tsv.gz',sep='\t')
 kk=draw[:,:,labels.index(kl)];epi=kk+draw[:,:,labels.index(ll)];total=draw.sum(2);bad=(epi<.1).mean(0)>.05;ratio=np.divide(kk,epi,out=np.full_like(kk,np.nan),where=epi>=.1);ratio[:,bad]=np.nan;der=pd.DataFrame({'spot_id':a.obs_names,'fraction_draws_epi_below_0.1':(epi<.1).mean(0),'ratio_unstable':bad,'ratio_of_posterior_means':kk.mean(0)/epi.mean(0)})
 for lab,x in {'Epi_total':epi,'I_state':ratio,'total_cells':total}.items():
  for nm,fn in [('mean',np.nanmean),('median',np.nanmedian),('q05',lambda x,axis:np.nanquantile(x,.05,axis=axis)),('q95',lambda x,axis:np.nanquantile(x,.95,axis=axis))]:der[lab+'_'+nm]=fn(x,axis=0)
 wt(der,O/'posterior_derived.tsv.gz');status={'execution_status':'COMPLETED' if conv.get('converged') or (final and epochs>=c['max_epochs']) else 'RUNNING' if not final else 'TIME_LIMITED','stop_reason':'NUMERICAL_STABLE' if conv.get('converged') else 'EPOCH_CAP' if epochs>=c['max_epochs'] else 'WALL_TIME' if final else 'INTERMEDIATE','numerical_status':'STABLE' if conv.get('converged') else 'PROVISIONAL_NOT_CONVERGED','biological_evidence':'EVALUATION_PENDING','epochs':epochs,'actual_n_draws':n_draws,'n_spots':a.n_obs,'n_labels':len(labels),'reference_status':st.get('status'),'reference_sha256':finger['signature_sha256'],'counts_sha256':finger['counts_sha256'],'convergence':conv,'updated_epoch':time.time()};js(O/'RESULT_STATUS.json',status)
 if final:
  status['draws_sha256']=sha(O/'w_sf_draws.npy');js(O/'RESULT_STATUS.json',status)
  scvi.settings.seed=seed('v7_tech_'+name);rp=m.sample_posterior(num_samples=100,batch_size=min(1000,a.n_obs),accelerator='gpu',device=[0],return_sites=['m_g','s_g_gene_add','detection_y_s','detection_mean_y_e','alpha_g_inverse'],return_samples=False,summary_frequency=1)['post_sample_means'];np.savez_compressed(O/'technical_parameter_means.npz',**rp);rp['w_sf']=draw.mean(0);ix=np.linspace(0,a.n_obs-1,min(256,a.n_obs),dtype=int);mu=m.module.model.compute_expected(rp,m.adata_manager,ind_x=ix)['mu'];obs=a.X[ix].toarray();assert np.isfinite(mu).all();wt(pd.DataFrame({'gene':a.var_names,'observed_mean':obs.mean(0),'expected_mean':mu.mean(0)}),O/'reconstruction_genes.tsv');js(O/'reconstruction_QC.json',{'log1p_gene_mean_pearson':float(pearsonr(np.log1p(obs.mean(0)),np.log1p(mu.mean(0))).statistic),'expected_observed_UMI_ratio':float(mu.sum()/obs.sum()),'finite':True,'interpretation':'plug-in posterior means, not joint-posterior predictive interval'});js(O/'MAPPING_DONE.json',status)
say('EXIT',name,epochs,conv.get('converged'),mode)
