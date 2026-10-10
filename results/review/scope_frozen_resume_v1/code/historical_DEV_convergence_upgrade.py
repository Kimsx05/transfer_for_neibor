# Historical one-time migration executed during DEV; retained for provenance.
# Do not rerun against active outputs. Use run_all.sh and checkpoint validation for continuation.
from pathlib import Path
import json,shutil
r=Path('/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_k32_reference_identifiability_v5/scope_frozen_resume_v1');p=r/'code/04_mapping.py';s=p.read_text();a=s.index('converged=False;t0=time.time()');b=s.index("scvi.settings.seed=seed('mapping_posterior_'+part);ns=",a)
new='''converged=False;t0=time.time();maxep=50 if smoke else cfg['mapping']['max_epochs']
baseline_path=O/'abundance_checkpoint_mean.tsv.gz';baseline=pd.read_csv(baseline_path,sep='\\t',index_col=0) if baseline_path.exists() else None
stab_path=O/'abundance_convergence_history.tsv';stability_history=pd.read_csv(stab_path,sep='\\t').to_dict('records') if stab_path.exists() else []
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
  current.to_csv(baseline_path,sep='\\t');baseline=current
 converged=bool(not smoke and epochs>=cfg['mapping'].get('minimum_epochs',5000) and loss_stable and abundance_stable);js(O/'convergence.json',{'epochs':epochs,'relative_window_change':rel,'window':w,'plateau_pass':loss_stable,'abundance_stable':abundance_stable,'converged':converged,'technical_smoke':smoke,'elapsed_this_process':time.time()-t0});say('MAP_CONVERGENCE',branch,part,epochs,rel,'abundance_stable',abundance_stable)
 if converged or smoke:break
'''
s=s[:a]+new+s[b:];p.write_text(s)
# Preserve the completed 3000-epoch stage as a diagnostic; its numerical plateau
# is explicitly superseded by the failed abundance stability check.
old=r/'05_deconvolution/MAIN_DEV';archive=r/'05_deconvolution/ARCHIVE_MAIN_DEV_3000';assert not archive.exists();old.rename(archive);old.mkdir();probe=r/'05_deconvolution/DEV_convergence_probe_5000';shutil.copytree(probe/'model',old/'model');shutil.copyfile(probe/'abundance_mean.tsv.gz',old/'abundance_checkpoint_mean.tsv.gz')
import csv
with (old/'training_history.tsv').open('w') as f:
 f.write((archive/'training_history.tsv').read_text());f.write('\n'.join((probe/'training_history.tsv').read_text().splitlines()[1:])+'\n')
review={'scope':'DEV-only technical convergence amendment before TEST','baseline_failure':'3000→5000 epochs: K32 median relative abundance change0.2473; Epi_other0.3266 despite ELBO plateau','decision':'Continue under unchanged priors/features with paired mean-abundance stability as well as ELBO convergence; no truth-based selection of best epoch','mapping_overrides':{'minimum_epochs':5000,'abundance_stability_samples':256,'abundance_median_relative_change_max':0.05,'abundance_q95_relative_change_max':0.20},'primary_abundances':['Epi_K32','Epi_other','total_epithelial'],'maximum_epochs':30000,'training_rule_for_TEST':'Same loss plus abundance-stability rule, without accessing TEST truth; no TEST-based setting changes','TEST_predictions_seen':False}
(r/'06_evaluation/DEV_CONVERGENCE_REVIEW.json').write_text(json.dumps(review,indent=2)+'\n')
import ast;ast.parse(s)
print('Adaptive abundance convergence enabled; MAIN resumes from preserved5000-epoch checkpoint.')
