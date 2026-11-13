from common import *
from expression import data
import torch,scvi,cellcharter as cc,anndata as ad
from sklearn.metrics import adjusted_rand_score
import resource
n=pd.read_parquet(V/'Cervilla/nodes.parquet');ids=n[n.analysis_qc].sample(1200,random_state=123).bin_id.to_numpy();x,ids,gn,sy=counts('Cervilla',ids=ids);tot=np.asarray(x.sum(0)).ravel();gi=np.argsort(tot)[-500:];a=data(x[:,gi],ids,gn[gi],'Cervilla');scvi.model.SCVI.setup_anndata(a,batch_key='section');m=scvi.model.SCVI(a,n_latent=10,n_layers=2,use_layer_norm='both',use_batch_norm='none',encode_covariates=True);start=time.time();m.train(max_epochs=2,batch_size=256,accelerator='gpu' if torch.cuda.is_available() else 'cpu',devices=1,enable_progress_bar=False);z=m.get_latent_representation();q=a[:100].copy();q.obs['section']='MPBC2';qm=scvi.model.SCVI.load_query_data(q,m);qm.train(max_epochs=1,batch_size=128,accelerator='gpu' if torch.cuda.is_available() else 'cpu',devices=1,enable_progress_bar=False,plan_kwargs={'weight_decay':0.0});assert qm.get_latent_representation().shape==(100,10)
f=cc.tl.GaussianMixture(n_clusters=2,random_state=123,trainer_params={'accelerator':'gpu' if torch.cuda.is_available() else 'cpu','devices':1,'max_epochs':10,'enable_progress_bar':False});f.fit(z);p=f.predict_proba(z);assert np.allclose(p.sum(1),1,atol=1e-5)
import harmonypy
h=harmonypy.run_harmony(z.T,pd.DataFrame({'section':['a']*600+['b']*600}),'section',theta=1,lamb=5,nclust=10,max_iter_harmony=2,random_state=123,ncores=4)
(R/'tables/expression_pilot.json').write_text(json.dumps(dict(actual_QC_bins=1200,genes=500,latent=10,epochs=2,seconds=time.time()-start,cuda=torch.cuda.is_available(),cuda_peak_MiB=torch.cuda.max_memory_allocated()/2**20 if torch.cuda.is_available() else None,maxRSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,query_new_batch_passed=True,GMM_predict_proba_passed=True,Harmony_fixed_parameters_passed=True,pilot_not_full_result=True),indent=2));print('EXPRESSION_PILOT_PASSED',flush=True)
