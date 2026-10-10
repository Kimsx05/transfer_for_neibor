from common import *
import sys,gc
from scipy.stats import spearmanr,pearsonr
from sklearn.metrics import roc_auc_score,average_precision_score,precision_recall_curve,roc_curve
O=R/'06_evaluation';threshold=json.loads((S/'06_evaluation/DEV_FROZEN.json').read_text())['threshold_on_posterior_mean'];labs=json.loads((R/'00_manifest/labels_S31.json').read_text());major=pd.read_csv(S/'01_scope/new_subtype_major_mapping.tsv',sep='\t').set_index('new_reference_label').reference_major.to_dict()
models=[(f'HIST_{b}_{s}',S/f'05_deconvolution/{b}_{s}',S/f'04_pseudospots/{s}','S31',s,True) for b in ['MAIN','NO_P07T'] for s in ['DEV','TEST']]+[(f'{l}_{c}',R/f'05_deconvolution/{l}_{c}',S/'04_pseudospots/DEV' if c=='FULL_DEV' else R/f'04_paired_depth_control/{c}',l,'DEV',False) for l,c in [('S30','FULL_DEV'),('S31','TARGET_UMI'),('S31','FIXED_CAPTURE'),('S30','TARGET_UMI'),('S30','FIXED_CAPTURE')]]
allspots=[];absorption=[];status=[];techcorr=[];curves=[]
for tag,d,inp,ls,split,hist in models:
 done=d/'MAPPING_DONE.json'
 if not done.exists():status.append(dict(model=tag,status='UNAVAILABLE',reason='MAPPING_DONE not yet present'));continue
 state=json.loads(done.read_text());ok=state['converged'];status.append(dict(model=tag,status='COMPLETED' if ok else 'NOT_CONVERGED',epochs=state['epochs'],historical=hist));manifest=pd.read_csv(inp/'manifest.tsv.gz',sep='\t').set_index('pseudospot_id');pred=pd.read_csv(d/'abundance_mean.tsv.gz',sep='\t',index_col=0);truth=pd.read_csv(inp/'truth_counts.tsv.gz',sep='\t',index_col=0);assert pred.index.equals(manifest.index) and truth.index.equals(manifest.index);labels=list(pred);n=len(pred);out=manifest.copy();out['model']=tag;out['labelset']=ls;out['historical']=hist;out['model_converged']=ok
 cache=O/(tag+'_derived.tsv.gz'); cache_fp={'draws_sha256':state.get('draws_sha256',state.get('w_sf_draws_sha256')),'spot_order_hash':key(out.index),'labels':labels,'ratio_min':0.1}; cf=O/(tag+'_derived_fingerprint.json')
 if cache.exists() and cf.exists() and json.loads(cf.read_text())==cache_fp:derived=pd.read_csv(cache,sep='\t',index_col=0);assert derived.index.equals(out.index)
 else:
  draws=np.load(d/'w_sf_draws.npy',mmap_mode='r');assert draws.shape==(1000,n,len(labels));chunks=[]
  for start in range(0,n,256):
   z=np.asarray(draws[:,start:start+256,:]);epi=z[:,:,labels.index('Epithelial')] if ls=='S30' else z[:,:,labels.index('Epi_K32')]+z[:,:,labels.index('Epi_other')];total=z.sum(2);near=(epi<.1).mean(0);dd=pd.DataFrame(index=out.index[start:start+256]);dd['epithelial_denominator_unstable']=near>.05;dd['fraction_draws_epi_below_min']=near;arrays={'epithelial':epi,'total_cells':total,'epithelial_fraction_all':epi/total}
   if ls=='S31':
    k=z[:,:,labels.index('Epi_K32')];rr=np.divide(k,epi,out=np.full_like(k,np.nan),where=epi>=.1);rr[:,near>.05]=np.nan;arrays.update(K32=k,other=z[:,:,labels.index('Epi_other')],K32_fraction_all=k/total,K32_fraction_epithelial=rr)
   for name,a in arrays.items():
    for stat,fn in [('mean',np.nanmean),('median',np.nanmedian),('q05',lambda x,axis:np.nanquantile(x,.05,axis=axis)),('q95',lambda x,axis:np.nanquantile(x,.95,axis=axis))]:dd[name+'_'+stat]=fn(a,axis=0)
   chunks.append(dd)
  derived=pd.concat(chunks);derived.to_csv(cache,sep='\t');js(cf,cache_fp);del draws;gc.collect()
 out=out.join(derived)
 targets={'epithelial':out.N_epithelial,'total_cells':out.N_total,'epithelial_fraction_all':out.N_epithelial/out.N_total}
 if ls=='S31':targets.update(K32=out.N_K32,other=out.N_other,K32_fraction_all=out.N_K32/out.N_total,K32_fraction_epithelial=out.N_K32/out.N_epithelial)
 for endpoint,actual in targets.items():
  out[endpoint+'_true']=actual;out[endpoint+'_error']=out[endpoint+'_mean']-actual;out[endpoint+'_abs_error']=abs(out[endpoint+'_error']);out[endpoint+'_coverage90']=(actual>=out[endpoint+'_q05'])&(actual<=out[endpoint+'_q95']);out[endpoint+'_width90']=out[endpoint+'_q95']-out[endpoint+'_q05']
 if ls=='S31':out['detected']=out.K32_mean>threshold;out['positive']=out.N_K32>0
 # Reconstruction technical summaries from the same frozen checkpoint means, never refit.
 if (d/'technical_terms.tsv.gz').exists():technical=pd.read_csv(d/'technical_terms.tsv.gz',sep='\t',index_col=0)
 else:
  import anndata as ad
  rp=np.load(d/'reconstruction_parameter_means.npz');a=ad.read_h5ad(d/'mapping_registered.h5ad',backed='r');batch=a.obs['_scvi_batch'].to_numpy();assert list(a.obs_names)==list(pred.index);a.file.close();branch='NO_P07T' if 'NO_P07T' in tag else 'MAIN';sig=pd.read_csv(S/f'03_reference/{branch}/signatures.tsv.gz',sep='\t',index_col=0)[labels];det=rp['detection_y_s'].reshape(-1);bg=rp['s_g_gene_add'].sum(1)[batch]*det;ct=(rp['w_sf']@(sig.to_numpy()*rp['m_g'].reshape(-1,1)).sum(0))*det;technical=pd.DataFrame({'detection_y_s_mean':det,'background_expected_UMI_from_means':bg,'cellular_expected_UMI_from_means':ct,'expected_UMI_from_means':bg+ct},index=pred.index)
 out=out.join(technical.drop(columns=[x for x in technical if x in out],errors='ignore'));out['background_fraction_from_means']=out.background_expected_UMI_from_means/out.expected_UMI_from_means
 for label in labels:
  if label in ['Epi_K32','Epi_other','Epithelial']:continue
  dd=out[['dataset','sample','source_group','scenario','depth_condition','model']].copy();dd['label']=label;dd['major']=major[label];dd['predicted']=pred[label];dd['true']=truth[label];dd['error']=dd.predicted-dd.true;dd['abs_error']=abs(dd.error);absorption.append(dd.reset_index())
 for sample,dg in out.groupby('sample'):
  for x in ['pre_full_UMI','post_full_UMI','pre_shared_UMI','post_shared_UMI','thinning_probability']:
   for y in ['detection_y_s_mean','background_fraction_from_means','total_cells_mean','epithelial_mean']:
    rho=float(spearmanr(dg[x],dg[y]).statistic) if dg[x].nunique()>1 and dg[y].nunique()>1 else np.nan;techcorr.append(dict(model=tag,sample=sample,source_group=dg.source_group.iloc[0],dataset=dg.dataset.iloc[0],x=x,y=y,Spearman=rho,n=len(dg)))
 if ls=='S31':
  for kind in ['PR','ROC']:
   if kind=='PR': yy,xx,tt=precision_recall_curve(out.positive,out.K32_mean);tt=np.r_[tt,np.nan]
   else:xx,yy,tt=roc_curve(out.positive,out.K32_mean)
   curves.append(pd.DataFrame({'model':tag,'kind':kind,'x':xx,'y':yy,'threshold':tt,'simulated_positive_fraction':out.positive.mean()}))
 out.reset_index().to_csv(O/(tag+'_per_spot.tsv.gz'),sep='\t',index=False);allspots.append(out.reset_index());say('SUMMARIZED',tag,n)
if not allspots:raise SystemExit('No completed model output')
D=pd.concat(allspots,ignore_index=True);wt(D,O/'per_spot_all.tsv.gz');wt(pd.DataFrame(status),O/'model_evaluation_status.tsv');wt(pd.DataFrame(techcorr),O/'technical_correlations_by_sample.tsv');wt(pd.concat(curves,ignore_index=True),O/'detection_curves.tsv.gz')
endpoints=['epithelial','total_cells','epithelial_fraction_all','K32','other','K32_fraction_all','K32_fraction_epithelial']
def metric(d,endpoint):
 y=d[endpoint+'_true'].to_numpy();z=d[endpoint+'_mean'].to_numpy();finite=np.isfinite(y)&np.isfinite(z);y=y[finite];z=z[finite];err=z-y;v=d.loc[finite];n=len(y);slope=np.polyfit(y,z,1) if n>1 and np.ptp(y)>0 else [np.nan,np.nan]
 return dict(endpoint=endpoint,n_spots=len(d),n_evaluable=n,n_samples=d['sample'].nunique(),n_groups=d.source_group.nunique(),MAE=np.mean(abs(err)) if n else np.nan,RMSE=np.sqrt(np.mean(err**2)) if n else np.nan,bias=np.mean(err) if n else np.nan,predicted_mean=np.mean(z) if n else np.nan,true_mean=np.mean(y) if n else np.nan,slope=slope[0],intercept=slope[1],Pearson=pearsonr(y,z).statistic if n>1 and np.ptp(y)>0 and np.ptp(z)>0 else np.nan,coverage90=v[endpoint+'_coverage90'].mean(),width90=v[endpoint+'_width90'].mean(),unstable_denominator_spots=int(d.epithelial_denominator_unstable.sum()))
rows=[];strata=[[],['dataset'],['sample','source_group','dataset'],['source_group'],['scenario','depth_condition'],['N_K32','depth_condition','N_epithelial'],['dataset','scenario','depth_condition'],['sample','source_group','dataset','scenario','depth_condition']]
for keys in strata:
 for vals,d in D.groupby(['model']+keys,dropna=False):
  vals=vals if isinstance(vals,tuple) else (vals,);base=dict(zip(['model']+keys,vals));base['stratification']='|'.join(keys) or 'pooled'
  for ep in endpoints:
   if ep+'_mean' not in d or d[ep+'_mean'].notna().sum()==0:continue
   rows.append({**base,**metric(d,ep)})
M=pd.DataFrame(rows);wt(M,O/'metrics_all_strata.tsv')
# Equal sampling weights apply before nonlinear RMSE and calibration calculations.
def weighted_metric(d,ep,nested):
 w=(1/d.groupby('sample')['sample'].transform('size')).to_numpy(float,copy=True)
 if nested:w/=d.groupby('source_group')['sample'].transform('nunique').to_numpy()
 y=d[ep+'_true'].to_numpy(float);z=d[ep+'_mean'].to_numpy(float);valid=np.isfinite(y)&np.isfinite(z);y=y[valid];z=z[valid];w=w[valid];w=w/w.sum();err=z-y;xm=np.sum(w*y);ym=np.sum(w*z);vx=np.sum(w*(y-xm)**2) if np.ptp(y)>0 else 0.;vy=np.sum(w*(z-ym)**2) if np.ptp(z)>0 else 0.;cov=np.sum(w*(y-xm)*(z-ym));sl=cov/vx if vx>0 else np.nan
 return dict(endpoint=ep,n_spots=len(d),n_evaluable=int(valid.sum()),n_samples=d['sample'].nunique(),n_groups=d.source_group.nunique(),MAE=np.sum(w*abs(err)),RMSE=np.sqrt(np.sum(w*err**2)),bias=np.sum(w*err),predicted_mean=ym,true_mean=xm,slope=sl,intercept=ym-sl*xm,Pearson=cov/np.sqrt(vx*vy) if vx>0 and vy>0 else np.nan,coverage90=np.sum(w*d.loc[valid,ep+'_coverage90']),width90=np.sum(w*d.loc[valid,ep+'_width90']))
weighted=[]
for tag,d in D.groupby('model'):
 for nested in [False,True]:
  for ep in endpoints:
   if ep+'_mean' in d and d[ep+'_mean'].notna().any():weighted.append(dict(model=tag,weighting='sample_equal_within_group_then_group_equal' if nested else 'sample_equal',**weighted_metric(d,ep,nested)))
W=pd.DataFrame(weighted);sample=W[W.weighting=='sample_equal'];wt(W,O/'metrics_equal_weight.tsv')
positive=[]
for tag,d in D[(D.labelset=='S31')&(D.N_K32>0)].groupby('model'):
 for nested in [False,True]:positive.append(dict(model=tag,weighting='sample_equal_within_group_then_group_equal' if nested else 'sample_equal',**weighted_metric(d,'K32',nested)))
wt(pd.DataFrame(positive),O/'positive_K32_metrics_equal_weight.tsv')
# Detection denominators and undefined precision are explicit; AUC is unavailable for one-class sources.
def detection(d):
 y=d.N_K32.to_numpy()>0;z=d.K32_mean.to_numpy();hit=z>threshold;tp=int((hit&y).sum());fp=int((hit&~y).sum());n0=int((~y).sum());n1=int(y.sum());nh=int(hit.sum());return dict(n=len(d),positive=n1,negative=n0,detected=nh,TP=tp,FP=fp,simulated_positive_fraction=y.mean(),FPR=fp/n0 if n0 else np.nan,recall=tp/n1 if n1 else np.nan,precision=tp/nh if nh else np.nan,precision_undefined=nh==0,AUC_undefined=not(n0 and n1),AUROC=roc_auc_score(y,z) if n0 and n1 else np.nan,AUPRC=average_precision_score(y,z) if n0 and n1 else np.nan,negative_mean=np.mean(z[~y]) if n0 else np.nan,negative_median=np.median(z[~y]) if n0 else np.nan,negative_q95=np.quantile(z[~y],.95) if n0 else np.nan,negative_max=np.max(z[~y]) if n0 else np.nan)
cls=[]
for keys in strata:
 for vals,d in D[D.labelset=='S31'].groupby(['model']+keys,dropna=False):
  vals=vals if isinstance(vals,tuple) else (vals,);cls.append({**dict(zip(['model']+keys,vals)),'stratification':'|'.join(keys) or 'pooled',**detection(d)})
C=pd.DataFrame(cls);wt(C,O/'detection_all_strata.tsv');c=C[C.stratification=='sample|source_group|dataset'];fields=['simulated_positive_fraction','FPR','recall','precision','AUROC','AUPRC','negative_mean'];cs=c.groupby('model')[fields].mean();cg=c.groupby(['model','source_group'])[fields].mean().groupby('model').mean();cs['weighting']='sample_equal_defined_sources';cg['weighting']='sample_equal_within_group_then_group_equal_defined_sources';cc=pd.concat([cs,cg]).reset_index();undef=c.groupby('model').agg(n_samples=('sample','size'),undefined_precision_samples=('precision_undefined','sum'),undefined_AUC_samples=('AUC_undefined','sum')).reset_index();wt(cc.merge(undef),O/'detection_equal_weight.tsv')
# Non-epithelial absorption, preserving label and existing major mapping.
A=pd.concat(absorption,ignore_index=True);wt(A,O/'non_epithelial_per_spot.tsv.gz');sm=A.groupby(['model','dataset','source_group','sample','label','major'])[['predicted','true','error','abs_error']].mean().reset_index();wt(sm,O/'non_epithelial_by_sample.tsv');gm=sm.groupby(['model','source_group','label','major'])[['predicted','true','error','abs_error']].mean().groupby(['model','label','major']).mean().reset_index();wt(gm,O/'non_epithelial_group_equal.tsv');am=A.groupby(['model','pseudospot_id','dataset','source_group','sample','major'])[['predicted','true']].sum().reset_index();am['error']=am.predicted-am.true;am['abs_error']=abs(am.error);wt(am.groupby(['model','dataset','source_group','sample','major'])[['predicted','true','error','abs_error']].mean().reset_index(),O/'non_epithelial_major_by_sample.tsv')
# Exact spot/bag paired differences before averaging.
pairs=[]
def pair(a,b,match,comparison):
 aa=D[D.model==a];bb=D[D.model==b]
 if aa.empty or bb.empty:return
 keep=['dataset','sample','source_group','scenario','depth_condition','N_K32','N_epithelial','pre_full_UMI']
 variables=[x for x in ['epithelial_mean','epithelial_abs_error','total_cells_mean','total_cells_abs_error','epithelial_fraction_all_mean','K32_mean','K32_abs_error','K32_fraction_epithelial_mean','detection_y_s_mean','background_fraction_from_means','post_full_UMI','thinning_probability'] if x in aa and x in bb and aa[x].notna().any() and bb[x].notna().any()]
 p=aa[list(dict.fromkeys(match+keep+variables))].merge(bb[match+variables],on=match,suffixes=('_A','_B'),validate='one_to_one');assert len(p)==len(aa)==len(bb);p['comparison']=comparison;p['model_A']=a;p['model_B']=b
 for x in variables:p['delta_B_minus_A__'+x]=p[x+'_B']-p[x+'_A']
 pairs.append(p)
pair('HIST_MAIN_DEV','S30_FULL_DEV',['pseudospot_id'],'full_DEV_S30_minus_S31')
pair('HIST_MAIN_TEST','HIST_NO_P07T_TEST',['pseudospot_id'],'historical_TEST_NO_P07T_minus_MAIN')
for rule in ['TARGET_UMI','FIXED_CAPTURE']:pair('S31_'+rule,'S30_'+rule,['bag_id','depth_condition'],'control_'+rule+'_S30_minus_S31')
for ls in ['S31','S30']:pair(ls+'_TARGET_UMI',ls+'_FIXED_CAPTURE',['bag_id','depth_condition'],'control_'+ls+'_FIXED_minus_TARGET')
if pairs:
 PP=pd.concat(pairs,ignore_index=True);wt(PP,O/'paired_per_spot.tsv.gz');cols=[x for x in PP if x.startswith('delta_')];ps=PP.groupby(['comparison','dataset','source_group','sample'])[cols].mean().reset_index();wt(ps,O/'paired_by_sample.tsv');pg=ps.groupby(['comparison','source_group'])[cols].mean().reset_index();wt(pg,O/'paired_by_source_group.tsv');wt(pg.groupby('comparison')[cols].mean().reset_index(),O/'paired_group_equal.tsv');wt(PP.groupby(['comparison','scenario','depth_condition','dataset','source_group','sample'])[cols].mean().reset_index(),O/'paired_stratified_samples.tsv')
# Historical reproduction compares exact sample equal endpoint definitions only.
old=pd.read_csv(S/'06_evaluation/metrics_overall.tsv',sep='\t');wt(old,R/'02_existing_output_audit/old_metrics_original.tsv');wt(sample[sample.model.str.startswith('HIST_')],R/'02_existing_output_audit/recomputed_sample_equal.tsv');wt(c[c.model.str.startswith('HIST_')],R/'02_existing_output_audit/historical_detection_by_sample.tsv')
js(O/'EVALUATION_DONE.json',{'time':datetime.datetime.now().isoformat(),'threshold':threshold,'threshold_refit':False,'models_summarized':list(D.model.unique()),'source_weighting':'Within each source_group average sample metrics, then equally average groups; macro precision averages defined sources and explicit undefined counts supplied','no_total_cell_calibration':'Truth constant20, slope/intercept/correlation undefined','S30_K32':'not applicable; NA, never zero or truth allocation','draws':'1000; epithelial sums and all ratios computed within draw','denominator':'Epi<0.1 in >5%draws -> epithelial K32 ratio NA; counts retained','historical_TEST':'Previously inspected; not new independent validation','per_spot_tests':False});say('EVALUATION_DONE')
