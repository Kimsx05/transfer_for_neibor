# FROZEN_VERSION_GUARD: reproduce in a new version directory; never overwrite v1.
from pathlib import Path as _FreezePath
if (_FreezePath(__file__).resolve().parents[1] / 'metadata/FROZEN.json').exists():
    raise RuntimeError('Frozen v1 is immutable; create v2 for changes.')
"""Reference calibration + observed-marker hierarchical annotation. No physical-neighbor voting."""
import os
from pathlib import Path
root=Path(__file__).resolve().parents[1]
os.environ.update(CELLTYPIST_FOLDER=str(root/'models/celltypist_cache'),OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='8')
if hasattr(os,'sched_setaffinity'):os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:8])
import numpy as np,pandas as pd,anndata as ad,scanpy as sc,celltypist,json,gc
from scipy import sparse
from sklearn.metrics import classification_report
from marker_rules import MAJOR,SUBTYPE,ANCHORS
np.random.seed(123);rng=np.random.default_rng(123)
config=json.loads((root/'metadata/classifier_config.json').read_text());genes=pd.Index(config['gene_list']);a=ad.read_h5ad(root/'cell_annotation.h5ad');held=ad.read_h5ad(root/'inputs/reference_holdout_temporary.h5ad')
gm=pd.read_csv(root/'metadata/gene_mapping_Cervilla.tsv',sep='\t');valid=gm.reference_symbol.notna();idx=genes.get_indexer(gm.loc[valid,'reference_symbol']);assert (idx>=0).all();transform=sparse.csr_matrix((np.ones(valid.sum(),dtype='int32'),(gm.loc[valid,'hd_feature_index'],idx)),shape=(a.n_vars,len(genes)));raw=(a.X@transform).tocsr();raw.sum_duplicates();raw.eliminate_zeros()
assert np.array_equal(held.var_names,genes)
def normalize(x):
 x=x.astype('float32').tocsr();tot=np.asarray(x.sum(1)).ravel();x=sparse.diags(1e4/np.maximum(tot,1)).dot(x).tocsr();x.data=np.log1p(x.data);return x

def predict(x,model):
 # Public CellTypist Model prediction with explicitly aligned sparse StandardScaler input.
 vi=genes.get_indexer(model.classifier.features);assert (vi>=0).all();out=[]
 for start in range(0,x.shape[0],5000):
  scaled=model.scaler.transform(x[start:start+5000,vi]).tocsr();scaled.data=np.minimum(scaled.data,10)
  decision,prob,label=model.predict_labels_and_prob(scaled)
  order=np.argsort(-decision,axis=1,kind='stable');rows=np.arange(len(label));classes=model.classifier.classes_
  out.append(pd.DataFrame({'prediction':classes[order[:,0]],'score':prob[rows,order[:,0]],'margin':decision[rows,order[:,0]]-decision[rows,order[:,1]],'second_candidate':classes[order[:,1]],'second_score':prob[rows,order[:,1]],'third_candidate':classes[order[:,2]] if len(classes)>2 else '', 'third_score':prob[rows,order[:,2]] if len(classes)>2 else np.nan}))
 return pd.concat(out,ignore_index=True)

def thin(x,target):
 x=x.copy().tocsr();total=np.asarray(x.sum(1)).ravel();p=np.minimum(1,target/np.maximum(total,1));x.data=rng.binomial(x.data.astype('int64'),np.repeat(p,np.diff(x.indptr))).astype('float32');x.eliminate_zeros();return x

def calibrate(p,true,modelname,condition):
 result=[]
 for label in sorted(set(true)|set(p.prediction)):
  q=p[p.prediction.eq(label)];truth=np.asarray(true)[q.index];margin=q.margin.to_numpy();correct=truth==label
  best=None
  # Empirical acceptance calibration, not a model parameter sweep or independent test accuracy.
  for threshold in np.unique(np.quantile(margin,np.linspace(0,.95,20))) if len(margin) else []:
   m=margin>=threshold;n=int(m.sum());precision=float(correct[m].mean()) if n else 0
   if n>=20 and precision>=.8:
    best=dict(margin_threshold=float(threshold),accepted_calibration_n=n,accepted_calibration_precision=precision);break
  row=dict(model=modelname,condition=condition,label=label,true_n=int((np.asarray(true)==label).sum()),predicted_n=len(q),unthresholded_precision=float(correct.mean()) if len(q) else None,calibrated=best is not None,margin_threshold=best['margin_threshold'] if best else None,accepted_calibration_n=best['accepted_calibration_n'] if best else 0,accepted_calibration_precision=best['accepted_calibration_precision'] if best else None)
  result.append(row)
 return result

def evidence(x,rules):
 ng={};tot={};present={}
 for label,markers in rules.items():
  ix=genes.get_indexer(markers);used=[g for g,i in zip(markers,ix) if i>=0];ix=ix[ix>=0];z=x[:,ix]
  ng[label]=z.getnnz(1);tot[label]=np.asarray(z.sum(1)).ravel();present[label]=used
 return ng,tot,present

counts=np.asarray(raw.sum(1)).ravel();low=a.obs.low_information.to_numpy();targets={'hd_median':max(10,int(np.median(counts[~low]))),'low_information_median':max(10,int(np.median(counts[low])))}
views={'full':held.X,'hd_median':thin(held.X,targets['hd_median']),'low_information_median':thin(held.X,targets['low_information_median'])}
cal=[];evalrows=[];hdpred={};models={}
logHD=normalize(raw)
for name in ['major','subtype','b_plasma']:
 model=celltypist.models.Model.load(str(root/f'models/celltypist_{name}.pkl'));models[name]=model
 subset=held.obs.major_label.isin(['B','Plasma']).to_numpy() if name=='b_plasma' else np.ones(held.n_obs,bool)
 truth=held.obs.loc[subset,'major_label' if name=='major' else 'reference_celltype_v2'].astype(str).to_numpy()
 for condition,x in views.items():
  pred=predict(normalize(x[subset]),model);cal+=calibrate(pred,truth,name,condition)
  report=classification_report(truth,pred.prediction,output_dict=True,zero_division=0)
  for label,metrics in report.items():
   if isinstance(metrics,dict):evalrows.append({'model':name,'condition':condition,'label':label,**metrics})
  pred['true_label']=truth;pred['source_cell_id']=held.obs_names[subset];pred.to_parquet(root/f'metadata/reference_predictions_{name}_{condition}.parquet',index=False,compression='zstd')
 print('PREDICT HD',name,flush=True);hdpred[name]=predict(logHD,model)
 hdpred[name].to_parquet(root/f'metadata/HD_predictions_{name}.parquet',index=False,compression='zstd')
caldf=pd.DataFrame(cal);caldf.to_csv(root/'metadata/reference_acceptance_calibration.tsv',sep='\t',index=False);pd.DataFrame(evalrows).to_csv(root/'metadata/reference_depth_validation.tsv',sep='\t',index=False)
lookup={(r['model'],r['condition'],r['label']):r for r in cal}
def accepted(pred,name):
 ok=np.zeros(len(pred),bool)
 for (label,is_low),inds in pred.assign(low=low).groupby(['prediction','low']).indices.items():
  full=lookup.get((name,'full',label),{});shallow=lookup.get((name,'low_information_median' if is_low else 'hd_median',label),{})
  if full.get('calibrated') and shallow.get('calibrated'):
   threshold=shallow['margin_threshold']
   ok[inds]=pred.margin.to_numpy()[inds]>=threshold
 return ok
ng,marker_counts,present=evidence(raw,MAJOR);sg,scount,spresent=evidence(raw,SUBTYPE);ag,_,apresent=evidence(raw,ANCHORS)
major=hdpred['major'];fine=hdpred['subtype'].copy();bselect=major.prediction.isin(['B','Plasma']).to_numpy();fine.loc[bselect,:]=hdpred['b_plasma'].loc[bselect,:].to_numpy()
major_ok=accepted(major,'major');fine_ok=accepted(hdpred['subtype'],'subtype');fine_ok[bselect]=accepted(hdpred['b_plasma'],'b_plasma')[bselect]
dictionary=pd.read_csv(root/'metadata/reference_label_dictionary.tsv',sep='\t').set_index('reference_celltype_v2').major_label.to_dict()
major_marker=np.array([ng[x][i] for i,x in enumerate(major.prediction)]);major_ok&=major_marker>=2
# Require observed B lineage anchor in addition to combined marker evidence.
b_anchor_ix=genes.get_indexer(['CD79A','MS4A1','CD19']);b_anchor=(raw[:,b_anchor_ix[b_anchor_ix>=0]].getnnz(1)>0);major_ok[major.prediction.eq('B')]&=b_anchor[major.prediction.eq('B')]
result=a.obs.copy();result['qc_flags']=result['qc_flags'].astype(object);result['individual_prediction']=hdpred['subtype'].prediction.to_numpy();result['individual_prediction_score']=hdpred['subtype'].score.to_numpy();result['major_prediction']=major.prediction.to_numpy();result['major_prediction_score']=major.score.to_numpy();result['major_prediction_margin']=major.margin.to_numpy();result['hierarchical_prediction']=fine.prediction.to_numpy();result['prediction_score']=fine.score.to_numpy();result['prediction_margin']=fine.margin.to_numpy();result['prediction_margin_scale']='top1_minus_top2_logit';result['second_candidate']=fine.second_candidate.to_numpy();result['second_score']=fine.second_score.to_numpy();result['third_candidate']=fine.third_candidate.to_numpy();result['third_score']=fine.third_score.to_numpy();result['major_label']='Unknown';result['subtype_label']='Unknown';result['annotation_status']='Unknown';result['unresolved_reason']='major_calibration_or_observed_marker_evidence_insufficient';result['major_marker_genes_detected']=major_marker;result['subtype_marker_genes_detected']=0;result['label_evidence']='insufficient_observed_expression';result['reference_version']=config['reference_version'];result['annotation_version']='frozen_subtypes_v1';result['malignant_status']='unknown_not_inferred'
for i,(m,f) in enumerate(zip(major.prediction,fine.prediction)):
 if not major_ok[i]:continue
 result.iloc[i,result.columns.get_loc('major_label')]=m
 result.iloc[i,result.columns.get_loc('subtype_label')]='unresolved_'+m.replace('/','_')
 result.iloc[i,result.columns.get_loc('annotation_status')]='major_only'
 n=int(sg.get(f,np.zeros(a.n_obs,dtype='int16'))[i]);result.iloc[i,result.columns.get_loc('subtype_marker_genes_detected')]=n
 reasons=[]
 if dictionary.get(f)!=m:reasons.append('major_subtype_model_disagreement')
 if not fine_ok[i]:reasons.append('subtype_depth_calibration_insufficient')
 if n<2:reasons.append('fewer_than_two_observed_subtype_markers')
 if f in ANCHORS and ag[f][i]<1:reasons.append('identity_anchor_not_detected')
 if f in ['Epithelial','Plasma','Mast']:reasons=['reference_broad_identity_only']
 if not reasons:
  result.iloc[i,result.columns.get_loc('subtype_label')]=f;result.iloc[i,result.columns.get_loc('annotation_status')]='subtype_supported';result.iloc[i,result.columns.get_loc('unresolved_reason')]='';result.iloc[i,result.columns.get_loc('label_evidence')]='major_and_subtype_reference_calibration;multiple_observed_identity_markers'
 else:
  result.iloc[i,result.columns.get_loc('unresolved_reason')]=';'.join(reasons);result.iloc[i,result.columns.get_loc('label_evidence')]='major_reference_calibration;multiple_observed_major_markers'
# Conservative potential mixture flag; no deletion and no assertion of true doublet.
epi=ng['Epithelial']>=3;mixed=np.zeros(a.n_obs,bool);conflict=np.repeat('',a.n_obs).astype(object)
for other in ['B','Plasma','Myeloid','T/NK','Endothelial','Fibroblast','Mast']:
 second=(ng[other]>=3)&(marker_counts[other]>=5);ratio=np.minimum(marker_counts['Epithelial'],marker_counts[other])/np.maximum(np.maximum(marker_counts['Epithelial'],marker_counts[other]),1)
 hit=epi&second&(ratio>=.5);mixed|=hit;conflict[hit]='observed_epithelial_and_'+other+'_programs;unconfirmed_mixture_or_contamination'
result.loc[mixed,'annotation_status']='mixed_or_ambiguous';result.loc[mixed,'subtype_label']='mixed_or_ambiguous';result.loc[mixed,'unresolved_reason']=conflict[mixed];result.loc[mixed,'label_evidence']='conflicting_observed_major_marker_combinations';result.loc[mixed,'qc_flags']=result.loc[mixed,'qc_flags'].map(lambda x:('' if x=='none' else x+';')+'marker_conflict_unconfirmed')
result['annotation_review_status']='computational_annotation_pending_user_review';result.to_parquet(root/'metadata/computed_cell_labels.parquet',index=True,compression='zstd')
# Marker evidence stays observed, with per-label detection, rather than imputed values.
rows=[]
for kind,rules in [('major',MAJOR),('subtype',SUBTYPE)]:
 for label,markers in rules.items():
  membership=result['major_label' if kind=='major' else 'subtype_label'].eq(label).to_numpy();n=int(membership.sum())
  for gene in markers:
   ix=genes.get_indexer([gene])[0];vals=raw[:,ix] if ix>=0 else None
   rows.append(dict(kind=kind,label=label,gene=gene,present_in_model_gene_universe=ix>=0,n_assigned=n,observed_detection_n=int(vals[membership].getnnz()) if ix>=0 and n else 0,observed_count_sum=int(vals[membership].sum()) if ix>=0 and n else 0))
pd.DataFrame(rows).to_csv(root/'metadata/marker_evidence.tsv',sep='\t',index=False)
summary=dict(n_cells=a.n_obs,annotation_status=result.annotation_status.value_counts().to_dict(),major_labels=result.major_label.value_counts().to_dict(),subtype_supported=result.loc[result.annotation_status.eq('subtype_supported'),'subtype_label'].value_counts().to_dict(),targets_raw_shared_counts=targets,n_reference_heldout=held.n_obs,calibration='Per-class logit margin threshold, >=20 accepted heldout predictions and empirical precision >=0.8 on matching HD thinning depth; full-depth class must also be calibratable, but its margin threshold is not transferred; same holdout used for calibration, optimistic and not independent HD validation',low_information='Kept independently; use separate median-depth thinning calibration',gc_b_rule='B major with >=2 observed B markers including pan-B anchor, GC B model calibration, and >=2 of BCL6/AICDA/RGS13 observed; no cycling or spatial evidence substitution',mixed='Potential epithelial/other lineage marker conflicts retained, not asserted doublets',malignancy='not inferred')
(root/'metadata/annotation_summary.json').write_text(json.dumps(summary,indent=2));(root/'metadata/label_rules.json').write_text(json.dumps(dict(major_markers=MAJOR,subtype_markers=SUBTYPE,identity_anchors=ANCHORS,calibration=summary['calibration'],low_information=summary['low_information'],scope='marker rules derived from frozen v2 frameworks; no physical-neighbor voting, no score-driven cell exclusion',score_interpretation='CellTypist one-vs-rest sigmoid scores are not validated accuracy; margin is logit top1-top2'),indent=2))
print(json.dumps(summary,indent=2),flush=True)
