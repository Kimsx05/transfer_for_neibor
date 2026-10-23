"""Same reference models and preprocessing; per-section depth calibration; candidates, not final labels."""
from common import *
os.environ['CELLTYPIST_FOLDER']=str(ROOT/'models/celltypist_cache')
import numpy as np,pandas as pd,anndata as ad,celltypist,sys
from scipy import sparse
from marker_rules import MAJOR,SUBTYPE,ANCHORS
name=sys.argv[1];meta,base,out=section(name);np.random.seed(123);rng=np.random.default_rng(123)
config=json.loads((ROOT/'metadata/classifier_config.json').read_text());genes=pd.Index(config['gene_list']);dictionary=pd.read_csv(ROOT/'metadata/reference_label_dictionary.tsv',sep='\t').set_index('reference_celltype_v2').major_label.to_dict()
a=ad.read_h5ad(out/'cells.h5ad');held=ad.read_h5ad(ROOT/'metadata/reference_holdout_temporary.h5ad');gm=pd.read_csv(ROOT/f'metadata/gene_mapping_{name}.tsv',sep='\t');valid=gm.reference_symbol.notna();trans=sparse.csr_matrix((np.ones(valid.sum(),dtype=np.int32),(gm.loc[valid,'hd_feature_index'],genes.get_indexer(gm.loc[valid,'reference_symbol']))),shape=(a.n_vars,len(genes)));raw=(a.X@trans).tocsr();raw.sum_duplicates();raw.eliminate_zeros();a.var['reference_symbol']=gm.reference_symbol.fillna('').to_numpy()
def norm(x):
 x=x.astype('float32').tocsr();x=sparse.diags(1e4/np.maximum(np.asarray(x.sum(1)).ravel(),1)).dot(x).tocsr();x.data=np.log1p(x.data);return x

def thin(x,target):
 x=x.copy().tocsr();p=np.minimum(1,target/np.maximum(np.asarray(x.sum(1)).ravel(),1));x.data=rng.binomial(x.data.astype('int64'),np.repeat(p,np.diff(x.indptr))).astype('float32');x.eliminate_zeros();return x

def predict(x,model):
 ix=genes.get_indexer(model.classifier.features);assert (ix>=0).all();scores=[]
 for start in range(0,len(x.indptr)-1,5000):
  q=model.scaler.transform(x[start:start+5000,ix]).tocsr();q.data=np.minimum(q.data,10);d,p,_=model.predict_labels_and_prob(q);scores.append(d.astype('float32'))
 return np.vstack(scores)

def top(d,classes):
 order=np.argsort(-d,axis=1,kind='stable');i=np.arange(len(d));return classes[order[:,0]],d[i,order[:,0]]-d[i,order[:,1]] if d.shape[1]>1 else np.full(len(d),np.nan),1/(1+np.exp(-np.clip(d[i,order[:,0]],-60,60))),classes[order[:,min(1,d.shape[1]-1)]]

def calibrate(d,classes,true,head,condition):
 pred,margin,_,_=top(d,classes);rows=[]
 for lab in classes:
  hit=pred==lab;m=margin[hit];correct=true[hit]==lab;best=None
  for th in np.unique(np.quantile(m,np.linspace(0,.95,20))) if len(m) and len(classes)>1 else []:
   yes=m>=th;n=yes.sum();precision=correct[yes].mean() if n else 0
   if n>=20 and precision>=.8:best=(float(th),int(n),float(precision));break
  rows.append(dict(head=head,condition=condition,label=lab,true_n=int((true==lab).sum()),predicted_n=int(hit.sum()),calibrated=best is not None,threshold=best[0] if best else np.nan,accepted_n=best[1] if best else 0,accepted_precision=best[2] if best else np.nan))
 return rows
low=a.obs.low_information.to_numpy();depth=np.asarray(raw.sum(1)).ravel();targets={k:max(10,int(np.median(depth[hit]))) for k,hit in [('ordinary',~low),('low_information',low)] if hit.any()};targets.setdefault('low_information',10)
views={'full':held.X,**{k:thin(held.X,v) for k,v in targets.items()}};lograw=norm(raw);hd={};cal=[];models={};hc={}
for head in ['major','subtype','b_plasma']:
 model=celltypist.models.Model.load(str(ROOT/f'models/celltypist_{head}.pkl'));models[head]=model;classes=np.asarray(model.classifier.classes_).astype(str);hc[head]=classes
 print('PREDICT',name,head,a.n_obs,flush=True);hd[head]=predict(lograw,model)
 for cond,x in views.items():
  h=held.obs.major_label.isin(['B','Plasma']).to_numpy() if head=='b_plasma' else np.ones(held.n_obs,bool);d=predict(norm(x[h]),model);maj=held.obs.major_label.to_numpy()[h].astype(str);fine=held.obs.reference_celltype_v2.to_numpy()[h].astype(str)
  cal+=calibrate(d,classes,maj if head=='major' else fine,head,cond)
  if head!='major':
   for major in sorted(set(dictionary.values())):
    ix=np.flatnonzero(np.array([dictionary[z] for z in classes])==major);hh=maj==major
    if len(ix)>=2 and hh.any():cal+=calibrate(d[hh][:,ix],classes[ix],fine[hh],head+'__'+major,cond)
cal=pd.DataFrame(cal);cal.to_csv(out/'metadata/depth_calibration.tsv',sep='\t',index=False);lookup={(z['head'],z['condition'],z['label']):z for z in cal.to_dict('records')}
def accept(pred,margin,heads):
 ok=np.zeros(len(pred),bool)
 for i,(label,head) in enumerate(zip(pred,heads)):
  full=lookup.get((head,'full',label),{});q=lookup.get((head,'low_information' if low[i] else 'ordinary',label),{})
  ok[i]=full.get('calibrated',False) and q.get('calibrated',False) and margin[i]>=q.get('threshold',np.inf)
 return ok
mp,mm,ms,msecond=top(hd['major'],hc['major']);sp,sm,ss,ssecond=top(hd['subtype'],hc['subtype']);fine=sp.copy();fm=sm.copy();fs=ss.copy();second=ssecond.copy();heads=np.repeat('subtype',a.n_obs).astype(object)
for major in hc['major']:
 hit=mp==major;head='b_plasma' if major in ['B','Plasma'] else 'subtype';ix=np.flatnonzero(np.array([dictionary[z] for z in hc[head]])==major)
 if not hit.any() or not len(ix):continue
 pred,margin,score,sec=top(hd[head][hit][:,ix],hc[head][ix]);fine[hit]=pred;fm[hit]=margin;fs[hit]=score;second[hit]=sec;heads[hit]=head+'__'+major if len(ix)>1 else head
major_ok=accept(mp,mm,np.repeat('major',a.n_obs));fine_ok=accept(fine,fm,heads)
o=a.obs;o['predicted_major']=mp;o['major_score']=ms;o['major_margin']=mm;o['major_calibrated']=major_ok;o['predicted_subtype_global']=sp;o['predicted_subtype']=fine;o['subtype_score']=fs;o['subtype_margin']=fm;o['subtype_second_candidate']=second;o['subtype_calibrated']=fine_ok;o['malignant_status']='unknown_not_inferred';o['annotation_version']='cell_annotation_v2_review'
# Marker presence from measured counts, with missing features recorded; no imputation.
evidence={};present={}
for label,markers in {**MAJOR,**SUBTYPE}.items():
 ix=genes.get_indexer(markers);ix=ix[ix>=0];q=raw[:,ix];key=label.replace('/','_').replace(' ','_');o['observed_markers__'+key]=q.getnnz(1).astype('int16');o['observed_marker_counts__'+key]=np.asarray(q.sum(1)).ravel().astype('int32');evidence[label]=(q.getnnz(1),np.asarray(q.sum(1)).ravel());present[label]=[g for g in markers if g in genes]
mixed=np.zeros(a.n_obs,bool);reason=np.repeat('',a.n_obs).astype(object);en,ec=evidence['Epithelial']
for lab in ['B','Plasma','Myeloid','T/NK','Endothelial','Fibroblast','Mast']:
 n,c=evidence[lab];ratio=np.minimum(ec,c)/np.maximum(np.maximum(ec,c),1);hit=(en>=3)&(n>=3)&(c>=5)&(ratio>=.5);mixed|=hit;reason[hit]='epithelial_and_'+lab+'_programs;possible_mixture_or_contamination'
o['mixed_or_ambiguous']=mixed;o['mixed_reason']=reason;o['model_eligible']=o.model_eligible&~mixed;o['model_role']=np.where(o.model_eligible,'structure_core','prediction_only');o['major_label']='Unknown';o['final_label']='Unknown';o['annotation_status']='pending_expression_cluster_review';o['lineage_unresolved']=True;o['subtype_unresolved']=True
# Master full gene raw counts unchanged; only obs/var metadata have been added.
a.write_h5ad(out/'cells.h5ad',compression='gzip');o.to_parquet(out/'metadata/reference_predictions.parquet',compression='zstd')
jwrite(out/'metadata/prediction_summary.json',dict(section=name,n_cells=a.n_obs,targets_actual_shared_counts=targets,major_candidates=pd.Series(mp).value_counts().to_dict(),subtype_candidates=pd.Series(fine).value_counts().to_dict(),major_calibrated=int(major_ok.sum()),subtype_calibrated=int(fine_ok.sum()),mixed_candidates=int(mixed.sum()),structure_core=int(o.model_eligible.sum()),reference_version=config['reference_version'],limitations='Same v1 heldout source groups reused for calibration, not independent HD truth; class conditional heads restrict global reference logits within predicted lineage; not spatial voting',marker_genes_available=present))
print('REFERENCE_PREDICTIONS_DONE',name,flush=True)
