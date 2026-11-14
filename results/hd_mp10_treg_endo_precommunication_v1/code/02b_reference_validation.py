from common import *
d=pd.read_csv(R/'tables/reference_holdout_predictions.tsv.gz',sep='\t');rows=[]
for dep,z in d.groupby('target_depth'):
 for t in sorted(z.truth.unique()):
  truth=z.truth.eq(t);pred=z.pruned.eq(t);tp=int((truth&pred).sum());rows.append(dict(depth=dep,type=t,truth_n=int(truth.sum()),predicted_n=int(pred.sum()),precision=tp/pred.sum()if pred.sum()else np.nan,recall=tp/truth.sum()))
out=pd.DataFrame(rows);old=R/'tables/annotation_validation.tsv'
if old.exists():
 prev=pd.read_csv(old,sep='\t');assert np.allclose(prev[['precision','recall']],out[['precision','recall']],equal_nan=True),'Frozen holdout results differ; do not silently reuse target rules'
put(out,'annotation_validation.tsv')
