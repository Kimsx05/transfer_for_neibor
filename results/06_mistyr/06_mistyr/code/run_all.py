import bootstrap
from bootstrap import ROOT
import json,hashlib,subprocess,os,sys
from pathlib import Path
import pandas as pd
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
p=json.loads((ROOT/'parameters.json').read_text()); man=pd.read_csv(ROOT/'section_input_summary.tsv',sep='\t')
man=pd.concat([man[man.section_id==p['pilot_section']],man[man.section_id!=p['pilot_section']]])
env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='')
events=[]
for i,row in enumerate(man.itertuples()):
 folder=ROOT/'by_dataset'/row.dataset/row.section_id
 sig=hashlib.sha256(''.join(sha(x) for x in [folder/'model_input.tsv',ROOT/'code/model_section.R',ROOT/'parameters.json',ROOT/'runtime/Rlib/mistyR/DESCRIPTION']).encode()).hexdigest()
 done=folder/'completion.json'
 if done.exists() and json.loads(done.read_text()).get('signature')==sig:
  print('CACHE_MATCH',row.section_id,flush=True);continue
 with open(ROOT/'logs'/f'{row.section_id}.log','w') as log:
  out=subprocess.run(['Rscript',str(ROOT/'code/model_section.R'),str(ROOT),str(folder)],env=env,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
 st=pd.read_csv(folder/'section_status.tsv',sep='\t') if (folder/'section_status.tsv').exists() else None
 ok=out.returncode==0 and st is not None and not (st.status=='BLOCKED').any()
 if i==0:
  assert ok, 'Pilot failed; remaining sections not started; inspect pilot log'
  au=pd.read_csv(folder/'predictor_exclusion_audit.tsv',sep='\t'); assert au.target_excluded.all() and au.deterministic_copy_excluded.all()
  assert not (au.target==au.predictor).any()
  cv=pd.read_csv(folder/'spatial_fold_audit.tsv',sep='\t');assert cv.fold_overlap.eq(0).all()
  pred=pd.read_csv(folder/'mp10_spatial_predictions.tsv',sep='\t')
  keys=[set(map(tuple,g[['section_id','barcode','spatial_fold']].to_numpy())) for _,g in pred.groupby('model')]
  assert len(keys)==2 and keys[0]==keys[1]
  (ROOT/'pilot_validation.json').write_text(json.dumps(dict(section_id=row.section_id,status='PASS',target_self_excluded=True,deterministic_copy_excluded=True,intra_only=True,same_complete_case_folds=True,train_test_overlap=0,signature=sig),indent=2))
 done.write_text(json.dumps(dict(signature=sig,exit_code=out.returncode,status='PASS_WITH_WARNING' if ok else 'BLOCKED'),indent=2))
 events.append(dict(section_id=row.section_id,dataset=row.dataset,exit_code=out.returncode,status='PASS_WITH_WARNING' if ok else 'BLOCKED'))
 pd.DataFrame(events).to_csv(ROOT/'execution_events.tsv',sep='\t',index=False)
 print('FINISHED',row.section_id,'exit',out.returncode,flush=True)
print('ALL_SECTIONS_ATTEMPTED',len(man),flush=True)
