import sys,json,traceback
from pathlib import Path
OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(OUT/'code'))
from core import *
from scipy.stats import false_discovery_control

def verify():
 from report import summarize
 assert summarize(), 'SEARCH_INCOMPLETE'
 c=pd.read_csv(OUT/'combination_summary.tsv',sep='\t',dtype={'combo_id':str});p=pd.read_parquet(OUT/'primary_relations.parquet');a=pd.read_parquet(OUT/'auxiliary_relations.parquet');checks=[]
 def check(name,ok):assert ok,name;checks.append(dict(check=name,passed=True))
 check('99_57_42_complete',len(c)==99 and c.status.eq('COMPLETE').all() and c.contains_inhouse.sum()==57)
 check('FINAL_only',p.global_status.eq('FINAL').all());check('primary_sum4K',len(p)==int((c.K*4).sum()))
 for cid,pp in p.groupby('combo_id'):
  check(cid+'_planned4K',len(pp)==int(c.set_index('combo_id').loc[cid,'K']*4))
  q=false_discovery_control(pp.p.fillna(1).to_numpy(float));ok=pp.p.notna().to_numpy();check(cid+'_within_BH',np.allclose(q[ok],pp.q_within_BH.to_numpy()[ok],rtol=1e-13,atol=1e-15))
 for meth in ['bh','by']:
  q=false_discovery_control(p.p.fillna(1).to_numpy(float),method=meth);ok=p.p.notna().to_numpy();check('global_'+meth,np.allclose(q[ok],p['q_global_'+meth.upper()].to_numpy()[ok],rtol=1e-13,atol=1e-15));check('NA_preserved_'+meth,p.loc[~ok,'q_global_'+meth.upper()].isna().all())
 frozen=a[(a.branch=='FROZEN_LABEL')&(a.track=='SOURCE_P05')&(a.target_kind=='NICHE')&a.metric.isin(['O','D'])&(a.subset=='ALL_POLICY_UNITS')];check('frozen_3564',len(frozen)==3564)
 master=pd.read_parquet(OUT/'shared/master_anchors.parquet')
 for r in c.itertuples():
  f=OUT/'combinations'/r.combo_id
  for stage in ['NICHE_FROZEN','RING_DONE','SPOTTED_MAIN_DONE','FROZEN_DONE']:check(r.combo_id+'_'+stage,valid(f,stage))
  lab=np.load(f/'labels.npz');ref=np.load(f/'reference.npz');check(r.combo_id+'_label_axis',len(lab['anchor_index'])==len(lab['accepted']) and master.iloc[lab['anchor_index']].dataset_id.isin(r.datasets.split('|')).all());check(r.combo_id+'_label_range',((lab['accepted']>=-1)&(lab['accepted']<r.K)).all());check(r.combo_id+'_reference_unique_member',len(np.unique(ref['reference_index']))==len(ref['reference_index']) and np.isin(ref['reference_index'],lab['anchor_index']).all())
 review=json.loads((OUT/'review/COMPLETE.json').read_text());check('review_at_most3',review['selected_combinations']<=3)
 logical=bytes_used();allocated=sum(f.stat().st_blocks*512 for f in OUT.rglob('*') if f.is_file() and not f.is_symlink());cfg=json.loads((OUT/'config.json').read_text());check('logical_and_allocated_budget',max(logical,allocated)<=cfg['budget_bytes']);check('review_budget',bytes_used(OUT/'review')<=cfg['figure_budget_bytes'])
 js(OUT/'audit/final_validation.json',dict(status='PASSED',checks=checks,primary_M=len(p),completed=99,global_BH_hits=int(p.q_global_BH.lt(.05).sum()),global_BY_hits=int(p.q_global_BY.lt(.05).sum()),within_BH_hits=int(p.q_within_BH.lt(.05).sum()),inhouse_global_supported=int((p.contains_inhouse&p.inhouse_actual_support&p.q_global_BH.lt(.05)).sum()),logical_bytes=logical,allocated_bytes=allocated,code_hash=sha(__file__)))
 with open(OUT/'RUN_SUMMARY.md','a') as f:f.write('\n最终独立校验：PASSED。99 个组合产物哈希、planned family 数量、SciPy BH/BY 数值、标签和参考索引、图包数量及逻辑/实际分配磁盘预算均通过。详见 audit/final_validation.json。\n')
 print('FINAL_VALIDATION_PASSED',flush=True)
try:verify()
except Exception:
 error=traceback.format_exc();js(OUT/'audit/final_validation.json',dict(status='FAILED',error=error));r=json.loads((OUT/'RESUME.json').read_text());r.update(status='SEARCH_INCOMPLETE',final_validation='FAILED');js(OUT/'RESUME.json',r)
 with open(OUT/'RUN_SUMMARY.md','a') as f:f.write('\n最终交付校验未通过，不能宣称执行完整完成。详见 audit/final_validation.json。\n')
 raise
