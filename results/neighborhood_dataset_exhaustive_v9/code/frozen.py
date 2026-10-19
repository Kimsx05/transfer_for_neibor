from core import *
from stats import summarize_spotted

def run(cid):
 folder=OUT/'combinations'/cid;folder.mkdir(exist_ok=True)
 if valid(folder,'FROZEN_DONE'):return
 guard(cid+':frozen',8*1024**2);c=pd.read_csv(OUT/'combinations.tsv',sep='\t',dtype={'combo_id':str}).set_index('combo_id').loc[cid];ds=c.datasets.split('|');s=pd.read_parquet(V8/'07_synthesis/all_section_effects.parquet',filters=[('k','==',10),('analysis_domain','==','MAIN'),('dataset_id','in',ds)])
 u,r=summarize_spotted(s,9,'FROZEN_LABEL');assert len(r[(r.scope=='ALL_DATASETS')&(r.target_kind=='NICHE')&(r.track=='SOURCE_P05')&(r.subset=='ALL_POLICY_UNITS')&r.metric.isin(['O','D'])])==36
 table(folder/'frozen_units.parquet',u);table(folder/'frozen_tests.parquet',r);js(folder/'frozen_source.json',dict(K=9,labels=str(V3/'03_main/SAMPLE_R12_LEIDEN/labels.parquet'),section_effects=str(V8/'07_synthesis/all_section_effects.parquet'),source='Original v8 E0; no source-domain equivalence claimed for REFIT',full_fixed_targets=list(range(9))))
 mark(folder,'FROZEN_DONE',[folder/'frozen_units.parquet',folder/'frozen_tests.parquet',folder/'frozen_source.json'])
if __name__=='__main__':run(sys.argv[1])
