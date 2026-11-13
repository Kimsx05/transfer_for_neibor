import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
import sys,pathlib,json,hashlib
import pandas as pd
from distributed import Client,LocalCluster
from arboreto.core import create_graph,SGBM_KWARGS,EARLY_STOP_WINDOW_LENGTH

def main():
 lin=sys.argv[1];B=pathlib.Path('results/mp10_immune_endothelial_mechanisms_v1');O=B/'objects/scenic';T=B/'tables';ex=pd.read_csv(O/f'{lin}_training.tsv.gz',sep='\t',index_col=0).T;targets=pd.read_csv(T/f'scenic_{lin}_seed456_targets.tsv',sep='\t').gene.tolist();targets=sorted(set(targets)&set(ex.columns));tfs=pathlib.Path('resources/mp10_regulatory_priors/human_tfs.txt').read_text().splitlines();out=O/f'{lin}_seed456_keytarget_adjacencies.tsv.gz'
 def sha(p):
  h=hashlib.sha256()
  with open(p,'rb') as f:
   for b in iter(lambda:f.read(1048576),b''):h.update(b)
  return h.hexdigest()
 signature=hashlib.sha256(json.dumps(dict(training=sha(O/f'{lin}_training.tsv.gz'),targets=targets,TFs=tfs,seed=456,main_adjacencies=sha(O/f'{lin}_seed123_adjacencies.tsv'),regulons=sha(T/f'scenic_{lin}_regulon_edges.tsv.gz')),sort_keys=True).encode()).hexdigest()
 cachemeta=O/f'{lin}_seed456_metadata.json'
 cache_valid=False
 if cachemeta.exists() and out.exists():
  old=json.loads(cachemeta.read_text());cache_valid=old.get('input_hash')==signature and old.get('adjacency_sha256')==sha(out)
 workers=min(2,len(os.sched_getaffinity(0)))
 if cache_valid:
  adj=pd.read_csv(out,sep='\t');workers=old.get('num_workers',workers);print('Verified second-seed GRN cache',lin)
 else:
  cluster=LocalCluster(n_workers=workers,threads_per_worker=1,processes=True,dashboard_address=None,memory_limit='16GB');client=Client(cluster)
  try:
   graph=create_graph(ex.values,ex.columns.tolist(),tfs,regressor_type='GBM',regressor_kwargs=SGBM_KWARGS,client=client,target_genes=targets,early_stop_window_length=EARLY_STOP_WINDOW_LENGTH,seed=456);adj=client.compute(graph,sync=True).sort_values('importance',ascending=False);adj.to_csv(out,sep='\t',index=False,compression='gzip')
  finally:client.close();cluster.close()
 mainadj=pd.read_csv(O/f'{lin}_seed123_adjacencies.tsv',sep='\t');key=pd.read_csv(T/f'scenic_{lin}_key_TF_selection.tsv',sep='\t');edges=pd.read_csv(T/f'scenic_{lin}_regulon_edges.tsv.gz',sep='\t');edges=edges[edges.TF.isin(key.source)];z=edges.merge(mainadj.rename(columns={'importance':'importance_seed123'}),on=['TF','target'],how='left').merge(adj.rename(columns={'importance':'importance_seed456'}),on=['TF','target'],how='left');z['recovered_seed456']=z.importance_seed456.notna()&(z.importance_seed456>0);z['self_target']=z.TF.eq(z.target);z['seed_comparison_evaluable']=(~z.self_target)&z.importance_seed123.notna()&z.importance_seed123.gt(0);z['stability_status']=['not_evaluable_self_target_signature_member' if a else 'not_evaluable_no_seed123_GRN_edge' if not b else 'recovered' if c else 'not_recovered' for a,b,c in zip(z.self_target,z.seed_comparison_evaluable,z.recovered_seed456)];z['second_seed_scope']='key targets, all original TF predictors, unchanged expression';z.to_csv(T/f'scenic_{lin}_connection_seed_stability.tsv.gz',sep='\t',index=False)
 meta=dict(edge_interpretation='self-target signature members excluded from GRN recovery denominator',lineage=lin,seed=456,num_workers=workers,input_hash=signature,adjacency_sha256=sha(out),n_targets=len(targets),key_TFs=sorted(set(key.source)),rule='top5 CollecTRI q ranks per primary model, independent of significance; union of their motif-pruned targets and named candidates/Top20',target_hash=hashlib.sha256('\n'.join(targets).encode()).hexdigest());(O/f'{lin}_seed456_metadata.json').write_text(json.dumps(meta,indent=2));print('SECOND SEED COMPLETE',lin,flush=True)
if __name__=='__main__':main()
