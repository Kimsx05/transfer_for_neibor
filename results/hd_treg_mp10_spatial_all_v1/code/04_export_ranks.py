from common import *
from scipy.io import mmwrite
cfg=yaml.safe_load((R/'frozen_config.yaml').read_text());genes=sorted(set(sum(cfg['signatures'].values(),[])));M=cfg['UCell_maxRank']
for s in SECTIONS:
 for kind in ['T32','T16','Epi64']:
  d=pd.read_parquet(R/s/f'{kind}_measure.parquet',columns=['region_id']+['rank_'+g for g in genes]);v=d.drop(columns='region_id').to_numpy(copy=True).T;v[v>=M]=0
  stem=R/'temporary'/f'{s}_{kind}';mmwrite(str(stem)+'_ucell_ranks.mtx',sparse.coo_matrix(v));pd.Series(genes).to_csv(str(stem)+'_ucell_genes.txt',index=False,header=False);d.region_id.to_csv(str(stem)+'_ucell_ids.txt',index=False,header=False)
