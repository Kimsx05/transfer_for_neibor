from common import *
from scipy.io import mmread
from scipy.stats import rankdata
cfg=yaml.safe_load((R/'frozen_config.yaml').read_text());sig=cfg['signatures'];checks=[]
for s in SECTIONS:
 for kind in ['T32','T16','Epi64']:
  d=pd.read_parquet(R/s/f'{kind}_regions.parquet').set_index('region_id');u=pd.read_csv(R/'tables'/f'{s}_{kind}_native_UCell.tsv.gz',sep='\t').set_index('region_id');assert d.index.equals(u.index)
  err=float(np.max(np.abs(d[list(sig)].to_numpy()-u[list(sig)].to_numpy())));assert err<1e-7,(s,kind,err);
  ranks=pd.read_parquet(R/s/f'{kind}_measure.parquet').set_index('region_id').filter(like='rank_').rename(columns=lambda c:c[5:]);computed=np.column_stack([ucell(ranks,g,cfg['UCell_maxRank'])for g in sig.values()]);assert np.max(np.abs(computed-u[list(sig)].to_numpy()))<1e-12
  d[list(sig)]=u[list(sig)];d['MP10']=mp10(ranks);d.reset_index().to_parquet(R/s/f'{kind}_regions.parquet',index=False)
  checks.append(dict(section=s,branch=kind,check='all_region_native_UCell_equivalence',max_error=err,n=len(d),passed=True))
q=pd.read_csv(R/'tables/UCell_native_raw_check.tsv',sep='\t').set_index('region_id');d=pd.read_parquet(R/'Cervilla/T32_regions.parquet').set_index('region_id').loc[q.index];err=np.max(np.abs(d[list(sig)].to_numpy()-q[list(sig)].to_numpy()));assert err<1e-12;checks.append(dict(section='Cervilla',branch='T32',check='public_UCell_full_raw_matrix_equivalence',max_error=err,n=len(q),passed=True))
x=mmread(R/'temporary/native_UCell_check.mtx').T.toarray();genes=np.loadtxt(R/'temporary/native_UCell_genes.txt',dtype=str);ti=np.flatnonzero(np.isin(genes,TOP));ranks=rankdata(-x,axis=1,method='average');expected=np.maximum(905-ranks[:,ti],0).sum(1)/sum(905-np.arange(1,21));err=float(np.max(np.abs(expected-d.MP10)));assert err<1e-12;checks.append(dict(section='Cervilla',branch='T32',check='historical_MP10_dense_average_rank_exact',max_error=err,n=len(d),passed=True))
put(pd.DataFrame(checks),'native_score_validation.tsv');print('NATIVE_CHECKS_OK',flush=True)
