from common import *
import importlib.util,anndata as ad
from scipy.stats import rankdata,spearmanr
sp=importlib.util.spec_from_file_location('measure',R/'code/05_measure.py');m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);rows=[]
for unit in ['cell','bin8','bin16']:
 if unit=='cell':
  a=ad.read_h5ad(A/'segmentation_assignment_v2/Cervilla/cells.h5ad');k=a.obs.qc_keep_analysis.to_numpy();x=sparse.csr_matrix(a.X[k]);ids=a.obs_names[k].to_numpy();sy=a.var.gene_symbol.to_numpy();del a
 elif unit=='bin8':x,d,_,sy,_=read8('Cervilla');ids=d.bin_id.to_numpy()
 else:x=sparse.load_npz(R/'models/Cervilla_bin16_counts.npz');d=pd.read_parquet(R/'tables/Cervilla_bin16_measurement.parquet');ids=d.bin_id.to_numpy();sy=pd.read_csv(R/'tables/hd_genes.tsv',sep='\t').symbol.to_numpy()
 ti=np.flatnonzero(np.isin(sy,TOP));out=m.auc_historical(x.indptr,x.indices,x.data,x.shape[1],ti);file=R/'tables'/f'Cervilla_{unit}_measured.parquet';d=pd.read_parquet(file);assert np.array_equal(ids,d.bin_id);d['mp10_historical_auc']=out;d.to_parquet(file,index=False)
 v=pd.read_csv(R/'tables'/f'{unit}_score_seed_validation.tsv.gz',sep='\t');ix=pd.Index(ids).get_indexer(v[v.seed==123].id);small=x[ix];K=int(np.ceil(.05*x.shape[1]));den=sum(K-i for i in range(1,21));ranks=rankdata(-small.toarray(),method='average',axis=1)[:,ti];expected=np.maximum(0,K-ranks).sum(1)/den;expected[np.asarray(small.sum(1)).ravel()==0]=np.nan;err=np.nanmax(np.abs(expected-out[ix]));assert err<1e-12
 v['auc_historical']=v.id.map(pd.Series(out,index=ids));v.to_csv(R/'tables'/f'{unit}_score_seed_validation.tsv.gz',sep='\t',index=False);rng=np.random.default_rng(123);z=small.copy();z.data=rng.binomial(z.data.astype(int),.5);z.eliminate_zeros();half=m.auc_historical(z.indptr,z.indices,z.data,z.shape[1],ti);put(pd.DataFrame(dict(id=ids[ix],historical_full=out[ix],historical_half=half)) ,f'{unit}_historical_depth_validation.tsv')
 rows.append(dict(unit=unit,algorithm='exact historical average-tie rank formula; raw all genes',max_abs_error_dense_scipy=err,n_checked=len(ix),half_depth_spearman=spearmanr(out[ix],half,nan_policy='omit').statistic,source='neighborhood_mp10_calibration_v5/code/common.py::auc_parts; reused by mp10_subtype_screen_round1_v2',version='HD_primary_historical_mean_rank_v1'));print('HISTORICAL_VALID',unit,err,flush=True)
put(pd.DataFrame(rows),'historical_score_validation.tsv');cfg=yaml.safe_load((R/'frozen_config.yaml').read_text());cfg['MP10']['primary_tie_method']='historical average rank, not native random zero ties';cfg['MP10']['implementation_source']='neighborhood_mp10_calibration_v5/code/common.py::auc_parts';cfg['MP10']['implementation_correction_before_association']=True;(R/'frozen_config.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False,allow_unicode=True))
