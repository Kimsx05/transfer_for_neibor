from common import *
import anndata as ad,pandas as pd,numpy as np,h5py
config=json.loads((ROOT/'metadata/classifier_config.json').read_text());genes=pd.Index(config['gene_list']);manifest=json.loads((ROOT/'metadata/reference_manifest.json').read_text());sel=pd.read_parquet(ROOT/'metadata/reference_holdout_selection.parquet')
a=ad.read_h5ad(manifest['object_path'],backed='r');ri=a.obs_names.get_indexer(sel.index);vi=a.var_names.get_indexer(genes);assert (ri>=0).all() and (vi>=0).all();order=np.argsort(ri);b=a[ri[order],vi].to_memory();b.obs=sel.iloc[order].copy();b.X=b.X.tocsr().astype('int32');b.write_h5ad(ROOT/'metadata/reference_holdout_temporary.h5ad',compression='gzip');a.file.close()
checks=[]
for name in ['Cervilla','D1_Johnson']:
 _,base,out=section(name);gm=pd.read_csv(ROOT/f'metadata/gene_mapping_{name}.tsv',sep='\t')
 with h5py.File(base/'filtered_feature_bc_matrix.h5') as f:ids=f['matrix/features/id'].asstr()[:];symbols=f['matrix/features/name'].asstr()[:]
 assert np.array_equal(ids,gm.feature_id) and np.array_equal(symbols,gm.source_symbol);actual=set(gm.reference_symbol.dropna());assert actual==set(genes)
 checks.append(dict(section=name,features=len(ids),mapped_genes=len(actual),exact_feature_order_matches_v1_mapping=True,compatible_models=['major','subtype','b_plasma'],normalization='sum duplicate mapped genes, restrict to same14327, normalize_total1e4 then log1p, same scaler and max10 clipping',reference_labels=32))
jwrite(ROOT/'metadata/reference_compatibility.json',checks);print('REFERENCE_READY',len(sel),checks,flush=True)
