from common import *
import numpy as np,pandas as pd,anndata as ad,h5py
from scipy import sparse
rng=np.random.default_rng(123);config=json.loads((V2/'metadata/classifier_config.json').read_text());genes=pd.Index(config['gene_list']);manifest=json.loads((V2/'metadata/reference_manifest.json').read_text());ref=ad.read_h5ad(manifest['object_path'],backed='r');dictionary=pd.read_csv(V2/'metadata/reference_label_dictionary.tsv',sep='\t').set_index('reference_celltype_v2').major_label.to_dict();dictionary.update({x:'Mural' for x in ['SMC','Pericyte']});selection=[]
for lab,q in ref.obs.groupby('reference_celltype_v2',observed=True):
 queues=[]
 for _,g in q.groupby(['dataset','orig.ident'],observed=True):queues.append(list(rng.permutation(ref.obs_names.get_indexer(g.index))))
 chosen=[]
 while len(chosen)<min(400,len(q)):
  for k in rng.permutation(len(queues)):
   if queues[k] and len(chosen)<min(400,len(q)):chosen.append(queues[k].pop())
 selection.extend(chosen)
ix=np.sort(selection);b=ref[ix,ref.var_names.get_indexer(genes)].to_memory();assert b.obs.reference_celltype_v2.nunique()==32;b.obs['major_v21']=b.obs.reference_celltype_v2.map(dictionary).astype(str);b.obs[['reference_celltype_v2','major_v21','dataset','orig.ident']].to_csv(ROOT/'metadata/SingleR_reference_selection.tsv',sep='\t',index_label='cell_id');ref.file.close()
def write_sparse(x,path,ids):
 x=x.tocsr();x.sum_duplicates();x.eliminate_zeros()
 with h5py.File(path,'w') as f:
  for k,v in [('data',x.data.astype('float64')),('indices',x.indices.astype('int32')),('indptr',x.indptr.astype('int32')),('shape',np.array(x.shape,dtype='int32'))]:f.create_dataset(k,data=v,compression='gzip')
  f.create_dataset('genes',data=genes.to_numpy().astype(object),dtype=h5py.string_dtype());f.create_dataset('cells',data=np.array(ids,dtype=object),dtype=h5py.string_dtype())
write_sparse(b.X,ROOT/'temporary/reference_shared_counts.h5',b.obs_names);meta=[]
for name,expected in [('Cervilla',97927),('D1_Johnson',210622)]:
 a=read_section(name);assert a.n_obs==expected and a.obs_names.is_unique;gm=pd.read_csv(V2/f'metadata/gene_mapping_{name}.tsv',sep='\t');valid=gm.reference_symbol.notna();t=sparse.csr_matrix((np.ones(valid.sum(),dtype='int32'),(gm.loc[valid,'hd_feature_index'],genes.get_indexer(gm.loc[valid,'reference_symbol']))),shape=(a.n_vars,len(genes)));x=(a.X@t).tocsr();write_sparse(x,ROOT/f'temporary/{name}_shared_counts.h5',a.obs_names)
 keep=['cell_id','section_id','n_counts','n_genes','low_information','mixed_or_ambiguous','mixed_reason','qc_flags','morphology_review_flag','model_role','model_eligible','scvi_training_core','umap_mapped','predicted_major','predicted_subtype','predicted_subtype_global','major_score','major_margin','major_calibrated','subtype_score','subtype_margin','subtype_calibrated','major_label','final_label','annotation_status','expression_cluster','lineage_cluster','x_um','y_um']
 o=a.obs[keep].copy();o.to_parquet(ROOT/f'metadata/{name}_v2_obs.parquet',compression='zstd');meta.append(dict(section=name,cells=a.n_obs,features=a.n_vars,shared_genes=len(genes),source=str(V2/name/'cells.h5ad'),source_stat=dict(size=(V2/name/'cells.h5ad').stat().st_size,mtime_ns=(V2/name/'cells.h5ad').stat().st_mtime_ns)))
 print('PREPARED',name,a.shape,flush=True)
js(ROOT/'metadata/inputs.json',meta);js(ROOT/'metadata/reference.json',dict(source=manifest['object_path'],reference_id=manifest['reference_version'],labels=32,selection='Up to400 per frozen subtype, round-robin dataset/specimen strata; all375 vCAF retained; seed123. Specimens are not assumed independent patients.',n_reference=len(b),shared_genes=len(genes),major_mapping=dictionary,reference_label_counts=b.obs.reference_celltype_v2.value_counts().to_dict()));print('REFERENCE',len(b),flush=True)
