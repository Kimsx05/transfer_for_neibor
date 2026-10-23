from common import *
from expression_evidence import load_all,normalize,group_evidence
import numpy as np,pandas as pd,scanpy as sc,anndata as ad
from scipy import sparse
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
import igraph as ig,leidenalg

def tag(s):return s.replace('/','_')
x,o,genes,z,u,physical=load_all();ob=pd.read_parquet(ROOT/'metadata/major_proposals.parquet');assert ob.index.equals(o.index);o=ob.copy();o['lineage_cluster_v21']='not_assigned_to_lineage';pcaout=np.full((len(o),30),np.nan,dtype='float32');info={}
for lineage,ii in o.loc[~o.proposed_major.eq('Unknown')].groupby('proposed_major').indices.items():
 # groupby indices above are positional in subset, so use label names for full-array positions.
 ix=np.flatnonzero(o.proposed_major.eq(lineage));core=ix[(o.n_genes.to_numpy()[ix]>30)&(~o.mixed_observed_programs.to_numpy()[ix])];small=np.setdiff1d(ix,core);prefix='lineage_'+tag(lineage)
 print('LINEAGE',lineage,len(ix),'informative',len(core),flush=True)
 if len(core)<40:
  o.loc[o.index[ix],'lineage_cluster_v21']=lineage+'::small_group';info[lineage]={'cells':len(ix),'informative':len(core),'method':'small_group_marker_review_no_forced_graph'}
 else:
  a=ad.AnnData(normalize(x[core]),obs=o.iloc[core][['section_id']].copy(),var=pd.DataFrame(index=genes));sc.pp.highly_variable_genes(a,flavor='seurat',n_top_genes=min(2000,a.n_vars),batch_key='section_id');hv=a.var.highly_variable.to_numpy()&~genes.str.match(r'^(MT-|RPL|RPS)').to_numpy() if hasattr(genes.str.match(r'^(MT-|RPL|RPS)'),'to_numpy') else a.var.highly_variable.to_numpy()&~genes.str.match(r'^(MT-|RPL|RPS)')
  # Preserve lineage marker genes in the feature space when detected, even if global HVG selection missed them.
  from programs import MAJOR
  from expression_evidence import SUBTYPE
  mapping=json.loads((ROOT/'metadata/reference.json').read_text())['major_mapping'];known=MAJOR.get(lineage,[])+sum([gs for lab,gs in SUBTYPE.items() if mapping.get(lab)==lineage],[]);hv|=genes.isin(known)&(np.asarray((a.X>0).sum(0)).ravel()>=3)
  a=a[:,hv].copy();sc.pp.scale(a,zero_center=False,max_value=10);ncomp=min(30,len(core)-2,a.n_vars-2);sc.tl.pca(a,n_comps=ncomp,svd_solver='arpack',random_state=123);pc=a.obsm['X_pca'];pcaout[core,:ncomp]=pc
  # kNN in measured expression PCA; no physical-coordinate graph or source-mixing requirement.
  nn=NearestNeighbors(n_neighbors=min(16,len(core)),n_jobs=8).fit(pc);dist,ind=nn.kneighbors(pc);row=np.repeat(np.arange(len(core)),ind.shape[1]-1);col=ind[:,1:].ravel();dd=dist[:,1:];scale=np.maximum(dd[:,-1],1e-6);wt=np.exp(-dd/scale[:,None]).ravel();g=ig.Graph(n=len(core),edges=list(zip(row.tolist(),col.tolist())),directed=False);g.es['weight']=wt;part=leidenalg.find_partition(g,leidenalg.RBConfigurationVertexPartition,weights='weight',resolution_parameter=.65,seed=123,n_iterations=2);labels=np.asarray(part.membership).astype(str);o.loc[o.index[core],'lineage_cluster_v21']=[lineage+'::'+k for k in labels];o.loc[o.index[small],'lineage_cluster_v21']=lineage+'::low_information_unclustered'
  pd.DataFrame({'gene':a.var_names}).to_csv(ROOT/f'metadata/{prefix}_features.tsv',sep='\t',index=False);info[lineage]={'cells':len(ix),'informative':len(core),'features':a.n_vars,'n_pcs':ncomp,'clusters':len(set(labels)),'method':'section-aware HVG+identity genes; sparse log1p10k variance scale clip10; centered PCA; kNN15; Leiden resolution0.65; no forced integration','source_limitation':'one section per source: technical and biological differences confounded','small_cells':len(small)}
 group_evidence(x,o,genes,'lineage_cluster_v21',prefix,universe=ix)
 o.to_parquet(ROOT/'metadata/lineage_structure.parquet',compression='zstd');js(ROOT/'metadata/lineage_structure_methods.json',info);print('DONE',lineage,flush=True)
np.savez_compressed(ROOT/'metadata/lineage_pca.npz',cell_id=o.index.to_numpy().astype(str),X_pca=pcaout);print('ALL_LINEAGE_STRUCTURE_COMPLETE',flush=True)
