from common import *
import anndata as ad
from scipy.stats import rankdata
import matplotlib.pyplot as plt

def auc(x,genes,top):
 K=int(np.ceil(.05*len(genes)));ix=pd.Index(genes).get_indexer([g for g in top if g in genes]);den=sum(K-r for r in range(1,min(len(ix),K)+1));rank=rankdata(-x,method='average');contribution=np.maximum(0,K-rank[ix])/den
 return contribution.sum(),contribution

def main():
 out=OUT/'07_source_validation';out.mkdir(exist_ok=True)
 if (out/'single_cell_audit.json').exists():return
 aroot=ROOT/'neighborhood_v1/00_inputs_and_audit/mp_completion';genes=pd.read_csv(aroot/'02_frozen_top20_gene_sets.csv');genes=genes[genes.MP.eq('EpiMP10')].sort_values('rank');top=genes.gene.tolist();assert len(top)==20;tab(OUT/'00_preflight/MP10_top20.tsv',genes)
 av=pd.read_csv(aroot/'05_gene_availability_by_section.csv');av=av[av.MP.eq('EpiMP10')];tab(OUT/'00_preflight/MP10_gene_availability.tsv',av)
 manifest=pd.read_csv(ROOT/'neighborhood_v1/00_inputs_and_audit/metadata/frozen_MP_scores_manifest.csv');tab(OUT/'00_preflight/MP10_score_assets.tsv',[dict(section=r.section,path=r.score_file,exists=Path(r.score_file).exists(),sha256=sha(r.score_file)) for r in manifest.itertuples()])
 path=ROOT/'subtype/v2/reference_training_v2/04_gene_harmonization/reference_v2_raw_counts_minimally_harmonized.h5ad';a=ad.read_h5ad(path,backed='r');md=a.obs.copy();md['cell_id']=md.index;md['position']=np.arange(len(md));rng=np.random.default_rng(123);sel=[]
 for typ,g in md.groupby('reference_celltype_v2',observed=True):
  groups=list(g.groupby('reference_batch',observed=True));caps=[len(v) for _,v in groups];alloc=allocate(caps,min(300,len(g)))
  for (_,v),n in zip(groups,alloc):sel.extend(rng.choice(v.position.to_numpy(),n,replace=False))
 sel=np.sort(sel);sub=md.iloc[sel].copy().reset_index(drop=True);x=a.X[sel].tocsr();universe=list(a.var_names);missing=[g for g in top if g not in universe];assert not missing;assert np.isfinite(x.data).all() and (x.data>=0).all();a.file.close()
 scores=[];contrib=[]
 for i in range(len(sub)):
  sc,c=auc(x[i].toarray().ravel(),universe,top);scores.append(sc);contrib.append(c)
 sub['MP10_Top20_AUCell']=scores;sub['raw_UMI']=np.asarray(x.sum(1)).ravel();tab(out/'single_cell_selected.tsv.gz',sub);sparse.save_npz(out/'selected_counts.npz',x);tab(out/'selected_genes.tsv',pd.DataFrame({'gene':universe}));np.save(out/'single_cell_gene_AUC_contribution.npy',contrib)
 summary=sub.groupby(['dataset','reference_batch','reference_celltype_v2'],observed=True).agg(n=('cell_id','size'),mean=('MP10_Top20_AUCell','mean'),median=('MP10_Top20_AUCell','median'),nonzero_fraction=('MP10_Top20_AUCell',lambda s:(s>0).mean())).reset_index();tab(out/'source_by_dataset_sample_type.tsv',summary)
 c=pd.DataFrame(contrib,columns=top);c['cell_type']=sub.reference_celltype_v2.to_numpy();tab(out/'gene_contribution_by_type.tsv',c.groupby('cell_type',observed=True).mean().reset_index())
 types=summary.groupby('reference_celltype_v2',observed=True)['mean'].mean().sort_values();fig,ax=plt.subplots(figsize=(9,9));ax.barh(types.index,types.values);ax.set_xlabel('Top20 AUCell: equal source-sample mean, full reference RNA universe');ax.set_title('Descriptive source audit; not spatial score calibration');fig.tight_layout();fig.savefig(out/'single_cell_source.png',dpi=160);plt.close(fig)
 # General center mixing: source-matched within reference_batch, counts-first.
 bags=[];members=[];valid_sources=[]
 for source,g in sub.groupby('reference_batch',observed=True):
  ep=g.index[g.reference_celltype_v2.eq('Epithelial')].to_numpy();ne=g.index[~g.reference_celltype_v2.eq('Epithelial')].to_numpy()
  if len(ep)<8 or len(ne)<8:continue
  valid_sources.append((source,ep,ne))
 for source,ep,ne in valid_sources:
  es=np.asarray(scores)[ep];lo=ep[es<=np.quantile(es,.25)];hi=ep[es>=np.quantile(es,.75)]
  for state,pool in [('fixed_all_epi',ep),('low_score',lo),('high_score',hi)]:
   for e in [0,5,10,15,20]:
    for rep in range(2):
     if len(pool)<e or len(ne)<20-e:continue
     ids=np.r_[rng.choice(pool,e,False),rng.choice(ne,20-e,False)];counts=np.asarray(x[ids].sum(0)).ravel().astype(np.int64);raw=int(counts.sum());prob=min(1,7005/raw) if raw else 0;thin=rng.binomial(counts,prob);sc,_=auc(thin,universe,top);bid=f'{source}::{state}::epi{e}::rep{rep}'
     bags.append(dict(bag_id=bid,source=source,state=state,n_epi=e,n_cells=20,raw_UMI=raw,attained_UMI=int(thin.sum()),target_UMI=7005,score=sc));members.extend(dict(bag_id=bid,cell_id=sub.cell_id.iloc[i]) for i in ids)
 tab(out/'general_pseudospot.tsv',bags);tab(out/'general_pseudospot_members.tsv.gz',members)
 js(out/'single_cell_audit.json',dict(status='PARTIAL_SOURCE_VALIDATION',n_cells=len(sub),n_types=sub.reference_celltype_v2.nunique(),n_source_samples=sub.reference_batch.nunique(),n_genes=len(universe),top20_missing=missing,score_method='Original top20; ranks average; K=ceil(.05*n_genes); original AUCell formula; full reference RNA universe',single_cell_input=str(path),input_size=path.stat().st_size,n_general_bags=len(bags),normal_vs_malignant='No independently verified cell malignancy annotation in this source; not inferred from sample suffix',source_patients='reference_batch is source sample, not verified independent patient',state_definition='MP10 itself, only scoring-fidelity check, not independent biological validation',pseudospot_scope='general source-matched mixing at target 7005 UMI; does not match all actual anchors/depth/gene universes; NOT environment-conditioned composition-only attribution',spatial_MP10_values_loaded=False,required_remaining=['environment-conditioned composition-only control with frozen labels and actual center composition','matched spatial gene-universe and depth coverage','independent malignant/normal cell labels','high-Epi and complete-ring real-data sensitivity after region freeze'],legacy_DestVI='v8 19-gene G0-UCell is auxiliary evidence; not validation of original Top20 AUCell'))
 log('Independent source audit complete')
if __name__=='__main__':main()
