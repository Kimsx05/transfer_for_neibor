from common import *
import anndata as ad
from scipy.stats import rankdata

def main():
 od=OUT/'07_source_validation';done=od/'general_pseudospot_status.json'
 if done.exists():return
 a=ad.read_h5ad(ROOT/'subtype/v2/reference_training_v2/04_gene_harmonization/reference_v2_raw_counts_minimally_harmonized.h5ad',backed='r');md=a.obs.copy();md['position']=np.arange(len(md));rng=np.random.default_rng(456);pools=[];selection=[]
 for ds,g in md.groupby('dataset',observed=True):
  eligible=[]
  for src,v in g.groupby('reference_batch',observed=True):
   ep=v[v.reference_celltype_v2.eq('Epithelial')];ne=v[~v.reference_celltype_v2.eq('Epithelial')]
   if len(ep)>=100 and len(ne)>=100:eligible.append((str(src),ep,ne))
  for src,ep,ne in sorted(eligible,key=lambda z:z[0])[:2]:
   ei=rng.choice(ep.position.to_numpy(),100,False);ni=rng.choice(ne.position.to_numpy(),100,False);selection.extend(ei);selection.extend(ni);pools.append((ds,src,ei,ni))
 ids=np.sort(selection);x=a.X[ids].tocsr();genes=list(a.var_names);a.file.close();ixmap={int(z):i for i,z in enumerate(ids)};top=pd.read_csv(OUT/'00_preflight/MP10_top20.tsv',sep='\t').gene.tolist();gi=pd.Index(genes).get_indexer(top);assert (gi>=0).all();K=int(np.ceil(.05*len(genes)));den=sum(K-r for r in range(1,21))
 def score(v):return np.maximum(0,K-rankdata(-v,method='average')[gi]).sum()/den
 bags=[];members=[]
 for ds,src,ei,ni in pools:
  ep=np.array([ixmap[int(i)] for i in ei]);ne=np.array([ixmap[int(i)] for i in ni]);es=np.array([score(x[i].toarray().ravel()) for i in ep]);order=np.argsort(es,kind='stable');low=ep[order[:25]];high=ep[order[-25:]]
  # Identical non-epithelial members across states at each epithelial number/repetition.
  for e in [0,5,10,15,20]:
   for rep in range(2):
    background=rng.choice(ne,20-e,False)
    for state,pool in [('fixed_all_epi',ep),('low_score',low),('high_score',high)]:
     if e==0 and state!='fixed_all_epi':continue
     chosen=np.r_[rng.choice(pool,e,False),background];v=np.asarray(x[chosen].sum(0)).ravel().astype(np.int64);raw=int(v.sum());thin=rng.binomial(v,min(1,7005/raw));sc=score(thin);bid=f'{src}::{state}::e{e}::r{rep}';bags.append(dict(bag_id=bid,dataset=ds,source_sample=src,state=state,n_epi=e,n_cells=20,raw_UMI=raw,attained_UMI=int(thin.sum()),target_UMI=7005,MP10_Top20_AUCell=sc));members.extend(dict(bag_id=bid,cell_id=str(md.index[ids[i]])) for i in chosen)
 tab(od/'general_pseudospot.tsv',bags);tab(od/'general_pseudospot_members.tsv.gz',members);tab(od/'pseudospot_source_selection.tsv',pd.DataFrame([dict(dataset=d,source_sample=s,n_epi=100,n_other=100) for d,s,e,n in pools]));status=dict(status='GENERAL_MIXING_COMPLETE_NOT_ASSOCIATION_ATTRIBUTION',n_bags=len(bags),n_sources=len(pools),n_datasets=len(set(d for d,s,e,n in pools)),genes=len(genes),seed=456,source_selection='up to two eligible source samples per dataset in lexical order; 100 Epithelial plus100 other cells per source; no cross-source mixing',state='MP10 itself; only fidelity, not independent biology',limitation='full reference RNA universe; not actual spatial universe/depth matched; not environment-conditioned composition-only; no true patient independence assumed');js(done,status)
 old=json.loads((od/'single_cell_audit.json').read_text());old['n_general_bags']=len(bags);old['pseudospot_details']=str(done);js(od/'single_cell_audit.json',old);log('General pseudo-spot actual bags '+str(len(bags)))
if __name__=='__main__':main()
