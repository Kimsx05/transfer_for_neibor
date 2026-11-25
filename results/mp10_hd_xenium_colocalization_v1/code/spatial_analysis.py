from base import *
from display import *
from unittest.mock import patch
from anndata import AnnData
import squidpy as sq
import squidpy.gr._nhood as nh
from scipy.stats import spearmanr
SQVER=sq.__version__
def enrich(s,d,g,labels,r,key):
    # Pull the moments from the very same Squidpy kernel invocation, so observed and null use identical graph, labels and RNG.
    v=pd.Series(np.asarray(labels),dtype='string');keep=v.notna().to_numpy();v=v[keep].astype('category');adj=g[keep][:,keep].tocsr()
    cats=v.cat.categories; obs=pd.DataFrame({'type':pd.Categorical(v.to_numpy(),categories=cats)},index=list(map(str,range(len(v)))))
    a=AnnData(X=sparse.csr_matrix((len(obs),0)),obs=obs);a.obsp['spatial_connectivities']=adj
    moments=[];orig=nh._permutation_moments_counts
    def kernel(*args):
        z=orig(*args);moments.extend(z);return z
    with patch.object(nh,'_permutation_moments_counts',kernel):
        result=sq.gr.nhood_enrichment(a,cluster_key='type',n_perms=1000,rng=123,n_jobs=1,normalization='none',copy=True,handle_nan='keep',show_progress_bar=False)
    ob=result.counts.astype('int64');delta=moments[0]/1000;exp=ob+delta;sd=np.sqrt(np.maximum((moments[1]-moments[0]*delta)/1000,0))
    zz=np.divide(-delta,sd,out=np.full_like(sd,np.nan),where=sd>0);np.testing.assert_allclose(result.zscore,zz,equal_nan=True)
    ratio=np.divide(ob,exp,out=np.full_like(exp,np.nan),where=exp>0)
    counts=v.value_counts();a1,a2=np.meshgrid(cats,cats,indexing='ij')
    t=pd.DataFrame({'from':a1.ravel(),'to':a2.ravel(),'observed':ob.ravel(),'expected':exp.ravel(),'null_sd':sd.ravel(),'z_score':zz.ravel(),'observed_expected':ratio.ravel()});t['n_from']=t['from'].map(counts);t['n_to']=t.to.map(counts);t['directed_edges']=adj.nnz;t['nodes']=adj.shape[0];t['seed']=123;t['permutations']=1000;t['squidpy']=SQVER
    put(t,d/f'{key}_{r}_matrix.tsv.gz')
    fig,ax=plt.subplots(figsize=(max(5.2,len(cats)*.22+1.8),max(4.7,len(cats)*.22+1)));z=zz.copy();np.fill_diagonal(z,np.nan)
    limit=max(1.,np.nanquantile(np.abs(z),.97)) if np.isfinite(z).any() else 1.;im=ax.imshow(z,cmap='RdBu_r',vmin=-limit,vmax=limit);ax.set_xticks(range(len(cats)),cats,rotation=90,fontsize=6);ax.set_yticks(range(len(cats)),cats,fontsize=6);fig.colorbar(im,ax=ax,shrink=.55,label='Neighborhood enrichment Z-score',extend='both');save(fig,s,f'{key}_{r}_matrix')
def composition(s,d,g,z,c,ids,coverage,r,high,state=None):
    # One row per evaluable center and neighbor category. All nodes remain in the denominator, including mixed and unresolved bins.
    lab=z.supported.astype(str);cats=sorted(set(lab)|{'Treg','Endothelial','Fibroblast','Myeloid'})
    feats={k:lab.eq(k).to_numpy() for k in cats}
    for maj in ['T/NK','Endothelial','Fibroblast','Myeloid']:feats[maj+' (all)']=z.major.eq(maj).to_numpy()
    n=np.asarray(g.sum(1)).ravel(); known=np.asarray(g@(np.logical_not(lab.isin(['Mixed','Unknown'])).astype('i4'))).ravel();fine=np.asarray(g@(z.flag.eq('subtype_supported').astype('i4'))).ravel();evaluable=(n>0)&(coverage>=.5)
    y=pd.DataFrame({'id':ids,'high':high,'coverage':coverage,'independent_neighbors':n,'known_neighbors':known,'supported_subtype_neighbors':fine,'evaluable':evaluable})
    effects=[];ev=np.flatnonzero(evaluable);summ=[]
    for k,flag in feats.items():
        ct=np.asarray(g@flag.astype('i4')).ravel();prop=np.divide(ct,n,out=np.full(len(n),np.nan),where=evaluable);density=np.divide(ct,coverage*np.pi*r**2/1e6,out=np.full(len(n),np.nan),where=evaluable)
        y['count_'+k]=ct;y['fraction_'+k]=prop
        # Unrepresented subtype in the section has no measurable identity-specific comparison.
        if not np.any(flag):prop[:]=np.nan;density[:]=np.nan;y['fraction_'+k]=prop
        groups=[ev[np.asarray(high)[ev]==b] for b in [True,False]]
        avg=lambda v,i:np.nanmean(v[i]) if len(i) and np.isfinite(v[i]).any() else np.nan
        h,o=groups;hp,op=avg(prop,h),avg(prop,o);hd,od=avg(density,h),avg(density,o)
        effects.append((s,r,k,len(h),len(o),int(flag.sum()),hp,op,hp-op,hd,od,hd-od,avg(n,h),avg(n,o),avg(coverage,h),avg(coverage,o)))
        if state is not None and len(ev)>2:
            for key,val in state.items():
                ok=ev[np.isfinite(val[ev])&np.isfinite(prop[ev])];rho=spearmanr(val[ok],prop[ok]).statistic if len(ok)>2 and np.ptp(val[ok])>0 and np.ptp(prop[ok])>0 else np.nan;summ.append((s,r,k,key,len(ok),rho))
    y.to_parquet(d/f'center_composition_{r}.parquet',index=False)
    put(pd.DataFrame(effects,columns=['sample','radius','neighbor','n_high','n_other','section_nodes','fraction_high','fraction_other','delta','density_high_mm2','density_other_mm2','delta_density_mm2','mean_neighbors_high','mean_neighbors_other','coverage_high','coverage_other']),d/f'effects_{r}.tsv')
    if summ:put(pd.DataFrame(summ,columns=['sample','radius','neighbor','state','n','Spearman']),d/f'continuous_neighbors_{r}.tsv')
def analyze(s):
    hd=s in HD;d=O/('hd' if hd else 'xenium')/s;z=pd.read_parquet(d/'labels.parquet');node=pd.read_parquet(d/'graph_nodes.parquet').set_index('id');assert node.index.equals(z.index)
    for r in ([50,200]):
        graph=d/f'graph_{r}.npz';g=sparse.load_npz(graph) if graph.exists() else None
        l=z.supported.copy()
        unresolved=z.major.isin(['T/NK','Fibroblast','Myeloid','Endothelial'])&z.flag.eq('major_only');l.loc[unresolved]=z.major[unresolved]+' (unresolved)'
        if hd:
            e=pd.read_parquet(d/'epi64_centers.parquet');g2=sparse.load_npz(d/f'epi64_independent_{r}.npz')
            composition(s,d,g2,z,e,e.region_id,e[f'coverage_{r}'].to_numpy(),r,e.H05.to_numpy())
        else:
            t=z.epi_eligible.to_numpy();cen=z[t];high=cen.chol_group.map({'Chol-high Epi':True,'other Epi':False}).to_numpy(dtype=float,na_value=np.nan)
            composition(s,d,g[t],z,cen,cen.index,node.loc[cen.index,f'coverage_{r}'].to_numpy(),r,high,{k:cen[k].to_numpy() for k in ['Chol_biosynthesis_8','MP10_covered','H05_transfer']})
        if g is None:continue
        for key,labels in [('major',z.major),('subtype',l)]:
            if not (d/f'{key}_{r}_matrix.tsv.gz').exists():enrich(s,d,g,labels,r,key)
        if not hd and z.chol_group.notna().any() and not (d/f'chol_{r}_matrix.tsv.gz').exists():
            # Full graph constructed before restricting to measured Epi and supported TME.
            lab=l.where(~z.major.isin(['Epithelial','Unknown','Mixed']) & ~unresolved)
            lab=lab.where(~z.epi_eligible,z.chol_group)
            enrich(s,d,g,lab,r,'chol')
    from coverage_tables import coverage
    coverage(s)
    if hd:
        from hd_mixture_review import review_weights
        review_weights(s)
    from summary_figures import plot_results
    plot_results(s)
if __name__=='__main__':
    for s in sys.argv[1:]:analyze(s)

