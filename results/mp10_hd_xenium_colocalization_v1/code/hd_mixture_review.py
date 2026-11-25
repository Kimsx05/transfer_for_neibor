from base import *
def review_weights(s):
    d=O/'hd'/s;z=pd.read_parquet(d/'labels.parquet');w=pd.read_parquet(d/'weights.parquet').reindex(z.index);fit=w.notna().all(axis=1).to_numpy();mat=w.fillna(0).to_numpy();e=pd.read_parquet(d/'epi64_centers.parquet')
    rows=[]
    for r in [200,50]:
        g=sparse.load_npz(d/f'epi64_independent_{r}.npz')
        n=np.asarray(g@fit.astype('i4')).ravel();ok=(n>0)&(e[f'coverage_{r}'].to_numpy()>=.5)
        v=np.divide(g@mat,n[:,None],out=np.full((len(n),mat.shape[1]),np.nan),where=ok[:,None])
        t=pd.DataFrame(v,columns=w.columns);t.insert(0,'fitted_neighbors',n);t.insert(0,'independent_neighbors',np.diff(g.indptr));t.insert(0,'H05',e.H05.to_numpy());t.insert(0,'region_id',e.region_id.to_numpy())
        t.to_parquet(d/f'epi64_neighbor_weights_{r}.parquet',index=False)
        for j,lab in enumerate(w.columns):
            h=ok&e.H05.to_numpy();other=ok&~e.H05.to_numpy()
            mean=lambda x:np.mean(x) if len(x) else np.nan
            yh,yo=mean(v[h,j]),mean(v[other,j]);rows.append((s,r,lab,h.sum(),other.sum(),yh,yo,yh-yo,mean(n[h]),mean(n[other])))
    put(pd.DataFrame(rows,columns=['sample','radius','component','H05_regions','other_regions','mean_weight_H05','mean_weight_other','difference','mean_fitted_neighbors_H05','mean_fitted_neighbors_other']),d/'mixture_weight_review.tsv')
if __name__=='__main__':
 for s in sys.argv[1:]:review_weights(s)

