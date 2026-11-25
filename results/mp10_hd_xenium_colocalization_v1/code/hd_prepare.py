from base import *
for s in HD:
    d=O/'hd'/s;d.mkdir(exist_ok=True); t=F/s
    meta=json.loads((t/'prepared.json').read_text())
    x,ids,g=tenx(meta['counts_path'])
    members=pd.read_parquet(t/'T16_members.parquet').set_index('bin_id').loc[ids]
    cols=['region_id','component','array_row_grid','array_col_grid','x_um','y_um','nCount','nFeature','n_bins']
    units=pd.read_parquet(t/'T16_measure.parquet',columns=cols)
    lab=pd.Index(units.region_id).get_indexer(members.region_id)
    assert min(lab)>=0 and members.index.is_unique and len(ids)==len(members) and units.region_id.is_unique
    m=sparse.csr_matrix((np.ones(len(ids),np.int32),(lab,np.arange(len(ids)))),shape=(len(units),len(ids)))
    y=(m@x).tocsr();y.eliminate_zeros()
    assert np.array_equal(np.asarray(x.sum(0)),np.asarray(y.sum(0))) and np.array_equal(np.asarray(y.sum(1)).ravel(),units.nCount)
    assert np.array_equal(np.bincount(lab),units.n_bins) and np.array_equal(np.diff(y.indptr),units.nFeature)
    nodes=pd.read_parquet(A/'hd_epi_niche_dual_methods_v1'/s/'nodes.parquet',columns=['bin_id','array_row','array_col','x_um','y_um','component','in_tissue_original']).set_index('bin_id').loc[ids]
    assert (nodes.in_tissue_original==1).all()
    for col in ['array_row','array_col']:
        assert np.array_equal(nodes[col].to_numpy()//2,units[col+'_grid'].to_numpy()[lab])
    assert np.array_equal(nodes.component,units.component.to_numpy()[lab])
    for col in ['x_um','y_um']:
        np.testing.assert_allclose(np.bincount(lab,weights=nodes[col])/units.n_bins,units[col])
    members.reset_index().to_parquet(d/'members.parquet',index=False)
    put(units,d/'units.tsv.gz')
    genes=pd.Index(g.name.unique()); collapse=sparse.csr_matrix((np.ones(len(g),np.int32),(np.arange(len(g)),genes.get_indexer(g.name))),shape=(len(g),len(genes)))
    cache=d/'counts';cache.mkdir(exist_ok=True); export(y@collapse,genes,units.region_id,cache)
    (d/'conservation.json').write_text(json.dumps(dict(section=s,genes=len(g),bins8=len(ids),bins16=len(units),umi=int(y.sum()),per_gene_conserved=True,grid_origin=[0,0],coordinate='mean physical centroid of actual members',input=meta['counts_path']),indent=2))
    print(s,len(units),int(y.sum()),flush=True)

