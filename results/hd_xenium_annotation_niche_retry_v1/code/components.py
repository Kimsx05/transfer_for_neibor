from base_data import *
sys.path.insert(0,str(O/'shared/deps'))
import tacco as tc
rng=np.random.default_rng(1729)
f=next((P/'subtype').rglob('reference_v2_final_30type_regression_input.h5ad'));ref=ad.read_h5ad(f,backed='r')
rs=ref.obs;fine=rs.reference_celltype_v2.astype(str);ma=fine.map(pd.read_csv(V/'shared/taxonomy.tsv',sep='\t').set_index('original').major.to_dict());idx=[]
for (_,lab),ii in rs.assign(coarse=ma).groupby(['orig.ident','reference_celltype_v2'],observed=True).indices.items():
 idx.extend(rng.choice(ii,min(8,len(ii)),replace=False))
idx=np.sort(idx);ref=ref[idx].to_memory();ref.obs['coarse']=pd.Categorical(ma.iloc[idx].to_numpy());assert ref.obs.coarse.notna().all();pd.Series(ref.obs_names).to_csv(O/'hd/tacco_reference_ids.tsv.gz',index=False)
ex=norm(ref.X);maj=list(ref.obs.coarse.cat.categories);pr=np.vstack([ex[ref.obs.coarse.eq(k).to_numpy()].mean(0).A1 for k in maj]);sel=set(np.concatenate([np.argsort(-(v-pr.mean(0)) / (.15+pr.mean(0)))[:100] for v in pr]));sel.update(np.flatnonzero(ref.var_names.isin(set(sum(M.values(),[])))));genes=ref.var_names[sorted(sel)]
for s in (sys.argv[1:] or HD):
 a=dataset(s);g=genes.intersection(a.var_names,sort=False);rr=ref[:,g].copy();q=a[:,g].copy();ii=np.flatnonzero((q.X.sum(1).A1>0)&(a.obs.nFeature>=5));ii=rng.permutation(ii)
 out=pd.DataFrame(np.nan,index=a.obs_names,columns=maj)
 for start in range(0,len(ii),10000):
  sub=q[ii[start:start+10000]].copy()
  res=tc.tl.annotate(sub,rr,'coarse',method='OT',bisections=2,platform_iterations=0,multi_center=None,verbose=0,assume_valid_counts=True)
  out.loc[res.index,res.columns]=res.values
 valid=out.notna().all(1);assert np.allclose(out[valid].sum(1),1,atol=1e-4)
 out.to_parquet(paths(s)/'tacco_weights.parquet');pd.Series(g).to_csv(paths(s)/'tacco_genes.tsv',index=False)

