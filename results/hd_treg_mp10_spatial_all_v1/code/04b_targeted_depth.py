from common import *
cfg=yaml.safe_load((R/'frozen_config.yaml').read_text());sig=cfg['signatures'];genes=sorted(set(sum(sig.values(),[])+TOP));M=cfg['UCell_maxRank'];gate();allrows=[]
for s in SECTIONS:
 x,b,ge,sy,path=load8(s);sel=np.array([np.flatnonzero(sy==g)[0]for g in genes]);ix=pd.Series(np.arange(len(b)),index=b.bin_id)
 for kind in ['T32','T16']:
  d=pd.read_parquet(R/s/f'{kind}_regions.parquet');v=d[d.measurement_eligible];rng=np.random.default_rng(123);chosen=np.concatenate([rng.choice(z,min(150,len(z)),False)for _,z in pd.Series(v.index).groupby(pd.cut(v.nFeature,[-1,500,1000,2000,4000,np.inf]).to_numpy(),observed=True)]);q=d.loc[chosen].copy();m=pd.read_parquet(R/s/f'{kind}_members.parquet');m=m[m.region_id.isin(q.region_id)];ri=pd.Series(np.arange(len(q)),index=q.region_id);agg=sparse.csr_matrix((np.ones(len(m),np.int64),(m.region_id.map(ri),m.bin_id.map(ix))),shape=(len(q),len(b)));z=(agg@x).tocsr();z.eliminate_zeros();assert np.array_equal(np.asarray(z.sum(1)).ravel(),q.nCount)
  z.data=rng.binomial(z.data.astype(int),.5);z.eliminate_zeros();rk=pd.DataFrame(selected_ranks(z.indptr,z.indices,z.data,z.shape[1],sel),columns=genes);q['half_Treg_main']=ucell(rk,sig['Treg_main'],M);q['half_T_lineage']=ucell(rk,sig['T_lineage'],M);q[['region_id','nCount','nFeature','Treg_main','T_lineage','half_Treg_main','half_T_lineage']].to_parquet(R/s/f'{kind}_targeted_depth.parquet',index=False)
  allrows.append(dict(section=s,branch=kind,n=len(q),sampling='fixed_seed_nFeature_stratified_among_measurement_eligible_only',Treg_half_rho=spearmanr(q.Treg_main,q.half_Treg_main).statistic,T_lineage_half_rho=spearmanr(q.T_lineage,q.half_T_lineage).statistic))
 print('TARGETED_DEPTH',s,flush=True)
put(pd.DataFrame(allrows),'targeted_T_depth_validation.tsv')
