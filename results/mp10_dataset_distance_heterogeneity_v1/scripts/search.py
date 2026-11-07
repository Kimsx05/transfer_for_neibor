from common import *

def main():
 reg=read('endpoints.tsv.gz').set_index('endpoint_id');u=read('unit_endpoints.tsv.gz');u=u[np.isfinite(u.value)].copy();u['z']=u.value;rr=u.endpoint_id.isin(reg.index[reg.kind.eq('r')]);u.loc[rr,'z']=transform(u.loc[rr,'value'],'r');u['section_only']=u.identity_level.eq('section_only').astype(int)
 allout=[];summaries=[]
 for threshold in [0,10,30]:
  isC=reg.module.isin(['C1','C2']);ids=reg.index[~isC if threshold==0 else isC];a=u[u.endpoint_id.isin(ids)].copy()
  if threshold:a=a[a.n_min.ge(threshold)&a.n_other_min.ge(threshold)]
  if threshold==30:a=a[a.endpoint_id.isin(reg.index[reg.range_type.eq('common_support_0_1000')&~reg.statistic.eq('partial_r')])]
  ds=a.groupby(['endpoint_id','dataset_id']).agg(mean_z=('z','mean'),var_z=('z','var'),n=('z','size'),n_sections=('n_sections','sum'),n_missing_unit=('section_only','sum'),n_positive=('value',lambda s:(s>0).sum()),n_spatial=('spatial_se','count')).reset_index();ds['n_min_rule']=threshold;summaries.append(ds)
  e=sorted(ds.endpoint_id.unique());pivot=lambda col:ds.pivot(index='endpoint_id',columns='dataset_id',values=col).reindex(index=e,columns=DS).to_numpy(float)
  mu=pivot('mean_z');va=pivot('var_z');n=pivot('n');ns=pivot('n_sections');mi=pivot('n_missing_unit');po=pivot('n_positive');sp=pivot('n_spatial');kind=reg.loc[e,'kind'].eq('r').to_numpy();descriptive=(reg.loc[e,'module'].eq('C2')&reg.loc[e,'range_type'].eq('full_range')).to_numpy()
  for comb in combos(DS):
   ix=[DS.index(d) for d in comb];k=len(ix);valid=np.isfinite(mu[:,ix]).all(1);validids=np.flatnonzero(valid)
   if not len(validids):continue
   mus=mu[valid][:,ix];vv=va[valid][:,ix];nn=n[valid][:,ix];missing=mi[valid][:,ix].sum(1);mean=mus.mean(1);cs=(vv/nn).sum(1)/k**2;df=(vv/nn).sum(1)**2/np.sum((vv/nn)**2/(nn-1),axis=1);se=np.sqrt(cs);eligible=(missing==0)&(nn>=2).all(1)&~descriptive[valid]&np.isfinite(se)&(se>0)
   pc=np.full(len(mean),np.nan);lo=pc.copy();hi=pc.copy();pc[eligible]=2*stats.t.sf(np.abs(mean[eligible]/se[eligible]),df[eligible]);crit=stats.t.ppf(.975,df[eligible]);lo[eligible]=mean[eligible]-crit*se[eligible];hi[eligible]=mean[eligible]+crit*se[eligible]
   pb=np.full(len(mean),np.nan)
   if k>=3:
    sd=mus.std(1,ddof=1)/np.sqrt(k);ok=(missing==0)&~descriptive[valid]&(sd>0);pb[ok]=2*stats.t.sf(np.abs(mean[ok]/sd[ok]),k-1)
   sign=mean[:,None]>=0;dsagree=((mus>=0)==sign).mean(1);unitpositive=po[valid][:,ix].sum(1);total=nn.sum(1);uagree=np.where(mean>=0,unitpositive/total,1-unitpositive/total)
   lomin=np.min((mus.sum(1)[:,None]-mus)/(k-1),axis=1) if k>1 else np.full(len(mean),np.nan);lomax=np.max((mus.sum(1)[:,None]-mus)/(k-1),axis=1) if k>1 else lomin.copy()
   kr=kind[valid];conv=lambda x:np.where(kr,np.tanh(x),x)
   allout.append(pd.DataFrame(dict(endpoint_id=np.array(e)[valid],combination=combo_id(comb),n_min_rule=threshold,n_datasets=k,n_records=total.astype(int),n_analysis_units=(total-missing).astype(int),n_sections=ns[valid][:,ix].sum(1).astype(int),n_section_only=missing.astype(int),estimate=conv(mean),conditional_low=conv(lo),conditional_high=conv(hi),p_conditional=pc,p_between_datasets=pb,dataset_direction_fraction=dsagree,unit_direction_fraction=uagree,leave_dataset_min=conv(lomin),leave_dataset_max=conv(lomax),n_spatial_records=sp[valid][:,ix].sum(1).astype(int))))
  print('search threshold',threshold,'endpoints',len(e),'rows',sum(len(x) for x in allout),flush=True)
 out=pd.concat(allout,ignore_index=True);ds=pd.concat(summaries,ignore_index=True);save(ds,'dataset_endpoint_summary.tsv.gz')
 fam=reg[['module','signature','score_method','statistic','range_type']].astype(str).agg('|'.join,axis=1);out['family']=out.endpoint_id.map(fam)+'|n'+out.n_min_rule.astype(str)
 for p,tag in [('p_conditional','conditional'),('p_between_datasets','between')]:
  out['q_'+tag+'_within_combination']=out.groupby(['family','combination'])[p].transform(lambda x:bh(x))
  out['q_'+tag+'_across_search']=out.groupby('family')[p].transform(lambda x:bh(x))
 out.drop(columns=['family','leave_dataset_min','leave_dataset_max']).to_csv(O/'tables/combination_associations.tsv.gz',sep='\t',index=False,na_rep='NA',float_format='%.4g');save([dict(combination=combo_id(c),datasets='|'.join(c),n_datasets=len(c)) for c in combos(DS)],'combinations.tsv')
 counts=out.groupby(['family','n_datasets']).agg(n_combinations=('combination','nunique'),n_rows=('endpoint_id','size'),n_conditional=('p_conditional','count'),n_between=('p_between_datasets','count')).reset_index();save(counts,'search_scope.tsv')
 print('search complete',len(out),flush=True)
if __name__=='__main__':main()
