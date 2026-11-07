from common import *

def fit(X,y,w=None):
 n=len(y);p=np.linalg.matrix_rank(X);df=n-p
 if df<6:return None
 sw=np.ones(n) if w is None else np.sqrt(w);xx=X*sw[:,None];yy=y*sw;inv=np.linalg.pinv(xx.T@xx);b=inv@xx.T@yy;res=yy-xx@b;h=np.sum(xx*(xx@inv),axis=1);bread=xx*(res/np.maximum(1-h,.02))[:,None];cov=inv@bread.T@bread@inv;se=np.sqrt(np.maximum(np.diag(cov),0));pv=2*stats.t.sf(np.abs(b/np.maximum(se,1e-15)),df);crit=stats.t.ppf(.975,df);return dict(beta=b,se=se,p=pv,low=b-crit*se,high=b+crit*se,df=df,residual=y-X@b,r2=1-np.sum((y-X@b)**2)/max(1e-15,np.sum((y-y.mean())**2)))
def design(x,ds,tech=None):
 dum=pd.get_dummies(pd.Series(ds),drop_first=True,dtype=float).to_numpy();a=np.column_stack([np.ones(len(x)),x,dum])
 if tech is not None:a=np.column_stack([a,tech])
 return a

def main():
 co=read('representative_coordinates.tsv');co=co[co.input_family.eq('niche_all_Epi_all')].set_index('record_id');meta=read('unit_metadata.tsv').set_index('record_id');reg=read('endpoints.tsv.gz').set_index('endpoint_id');u=read('unit_endpoints.tsv.gz');u=u[u.identity_level.eq('analysis_unit')&u.record_id.isin(co.index)&u.value.notna()]
 bg=read('unit_background.tsv.gz');epibg=bg[bg.input_family.eq('cell_composition_whole_tissue')&bg.feature.eq('Epithelial')].set_index('record_id').value
 keep=reg.module.str.startswith('B3')|reg.module.eq('C4')|((reg.module.isin(['C1','C2']))&reg.range_type.eq('common_support_0_1000'));all_u=u.copy();u=u[u.endpoint_id.isin(reg.index[keep])];out=[];objects={};PC=(co[['PC1','PC2']]-co[['PC1','PC2']].mean())/co[['PC1','PC2']].std()
 for eid,z in u.groupby('endpoint_id',sort=False):
  r=reg.loc[eid]
  if r.module in ['C1','C2']:z=z[z.n_min.ge(10)&z.n_other_min.ge(10)]
  if len(z)<12:continue
  x=PC.loc[z.record_id].to_numpy();y=transform(z.value.to_numpy(),r.kind);ds=z.dataset_id.to_numpy()
  if len(np.unique(ds))<2:continue
  X=design(x,ds);res=fit(X,y)
  if res is None:continue
  base=fit(design(np.empty((len(y),0)),ds),y);gain=res['r2']-base['r2'];m=meta.loc[z.record_id];tech=np.column_stack([np.log1p(m.n_P95 if r.module in ['C1','C2'] else m.n_Epi),m.niche_coverage]);validtech=np.isfinite(tech).all(1);tr=fit(design(x[validtech],ds[validtech],tech[validtech]),y[validtech]) if validtech.sum()>=12 else None
  ev=epibg.reindex(z.record_id).to_numpy();ve=np.isfinite(ev);cr=fit(design(x[ve],ds[ve],ev[ve,None]),y[ve]) if ve.sum()>=12 else None
  se=z.spatial_se.to_numpy();ok=np.isfinite(se)&(se>0);wr=sr=None
  if ok.sum()>=12:
   floor=max(np.var(res['residual'],ddof=1),1e-6);weights=1/(se[ok]**2+floor);wr=fit(design(x[ok],ds[ok]),y[ok],weights);sr=fit(design(x[ok],ds[ok]),y[ok])
  dsmean=pd.DataFrame(dict(ds=ds,PC1=x[:,0],PC2=x[:,1],y=y)).groupby('ds').mean()
  for j in range(2):
   between=stats.pearsonr(dsmean.iloc[:,j],dsmean.y) if len(dsmean)>=3 and dsmean.iloc[:,j].std()>0 and dsmean.y.std()>0 else (np.nan,np.nan)
   out.append(dict(endpoint_id=eid,ecology_axis='PC'+str(j+1),n_units=len(y),n_datasets=len(dsmean),beta=res['beta'][j+1],low=res['low'][j+1],high=res['high'][j+1],p_HC3=res['p'][j+1],ecology_incremental_R2=gain,between_dataset_r=between[0],between_dataset_p=between[1],composition_adjusted_beta=cr['beta'][j+1] if cr else np.nan,technical_adjusted_beta=tr['beta'][j+1] if tr else np.nan,technical_adjusted_p=tr['p'][j+1] if tr else np.nan,n_with_spatial_se=int(ok.sum()),spatial_weighted_beta=wr['beta'][j+1] if wr else np.nan,OLS_same_spatial_subset_beta=sr['beta'][j+1] if sr else np.nan))
  objects[eid]=(z,x,y,ds)
 out=pd.DataFrame(out);keys=reg[['module','signature','score_method','statistic']].astype(str).agg('|'.join,axis=1);out['family']=out.endpoint_id.map(keys);out['q_within_endpoint_family']=out.groupby('family').p_HC3.transform(lambda p:bh(p));out['q_across_ecology_search']=bh(out.p_HC3);out['q_between_dataset_search']=bh(out.between_dataset_p);save(out.drop(columns='family'),'ecology_associations.tsv.gz')
 cand=out.merge(reg.reset_index(),on='endpoint_id');cand=cand[cand.signature.ne('Top100')&cand.statistic.isin(['r','paired_delta_r'])&cand.n_units.ge(25)&cand.n_datasets.ge(4)&~cand.score_method.eq('UCell_RNA')].copy();cand['rank_score']=cand.ecology_incremental_R2.clip(0)*cand.n_units/55/(1+np.abs(cand.beta-cand.technical_adjusted_beta));sel=[]
 for module,z in cand.groupby('module'):
  z=z.sort_values('rank_score',ascending=False).drop_duplicates('canonical_program').head(5);sel.extend(z.endpoint_id.tolist())
 extra=out.merge(reg.reset_index(),on='endpoint_id');extra=extra[extra.module.eq('C1')&extra.statistic.eq('partial_r')&extra.score_method.eq('AddModuleScore_SCT_data')&extra.n_units.ge(25)].sort_values('p_HC3').head(1);sel.extend(extra.endpoint_id.tolist());sel=list(dict.fromkeys(sel));loo=[];within=[];sensitivity=[]
 for eid in sel:
  if eid not in objects:continue
  z,x,y,ds=objects[eid];r=reg.loc[eid]
  for d in sorted(set(ds)):
   ok=ds!=d;rr=fit(design(x[ok],ds[ok]),y[ok])
   if rr:
    for j in range(2):loo.append(dict(endpoint_id=eid,left_out=d,ecology_axis='PC'+str(j+1),n_units=ok.sum(),beta=rr['beta'][j+1],low=rr['low'][j+1],high=rr['high'][j+1],p_HC3=rr['p'][j+1]))
   ok=ds==d
   if ok.sum()>=4:
    for j in range(2):
     v=stats.pearsonr(x[ok,j],y[ok]) if np.std(x[ok,j])>0 and np.std(y[ok])>0 else (np.nan,np.nan);within.append(dict(endpoint_id=eid,dataset_id=d,ecology_axis='PC'+str(j+1),n_units=ok.sum(),r=v[0],p_descriptive=v[1]))
  # All prespecified counterpart rows, not inferred from the original correlation.
  alt=reg[reg.module.eq(r.module)&reg.canonical_program.eq(r.canonical_program)&reg.candidate_type.eq(r.candidate_type)]
  for aid,ar in alt.iterrows():
   if aid==eid:continue
   az=all_u[all_u.endpoint_id.eq(aid)]
   if ar.module in ['C1','C2']:az=az[az.n_min.ge(10)&az.n_other_min.ge(10)]
   ids=z.record_id[z.record_id.isin(az.record_id)];look=z.set_index('record_id');alook=az.set_index('record_id')
   if len(ids)<12:continue
   xx=PC.loc[ids].to_numpy();dd=look.loc[ids,'dataset_id'].to_numpy();X=design(xx,dd);aa=fit(X,transform(alook.loc[ids,'value'].to_numpy(),ar.kind));bb=fit(X,transform(look.loc[ids,'value'].to_numpy(),r.kind))
   if aa and bb:
    for j in range(2):sensitivity.append(dict(endpoint_id=eid,alternate_endpoint_id=aid,ecology_axis='PC'+str(j+1),n_paired_units=len(ids),base_beta_same_units=bb['beta'][j+1],alternate_beta=aa['beta'][j+1],alternate_p=aa['p'][j+1]))
 save(loo,'ecology_leave_dataset_out.tsv');save(within,'ecology_within_dataset.tsv');save(sensitivity,'ecology_sensitivity.tsv.gz');save(cand[cand.endpoint_id.isin(sel)].sort_values('rank_score',ascending=False),'ecology_discussion_candidates.tsv')
 print('ecology complete',len(out),'associations',len(sel),'leaveout endpoints',flush=True)
if __name__=='__main__':main()
