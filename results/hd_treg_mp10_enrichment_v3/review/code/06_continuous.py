from common import *
rows=[];diagnostics=[];robust=[];roi=pd.read_csv(R/'tables/roi_manifest.tsv',sep='\t')
for s in SECTIONS:
 n=pd.read_parquet(R/s/'neighborhoods.parquet')
 for r in RADII:
  d=n[n.radius_um.eq(r)&n.U_r].reset_index(drop=True);N=len(d);covs=['T_lineage','nCount','Epi_nCount','Epi_n','coverage','T_background_balance'];C=np.ones((N,1));names=['intercept']
  for col in covs:
   v=rankz(d[col]);rr=v-C@np.linalg.lstsq(C,v,rcond=None)[0]
   if np.var(rr)>.1:C=np.column_stack([C,v]);names.append(col)
  x=rankz(d.Treg_main);X=np.column_stack([C,x]);Y=np.column_stack([rankz(d.Epi_MP10),rankz(d.fraction_H05),rankz(d.fraction_H10)]);beta=np.linalg.lstsq(X,Y,rcond=None)[0][-1];raw=[spearmanr(d.Treg_main,d[col]).statistic if d[col].nunique()>1 else np.nan for col in ['Epi_MP10','fraction_H05','fraction_H10']]
  designs=pd.DataFrame(X,columns=names+['Treg_main']);designs.insert(0,'region_id',d.region_id);designs.to_parquet(R/s/f'continuous_design_r{r}.parquet',index=False)
  for col in ['nCount','nFeature','T_lineage','T_UMI','Epi_nCount','Epi_n','coverage','Tfh','Activated_T']:
   diagnostics.append(dict(section=s,radius_um=r,feature=col,N=N,Treg_spearman=spearmanr(d.Treg_main,d[col]).statistic if d[col].nunique()>1 else np.nan,MP10_spearman=spearmanr(d.Epi_MP10,d[col]).statistic if d[col].nunique()>1 else np.nan))
  for size in [400,800]:
   block=block_ids(d,size);nb=block.max()+1;xx=np.zeros((nb,X.shape[1],X.shape[1]));xy=np.zeros((nb,X.shape[1],Y.shape[1]));np.add.at(xx,block,X[:,:,None]*X[:,None,:]);np.add.at(xy,block,X[:,:,None]*Y[:,None,:]);rng=np.random.default_rng(123);w=rng.multinomial(nb,np.ones(nb)/nb,size=999);ax=(w@xx.reshape(nb,-1)).reshape(-1,X.shape[1],X.shape[1]);ay=(w@xy.reshape(nb,-1)).reshape(-1,X.shape[1],Y.shape[1]);bs=[]
   for aa,bb in zip(ax,ay):
    bs.append(np.linalg.lstsq(aa,bb,rcond=None)[0][-1]if np.linalg.matrix_rank(aa)==aa.shape[0]else np.full(Y.shape[1],np.nan))
   bs=np.array(bs);rr=x-C@np.linalg.lstsq(C,x,rcond=None)[0]
   for j,col in enumerate(['Epi_MP10','fraction_H05','fraction_H10']):
    good=np.isfinite(bs[:,j]);ok=N>=50 and nb>=12 and good.mean()>=.95 and np.var(rr)>.1 and d[col].nunique()>1;lo,hi=np.quantile(bs[good,j],[.025,.975])if ok else[np.nan,np.nan]
    rows.append(dict(section=s,radius_um=r,readout=col,N=N,block_um=size,n_blocks=nb,spearman=raw[j],adjusted_rank_beta=float(beta[j]),low=lo,high=hi,bootstrap_valid_fraction=good.mean(),P=np.nan,covariates=','.join(names[1:]),reason='descriptive_with_block_CI'if ok else'insufficient_support_or_collinearity',no_new_formal_P=True))
   # All blocks, no chosen effect direction.
   for bid in range(nb):
    keep=block!=bid;b=np.linalg.lstsq(X[keep],Y[keep,0],rcond=None)[0][-1];robust.append(dict(section=s,radius_um=r,diagnostic='continuous_leave_one_block_out',block_um=size,excluded_block=bid,removed=int((~keep).sum()),remaining=int(keep.sum()),adjusted_rank_beta=b))
  gone=np.zeros(N,bool)
  for q in roi[roi.section.eq(s)&roi.roi.str.startswith('candidate')].itertuples():gone|=d.x_um.between(q.x0,q.x1)&d.y_um.between(q.y0,q.y1)
  keep=~gone;b=np.linalg.lstsq(X[keep],Y[keep,0],rcond=None)[0][-1]if keep.sum()>X.shape[1]else np.nan;robust.append(dict(section=s,radius_um=r,diagnostic='exclude_original_candidate_ROI_union',removed=int(gone.sum()),remaining=int(keep.sum()),adjusted_rank_beta=b))
  # Same fixed ROI exclusion for binary summaries, without changing remaining exposures.
  for label in LABELS:
   for h in HIGHS:
    a=int((d.loc[keep,label]&d.loc[keep,'near_'+h].astype(bool)).sum());b=int(d.loc[keep,label].sum()-a);c=int(d.loc[keep,'near_'+h].sum()-a);dd=int(keep.sum()-a-b-c);fx=effects(a,b,c,dd);robust.append(dict(section=s,radius_um=r,label_definition=label,high_definition=h,diagnostic='binary_exclude_original_candidate_ROI_union',removed=int(gone.sum()),remaining=int(keep.sum()),FE=float(fx['FE']),PR=float(fx['PR']),OR=float(fx['OR'])))
 print('CONTINUOUS',s,flush=True)
put(pd.DataFrame(rows),'continuous_associations.tsv');put(pd.DataFrame(diagnostics),'continuous_quality_diagnostics.tsv');put(pd.DataFrame(robust),'continuous_and_ROI_robustness.tsv')
