from common import *
cfg=yaml.safe_load((R/'frozen_config.yaml').read_text());B=cfg['permutation']['B'];allrows=[];edgerows=[];robust=[];curves=[];hist=[];exch=[];bootdetails=[]
hot=pd.read_parquet(R/'tables/hotspot_membership.parquet')
for sn,s in enumerate(SECTIONS):
 path=R/s/'enrichment_effects.tsv'
 if path.exists()and os.environ.get('FORCE')!='1':print('CACHE',s,flush=True);continue
 n=pd.read_parquet(R/s/'neighborhoods.parquet');e=pd.read_parquet(R/s/'Epi_regions.parquet').set_index('region_id');pp=pd.read_parquet(V/s/'region_pairs.parquet');pp=pp[pp.T_grid_um.eq(32)&pp.radius_um.eq(200)&~pp.shared_original_bins];rows=[];edges=[];rb=[];cv=[];hs=[];ee=[];bd=[]
 for radius in RADII:
  u=n[n.radius_um.eq(radius)&n.U_r].reset_index(drop=True);N=len(u);group=u.stratum.to_numpy();L=u[LABELS].to_numpy(np.float32);H=u[['near_'+h for h in HIGHS]].to_numpy(np.float32);D=u[['n_'+h for h in HIGHS]].to_numpy(np.float32);pos=L.sum(0);NH=H.sum(0);aobs=L.T@H;Oobs=L.T@D
  permA=np.zeros((B,2,2),np.float32);permO=np.zeros_like(permA);curve_counts={};curve_null={}
  if radius==200:
   lookup=pd.Series(np.arange(N),index=u.region_id);q=pp[pp.T_region_id.isin(lookup.index)].copy();q['i']=q.T_region_id.map(lookup);q['step']=np.ceil((q.centroid_distance_um-1e-9)/10).clip(0,20).astype(int)
   for h in HIGHS:
    z=q[e.loc[q.Epi_region_id,h].to_numpy()];cc=np.zeros((N,21),np.float32);np.add.at(cc,(z.i.to_numpy(),z.step.to_numpy()),1);curve_counts[h]=np.cumsum(cc,axis=1);curve_null[h]=np.zeros((B,2,21),np.float32)
  for start,lp in perm_batches(L,group,B,123+sn*10+radius):
   for j in range(2):
    permA[start:start+len(lp),j,:]=lp[:,:,j]@H;permO[start:start+len(lp),j,:]=lp[:,:,j]@D
    if radius==200:
     for h in HIGHS:curve_null[h][start:start+len(lp),j,:]=lp[:,:,j]@curve_counts[h]
  # Same permutation indices for both labels/high definitions; only T labels move.
  for j,label in enumerate(LABELS):
   groups=[ix for _,ix in u.groupby('stratum').groups.items()];mixed=[np.asarray(ix)for ix in groups if 0<L[ix,j].sum()<len(ix)];moving=sum(len(ix)for ix in mixed);movpos=sum(L[ix,j].sum()for ix in mixed);validnull=moving>=20 and movpos>=5 and len(mixed)>=2;reason='ok'if validnull else'insufficient_label_exchangeability';ee.append(dict(section=s,radius_um=radius,label=label,N=N,positive=int(pos[j]),strata=len(groups),mixed_strata=len(mixed),exchangeable_objects=moving,movable_positive=int(movpos),nonexchangeable_fraction=1-moving/N if N else np.nan,null_eligible=validnull,reason=reason))
   for k,h in enumerate(HIGHS):
    a=float(aobs[j,k]);b=float(pos[j]-a);c=float(NH[k]-a);d=float(N-pos[j]-c);fx=effects(a,b,c,d);na=permA[:,j,k];nfx=effects(na,pos[j]-na,NH[k]-na,N-pos[j]-NH[k]+na);P,Z,mu,sd=empirical(float(fx['FE']),nfx['FE'],validnull)
    base=dict(section=s,radius_um=radius,high_definition=h,label_definition=label);rec=dict(**base,a=int(a),b=int(b),c=int(c),d=int(d),N=N,label_positive=int(pos[j]),near_high_denominator=int(NH[k]),comparison_denominator=int(N-NH[k]),positive_fraction_high=float(divide(a,a+c)),positive_fraction_comparison=float(divide(b,b+d)),overall_positive_fraction=float(divide(a+b,N)),P=P,null_mean_FE=mu,null_sd_FE=sd,permutations=B,seed=123+sn*10+radius,estimability_reason=reason if not validnull else'constant_or_undefined_null'if not np.isfinite(P)else'estimable_restricted_random_label',family='primary_FE'if label=='FOXP3_supported'and h=='H05'else'auxiliary_FE')
    for m,v in fx.items():rec[m]=float(v);rec[m+'_status']=status(float(v))
    rows.append(rec);po,zo,mo,so=empirical(Oobs[j,k],permO[:,j,k],validnull);edges.append(dict(**base,O=int(Oobs[j,k]),null_mean=mo,null_sd=so,Z=zo,P=po,permutations=B,seed=123+sn*10+radius,N_T=N,n_label=int(pos[j]),total_valid_high_edges=int(D[:,k].sum()),reason=reason if not validnull else'zero_null_sd'if so<1e-12 else'estimable',family='auxiliary_edge'))
    for stat,vals in [('a',permA[:,j,k]),('edge_count',permO[:,j,k])]:
     values,counts=np.unique(vals,return_counts=True)
     for value,count in zip(values,counts):hs.append(dict(**base,statistic=stat,null_value=value,frequency=count))
    # Center-based spatial block bootstrap. Degenerate draws remain in accounting.
    for size in [400,800]:
     block=block_ids(u,size);nb=block.max()+1 if N else 0;arr=np.column_stack([L[:,j]*H[:,k],L[:,j]*(1-H[:,k]),(1-L[:,j])*H[:,k],(1-L[:,j])*(1-H[:,k])]);bc=np.column_stack([np.bincount(block,weights=arr[:,z],minlength=nb)for z in range(4)]);rng=np.random.default_rng(123);weights=rng.multinomial(nb,np.ones(nb)/nb,size=999)if nb else np.zeros((999,0));boot=weights@bc;bf=effects(*boot.T);pb=len(np.unique(block[L[:,j]>0]));record=dict(**base,block_um=size,n_blocks=nb,positive_blocks=pb,B=999,zero_positive_draw_fraction=float(((boot[:,0]+boot[:,1])==0).mean()),no_high_draw_fraction=float(((boot[:,0]+boot[:,2])==0).mean()))
     for metric,vals in bf.items():
      finite=np.isfinite(vals);rel=nb>=12 and pb>=3 and finite.mean()>=.95 and np.isfinite(float(fx[metric]));lo,hi=np.quantile(vals[finite],[.025,.975],method='nearest')if rel else[np.nan,np.nan];record.update({metric:float(fx[metric]),metric+'_low':lo,metric+'_high':hi,metric+'_NA_fraction':float(np.isnan(vals).mean()),metric+'_Inf_fraction':float(np.isinf(vals).mean()),metric+'_CI_status':'estimable'if rel else'insufficient_blocks_or_degenerate_draws'})
      if size==400:rec.update({metric+'_low':lo,metric+'_high':hi,metric+'_bootstrap_NA_fraction':float(np.isnan(vals).mean()),metric+'_bootstrap_Inf_fraction':float(np.isinf(vals).mean()),metric+'_CI_status':record[metric+'_CI_status']})
     bd.append(record)
    # Prespecified positive spatial components; no MP10 used in grouping or exclusion.
    hh=hot[hot.section.eq(s)&hot.label.eq(label)]
    for hotspot,qq in hh.groupby('hotspot'):
     gone=u.region_id.isin(qq.region_id);keep=~gone;aa=int((u.loc[keep,label]&u.loc[keep,'near_'+h].astype(bool)).sum());bb=int(u.loc[keep,label].sum()-aa);cc=int(u.loc[keep,'near_'+h].sum()-aa);dd=int(keep.sum()-aa-bb-cc);fe=effects(aa,bb,cc,dd)
     rb.append(dict(**base,diagnostic='leave_one_positive_hotspot_out',hotspot=int(hotspot),removed_centers=int(gone.sum()),remaining_N=int(keep.sum()),remaining_positive=aa+bb,FE=float(fe['FE']),PR=float(fe['PR']),OR=float(fe['OR']),P=np.nan,exclusion='members of fixed64um linked positive component; unchanged background and remaining Epi neighborhoods'))
    if radius==200:
     obs=L[:,j]@curve_counts[h];null=curve_null[h][:,j,:];means=null.mean(0);low,high=np.quantile(null,[.025,.975],axis=0)
     for st in range(21):cv.append(dict(**base,distance_um=st*10,N_fixed_U200=N,n_positive=int(pos[j]),observed_edges=float(obs[st]),mean_neighbors_per_positive=float(divide(obs[st],pos[j])),null_mean_edges=float(means[st]),null_low_edges=float(low[st]),null_high_edges=float(high[st]),observed_over_null=float(divide(obs[st],means[st])),null_low_ratio=float(divide(low[st],means[st])),null_high_ratio=float(divide(high[st],means[st])),null_eligible=validnull,permutations=B,name='source_excluded_cross_proximity_curve',P=np.nan,envelope='pointwise restricted-label, not simultaneous and not standard Kcross'))
  print('ENRICHMENT',s,radius,N,flush=True)
 for name,frame in [('enrichment_effects',rows),('neighborhood_enrichment',edges),('spatial_robustness',rb),('cross_distance_results',cv),('null_histograms',hs),('exchangeability_summary',ee),('bootstrap_diagnostics',bd)]:pd.DataFrame(frame).to_csv(R/s/(name+'.tsv'),sep='\t',index=False,na_rep='NA')
print('ENRICHMENT_ALL_DONE',flush=True)
