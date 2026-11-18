from fu_common import *
from statsmodels.stats.multitest import multipletests
F=yaml.safe_load((R/'frozen_config.yaml').read_text());d=pd.read_parquet(R/'tables/Cervilla_common_centers_r200.parquet').set_index('T_region_id');y=np.column_stack([zrank(d.MP10_bin),zrank(d.MP10_cell)]);roi=pd.read_csv(R/'tables/roi_manifest_original.tsv',sep='\t');roi=roi[roi.section.eq('Cervilla')&roi.roi.str.startswith('candidate')];hot=np.zeros(len(d),bool)
for q in roi.itertuples():hot|=d.x_um.between(q.x0,q.x1)&d.y_um.between(q.y0,q.y1)
d['exclude_hotspot']=hot;d['zrank_MP10_bin']=y[:,0];d['zrank_MP10_cell']=y[:,1];d.reset_index().to_parquet(R/'tables/Cervilla_paired_analysis.parquet',index=False)
rows=[];boots=[];lobo=[]
def fit_pair(X,Y,dat,block,variant,do_p=False):
 n=len(dat);blocks=pd.factorize(pd.MultiIndex.from_arrays([dat.component,np.floor(dat.x_um/block),np.floor(dat.y_um/block)]))[0];nb=blocks.max()+1;beta=np.linalg.lstsq(X,Y,rcond=None)[0][-1];delta=float(beta[1]-beta[0]);assert np.isclose(delta,np.linalg.lstsq(X,Y[:,1]-Y[:,0],rcond=None)[0][-1]);rec=dict(variant=variant,block_um=block,n=n,n_blocks=nb,beta_bin=float(beta[0]),beta_cell=float(beta[1]),delta=delta,P_delta=np.nan,status='descriptive')
 for name in ['beta_bin','beta_cell','delta']:rec[name+'_low']=np.nan;rec[name+'_high']=np.nan
 if n>=50 and nb>=12:
  rng=np.random.default_rng(123);ix=[np.flatnonzero(blocks==j)for j in range(nb)];B=F['bootstrap_reps'];bs=[]
  for rep in range(B):
   z=np.concatenate([ix[j]for j in rng.integers(0,nb,nb)]);bb=np.linalg.lstsq(X[z],Y[z],rcond=None)[0][-1];bs.append([bb[0],bb[1],bb[1]-bb[0]])
  bs=np.array(bs)
  for j,name in enumerate(['beta_bin','beta_cell','delta']):rec[name+'_low'],rec[name+'_high']=np.quantile(bs[:,j],[.025,.975])
  for j,b in enumerate(bs):boots.append(dict(variant=variant,block_um=block,rep=j,beta_bin=b[0],beta_cell=b[1],delta=b[2]))
  rec['status']='paired_within_section_exploratory'
  Z=X[:,:-1];resx=X[:,-1]-Z@np.linalg.lstsq(Z,X[:,-1],rcond=None)[0];rec['Treg_residual_variance']=float(np.var(resx))
  if do_p and np.var(resx)>.1:
   yy=Y[:,1]-Y[:,0];fit0=Z@np.linalg.lstsq(Z,yy,rcond=None)[0];res=yy-fit0;h=np.sum(Z*np.linalg.pinv(Z).T,axis=1);res/=np.sqrt(np.maximum(1-h,1e-6));weight=np.linalg.pinv(X)[-1];gg=np.bincount(blocks,weights=weight*res,minlength=nb);null=rng.choice([-1.,1.],size=(B,nb))@gg;rec['P_delta']=(1+(np.abs(null)>=abs(delta)).sum())/(B+1)
  elif do_p:rec['status']='descriptive_Treg_not_separable'
  if variant=='common':
   for j in range(nb):
    keep=blocks!=j;bb=np.linalg.lstsq(X[keep],Y[keep],rcond=None)[0][-1];v=dat.iloc[ix[j]];lobo.append(dict(block_um=block,left_out_block=j,removed_n=len(v),center_x=v.x_um.mean(),center_y=v.y_um.mean(),beta_bin=bb[0],beta_cell=bb[1],delta=bb[1]-bb[0]))
 return rec
for version in ['common','competition']:
 dd=pd.read_parquet(R/'tables'/f'design_r200_{version}.parquet').set_index('T_region_id');assert dd.index.equals(d.index);X=dd.to_numpy(float)
 for block in ([400,800]if version=='common'else[400]):rows.append(fit_pair(X,y,d,block,version,True))
 if version=='common':
  rows.append(fit_pair(X[~hot],y[~hot],d.loc[~hot],400,'exclude_original_hotspot',True));X0=np.column_stack([np.ones(len(d)),X[:,-1]]);rows.append(fit_pair(X0,y,d,400,'unadjusted',False))
res=pd.DataFrame(rows);res['test_family']=np.where(res.variant.eq('common')&res.block_um.eq(400),'single_primary_delta',np.where(res.variant.eq('unadjusted'),'descriptive_only','unified_auxiliary_delta'));res['BH_q_auxiliary']=np.nan;ix=res.index[res.test_family.eq('unified_auxiliary_delta')&res.P_delta.notna()];res.loc[ix,'BH_q_auxiliary']=multipletests(res.loc[ix,'P_delta'],method='fdr_bh')[1];put(res,'paired_effects.tsv');pd.DataFrame(boots).to_parquet(R/'tables/paired_bootstrap.parquet',index=False);put(pd.DataFrame(lobo),'leave_one_block_out.tsv')
put(pd.DataFrame([dict(rule='union original candidate rectangles256um; not512um display',excluded_T_centers=int(hot.sum()),remaining_T_centers=int((~hot).sum()),candidate_center_distance_um=float(np.hypot(roi.x0.iloc[0]-roi.x0.iloc[1],roi.y0.iloc[0]-roi.y0.iloc[1])))]),'hotspot_exclusion.tsv')
# Expose actual source-depth differences; identical names in v1 were not identical numbers.
put(pd.DataFrame([dict(n=len(d),depth_bin_median=d.depth_bin.median(),depth_cell_median=d.depth_cell.median(),depth_rank_rho=spearmanr(d.depth_bin,d.depth_cell).statistic,depth_values_identical=bool(np.array_equal(d.depth_bin,d.depth_cell)),raw_MP10_rho=spearmanr(d.MP10_bin,d.MP10_cell).statistic,same_common_design=True)]),'source_measurement_correspondence.tsv')
from analysis_stats import fit_effect
b=pd.read_parquet(V/'Cervilla/T32_neighborhoods.parquet');b=b[b.radius_um.eq(200)&b.source_rule.eq('exclude_shared')&b.association_eligible].set_index('T_region_id');ab=[]
for name,v in [('A_main_full',b),('B_common_centers_full',b.loc[d.index])]:ab.append(dict(step=name,unit='bin',**fit_effect(v,'Treg_main',block_um=400)))
for unit in ['bin','cell']:
 z=res[res.variant.eq('common')&res.block_um.eq(400)].iloc[0];ab.append(dict(step='C_identical_common_design',unit=unit,n=len(d),effect=z['beta_'+unit],low=z['beta_'+unit+'_low'],high=z['beta_'+unit+'_high'],P=np.nan,status='source_beta_CI_only; primary_test_is_delta'))
put(pd.DataFrame(ab),'stepwise_decomposition.tsv');print('PAIRED_MODELS_DONE',flush=True)
