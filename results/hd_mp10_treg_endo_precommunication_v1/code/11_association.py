from common import *
from scipy.stats import rankdata,spearmanr
from statsmodels.stats.multitest import multipletests
from datetime import datetime
assert (R/'niche_frozen.json').exists(),'Freeze niche before inspecting relationships'
reg=pd.read_parquet(R/'tables/Cervilla_epithelial_regions.parquet');mem=pd.read_csv(R/'tables/epithelial_region_membership.tsv.gz',sep='\t');formal=pd.read_parquet(R/'tables/Cervilla_cell_formal_metadata.parquet');audit=[];assoc=[];niche_res=[];high_res=[];relations=[]

def zscore(x):return (x-np.mean(x))/max(np.std(x),1e-12)
def infer(d,x,y='MP10',cov=None,blocksize=400):
 cov=cov or [];cols=[x,y]+cov;sel=d[cols].notna().all(axis=1)&np.isfinite(d[cols]).all(axis=1);v=d[sel].copy();block=(v.x_um//blocksize).astype(int).astype(str)+'_'+(v.y_um//blocksize).astype(int).astype(str);labs=sorted(block.unique());idx=[np.flatnonzero(block.to_numpy()==b)for b in labs];base=dict(n=len(v),n_blocks=len(labs),block_um=blocksize,nonzero_predictor=int((v[x]!=0).sum()),predictor=x,outcome=y,covariates=';'.join(cov))
 binary_support_ok=True
 if v[x].nunique()==2:
  by=v.groupby(x).size();byblocks=v.assign(_block=block).groupby(x)._block.nunique();binary_support_ok=bool(by.min()>=20 and byblocks.min()>=8);base['min_group_n']=int(by.min());base['min_group_blocks']=int(byblocks.min())
 if len(v)<50 or len(labs)<12 or v[x].nunique()<2 or v[y].nunique()<3 or not binary_support_ok:return {**base,'status':'descriptive_not_enough_independent_spatial_support','effect':spearmanr(v[x],v[y]).statistic if len(v)>3 else np.nan,'low':np.nan,'high':np.nan,'P':np.nan}
 X=np.column_stack([np.ones(len(v)),zscore(rankdata(v[x]))]+[zscore(v[c].to_numpy(float))for c in cov]);Y=zscore(rankdata(v[y]));beta=np.linalg.lstsq(X,Y,rcond=None)[0][1];nullX=np.delete(X,1,axis=1);nullfit=nullX@np.linalg.lstsq(nullX,Y,rcond=None)[0];residual=Y-nullfit;pinv=np.linalg.pinv(X);rng=np.random.default_rng(123);boots=[];nulls=[];bmap=pd.Index(labs).get_indexer(block)
 for b in range(499):
  take=np.concatenate([idx[k]for k in rng.integers(0,len(idx),len(idx))]);boots.append(np.linalg.lstsq(X[take],Y[take],rcond=None)[0][1]);sign=rng.choice([-1,1],len(idx));nulls.append((pinv@(nullfit+residual*sign[bmap]))[1])
 return {**base,'status':'within_section_spatial_block_inference','effect':beta,'low':np.quantile(boots,.025),'high':np.quantile(boots,.975),'P':(1+np.sum(np.abs(nulls)>=abs(beta)))/500,'effect_definition':'Spearman rho if unadjusted; standardized rank slope if adjusted','inference':'499 spatial-block pair bootstrap CI; restricted-null block sign wild bootstrap P; not patient replication'}

for radius in [50,200]:
 n=pd.read_parquet(R/'tables'/f'Cervilla_cell_neighborhood_r{radius}.parquet');m=mem.merge(n,on='cell_id',validate='one_to_one',suffixes=('','_center'));agg=m.groupby('region_id').agg(n_centers=('cell_id','size'),known_fraction=('known_fraction','mean'),observed_area_fraction=('observed_area_fraction','mean'),Endothelial_count=('count_Endothelial','mean'),observed_area_um2=('observed_area_um2','mean'),Endothelial_distance=('distance_Endothelial_observed_um','mean'),distance_censored_fraction=('distance_Endothelial_censored','mean'),eligible_center_fraction=('neighborhood_eligible','mean'))
 pcols=[c for c in n if c.startswith('P_')];ps=m.groupby('region_id')[pcols].mean();agg=agg.join(ps);agg['Endothelial_density']=agg.Endothelial_count/agg.observed_area_um2*1e6;v=reg.merge(agg,on='region_id',validate='one_to_one');v['radius_um']=radius;v['association_eligible']=v.score_eligible&(v.observed_area_fraction>=.5)&(v.eligible_center_fraction>=.8);v['log_nCount']=np.log1p(v.nCount_RNA);v['epi_marker_balance']=v.Epi_UMI/(v.Epi_UMI+v.Fib_UMI+v.T_UMI+v.Endo_UMI).clip(lower=1);v.loc[v.distance_censored_fraction>.2,'Endothelial_distance']=np.nan
 for mode,cov in [('unadjusted',[]),('conditional',['log_nCount','epi_marker_balance','observed_area_fraction'])]:
  for target in ['Endothelial_density','Endothelial_distance']:
   result=infer(v[v.association_eligible],target,cov=cov);result.update(radius_um=radius,model=mode,family='continuous_primary');assoc.append(result)
 # Additional vessel-distance conditioning for endothelial abundance, and conservative block-scale sensitivity.
 extra=infer(v[v.association_eligible],'Endothelial_density',cov=['log_nCount','epi_marker_balance','observed_area_fraction','Endothelial_distance']);extra.update(radius_um=radius,model='conditional_plus_vessel_distance',family='continuous_sensitivity');assoc.append(extra)
 ext=infer(v[v.association_eligible],'Endothelial_density',cov=['log_nCount','epi_marker_balance','observed_area_fraction'],blocksize=800);ext.update(radius_um=radius,model='800um_block_sensitivity',family='continuous_sensitivity');assoc.append(ext)
 ext2=infer(v[v.association_eligible],'Endothelial_density',cov=['log_nCount','epi_marker_balance','observed_area_fraction','known_fraction']);ext2.update(radius_um=radius,model='conditional_plus_annotation_support',family='continuous_sensitivity');assoc.append(ext2)
 ext3=infer(v[v.association_eligible],'Endothelial_distance',cov=['log_nCount','epi_marker_balance','observed_area_fraction'],blocksize=800);ext3.update(radius_um=radius,model='800um_block_sensitivity',family='continuous_sensitivity');assoc.append(ext3)
 # Prespecified MP10-tail contrasts are descriptive sensitivities, never replace continuous primary results.
 for group in ['group_H50','group_H10','group_H05']:
  for state,part in v[v.association_eligible].groupby(group):
   high_res.append(dict(radius_um=radius,group_definition=group,state=state,n=len(part),median_MP10=float(part.MP10.median()),median_Endothelial_density=float(part.Endothelial_density.median()),median_Endothelial_distance=float(part.Endothelial_distance.median())))
 # Region-level labels only summarize its pre-frozen Epi-center niche membership; mixed regions remain separate.
 for level in ['niche','niche_K2']:
  composition=pd.crosstab(m.region_id,m[level],normalize='index');composition=composition.drop(columns=[0],errors='ignore');best=composition.idxmax(axis=1);frac=composition.max(axis=1);v[level+'_region_summary']=v.region_id.map(best).where(v.region_id.map(frac)>=.6,0).fillna(0).astype(int)
  for k in sorted(v[level+'_region_summary'].unique()):
   if k==0:continue
   vv=v[v.association_eligible&v[level+'_region_summary'].gt(0)].copy();vv['niche_indicator']=vv[level+'_region_summary'].eq(k).astype(int)
   # Binary predictor: infer() accepts binary after explicit variance check below.
   result=infer(vv,'niche_indicator',cov=['log_nCount','epi_marker_balance','observed_area_fraction']);result.update(radius_um=radius,niche=int(k),level=level,model='conditional',family='niche_primary' if level=='niche' else 'niche_coarse_sensitivity',n_in=int(vv.niche_indicator.sum()),n_out=int((1-vv.niche_indicator).sum()));niche_res.append(result)
 v.to_parquet(R/'tables'/f'Cervilla_region_association_r{radius}.parquet',index=False);relations.append(v)
 # Explicit missing endpoints, not zero Treg abundance.
 for endpoint in ['Treg_density','Treg_over_CD4','Treg_over_T','Treg_distance','reverse_Treg_vs_Tconv_Tfh']:
  assoc.append(dict(radius_um=radius,predictor=endpoint,outcome='MP10',family='continuous_primary',model='not_tested',status='not_evaluable_Treg_identity_not_established',n=0,P=np.nan,effect=np.nan,low=np.nan,high=np.nan))
# FDR only within predefined families/model purpose. Unadjusted and conditional retained together.
a=pd.DataFrame(assoc);a['BH_q']=np.nan
for fam,g in a[a.P.notna()].groupby('family'):a.loc[g.index,'BH_q']=multipletests(g.P,method='fdr_bh')[1]
b=pd.DataFrame(niche_res);b['BH_q']=np.nan
for fam,g in b[b.P.notna()].groupby('family'):b.loc[g.index,'BH_q']=multipletests(g.P,method='fdr_bh')[1]
put(a,'association_results.tsv');put(b,'niche_MP10_association.tsv');put(pd.DataFrame(high_res),'MP10_group_sensitivity.tsv');(R/'models/association_execution.json').write_text(json.dumps(dict(started_after_niche_freeze=True,completed=datetime.now().isoformat(),MP10_source='historical mean-tie rank AUC on raw counts pooled only over credible Epi in fixed 64um connected regions',Treg='not_evaluable; candidate scores/markers retained; no Treg abundance imputation',spatial_unit='region, not copied score per cell'),indent=2));print(a.to_string(index=False));print(b.to_string(index=False))
