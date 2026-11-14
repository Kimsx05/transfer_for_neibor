"""One bounded target-purity diagnostic. Original annotation/niche/primary tests unchanged."""
from common import *
from scipy.spatial import cKDTree
from scipy.stats import rankdata,spearmanr
from statsmodels.stats.multitest import multipletests
import ast
tr=ast.parse((R/'code/11_association.py').read_text());exec(compile(ast.Module(body=[n for n in tr.body if isinstance(n,ast.FunctionDef)],type_ignores=[]),'association_functions','exec'))
d=pd.read_parquet(R/'tables/Cervilla_cell_formal_metadata.parquet');endo=d[d.type_supported.eq('Endothelial')].copy();endo['raw_Epi_mixed']=endo.Epi_detected>=2;endo['raw_Fib_mixed']=endo.Fib_detected>=2;endo['sensitivity_keep']=~(endo.raw_Epi_mixed|endo.raw_Fib_mixed);put(endo[['cell_id','x_um','y_um','Endo_detected','Epi_detected','Fib_detected','raw_Epi_mixed','raw_Fib_mixed','sensitivity_keep']],'Endothelial_mixture_diagnostic.tsv.gz');targets=endo[endo.sensitivity_keep];mem=pd.read_csv(R/'tables/epithelial_region_membership.tsv.gz',sep='\t');out=[]
for radius in [50,200]:
 n=pd.read_parquet(R/'tables'/f'Cervilla_cell_neighborhood_r{radius}.parquet');ct=np.zeros(len(n));near=np.full(len(n),np.nan)
 for c in n.component.unique():
  ti=targets[targets.component==c];ix=np.flatnonzero(n.component==c)
  if len(ti):
   tree=cKDTree(ti[['x_um','y_um']]);q=n.iloc[ix][['x_um','y_um']];ct[ix]=tree.query_ball_point(q,radius+1e-7,return_length=True,workers=1);near[ix]=tree.query(q,workers=1)[0]
 n['unmixed_count']=ct;n['unmixed_distance_censored']=~np.isfinite(near)|(near>n.distance_to_observed_boundary_um);n['unmixed_distance']=np.where(n.unmixed_distance_censored,np.nan,near);a=mem[['region_id','cell_id']].merge(n[['cell_id','unmixed_count','unmixed_distance','unmixed_distance_censored']],on='cell_id',validate='one_to_one').groupby('region_id').agg(unmixed_count=('unmixed_count','mean'),unmixed_distance=('unmixed_distance','mean'),unmixed_censored_fraction=('unmixed_distance_censored','mean'));v=pd.read_parquet(R/'tables'/f'Cervilla_region_association_r{radius}.parquet').merge(a,on='region_id');v['Endothelial_density_unmixed']=v.unmixed_count/v.observed_area_um2*1e6;v['Endothelial_distance_unmixed']=v.unmixed_distance.where(v.unmixed_censored_fraction<=.2)
 for target in ['Endothelial_density_unmixed','Endothelial_distance_unmixed']:
  z=infer(v[v.association_eligible],target,cov=['log_nCount','epi_marker_balance','observed_area_fraction']);z.update(radius_um=radius,model='exclude_Endo_with_two_Epi_or_Fib_markers',family='target_mixture_sensitivity',n_target_original=len(endo),n_target_sensitivity=len(targets));out.append(z)
 put(v[['region_id','association_eligible','Endothelial_density_unmixed','Endothelial_distance_unmixed','unmixed_censored_fraction']],f'Endo_mixture_region_sensitivity_r{radius}.tsv.gz')
z=pd.DataFrame(out);z['BH_q']=np.nan;ix=z.P.notna();z.loc[ix,'BH_q']=multipletests(z.loc[ix,'P'],method='fdr_bh')[1];put(z,'association_Endo_mixture_sensitivity.tsv');print(z.to_string(index=False))
