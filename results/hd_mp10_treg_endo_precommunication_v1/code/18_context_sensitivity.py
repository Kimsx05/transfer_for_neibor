"""One explicitly requested tissue-context sensitivity; no refitting of units, scores or niches."""
from common import *
from scipy.spatial import cKDTree
from scipy.stats import rankdata,spearmanr
from statsmodels.stats.multitest import multipletests
import ast
# Reuse exactly the tested statistical implementation without running the primary analysis again.
tree=ast.parse((R/'code/11_association.py').read_text());defs=ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef)],type_ignores=[]);exec(compile(defs,'11_association_functions','exec'))
b=pd.read_csv(A/'hd_bin8_banksy_v1/Cervilla/metadata.tsv.gz',sep='\t');reg=pd.read_parquet(R/'tables/Cervilla_epithelial_regions.parquet');dist,ix=cKDTree(b[['x_um','y_um']]).query(reg[['x_um','y_um']],workers=1);context=reg[['region_id']].copy();context['context_bin_id']=b.bin_id.to_numpy()[ix];context['context_distance_um']=dist;context['existing_l08_domain']=np.where(dist<=8,b.banksy_l08_cluster.to_numpy()[ix].astype(str),'unmapped');context['mapping_purpose']='spatial-context covariate at regional centroid; not cell/bin identity matching or cell-type transfer'
mem=pd.read_csv(R/'tables/epithelial_region_membership.tsv.gz',sep='\t');f=pd.read_parquet(R/'tables/Cervilla_cell_formal_metadata.parquet');edge=mem[['region_id','cell_id']].merge(f[['cell_id','distance_to_observed_boundary_um']],on='cell_id').groupby('region_id').distance_to_observed_boundary_um.mean();context['mean_boundary_distance_um']=context.region_id.map(edge);put(context,'regional_existing_domain_context.tsv')
res=[]
for radius in [50,200]:
 v=pd.read_parquet(R/'tables'/f'Cervilla_region_association_r{radius}.parquet').merge(context,on='region_id',validate='one_to_one');v=v[v.association_eligible].copy();counts=v.existing_l08_domain.value_counts();v['domain_covariate']=v.existing_l08_domain.where(v.existing_l08_domain.map(counts)>=30,'small_or_unmapped');dm=pd.get_dummies(v.domain_covariate,prefix='domain',drop_first=True,dtype=float);v=pd.concat([v,dm],axis=1);v['log_boundary_distance']=np.log1p(v.mean_boundary_distance_um);cov=['log_nCount','epi_marker_balance','known_fraction','observed_area_fraction','log_boundary_distance']+dm.columns.tolist()
 for target in ['Endothelial_density','Endothelial_distance']:
  z=infer(v,target,cov=cov);z.update(radius_um=radius,model='existing_BANKSY_context_and_edge_sensitivity',family='tissue_context_sensitivity',domain_terms=len(dm.columns));res.append(z)
put(pd.DataFrame(res).assign(BH_q=multipletests([r['P'] for r in res],method='fdr_bh')[1]),'association_tissue_context.tsv');print(pd.DataFrame(res).to_string(index=False))
