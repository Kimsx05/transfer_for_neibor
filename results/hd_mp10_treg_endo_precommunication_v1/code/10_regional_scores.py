from common import *
import importlib.util,anndata as ad
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr
sp=importlib.util.spec_from_file_location('measure',R/'code/05_measure.py');mod=importlib.util.module_from_spec(sp);sp.loader.exec_module(mod);auc=mod.auc_sparse;historical=mod.auc_historical
gate();d=pd.read_parquet(R/'tables/Cervilla_cell_annotation.parquet');nodes=pd.read_parquet(A/'hd_epi_niche_dual_methods_v1/Cervilla/nodes.parquet',columns=['x_um','y_um','component']);_,near=cKDTree(nodes[['x_um','y_um']]).query(d[['x_um','y_um']]);d['component']=nodes.component.to_numpy()[near];e=d[d.mp10_epithelial_eligible].copy();e['gx']=np.floor(e.x_um/64).astype(int);e['gy']=np.floor(e.y_um/64).astype(int)
# Within-grid connected Epi-only groups prevent pooling across gaps >16um; boundaries are input/geometry supported, not a new histology segmentation.
mem=[]
for (comp,gx,gy),g in e.groupby(['component','gx','gy']):
 pairs=cKDTree(g[['x_um','y_um']]).query_pairs(16,output_type='ndarray');G=sparse.coo_matrix((np.ones(len(pairs)),(pairs[:,0],pairs[:,1])),shape=(len(g),len(g))).tocsr();_,cc=connected_components(G,directed=False)
 for k,grp in g.assign(cc=cc).groupby('cc'):
  rid=f'Cervilla_C{comp}_G{gx}_{gy}_{k}';mem.extend(dict(region_id=rid,cell_id=c,n_epi=len(grp),x_um=x,y_um=y,component=comp)for c,x,y in grp[['cell_id','x_um','y_um']].itertuples(index=False,name=None))
m=pd.DataFrame(mem);put(m,'epithelial_region_membership.tsv.gz');a=ad.read_h5ad(A/'segmentation_assignment_v2/Cervilla/cells.h5ad');ids=a.obs_names;cols=pd.Index(ids).get_indexer(m.cell_id);assert(cols>=0).all();regions=sorted(m.region_id.unique());ridx=pd.Index(regions).get_indexer(m.region_id);G=sparse.csr_matrix((np.ones(len(m),np.int64),(ridx,cols)),shape=(len(regions),a.n_obs));y=G@a.X;sy=a.var.gene_symbol.to_numpy();del a;assert int(y.sum())==int(d.set_index('cell_id').loc[m.cell_id,'nCount_RNA'].sum());reg=m.groupby('region_id').agg(n_epi=('cell_id','size'),x_um=('x_um','mean'),y_um=('y_um','mean'),component=('component','first')).loc[regions].reset_index();reg['nCount_RNA']=np.asarray(y.sum(1)).ravel();reg['nFeature_RNA']=np.asarray((y>0).sum(1)).ravel();reg=panel(y,reg,sy);ti=np.flatnonzero(np.isin(sy,TOP));reg['MP10_native_observed_seed123']=auc(y.indptr,y.indices,y.data,y.shape[1],ti,123,True);reg['MP10']=historical(y.indptr,y.indices,y.data,y.shape[1],ti);reg['score_eligible']=(reg.n_epi>=5)&(reg.MP10_detected>=5)&(reg.nFeature_RNA>=905);reg['score_reason']=np.where(reg.n_epi<5,'fewer_than_5_credible_Epi',np.where(reg.MP10_detected<5,'fewer_than_5_Top20_detected',np.where(reg.nFeature_RNA<905,'positive_gene_ranking_does_not_cover_auc_window','eligible')))
# Fixed-depth binomial thinning is a measurement sensitivity, not a new selected score.
rng=np.random.default_rng(123)
for name,prob in [('half',np.repeat(.5,len(reg))),('depth3000',np.minimum(1,3000/np.maximum(1,reg.nCount_RNA)))]:
 z=y.copy();z.data=rng.binomial(z.data.astype(int),np.repeat(prob,np.diff(z.indptr)));z.eliminate_zeros();reg['MP10_'+name]=historical(z.indptr,z.indices,z.data,z.shape[1],ti)
for seed in [124,125]:reg['MP10_seed'+str(seed)]=auc(y.indptr,y.indices,y.data,y.shape[1],ti,seed,True)
el=reg.score_eligible;quant=np.quantile(reg.loc[el,'MP10'],[.5,.9,.95],method='linear');reg['group_H50']=np.where(el,np.where(reg.MP10>quant[0],'H50','L50'),'not_scorable');reg['group_H10']=np.where(el,np.where(reg.MP10>quant[1],'H10',np.where(reg.MP10<=quant[0],'L50','middle')),'not_scorable');reg['group_H05']=np.where(el,np.where(reg.MP10>quant[2],'H05',np.where(reg.MP10<=quant[0],'L50','middle')),'not_scorable');reg.to_parquet(R/'tables/Cervilla_epithelial_regions.parquet',index=False)
val=[]
for k in ['MP10_half','MP10_depth3000','MP10_seed124','MP10_seed125','nCount_RNA','nFeature_RNA']:
 val.append(dict(score='region64um_connected_Epi_only',comparison=k,n=int(el.sum()),spearman=spearmanr(reg.loc[el,'MP10'],reg.loc[el,k]).statistic))
put(pd.DataFrame(val),'regional_MP10_validation.tsv');(R/'models/regional_score_definition.json').write_text(json.dumps(dict(grid_um=64,max_Epi_graph_link_um=16,genes=TOP,aucMaxRank=905,historical_tie_method='average ranks; exactly neighborhood_mp10_calibration_v5/code/common.py auc_parts',native_sensitivity_keepZeroesAsNA=True,min_Epi=5,min_detected_Top20=5,min_positive_genes=905,quantile_type=7,quantiles=quant.tolist(),n_regions=len(reg),n_score_eligible=int(el.sum()),raw_UMI_conservation=True,score_unit='region; never duplicate as independent single-cell observations'),indent=2));print(reg.groupby('score_reason').size());print(pd.DataFrame(val).to_string(index=False))
