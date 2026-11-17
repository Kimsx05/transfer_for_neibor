from common import *
import ast,anndata as ad
from scipy.spatial import distance_matrix
cfg=yaml.safe_load((R/'frozen_config.yaml').read_text());gate()
axes=[('TIGIT_NECTIN2','TIGIT','NECTIN2'),('CD96_NECTIN1','CD96','NECTIN1')]
put(pd.DataFrame([dict(axis=a,DB_version='CellChatDB.human_1.6.1',DB_ligand=t,DB_receptor=e,source='mp10_cellchat_high_definition_sensitivity_v1/RESULTS.md H05 Treg input table',DB_sender='Treg',DB_receiver='Epi_High',interpretation='database encoding only; receptor-bearing T and epithelial binding partner do not prove direction; no LR probability calculated',selected_before_HD_association=True)for a,t,e in axes]),'candidate_axes_frozen.tsv')
# Reuse the authors native distance helper only. This avoids pinned obsolete imports and unused Gi* machinery.
source=(R/'sources/spottedpy_hotspot_helper.py').read_text();tree=ast.parse(source);f=next(n for n in tree.body if isinstance(n,ast.FunctionDef)and n.name=='calculateDistancesHelper');code=ast.Module(body=[f],type_ignores=[]);ns={'np':np,'pd':pd,'distance_matrix':distance_matrix};exec(compile(code,'SpottedPy_original_calculateDistancesHelper','exec'),ns);distancefn=ns['calculateDistancesHelper']
summary=[];hot=[];readiness=[]
for s in SECTIONS:
 t=pd.read_parquet(R/s/'T32_regions.parquet').set_index('region_id');e=pd.read_parquet(R/s/'Epi64_regions.parquet').set_index('region_id');pairs=pd.read_parquet(R/s/'region_pairs.parquet');pp=pairs[(pairs.T_grid_um==32)&~pairs.shared_original_bins]
 for axis,tg,eg in axes:
  for radius in [50,200]:
   p=pp[pp.radius_um==radius].copy();p['T_gene_detected']=t.loc[p.T_region_id,'raw_'+tg].to_numpy()>0;p['Epi_gene_detected']=e.loc[p.Epi_region_id,'raw_'+eg].to_numpy()>0;p['Treg_like']=t.loc[p.T_region_id,'evidence'].eq('Treg_like_multigene').to_numpy();p['Epi_H50']=e.loc[p.Epi_region_id,'H50'].to_numpy();p['Epi_L50']=e.loc[p.Epi_region_id,'L50'].to_numpy();p['Epi_H10']=e.loc[p.Epi_region_id,'H10'].to_numpy();p['Epi_H05']=e.loc[p.Epi_region_id,'H05'].to_numpy()
   p['axis']=axis;p['section']=s;p.to_parquet(R/s/f'candidate_{axis}_r{radius}.parquet',index=False)
   for grp in ['all','Epi_L50','Epi_H50','Epi_H10','Epi_H05']:
    v=p if grp=='all' else p[p[grp]];supported=v[v.T_gene_detected&v.Epi_gene_detected];like=supported[supported.Treg_like]
    summary.append(dict(section=s,axis=axis,radius_um=radius,Epi_group=grp,source_rule='no_shared_8um_bins',region_pairs=len(v),both_genes_detected_pairs=len(supported),T_regions_with_neighbor_expression=supported.T_region_id.nunique(),Epi_regions_with_neighbor_expression=supported.Epi_region_id.nunique(),Treg_like_supported_T_regions=like.T_region_id.nunique(),pair_fraction=len(supported)/len(v)if len(v)else np.nan))
  a=[v for v in summary if v['section']==s and v['axis']==axis and v['radius_um']==200 and v['Epi_group']=='all'][0]
  readiness.append(dict(section=s,axis=axis,category='regional_candidate_axis_support'if a['Treg_like_supported_T_regions']>=3 else 'general_T_lineage_only'if a['T_regions_with_neighbor_expression']>=3 else 'insufficient_target_signal',credible_Treg_cell_identity=False,formal_communication_ready=False,reason='regional signals do not establish sender identity or signaling; no ligand-receptor score fitted'))
 # SpottedPy nearest-distance helper: predefined evidence hotspots and MP10 groups; no Gi* claimed.
 # Batch is an independently verified input component. Coordinates supplied in microns.
 b=pd.read_parquet(R/s/'bin_evidence.parquet');p=b[b.eligible_pool];grid=np.zeros((int(b.array_row.max())+3,int(b.array_col.max())+3),bool);grid[p.array_row.to_numpy(int)+1,p.array_col.to_numpy(int)+1]=True;edge=ndimage.distance_transform_edt(grid)*8
 edgeby=p.copy();edgeby['edge_um']=edge[p.array_row.to_numpy(int)+1,p.array_col.to_numpy(int)+1];edge_tree=cKDTree(edgeby[['x_um','y_um']]);edgedist=edgeby.edge_um.to_numpy()
 for comp,tt in t[t.measurement_eligible].groupby('component'):
  ee=e[e.measurement_eligible&(e.component==comp)];tt=tt[tt.evidence.eq('Treg_like_multigene')]
  if not len(tt):continue
  # Bounded exact evaluation: at most64 primary positions per call, no global N x N matrix.
  for start in range(len(tt)):
   u=tt.iloc[start:start+1];excluded=pairs.loc[(pairs.T_grid_um==32)&pairs.T_region_id.eq(u.index[0])&pairs.shared_original_bins,'Epi_region_id'];ee0=ee;ee=ee0.loc[~ee0.index.isin(excluded)];obs=pd.concat([u[['x_um','y_um']],ee[['x_um','y_um']]]);obs['array_row']=obs.y_um;obs['array_col']=obs.x_um;obs['Treg_evidence_hotspot']=np.nan;obs.loc[u.index,'Treg_evidence_hotspot']=1.
   for g in ['H50','H10','L50']:obs['Epi_'+g]=np.nan;obs.loc[ee.index[ee[g]],'Epi_'+g]=1.
   aa=ad.AnnData(obs=obs);z=distancefn(aa,['Treg_evidence_hotspot'],['Epi_H50','Epi_H10','Epi_L50'],str(comp),False,False)
   # Missing target produces no native rows; add explicit unavailable rows, never maximum-distance fill.
   for rid in u.index:
    for g in ['H50','H10','L50']:
     v=z[(z.primary_index==rid)&z.comparison_variable.eq('Epi_'+g)]if len(z)else pd.DataFrame();dist=float(v.min_distance.iloc[0])if len(v)else np.nan;pos=u.loc[rid,['x_um','y_um']].to_numpy(float);bound=float(edgedist[edge_tree.query(pos)[1]]);hot.append(dict(section=s,component=comp,T_region_id=rid,comparison='Epi_'+g,native_distance_um=dist,boundary_distance_um=bound,observable=bool(np.isfinite(dist)and dist<=bound),distance_um=dist if np.isfinite(dist)and dist<=bound else np.nan,status='observed'if np.isfinite(dist)and dist<=bound else 'edge_censored'if np.isfinite(dist)else 'no_target_in_component',source_rule='exclude_any_shared_original_8um_bins',method='SpottedPy0.1.7 native calculateDistancesHelper one source at a time; no Gi*; input-component isolated'))
   ee=ee0
 print('CANDIDATES_SPOTTED',s,flush=True)
put(pd.DataFrame(summary),'candidate_expression_support.tsv');put(pd.DataFrame(readiness),'communication_readiness.tsv');put(pd.DataFrame(hot),'spottedpy_hotspot_distances.tsv.gz')
(R/'sources/SpottedPy_execution.json').write_text(json.dumps({'version_from_setup':'0.1.7','tree_sha':json.loads((R/'sources/spottedpy_tree.json').read_text()).get('sha'),'function':'calculateDistancesHelper','implementation':'unaltered function AST from source; only numpy/pandas/scipy dependencies; one primary per call with shared-source targets removed; no native Gi* hotspot fit; frozen evidence/MP10 indicators as positions','coordinates':'physical um; array_row and array_col argument fields repurposed as y/x um explicitly','empty_hotspot_default_to_max_distance':False,'split_by_input_component':True},indent=2))
