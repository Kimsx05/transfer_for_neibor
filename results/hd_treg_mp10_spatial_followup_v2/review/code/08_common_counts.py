from fu_common import *
import anndata as ad
# Geometry is fixed before evaluating new scores or effects. This is a bounded reconstruction,
# following the documented failure to preserve both sets of entire historical score regions.
config=dict(frozen_at=stamp(),scope='Cervilla 200um common T centers; exploratory common-support sensitivity',cell_rule='retain entire trusted source cells only when every assigned 2um position maps inside the common8um source footprint, one global64um grid; merge touching cell footprints into connected scoring units',bin_rule='sum each supplied8um bin in retained cell footprints once per scoring unit',score_rule='historical18085 background, average ranks K905; both sources nFeature>=905, MP10 detected>=5, epithelial detected>=2',neighborhood='preexisting200um region neighbors; exclude all shared T32 source bins',statistics='center averages across paired qualifying score units; no new formal P value; same covariate framework including both reconstructed depths',limitations='equal8um footprints do not imply equal RNA sources; whole cells are not divided; no new segmentation')
p=R/'common_support_reconstruction.yaml'
if not p.exists():p.write_text(yaml.safe_dump(config,sort_keys=False))
b=pd.read_parquet(V/'Cervilla/bin_evidence.parquet');bc=b.set_index('bin_id');em=pd.read_parquet(V/'Cervilla/Epi64_members.parquet');eb={k:set(v.bin_id)for k,v in em.groupby('region_id')};cb=pd.read_parquet(R/'tables/Cervilla_cell_region_source8.parquet');cb={k:set(v.bin_id)for k,v in cb.groupby('region_id')}
cm=pd.read_parquet(R/'tables/Cervilla_cell_region_source_cells.parquet');cellregion=cm.set_index('cell_label').region_id.to_dict();cellids=cm.set_index('cell_label').cell_id.to_dict();q=pd.read_parquet(A/'segmentation_assignment_v2/Cervilla/bin_to_cell.parquet',columns=['array_row','array_col','cell_label']);q=q[q.cell_label.isin(cm.cell_label)].copy();q['array_row']=q.array_row//4;q['array_col']=q.array_col//4;q=q.drop_duplicates();q=q.merge(b[['bin_id','array_row','array_col']],on=['array_row','array_col'],how='left',validate='many_to_one');missing=set(q.loc[q.bin_id.isna(),'cell_label']);foot={k:set(v.bin_id)for k,v in q[~q.cell_label.isin(missing)].groupby('cell_label')};grid={}
for k,v in q[~q.cell_label.isin(missing)].groupby('cell_label'):
 g=set(zip(v.array_row//8,v.array_col//8))
 if len(g)==1:grid[k]=next(iter(g))
perregion={k:list(v.cell_label)for k,v in cm.groupby('region_id')};bp=pd.read_parquet(V/'Cervilla/region_pairs.parquet');bp=bp[bp.T_grid_um.eq(32)&bp.radius_um.eq(200)&~bp.shared_original_bins];cp=pd.read_parquet(R/'tables/Cervilla_cell_region_pairs.parquet');cp=cp[cp.radius_um.eq(200)&~cp.shared_original_bins];cent=pd.read_parquet(R/'tables/Cervilla_common_centers_r200.parquet');tm=pd.read_parquet(V/'Cervilla/T32_members.parquet');tm={k:set(v.bin_id)for k,v in tm.groupby('region_id')};units=[];membership=[]
for rr in cent.itertuples():
 ei=bp.loc[bp.T_region_id.eq(rr.T_region_id),'Epi_region_id'];ci=cp.loc[cp.T_region_id.eq(rr.T_region_id),'cell_region_id'];bs=set().union(*(eb[i]for i in ei));cs=set().union(*(cb.get(i,set())for i in ci));common=(bs&cs)-tm.get(rr.T_region_id,set());groups={}
 for cr in ci:
  for k in perregion[cr]:
   if k in grid and foot[k]<=common:groups.setdefault(grid[k],[]).append(k)
 for g,allks in groups.items():
  # Distinct historical cell score regions can share an8um source bin. Merge touching
  # whole-cell footprints within this fixed grid before summing bin counts exactly once.
  parent={k:k for k in allks};owner={}
  def root(k):
   while parent[k]!=k:parent[k]=parent[parent[k]];k=parent[k]
   return k
  for k in allks:
   for bid in foot[k]:
    br,co=bc.loc[bid,['array_row','array_col']].astype(int)
    for dr in [-1,0,1]:
     for dc in [-1,0,1]:
      other=owner.get((br+dr,co+dc))
      if other is not None:parent[root(k)]=root(other)
    owner[(br,co)]=k
  connected={}
  for k in allks:connected.setdefault(root(k),[]).append(k)
  for j,ks in enumerate(connected.values()):
   bins=set().union(*(foot[k]for k in ks));uid=f'{rr.T_region_id}|{g[0]}_{g[1]}|{j}'
   units.append(dict(unit_id=uid,T_region_id=rr.T_region_id,source_cell_region=';'.join(sorted(set(cellregion[k]for k in ks))),cells=ks,bins=sorted(bins),n_cells=len(ks),n_bins=len(bins)))
   membership.extend([dict(unit_id=uid,source='whole_cell',member_id=cellids[k])for k in ks]);membership.extend([dict(unit_id=uid,source='supplied_bin8',member_id=x)for x in bins])
print('GEOMETRIC_UNITS',len(units),len(set(u['T_region_id']for u in units)),flush=True)
members=pd.DataFrame(membership);members['T_region_id']=members.unit_id.str.split('|').str[0];assert not members.duplicated(['T_region_id','source','member_id']).any()
members.to_parquet(R/'tables/common_support_reconstruction_members.parquet',index=False)
x,d,ge,sy,_=v1.load8('Cervilla');di=pd.Series(np.arange(len(d)),index=d.bin_id);a=ad.read_h5ad(A/'segmentation_assignment_v2/Cervilla/cells.h5ad',backed='r');assert list(a.var_names)==list(ge);labels=sorted(set(k for u in units for k in u['cells']));ix=np.array([a.obs_names.get_loc(cellids[k])for k in labels]);order=np.argsort(ix);cx=a.X[ix[order]].tocsr();labelrow={labels[o]:j for j,o in enumerate(order)};assert np.all(cx.data==np.floor(cx.data));a.file.close()
selgenes=v1.TOP+v1.PANELS['Epi'];sel=np.array([np.flatnonzero(sy==g)[0]for g in selgenes]);rows=[]
for u in units:
 bi=di.loc[u['bins']].to_numpy();ci=[labelrow[k]for k in u['cells']];ys=[sparse.csr_matrix(x[bi].sum(0)),sparse.csr_matrix(cx[ci].sum(0))];z={k:v for k,v in u.items()if k not in ['cells','bins']}
 for source,y,expected in zip(['bin','cell'],ys,[x[bi].sum(),cx[ci].sum()]):
  y.eliminate_zeros();assert y.sum()==expected;r=v1.selected_ranks(y.indptr,y.indices,y.data,18085,sel)[0];raw=y[:,sel].toarray()[0];mp=float(np.maximum(905-r[:20],0).sum()/sum(905-np.arange(1,21)));nf=y.nnz;nd=int((raw[:20]>0).sum());ne=int((raw[20:]>0).sum());z.update({source+'_MP10':mp,source+'_nCount':int(y.sum()),source+'_nFeature':nf,source+'_MP10_detected':nd,source+'_Epi_detected':ne,source+'_eligible':nf>=905 and nd>=5 and ne>=2})
 z['paired_eligible']=z['bin_eligible']and z['cell_eligible'];rows.append(z)
z=pd.DataFrame(rows);put(z,'common_support_reconstructed_units.tsv');good=z[z.paired_eligible];out=good.groupby('T_region_id').agg(MP10_bin=('bin_MP10','mean'),MP10_cell=('cell_MP10','mean'),depth_bin=('bin_nCount','mean'),depth_cell=('cell_nCount','mean'),n_paired_units=('unit_id','size')).reset_index();out=cent.drop(columns=['MP10_bin','MP10_cell','depth_bin','depth_cell']).merge(out,on='T_region_id');out.to_parquet(R/'tables/common_support_reconstructed_centers.parquet',index=False)
summary=dict(geometric_units=len(z),paired_eligible_units=len(good),common_centers=len(cent),paired_eligible_centers=len(out),excluded_geometry_cells=len(cm)-len(labels),input_genes=18085,integer_counts=True,subset_sums_conserved=True,formal_P_computed=False,rule='all original cell counts retained, equal8um source footprint, no area splitting')
if len(out)>=20:
 nuis=['T_lineage','T_nCount','observed_circle_fraction','T_background_balance','depth_bin','depth_cell'];X=np.ones((len(out),1));kept=[]
 for col in nuis:
  v=zrank(out[col]);res=v-X@np.linalg.lstsq(X,v,rcond=None)[0]
  if np.var(res)>.1:X=np.column_stack([X,v]);kept.append(col)
 X=np.column_stack([X,zrank(out.Treg_main)]);Y=np.column_stack([zrank(out.MP10_bin),zrank(out.MP10_cell)]);beta=np.linalg.lstsq(X,Y,rcond=None)[0][-1];summary.update(beta_bin=beta[0],beta_cell=beta[1],delta_beta=beta[1]-beta[0],covariates=','.join(kept))
else:summary['association_status']='insufficient jointly scoreable reconstructed centers; descriptive only'
put(pd.DataFrame([summary]),'common_support_reconstruction_summary.tsv');print(summary,flush=True)
# Same fixed rank/paired block bootstrap used in the main source comparison; descriptive CIs only.
if len(out)>=50:
 pd.DataFrame(X,columns=['intercept']+kept+['Treg_main'],index=out.T_region_id).reset_index().to_parquet(R/'tables/common_support_design.parquet',index=False)
 res=[]
 for block in [400,800]:
  bb=pd.factorize(pd.MultiIndex.from_arrays([out.component,np.floor(out.x_um/block),np.floor(out.y_um/block)]))[0];nb=bb.max()+1;ix=[np.flatnonzero(bb==i)for i in range(nb)];rng=np.random.default_rng(123);vals=[]
  for _ in range(999):
   jj=np.concatenate([ix[i]for i in rng.integers(0,nb,nb)]);be=np.linalg.lstsq(X[jj],Y[jj],rcond=None)[0][-1];vals.append([be[0],be[1],be[1]-be[0]])
  vals=np.array(vals);rec=dict(block_um=block,n=len(out),n_blocks=nb,P=np.nan,scope='common8um support whole-cell reconstruction; descriptive paired CI')
  for j,key in enumerate(['beta_bin','beta_cell','delta']):rec[key]=[beta[0],beta[1],beta[1]-beta[0]][j];rec[key+'_low'],rec[key+'_high']=np.quantile(vals[:,j],[.025,.975])
  res.append(rec)
 put(pd.DataFrame(res),'common_support_paired_effects.tsv')
