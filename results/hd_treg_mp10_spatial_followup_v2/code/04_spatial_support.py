from fu_common import *
b=pd.read_parquet(V/'Cervilla/bin_evidence.parquet');t=pd.read_parquet(V/'Cervilla/T32_regions.parquet');e=pd.read_parquet(V/'Cervilla/Epi64_regions.parquet');e=e[e.measurement_eligible].set_index('region_id');em=pd.read_parquet(V/'Cervilla/Epi64_members.parquet');eb={k:set(v.bin_id)for k,v in em[em.region_id.isin(e.index)].groupby('region_id')}
c=pd.read_parquet(OLD/'tables/Cervilla_epithelial_regions.parquet');c=c[c.score_eligible].copy();cm=pd.read_csv(OLD/'tables/epithelial_region_membership.tsv.gz',sep='\t');cm=cm[cm.region_id.isin(c.region_id)].copy();cm['cell_label']=cm.cell_id.str.rsplit('-',n=1).str[-1].astype(int)
q=pd.read_parquet(A/'segmentation_assignment_v2/Cervilla/bin_to_cell.parquet',columns=['array_row','array_col','cell_label']);q=q[q.cell_label.isin(cm.cell_label)].merge(cm[['cell_label','region_id']],on='cell_label',validate='many_to_one');q['array_row']=q.array_row//4;q['array_col']=q.array_col//4
src=q[['region_id','array_row','array_col']].drop_duplicates().merge(b[['bin_id','array_row','array_col']],on=['array_row','array_col'],validate='many_to_one');src[['region_id','bin_id']].drop_duplicates().to_parquet(R/'tables/Cervilla_cell_region_source8.parquet',index=False);cm.to_parquet(R/'tables/Cervilla_cell_region_source_cells.parquet',index=False);cb={k:set(v.bin_id)for k,v in src.groupby('region_id')}
dist,j=cKDTree(b[['x_um','y_um']]).query(c[['x_um','y_um']]);c['component_current']=b.component.iloc[j].to_numpy();c['component_match_distance_um']=dist;c['coordinate_supported']=dist<=16;c.to_parquet(R/'tables/Cervilla_cell_regions.parquet',index=False);c=c.set_index('region_id');shared=pd.read_parquet(V/'Cervilla/cell_region_shared_source.parquet');shared=shared.groupby('T_region_id').region_id.apply(set).to_dict();binpairs=pd.read_parquet(V/'Cervilla/region_pairs.parquet');binpairs=binpairs[binpairs.T_grid_um.eq(32)&~binpairs.shared_original_bins];rows=[];pairs=[];support=[];bc=b.set_index('bin_id');common_attempt=[]
for rad in [50,200]:
 tree=cKDTree(c[['x_um','y_um']]);near=tree.query_ball_point(t.loc[t.measurement_eligible,['x_um','y_um']],rad)
 for rr,js in zip(t[t.measurement_eligible].itertuples(),near):
  cc=c.iloc[js];cc=cc[cc.coordinate_supported&cc.component_current.eq(rr.component)];raw_ids=list(cc.index);cc=cc[~cc.index.isin(shared.get(rr.region_id,set()))];ei=list(binpairs.loc[binpairs.radius_um.eq(rad)&binpairs.T_region_id.eq(rr.region_id),'Epi_region_id']);ci=list(cc.index)
  for eid in raw_ids:pairs.append(dict(T_region_id=rr.region_id,radius_um=rad,cell_region_id=eid,distance_um=float(np.hypot(c.loc[eid,'x_um']-rr.x_um,c.loc[eid,'y_um']-rr.y_um)),shared_original_bins=eid in shared.get(rr.region_id,set())))
  bs=set().union(*(eb[i]for i in ei))if ei else set();cs=set().union(*(cb.get(i,set())for i in ci))if ci else set();both=bs&cs
  z=dict(T_region_id=rr.region_id,radius_um=rad,bin_regions=len(ei),cell_regions=len(ci),bin_source_positions=len(bs),cell_source_positions=len(cs),common_positions=len(both),bin_only_positions=len(bs-cs),cell_only_positions=len(cs-bs),Jaccard=len(both)/len(bs|cs)if bs|cs else np.nan,bin_MP10=float(e.loc[ei,'MP10'].mean())if ei else np.nan,cell_MP10=float(c.loc[ci,'MP10'].mean())if ci else np.nan,bin_depth=float(e.loc[ei,'nCount'].mean())if ei else np.nan,cell_depth=float(c.loc[ci,'nCount_RNA'].mean())if ci else np.nan)
  for name,ss in [('bin',bs),('cell',cs),('common',both)]:
   v=bc.loc[list(ss)]if ss else bc.iloc[:0]
   for field,fun in [('x_um','min'),('x_um','max'),('y_um','min'),('y_um','max')]:z[f'{name}_{field}_{fun}']=getattr(v[field],fun)()if len(v)else np.nan
  rows.append(z)
  if rad==200:
   # Strict common8um footprint using only intact preexisting score regions. No cell count slicing.
   available=both;bi=ei;ci2=ci
   while True:
    bi=[i for i in bi if eb[i] and eb[i]<=available];ci2=[i for i in ci2 if cb.get(i,set()) and cb[i]<=available];nb=set().union(*(eb[i]for i in bi))if bi else set();nc=set().union(*(cb[i]for i in ci2))if ci2 else set();new=nb&nc
    if new==available:break
    available=new
   common_attempt.append(dict(T_region_id=rr.region_id,common_intact_bin_regions=len(bi),common_intact_cell_regions=len(ci2),common_positions=len(available),eligible=bool(bi and ci2),bin_MP10=float(e.loc[bi,'MP10'].mean())if bi else np.nan,cell_MP10=float(c.loc[ci2,'MP10'].mean())if ci2 else np.nan,rule='iterated intersection footprint; retain only whole score regions, no cell slicing or score replication'))
pd.DataFrame(pairs).to_parquet(R/'tables/Cervilla_cell_region_pairs.parquet',index=False);put(pd.DataFrame(rows),'Cervilla_source_coverage.tsv');put(pd.DataFrame(common_attempt),'common_support_feasibility.tsv')
common=pd.read_parquet(R/'tables/Cervilla_common_centers_r200.parquet');q=pd.DataFrame(rows);q=q[q.radius_um.eq(200)].set_index('T_region_id').loc[common.T_region_id];assert np.allclose(q.bin_MP10,common.MP10_bin)and np.allclose(q.cell_MP10,common.MP10_cell)
put(pd.DataFrame([dict(n_common_T_centers=len(q),median_Jaccard=q.Jaccard.median(),q10_Jaccard=q.Jaccard.quantile(.1),q90_Jaccard=q.Jaccard.quantile(.9),exact_equal_footprints=int((q.bin_only_positions.eq(0)&q.cell_only_positions.eq(0)).sum()),strict_common_whole_region_centers=int(pd.DataFrame(common_attempt).eligible.sum()),scope='same T centers is not same epithelial support; footprints on8um source grid, not histological area')]),'source_coverage_summary.tsv');print('SUPPORT_DONE',flush=True)
