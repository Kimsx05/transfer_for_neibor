from base6 import *
def run(s):
 o=R/s; d=pd.read_parquet(o/'cells_geometry.parquet');assert d[['x','y']].notna().all().all()
 ov=pd.read_parquet(o/'nuclear_overlap.parquet');old=pd.read_parquet(V5/s/'cells_scored.parquet',columns=['id','pred','mixed','pTreg','group'])
 parent=ov.sort_values('intersection').drop_duplicates('new_id',keep='last').merge(old,left_on='old_id',right_on='id').set_index('new_id')
 for k in ['old_id','pred','mixed','pTreg','group']:d['v5_'+k]=d.id.map(parent[k])
 # Only flag coalesced instances where two substantial old cores each predominantly occupy this mask.
 combined=ov[(ov.overlap_old>=.65)&(ov.overlap_new>=.2)].groupby('new_id').size()
 multi=d.id.map(combined).fillna(0)>=2
 ap=pd.read_parquet(o/'assignment.parquet');raw=pd.read_parquet(o/'raw_signal.parquet')
 amb=raw[raw.ambiguous & raw[direct+core].sum(axis=1).gt(0)].merge(ap[['source_bin_index','nearest_id','runner_id']],on='source_bin_index')
 d['adjacent_ambiguous_key_bins']=d.id.map(pd.concat([amb.nearest_id,amb.runner_id]).value_counts()).fillna(0).astype(int)
 # Retained key observations can have independent in-nucleus support. Otherwise adjacent contested signal flagged for review.
 inkey=raw[raw.id.gt(0)&raw.id.eq(raw.nucleus_id)&raw[direct+core].sum(axis=1).gt(0)].groupby('id').size()
 disputed=(d.adjacent_ambiguous_key_bins>0)&~d.id.isin(inkey.index)&(d.T_width>0)
 if "foreign_key_bins" in d:disputed=disputed|d.foreign_key_bins.gt(0)|~d.proseg_modelled
 d['competing_RNA']=d[epi+fib+mye].sum(axis=1)>0
 d['mixed_v6']=np.select([multi,disputed,(d.major!='T')&(d.major_support>=.70)&(d.p_T<.3)&(d.T_width>0),d.competing_RNA&(d.T_width>0)],['unresolved_doublet','assignment_uncertain','nonT_dominant','minor_background'],default='none')
 newhint=(d.nCount>0)&(d[T+reg].sum(axis=1).gt(0)|d.p_T.ge(.3)|d.major.eq('T'))
 cand=pd.read_parquet(V5/s/'candidate_ledger.parquet')[['unit','sources','x','y']];new=d.loc[newhint,['id','x','y']].copy();new['unit']='v6:'+new.id.astype(str);new['sources']='v6_cell_RNA_or_model_hint';cand=pd.concat([cand,new[cand.columns]],ignore_index=True);assert cand.unit.is_unique
 tree=cKDTree(cand[['x','y']]);xy=d[['x','y']].to_numpy();dst,near=tree.query(xy,workers=3);d['in_D']=dst<=200
 origin=np.floor((cand[['x','y']].min().to_numpy()-216)/8)*8;nx,ny=np.ceil((cand[['x','y']].max().to_numpy()+216-origin)/8).astype(int)+1
 gx,gy=np.meshgrid(np.arange(nx),np.arange(ny));mask=tree.query(np.c_[gx.ravel(),gy.ravel()]*8+origin,workers=3)[0].reshape(ny,nx)<=200
 lab,n=ndimage.label(mask,structure=np.ones((3,3)));np.savez_compressed(o/'D.npz',roi=lab.astype('uint16'),origin=origin,step_um=8)
 idx=np.rint((cand[['x','y']].to_numpy()-origin)/8).astype(int);cand['roi']=lab[idx[:,1],idx[:,0]];d['roi']=np.where(d.in_D,cand.roi.to_numpy()[near],0)
 nat=pd.read_parquet(V/s/'bin_evidence.parquet',columns=['array_row','array_col','component','bin_id']);au=js(N/s/'metadata/input_audit.json');aff=np.asarray(au['array_col_row_affine_fullres']);rc=((xy/au['native_mpp']-aff[2])@np.linalg.inv(aff[:2]));cc,rr=np.floor((rc+.5)/4).astype(int).T
 nat['key']=nat.array_row.astype('int64')*100000+nat.array_col;d['component']=pd.Series(rr*100000+cc).map(nat.set_index('key').component).fillna(-1).astype(int)
 d['geometry_assignable']=~multi & ~disputed;d['UMI_ok']=d.nCount>=25;d['direct_ok']=d.direct_width>=1;d['major_T_ok']=d.p_T>=.5;d['sub_Treg_ok']=d.p_Treg_given_T>=.5
 shared=d.in_D&d.geometry_assignable&d.UMI_ok&d.direct_ok
 d['A']=shared&d.major_T_ok&d.FOXP3.gt(0);d['B']=shared&d.major_T_ok&d.FOXP3.eq(0)&d.other_core_width.gt(0)&d.sub_Treg_ok
 d['main']=d.A|d.B;d['strict']=d.A&d.other_core_width.gt(0)&d.adjacent_ambiguous_key_bins.eq(0)
 d['extended']=shared&d.p_T.ge(.3)&(d.FOXP3.gt(0)|(d.other_core_width.gt(0)&d.p_Treg_given_T.ge(.3)));assert (d.extended|~d.main).all()
 d['control']=shared&d.major_T_ok&d.FOXP3.eq(0)&d.nonTreg_support.ge(.5)&d.p_Treg_given_T.lt(.2)
 d['epi']=d.in_D&~multi&~disputed&d.major.eq('Epithelial')&d.p_epi.ge(.5)&d.Epi_width.ge(2)
 d['outcome']=np.select([d.A,d.B,d.control,d.extended,multi,disputed,d.nCount.lt(25),d.direct_width.eq(0),~d.major_T_ok],['A','B','other_T','extended_only','unresolved_doublet','assignment_uncertain','low_UMI','no_direct_T','weak_major_T'],default='insufficient_subtype')
 d.to_parquet(o/'identity_frozen.parquet',index=False)
 # Many-to-many source mapping by measured member or overlapping nucleus.
 links=[];nuc=ov[ov.substantial].copy();nuc['unit']='n:'+nuc.old_id.astype(str);links.append(nuc.rename(columns={'new_id':'id'})[['unit','id']]);links.append(new[['unit','id']])
 rr=raw[raw.id.gt(0)].copy();rr['unit']='2:'+rr.source_bin_index.astype(str);links.append(rr[['unit','id']])
 ap['key']=(ap.array_row.astype('int64')//4)*100000+ap.array_col//4;par=ap.loc[ap.id.gt(0),['key','id']].drop_duplicates();grid=nat.merge(par,on='key');grid['unit']='8:'+grid.bin_id;links.append(grid[['unit','id']])
 for k in ['T16','T32']:
  mm=pd.read_parquet(V/s/(k+'_members.parquet'));mm['unit']='r:'+mm.region_id;links.append(mm[mm.unit.isin(cand.unit)].merge(grid[['bin_id','id']],on='bin_id')[['unit','id']])
 for k,pre in [('T8_from2','8grid:'),('T16_from2','16grid:')]:
  m=pd.read_parquet(U/s/(k+'_members.parquet'));m['unit']=pre+m.id;m['id']=ap.id.to_numpy()[m.source_bin_index];links.append(m[m.id.gt(0)][['unit','id']])
 ll=pd.concat(links,ignore_index=True);ll=ll[ll.unit.isin(cand.unit)].drop_duplicates();ll.to_parquet(o/'source_links.parquet',index=False)
 cand['n_final_cells']=cand.unit.map(ll.unit.value_counts()).fillna(0).astype(int);cand.to_parquet(o/'sources.parquet',index=False)
 q=d[d.id.isin(ll.id)].copy();q['review_id']='cell:'+q.id.astype(str);u=cand[cand.n_final_cells.eq(0)].copy();u['review_id']=u.unit
 queue=pd.concat([q[['review_id','id','x','y','roi','outcome']],u[['review_id','x','y','roi']].assign(outcome='unlocalized')],ignore_index=True);assert queue.review_id.is_unique
 queue['review_state']='pending';queue.to_parquet(o/'review_queue.parquet',index=False)
 put(cand.groupby('roi').agg(sources=('unit','size'),unlocalized=('n_final_cells',lambda x:(x==0).sum()),min_x=('x','min'),max_x=('x','max'),min_y=('y','min'),max_y=('y','max')).reset_index(),o/'roi_manifest.tsv')
 pr=pd.read_parquet(o/'priority_old.parquet').rename(columns={'id':'old_id'}).merge(ov[ov.substantial],on='old_id',how='left').merge(d[['id','outcome','mixed_v6','p_T','p_Treg_given_T','FOXP3','direct_width','nCount']],left_on='new_id',right_on='id',how='left',suffixes=('_old','_new'))
 put(pr,o/'priority_trace.tsv');print(s,len(cand),len(queue),d[['in_D','A','B','extended','strict','control','epi']].sum().to_dict(),flush=True)
if __name__=='__main__':
 for s in sys.argv[1:] or sections:run(s)
