from common import *
rows=[]
for sn,s in enumerate(SECTIONS):
 n=pd.read_parquet(R/s/'neighborhoods.parquet');e0=pd.read_parquet(R/s/'Epi_regions.parquet');pp=pd.read_parquet(V/s/'region_pairs.parquet');pp=pp[pp.T_grid_um.eq(32)&pp.radius_um.eq(200)&~pp.shared_original_bins]
 for radius in RADII:
  t=n[n.radius_um.eq(radius)&n.coverage.ge(.5)].reset_index(drop=True);e=e0[e0.component.isin(t.component)].reset_index(drop=True);ti=pd.Series(np.arange(len(t)),index=t.region_id);ei=pd.Series(np.arange(len(e)),index=e.region_id);q=pp[pp.centroid_distance_um.le(radius+1e-8)&pp.T_region_id.isin(ti.index)&pp.Epi_region_id.isin(ei.index)];it=ti.loc[q.T_region_id].to_numpy();ie=ei.loc[q.Epi_region_id].to_numpy();tt=pairlist(t,radius);ep=pairlist(e,radius);S=sparse.csr_matrix((np.ones(len(tt)*2,np.float32),(np.r_[tt[:,0],tt[:,1]],np.r_[tt[:,1],tt[:,0]])),shape=(len(t),len(t)));L=t[LABELS].to_numpy(np.float32);D=np.column_stack([np.bincount(it,weights=e[h].to_numpy()[ie],minlength=len(t))for h in HIGHS]).astype(np.float32);group=strata(t,t.Epi_n);nullAB=np.zeros((999,2,2));nullAA=np.zeros((999,2))
  for start,lp in perm_batches(L,group,999,123+sn*10+radius):
   for j in range(2):
    ll=lp[:,:,j];nullAB[start:start+len(lp),j,:]=ll@D;nullAA[start:start+len(lp),j]=np.sum(ll*(S@ll.T).T,axis=1)/2
  for j,label in enumerate(LABELS):
   npT=int(L[:,j].sum());taa=int((L[tt[:,0],j]*L[tt[:,1],j]).sum());O=L[:,j]@D
   for k,h in enumerate(HIGHS):
    nh=int(e[h].sum());ee=e[h].to_numpy();eaa=int((ee[ep[:,0]]&ee[ep[:,1]]).sum())
    for direction,nA,nB,aa,naa in [('label_to_high_Epi',npT,nh,taa,nullAA[:,j]),('high_Epi_to_label',nh,npT,eaa,np.full(999,eaa))]:
     ok=nA>=2 and nB>0 and aa>0;ms=float(O[k]/aa)if ok else np.nan;nms=float((O[k]/(nA*nB))/(aa/(nA*(nA-1)/2)))if ok else np.nan
     good=(naa>0)&(nA>=2)&(nB>0);nm=np.full(999,np.nan);nn=np.full(999,np.nan);nm[good]=nullAB[good,j,k]/naa[good]
     if nA>=2 and nB>0:nn[good]=(nullAB[good,j,k]/(nA*nB))/(naa[good]/(nA*(nA-1)/2))
     rec=dict(section=s,radius_um=radius,high_definition=h,label_definition=label,direction=direction,n_A=nA,n_B=nB,I_AB=int(O[k]),I_AA=aa,MS=ms,NMS=nms,MS_status=status(ms),NMS_status=status(nms),N_T_pool=len(t),N_Epi_pool=len(e),isolated_T_included=int(t.Epi_n.eq(0).sum()),reason='estimable_descriptive'if ok else'absent_type_or_nA_lt2_or_no_AA_pairs',pair_convention='AB once; AA unordered; no diagonal; whole shared-source AB excluded',null_B=999,null_valid_fraction=float(good.mean()),P=np.nan)
     for name,val in [('MS',nm),('NMS',nn)]:
      usable=good.mean()>=.95;rec[name+'_null_mean']=float(np.nanmean(val))if usable else np.nan;lo,hi=np.nanquantile(val,[.025,.975])if usable else[np.nan,np.nan];rec[name+'_null_low']=lo;rec[name+'_null_high']=hi
     rows.append(rec)
  print('MIXING',s,radius,flush=True)
put(pd.DataFrame(rows),'mixing_scores.tsv')
