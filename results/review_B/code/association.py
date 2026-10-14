from common import *
from scipy.linalg import qr
METHODS=['RAW','ADJUSTED_MAIN','ADJUSTED_CENTER_TYPE']
def residualize(X,C):
 # Center/normalize rank columns for stable QR; remove constants explicitly.
 cv=np.asarray(C,float);sd=cv.std(0);keep=sd>1e-12;C=np.column_stack([np.ones(len(cv)),(cv[:,keep]-cv[:,keep].mean(0))/sd[keep]])
 Q,R,piv=qr(C,mode='economic',pivoting=True);diag=abs(np.diag(R));rank=int((diag>diag.max()*1e-10).sum());note='removed_constant='+','.join(map(str,np.flatnonzero(~keep)))+';QR_rank='+str(rank)+'/'+str(C.shape[1])
 if rank<C.shape[1]:return None,rank,note+';COLLINEAR_DESIGN'
 return X-Q[:,:rank]@(Q[:,:rank].T@X),rank,note

def correlate(X,Y):
 X=X-X.mean(0);Y=Y-Y.mean(0);den=np.sqrt(np.sum(X*X,0)[:,None]*np.sum(Y*Y,0)[None,:]);return np.divide(X.T@Y,den,out=np.full((X.shape[1],Y.shape[1]),np.nan),where=den>1e-12)

def section_associations():
 if done('M04_SECTION_ASSOCIATIONS_READY'):return
 m,q,center=load();sets=masks(m,q);score=m.score_eligible.to_numpy();cov=pd.read_csv(V2/'00_preflight/coverage63.tsv',sep='\t');qq=queries();groups={k:[v for v in qq if v['accepted_label']==k] for k in [-99]+list(range(1,9))}
 # Shared draw for all queried types/niches/rings/layers within the main scored pool.
 selected=np.zeros(len(m),bool);blockrows=[]
 for sec,g in m[m.score_eligible].groupby('section_id'):
  z=g.copy();z['block_row']=z.array_row//8;z['block_col']=z.array_col//8;z['stable_draw']=[hashlib.sha256(('123|'+uid).encode()).hexdigest() for uid in z.spot_uid];keep=z.sort_values('stable_draw').drop_duplicates(['block_row','block_col']);selected[keep.index]=True
  blockrows.extend(keep[['spot_uid','section_id','block_row','block_col','stable_draw']].to_dict('records'))
 tab('05_sensitivity/spatial_block_selected_anchors.tsv',blockrows);sets['SPARSE_BLOCK']=selected
 for ci,cr in enumerate(cov.itertuples()):
  dest=OUT/'04_MP10_association/section_cache'/f'{cr.section_id}.parquet';dest.parent.mkdir(exist_ok=True)
  if dest.exists():continue
  rows=[];secix=m.index[m.section_id.eq(cr.section_id)].to_numpy();meta=m.loc[secix].iloc[0] if len(secix) else None
  for aset,mask in sets.items():
   arr=q if aset!='TME_CONDITIONAL' else q/q[:,:,NON].sum(2,keepdims=True)
   for k,qs in groups.items():
    ix=secix[mask[secix]&score[secix]&(True if k==-99 else m.loc[secix,'accepted_label'].to_numpy()==k)];n=len(ix);cells=[t['cell_index'] for t in qs];X=arr[ix][:,:,cells].reshape(n,3*len(cells));Y=m.loc[ix,LAYERS].to_numpy(float);C=np.column_stack([center[ix,EPI],np.log1p(m.loc[ix,'center_nCount'].to_numpy(float))]);complete=np.isfinite(X).all(1)&np.isfinite(Y).all(1)&np.isfinite(C).all(1);X=X[complete];Y=Y[complete];C=C[complete];ii=ix[complete];n=len(ii);raw=np.full((3*len(cells),2),np.nan);main=raw.copy();extra=raw.copy();notes='';status=np.full(3*len(cells),'INSUFFICIENT_CENTERS',dtype=object);extra_status=status.copy();cov_rank=0
    if n>=30:
     xvar=np.array([len(np.unique(x))>=3 for x in X.T]);yvar=all(len(np.unique(y))>=3 for y in Y.T)
     XR=rankdata(X,axis=0,method='average');YR=rankdata(Y,axis=0,method='average');CR=rankdata(C,axis=0,method='average');raw=correlate(XR,YR);both=np.column_stack([XR,YR]);res,cov_rank,notes=residualize(both,CR);status[:]='CALCULABLE'
     if not yvar:status[:]='SCORE_TOO_FEW_DISTINCT'
     status[~xvar]='COMPOSITION_TOO_FEW_DISTINCT';bad=(status!='CALCULABLE');raw[bad]=np.nan
     if res is None:status[status=='CALCULABLE']='COLLINEAR_DESIGN'
     elif n-cov_rank<3:status[status=='CALCULABLE']='INSUFFICIENT_RESIDUAL_DF'
     else:
      main=correlate(res[:,:X.shape[1]],res[:,X.shape[1]:]);main[bad]=np.nan;status[(~np.isfinite(main).all(1))&(status=='CALCULABLE')]='ZERO_RESIDUAL_VARIANCE'
     extra_status=status.copy()
     for j,c in enumerate(cells):
      cx=center[ii,c];jj=np.array([j,len(cells)+j,2*len(cells)+j]);ok=np.isfinite(cx)
      if not ok.all():extra_status[jj]='MISSING_CENTER_TYPE_COVARIATE';continue
      E=np.column_stack([CR,rankdata(cx,method='average')]);ER,rank,note=residualize(np.column_stack([XR[:,jj],YR]),E)
      if ER is None:extra_status[jj]='COLLINEAR_DESIGN';continue
      if n-rank<3:extra_status[jj]='INSUFFICIENT_RESIDUAL_DF';continue
      extra[jj]=correlate(ER[:,:3],ER[:,3:]);extra[jj[bad[jj]]]=np.nan;extra_status[jj[~np.isfinite(extra[jj]).all(1)]]='ZERO_RESIDUAL_VARIANCE'
     # Same usable section per method/ring/type across both score layers.
     for mat in [raw,main,extra]:mat[~np.isfinite(mat).all(1)]=np.nan
    for j,t in enumerate(qs):
     for rn,ring in enumerate(RINGS):
      col=rn*len(cells)+j
      for method,mat in [('RAW',raw),('ADJUSTED_MAIN',main),('ADJUSTED_CENTER_TYPE',extra)]:
       for li,layer in enumerate(LAYERS):
        value=mat[col,li];st='CALCULABLE' if np.isfinite(value) else (extra_status[col] if method=='ADJUSTED_CENTER_TYPE' else status[col]);st=st if st!='CALCULABLE' else 'ZERO_VARIANCE_OR_UNAVAILABLE'
        rows.append(dict(analysis_set=aset,query_id=t['query_id'],range=t['range'],accepted_label=k,celltype=t['celltype'],ring=ring,score_layer=layer,method=method,dataset_id=cr.dataset_id,section_id=cr.section_id,analysis_unit_uid=meta.analysis_unit_uid if meta is not None else '',analysis_unit_type=meta.analysis_unit_type if meta is not None else '',known_patient=bool(meta.source_patient_identity_known) if meta is not None else False,n_centers=n,rho=value,rho_for_fisher=np.clip(value,-1+1e-7,1-1e-7),fisher_clipped=bool(np.isfinite(value) and abs(value)>1-1e-7),status='CALCULABLE' if np.isfinite(value) else st,covariate_audit=notes,main_QR_rank=cov_rank))
  pd.DataFrame(rows).to_parquet(dest,index=False);log(f'M04 associations section {ci+1}/63 {cr.section_id}')
 checkpoint('M04_SECTION_ASSOCIATIONS_READY')

def aggregate():
 if done('M04_MP10_ASSOCIATIONS_DONE'):return
 paths=sorted((OUT/'04_MP10_association/section_cache').glob('*.parquet'));assert len(paths)==63;s=pd.concat([pd.read_parquet(p) for p in paths],ignore_index=True);tab('04_MP10_association/MP10_section_associations.tsv.gz',s);s['fisher_z']=np.arctanh(s.rho_for_fisher);units=[]
 for scope in SCOPES:
  z=s[s.rho.notna()&(True if scope=='ALL_DATASETS' else s.dataset_id.eq(scope))];key=['analysis_set','query_id','range','accepted_label','celltype','ring','score_layer','method','analysis_unit_uid','analysis_unit_type','known_patient'];u=z.groupby(key,dropna=False).agg(mean_fisher_z=('fisher_z','mean'),n_sections=('section_id','nunique'),n_centers=('n_centers','sum'),n_datasets=('dataset_id','nunique'),n_clipped_sections=('fisher_clipped','sum')).reset_index();u['rho']=np.tanh(u.mean_fisher_z);u['scope']=scope;units.append(u)
 u=pd.concat(units,ignore_index=True);tab('04_MP10_association/MP10_unit_associations.tsv.gz',u)
 # Actual unit coverage and effect for all planned sensitivity cells, not only successful groups.
 key=['analysis_set','scope','query_id','range','accepted_label','celltype','ring','score_layer','method'];desc=u.groupby(key).agg(n_units=('rho','count'),median_rho=('rho','median'),mean_rho=('rho','mean'),n_positive=('rho',lambda x:int((x>0).sum())),n_negative=('rho',lambda x:int((x<0).sum()))).reset_index()
 plan=pd.DataFrame([dict(analysis_set=aset,scope=sc,**q,ring=r,score_layer=layer,method=method) for aset in SETS+['SPARSE_BLOCK'] for sc in SCOPES for q in queries() for r in RINGS for layer in LAYERS for method in METHODS]).drop(columns='cell_index');desc=plan.merge(desc,on=key,how='left',validate='one_to_one');desc['n_units']=desc.n_units.fillna(0).astype(int);desc['status']=np.where(desc.n_units.ge(5),'DESCRIPTIVE_EFFECT','INSUFFICIENT_UNITS');tab('05_sensitivity/MP10_all_sensitivity_effects.tsv.gz',desc)
 main=u[u.analysis_set.eq('MAIN')];groups={k:g for k,g in main.groupby(['scope','query_id','ring','score_layer','method'])};tests=[]
 for subset in ['ALL_POLICY_UNITS','VERIFIED_PATIENTS_ONLY']:
  for sc in SCOPES:
   for query in queries():
    for r in RINGS:
     for layer in LAYERS:
      for method in METHODS:
       g=groups.get((sc,query['query_id'],r,layer,method),main.iloc[:0]);g=g[g.known_patient] if subset=='VERIFIED_PATIENTS_ONLY' else g;tests.append(dict(unit_subset=subset,scope=sc,**query,ring=r,score_layer=layer,method=method,**signed(g.rho),n_verified_units=int(g.known_patient.sum()),n_other_units=int((~g.known_patient).sum())))
 test=add_bh(pd.DataFrame(tests),['unit_subset','method'],1008,'association_q');tab('04_MP10_association/MP10_scope_associations.tsv',test);checkpoint('M04_MP10_ASSOCIATIONS_DONE');log('M04 all ranges/layers/methods complete; primary family 1008 per scope')
if __name__=='__main__':section_associations();aggregate()
