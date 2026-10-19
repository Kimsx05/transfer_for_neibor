from core import *
from scipy.stats import rankdata,wilcoxon
from scipy.linalg import qr
extract(V7/'code/common.py',['signed'],globals())
extract(V7/'code/association.py',['residualize','correlate'],globals())

def section_rings(m,q,center,K,meta,aset='MAIN'):
 rows=[];ar=[];base=(m[['R1_valid_n','R2_valid_n','R3_valid_n']].to_numpy()>=[4,8,12]).all(1)&m.accepted_label.ge(0).to_numpy();sc=base&m.epi_score_valid.to_numpy()&m.n_signature_available_model.eq(20).to_numpy()&m.n_signature_available_raw.eq(20).to_numpy()&np.isfinite(m[LAYERS]).all(1).to_numpy()
 if aset=='COMPLETE_RINGS':base&=m.complete123.to_numpy();sc&=base
 if aset in ['TME_MATCHED_RAW','TME_CONDITIONAL']:
  base&=(q[:,:,NON].sum(2)>=.05).all(1);sc&=base
  if aset=='TME_CONDITIONAL':q=q/q[:,:,NON].sum(2,keepdims=True)
 if aset=='SPARSE_BLOCK':
  ii=np.flatnonzero(sc);d=m.iloc[ii].copy();d['draw']=[hashlib.sha256(('123|'+uid).encode()).hexdigest() for uid in d.spot_uid];d['br']=d.array_row//8;d['bc']=d.array_col//8;ix=d.sort_values('draw').drop_duplicates(['br','bc']).index.to_numpy();sc[:]=False;sc[ix]=True
 bg=q[base].mean(0) if base.any() else np.full((3,30),np.nan)
 for k in range(K):
  use=base&m.accepted_label.eq(k).to_numpy();rest=base&~use;nt=int(use.sum());nr=int(rest.sum());valid=nt>=20 and nr>=20;target=q[use].mean(0) if nt else np.full((3,30),np.nan);ex=target-bg;other=q[rest].mean(0) if nr else np.full((3,30),np.nan)
  for j in NON:
   for r in range(3):rows.append(dict(**meta,analysis_set=aset,niche=k,celltype=TYPES[j],endpoint=f'R{r+1}',kind='BACKGROUND',raw=target[r,j],background=bg[r,j],effect=ex[r,j],rest=other[r,j],n_target=nt,n_rest=nr,eligible=valid,status='CALCULABLE' if valid else 'INSUFFICIENT_TARGET_OR_REST'))
   for a,b in [(0,1),(0,2),(1,2)]:rows.append(dict(**meta,analysis_set=aset,niche=k,celltype=TYPES[j],endpoint=f'R{a+1}-R{b+1}',kind='RINGPAIR',raw=target[a,j]-target[b,j],background=bg[a,j]-bg[b,j],effect=ex[a,j]-ex[b,j],rest=other[a,j]-other[b,j],n_target=nt,n_rest=nr,eligible=valid,status='CALCULABLE' if valid else 'INSUFFICIENT_TARGET_OR_REST'))
 for k in [-99]+list(range(K)):
  ii=np.flatnonzero(sc&(True if k==-99 else m.accepted_label.eq(k).to_numpy()));n=len(ii);X=q[ii][:,:,NON].reshape(n,87);Y=m.iloc[ii][LAYERS].to_numpy(float);C=np.column_stack([center[ii,EPI],np.log1p(m.iloc[ii].center_nCount.to_numpy(float))]);complete=np.isfinite(X).all(1)&np.isfinite(Y).all(1)&np.isfinite(C).all(1);X=X[complete];Y=Y[complete];C=C[complete];ii=ii[complete];n=len(ii);mats={key:np.full((87,2),np.nan) for key in ['RAW','ADJUSTED_MAIN','ADJUSTED_CENTER_TYPE']};status=np.full(87,'INSUFFICIENT_CENTERS',object);extra_status=status.copy();notes='';rank=0
  if n>=30:
   xvar=np.array([len(np.unique(x))>=3 for x in X.T]);yvar=all(len(np.unique(y))>=3 for y in Y.T);XR=rankdata(X,axis=0,method='average');YR=rankdata(Y,axis=0,method='average');CR=rankdata(C,axis=0,method='average');mats['RAW']=correlate(XR,YR);res,rank,notes=residualize(np.column_stack([XR,YR]),CR);status[:]='CALCULABLE'
   if not yvar:status[:]='SCORE_TOO_FEW_DISTINCT'
   status[~xvar]='COMPOSITION_TOO_FEW_DISTINCT';bad=status!='CALCULABLE';mats['RAW'][bad]=np.nan
   if res is None:status[status=='CALCULABLE']='COLLINEAR_DESIGN'
   elif n-rank<3:status[status=='CALCULABLE']='INSUFFICIENT_RESIDUAL_DF'
   else:mats['ADJUSTED_MAIN']=correlate(res[:,:87],res[:,87:]);mats['ADJUSTED_MAIN'][bad]=np.nan;status[(~np.isfinite(mats['ADJUSTED_MAIN']).all(1))&(status=='CALCULABLE')]='ZERO_RESIDUAL_VARIANCE'
   extra_status=status.copy()
   for ci,j in enumerate(NON):
    jj=np.array([ci,29+ci,58+ci]);cx=center[ii,j]
    if not np.isfinite(cx).all():extra_status[jj]='MISSING_CENTER_TYPE_COVARIATE';continue
    ER,rnk,note=residualize(np.column_stack([XR[:,jj],YR]),np.column_stack([CR,rankdata(cx,method='average')]))
    if ER is None:extra_status[jj]='COLLINEAR_DESIGN';continue
    if n-rnk<3:extra_status[jj]='INSUFFICIENT_RESIDUAL_DF';continue
    mats['ADJUSTED_CENTER_TYPE'][jj]=correlate(ER[:,:3],ER[:,3:]);mats['ADJUSTED_CENTER_TYPE'][jj[bad[jj]]]=np.nan
   for mat in mats.values():mat[~np.isfinite(mat).all(1)]=np.nan
  for method,mat in mats.items():
   for rn in range(3):
    for ci,j in enumerate(NON):
     for li,layer in enumerate(LAYERS):
      rho=mat[rn*29+ci,li];why=(extra_status if method=='ADJUSTED_CENTER_TYPE' else status)[rn*29+ci];ar.append(dict(**meta,analysis_set=aset,niche=k,celltype=TYPES[j],endpoint=f'R{rn+1}',layer=layer,method=method,n_centers=n,rho=rho,status='CALCULABLE' if np.isfinite(rho) else why if why!='CALCULABLE' else 'ZERO_VARIANCE_OR_UNAVAILABLE',covariate_audit=notes))
 return (pd.DataFrame(rows).iloc[:0] if aset=='SPARSE_BLOCK' else pd.DataFrame(rows)),pd.DataFrame(ar)

def aggregate(s,a,K):
 units=[];tests=[]
 for scope in ['ALL_DATASETS']+sorted(s.dataset_id.unique()):
  ss=s if scope=='ALL_DATASETS' else s[s.dataset_id.eq(scope)];aa=a if scope=='ALL_DATASETS' else a[a.dataset_id.eq(scope)]
  for kind,d,keys in [('COMPOSITION',ss,['analysis_set','niche','celltype','kind','endpoint']),('ASSOCIATION',aa,['analysis_set','niche','celltype','method','endpoint','layer'])]:
   for key,g in d.groupby(keys,sort=False):
    base=dict(zip(keys,key));gg=g[g.eligible&g.effect.notna()] if kind=='COMPOSITION' else g[g.rho.notna()];u=[]
    for uid,v in gg.groupby('analysis_unit_uid',sort=False):
     effect=v.effect.mean() if kind=='COMPOSITION' else np.tanh(np.arctanh(np.clip(v.rho.to_numpy(float),-1+1e-7,1-1e-7)).mean());r=dict(scope=scope,category=kind,**base,analysis_unit_uid=uid,effect=effect,known_patient=bool(v.known_patient.all()),n_sections=len(v),datasets='|'.join(sorted(v.dataset_id.unique())),raw=v.raw.mean() if kind=='COMPOSITION' else np.nan);u.append(r);units.append(r)
    for subset in ['ALL_POLICY_UNITS','VERIFIED_PATIENTS_ONLY']:
     uu=[v for v in u if subset=='ALL_POLICY_UNITS' or v['known_patient']];family='RING|'+(base['kind'] if kind=='COMPOSITION' else 'ASSOCIATION|'+base['method'])+'|'+subset;stats=signed([v['effect'] for v in uu]);
     if base['analysis_set']!='MAIN':stats.update(p=np.nan,W=np.nan,actual_method='DESCRIPTIVE_SENSITIVITY_NO_TEST')
     tests.append(dict(scope=scope,category=kind,subset=subset,**base,family_id=family,planned_n=3*K*29 if kind=='COMPOSITION' else 6*(K+1)*29,**stats))
 r=pd.DataFrame(tests);r['q_within_BH']=np.nan
 for key,idx in r.groupby(['scope','family_id','analysis_set']).groups.items():
  assert len(idx)==r.loc[idx,'planned_n'].iloc[0],(key,len(idx));r.loc[idx,'q_within_BH']=correction(r.loc[idx,'p'])
 return pd.DataFrame(units),r

def run(cid,sets=['MAIN']):
 folder=OUT/'combinations'/cid
 if sets==['MAIN'] and valid(folder,'RING_DONE'):return
 assert valid(folder,'NICHE_FROZEN');guard(cid+':rings',20*1024**2);labels=np.load(folder/'labels.npz');idx=labels['anchor_index'];K=json.loads((folder/'selection.json').read_text())['K'];allm=pd.read_parquet(V7/'01_ring123/anchor_metadata.parquet');m=allm.iloc[idx].reset_index(drop=True);m['accepted_label']=labels['accepted'];arr=np.load(V7/'01_ring123/ring123_arrays.npz');assert np.array_equal(arr['spot_uid'],allm.spot_uid);q=arr['q'][idx];center=arr['center_mean'][idx];cov=pd.read_csv(V2/'00_preflight/coverage63.tsv',sep='\t');members=pd.read_csv(OUT/'combinations.tsv',sep='\t',dtype={'combo_id':str}).set_index('combo_id').loc[cid,'datasets'].split('|');cs=[];ass=[]
 for sec in cov[cov.dataset_id.isin(members)].section_id:
  take=np.flatnonzero(m.section_id.eq(sec));p=V8/'01_geometry'/sec/'metadata.json';meta=json.loads(p.read_text());meta={k:meta[k] for k in ['dataset_id','section_id','analysis_unit_uid','known_patient']}
  for aset in sets:
   dest=folder/'ring_sections'/f'{sec}_{aset}_composition.parquet';adest=folder/'ring_sections'/f'{sec}_{aset}_association.parquet';state=dest.with_suffix('.json')
   if state.exists():
    d=json.loads(state.read_text());assert d['label_hash']==sha(folder/'labels.npz') and d['composition']==sha(dest) and d['association']==sha(adest);s=pd.read_parquet(dest);a=pd.read_parquet(adest)
   else:
    s,a=section_rings(m.iloc[take].reset_index(drop=True),q[take],center[take],K,meta,aset);table(dest,s);table(adest,a);js(state,dict(label_hash=sha(folder/'labels.npz'),composition=sha(dest),association=sha(adest)))
   cs.append(s);ass.append(a)
 s=pd.concat(cs,ignore_index=True);a=pd.concat(ass,ignore_index=True);u,r=aggregate(s,a,K);suffix='' if sets==['MAIN'] else '_sensitivity';table(folder/f'ring_units{suffix}.parquet',u);table(folder/f'ring_tests{suffix}.parquet',r)
 if sets==['MAIN']:mark(folder,'RING_DONE',[folder/'ring_units.parquet',folder/'ring_tests.parquet']+list((folder/'ring_sections').glob('*')))
 log(cid+' ring completed '+','.join(sets))
if __name__=='__main__':run(sys.argv[1])
