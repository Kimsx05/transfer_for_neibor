from common import *

def summarize():
 if done('M02_BACKGROUND_READY'):log('M02 checkpoint reused');return
 m,q,center=load();sets=masks(m,q);cov=pd.read_csv(V2/'00_preflight/coverage63.tsv',sep='\t');p=pairs();rows=[];background=[];closure=[];ident=[]
 for aset,mask in sets.items():
  arr=q if aset!='TME_CONDITIONAL' else q/np.sum(q[:,:,NON],axis=2,keepdims=True)
  for cr in cov.itertuples():
   idx=m.index[mask&m.section_id.eq(cr.section_id)].to_numpy();g=m.loc[idx];bg=arr[idx].mean(0) if len(idx) else np.full((3,30),np.nan)
   for r in range(3):
    for j in ([v for v in range(30) if v!=EPI] if aset=='TME_CONDITIONAL' else range(30)):background.append(dict(analysis_set=aset,dataset_id=cr.dataset_id,section_id=cr.section_id,ring=RINGS[r],celltype=TYPES[j],n_background=len(idx),background=bg[r,j]))
   for k in range(9):
    ai=idx[g.accepted_label.to_numpy()==k];bi=idx[g.accepted_label.to_numpy()!=k];valid=len(ai)>=20 and len(bi)>=20;tm=arr[ai].mean(0) if len(ai) else np.full((3,30),np.nan);rest=arr[bi].mean(0) if len(bi) else np.full((3,30),np.nan);lam=len(ai)/len(idx) if len(idx) else np.nan;e=tm-bg;dr=tm-rest
    if len(ai) and len(bi):err=np.max(abs(e-(1-lam)*dr));assert err<1e-12;ident.append(dict(analysis_set=aset,section_id=cr.section_id,niche=f'niche{k}',max_error=err))
    if aset=='MAIN':
     for r in range(3):closure.append(dict(dataset_id=cr.dataset_id,section_id=cr.section_id,analysis_unit_uid=g.analysis_unit_uid.iloc[0] if len(g) else '',known_patient=bool(g.source_patient_identity_known.iloc[0]) if len(g) else False,niche=f'niche{k}',ring=RINGS[r],n_target=len(ai),Epi_ring=tm[r,EPI],TME_ring=np.sum(tm[r,NON]),sum_Epi_TME=tm[r,EPI]+np.sum(tm[r,NON]),status='AVAILABLE' if len(ai) else 'NO_TARGET'))
    targets=p[p.accepted_label.eq(k)]
    for pr in targets.itertuples():
     j=pr.cell_index
     for r in range(3):rows.append(dict(analysis_set=aset,pair_id=pr.pair_id,review_tier=pr.review_tier,niche=pr.niche,celltype=pr.celltype,ring=RINGS[r],dataset_id=cr.dataset_id,section_id=cr.section_id,analysis_unit_uid=g.analysis_unit_uid.iloc[0] if len(g) else '',analysis_unit_type=g.analysis_unit_type.iloc[0] if len(g) else '',known_patient=bool(g.source_patient_identity_known.iloc[0]) if len(g) else False,n_target=len(ai),n_rest=len(bi),n_background=len(idx),target_weight=lam,target=tm[r,j],background=bg[r,j],excess=e[r,j],rest=rest[r,j],target_rest=dr[r,j],eligible=valid,status='CALCULABLE' if valid else 'NO_TARGET' if not len(ai) else 'INSUFFICIENT_TARGET_OR_REST'))
  log('M02 section summaries '+aset)
 s=pd.DataFrame(rows);tab('02_background_composition/section_ring123_composition.tsv.gz',s);tab('02_background_composition/background_by_section_ring.tsv.gz',background);tab('02_background_composition/background_identity_audit.tsv.gz',ident);tab('05_sensitivity/composition_closure_diagnostics.tsv',closure)
 u=[]
 for scope in SCOPES:
  z=s if scope=='ALL_DATASETS' else s[s.dataset_id.eq(scope)]
  for summary,zz in [('INFERENTIAL',z[z.eligible]),('DESCRIPTIVE',z[z.n_target.gt(0)])]:
   key=['analysis_set','pair_id','review_tier','niche','celltype','ring','analysis_unit_uid','analysis_unit_type','known_patient'];x=zz.groupby(key,dropna=False).agg(target=('target','mean'),background=('background','mean'),excess=('excess','mean'),rest=('rest','mean'),target_rest=('target_rest','mean'),target_weight=('target_weight','mean'),n_sections=('section_id','nunique'),n_anchors=('n_target','sum'),n_datasets=('dataset_id','nunique')).reset_index();x['scope']=scope;x['summary_type']=summary;u.append(x)
 u=pd.concat(u,ignore_index=True);tab('02_background_composition/unit_ring123_composition.tsv.gz',u)
 cl=pd.DataFrame(closure);cunits=[]
 for scope in SCOPES:
  z=cl[cl.n_target.gt(0)&(True if scope=='ALL_DATASETS' else cl.dataset_id.eq(scope))];z=z.groupby(['analysis_unit_uid','known_patient','niche','ring']).agg(Epi_ring=('Epi_ring','mean'),TME_ring=('TME_ring','mean'),n_sections=('section_id','nunique')).reset_index();z['scope']=scope;cunits.append(z)
 tab('05_sensitivity/closure_unit_summary.tsv',pd.concat(cunits,ignore_index=True))
 # Direction/coverage summaries for all finite sensitivity sets, without reusing primary q.
 su=u[u.summary_type.eq('INFERENTIAL')].groupby(['analysis_set','scope','pair_id','niche','celltype','ring']).agg(n_units=('excess','count'),median_excess=('excess','median'),mean_excess=('excess','mean'),median_target=('target','median'),median_rest_difference=('target_rest','median'),n_positive=('excess',lambda x:int((x>0).sum())),n_negative=('excess',lambda x:int((x<0).sum()))).reset_index()
 # Fill planned cells including completely empty scope/combination.
 planned=pd.DataFrame([dict(analysis_set=aset,scope=scope,pair_id=pr.pair_id,niche=pr.niche,celltype=pr.celltype,ring=r) for aset in SETS for scope in SCOPES for pr in p.itertuples() for r in RINGS]);su=planned.merge(su,how='left',validate='one_to_one');su['n_units']=su.n_units.fillna(0).astype(int);su['status']=np.where(su.n_units.ge(5),'DESCRIPTIVE_SENSITIVITY','INSUFFICIENT_UNITS')
 tab('05_sensitivity/TME_conditional_sensitivity.tsv.gz',su[su.analysis_set.isin(['MAIN','TME_MATCHED_RAW','TME_CONDITIONAL'])]);tab('05_sensitivity/complete_rings_composition.tsv.gz',su[su.analysis_set.isin(['MAIN','COMPLETE_RINGS'])]);checkpoint('M02_BACKGROUND_READY')

def statistics():
 if done('M03_COMPOSITION_STATISTICS_READY'):return
 u=pd.read_csv(OUT/'02_background_composition/unit_ring123_composition.tsv.gz',sep='\t');s=pd.read_csv(OUT/'02_background_composition/section_ring123_composition.tsv.gz',sep='\t');p=pairs();zero=[];pair=[]
 u=u[u.analysis_set.eq('MAIN')&u.summary_type.eq('INFERENTIAL')];groups={(sc,pid):g for (sc,pid),g in u.groupby(['scope','pair_id'])}
 for subset in ['ALL_POLICY_UNITS','VERIFIED_PATIENTS_ONLY']:
  for scope in SCOPES:
   for pr in p.itertuples():
    g=groups.get((scope,pr.pair_id),u.iloc[:0]);g=g[g.known_patient] if subset=='VERIFIED_PATIENTS_ONLY' else g
    ids=['analysis_unit_uid','analysis_unit_type','known_patient'];ex=g.pivot(index=ids,columns='ring',values='excess').reindex(columns=RINGS);raw=g.pivot(index=ids,columns='ring',values='target').reindex(columns=RINGS);assert ex.notna().all().all();n=len(ex);known=int(np.sum(ex.index.get_level_values('known_patient'))) if n else 0
    base=dict(unit_subset=subset,scope=scope,pair_id=pr.pair_id,niche=pr.niche,celltype=pr.celltype,review_tier=pr.review_tier,n_verified_units=known,n_other_units=n-known)
    for r in RINGS:zero.append(dict(**base,ring=r,**signed(ex[r]),raw_target_median=raw[r].median()))
    for a,b in [('R1','R2'),('R1','R3'),('R2','R3')]:pair.append(dict(**base,ring_pair=a+'-'+b,**signed(ex[a]-ex[b]),raw_pair_mean=(raw[a]-raw[b]).mean(),raw_pair_median=(raw[a]-raw[b]).median()))
 z=add_bh(pd.DataFrame(zero),['unit_subset'],417,'background_zero_q');r=add_bh(pd.DataFrame(pair),['unit_subset'],417,'ring_pairwise_q');tab('03_ring_statistics/ring123_background_tests.tsv',z);tab('03_ring_statistics/ring123_pairwise_tests.tsv',r)
 checkpoint('M03_COMPOSITION_STATISTICS_READY');log('M03 417 planned background and 417 ring-pair tests per scope; separate verified sensitivity')
if __name__=='__main__':summarize();statistics()
