from common import *
TARGETS=[('NICHE',f'niche{i}',i) for i in range(9)]+[('CELL',TYPES[j],j) for j in NON]

def stratify(t,ids):
 if len(ids)<40:return None,'STRATIFIED_BACKGROUND_UNAVAILABLE'
 cov=np.column_stack([t.epithelial_fraction.to_numpy(float)[ids],np.log1p(t.center_nCount.to_numpy(float)[ids])])
 if not np.isfinite(cov).all():return None,'STRATIFIED_COVARIATE_MISSING'
 for bins in [3,2]:
  cols=[np.zeros(len(ids)) if len(np.unique(cov[:,i]))==1 else np.asarray(pd.qcut(cov[:,i],bins,labels=False,duplicates='drop'),float) for i in range(2)]
  if any(not np.isfinite(c).all() for c in cols):continue
  codes=cols[0].astype(int)*bins+cols[1].astype(int);labs,ct=np.unique(codes,return_counts=True)
  if len(labs)>1 and ct.min()>=20:return codes,f'{bins}x{bins}_ties_merged_actual_{len(np.unique(cols[0]))}x{len(np.unique(cols[1]))}_nonempty{len(labs)}'
 return None,'STRATIFIED_BACKGROUND_UNAVAILABLE'

def metric_record(metric,ids,Hmask,values,targetmask,t,component,ntarget,base_status):
 hm=Hmask[ids];nh=int(hm.sum());ne=len(ids);non=~hm;v=np.asarray(values,float);meanH=v[hm].mean() if nh else np.nan;meanE=v.mean() if ne else np.nan;meanN=v[non].mean() if non.any() else np.nan;effect=meanH-meanE;deltaNon=meanH-meanN;identity=(1-nh/ne)*deltaNon if ne else np.nan
 if np.isfinite(effect) and np.isfinite(identity):assert abs(effect-identity)<1e-9*max(1,abs(effect)),('SHRINKAGE_FAIL',metric,effect,identity)
 member=int(targetmask[ids].sum());eligible=ne>=30 and nh>=20 and non.sum()>=20 and (member>=20 and ne-member>=20 if metric=='O' else ntarget>=5)
 status=base_status
 if status=='READY':status='EMPTY_METRIC_DOMAIN' if not ne else 'NO_HOT_IN_METRIC_DOMAIN' if not nh else 'CALCULABLE' if eligible else 'DESCRIPTIVE_ONLY_INSUFFICIENT_COVERAGE'
 if base_status!='READY':eligible=False;effect=np.nan
 codes,strat_status=stratify(t,ids);strat_effect=np.nan;strat_bg=np.nan
 if codes is not None and nh and base_status=='READY':
  strat_bg=sum((hm&(codes==b)).sum()/nh*v[codes==b].mean() for b in np.unique(codes));strat_effect=meanH-strat_bg
 comps=np.unique(component[ids[hm]]);cm=[v[hm&(component[ids]==c)].mean() for c in comps];component_effect=np.mean(cm)-meanE if len(cm) and base_status=='READY' else np.nan
 def cv(col,mask):
  z=t[col].to_numpy(float)[ids][mask];return float(np.nanmedian(z)) if len(z) else np.nan
 failures=[]
 for label,value,minimum in [('E',ne,30),('H',nh,20),('NONHOT',int(non.sum()),20)]+([('TARGET_IN_E',member,20),('REST_IN_E',ne-member,20)] if metric=='O' else [('TARGET_FULL',ntarget,5)]):
  if value<minimum:failures.append(label+'_LT_'+str(minimum))
 return dict(coverage_fail_reasons=';'.join(failures),metric=metric,effect=effect,effect_pp=100*effect if metric=='O' else np.nan,mean_hot=meanH,mean_background=meanE,median_hot=np.median(v[hm]) if nh else np.nan,median_background=np.median(v) if ne else np.nan,mean_nonhot=meanN,effect_vs_nonhot=deltaNon,shrinkage_error=effect-identity,n_E=ne,n_H=nh,n_nonhot=int(non.sum()),n_target_in_E=member,n_rest_in_E=ne-member,n_target_full=ntarget,n_hot_target=int(targetmask[ids[hm]].sum()),hot_fraction=nh/ne if ne else np.nan,eligible=bool(eligible and base_status=='READY'),status=status,stratified_effect=strat_effect,stratified_background=strat_bg,stratification_status=strat_status,component_equal_effect=component_effect,n_hot_components=len(comps),Epi_median_H=cv('epithelial_fraction',hm),Epi_median_E=cv('epithelial_fraction',np.ones(ne,bool)),nCount_median_H=cv('center_nCount',hm),nCount_median_E=cv('center_nCount',np.ones(ne,bool)))

def section_effects(sec):
 dest=OUT/'03_niche_geometry'/f'{sec}__all_effects.parquet'
 if dest.exists() and done('relations_'+sec):return
 t,a,g=load_section(sec);meta=json.loads((OUT/'01_geometry'/sec/'metadata.json').read_text());n=len(t);xy=a['xy_pitch'];physical=t.physical_component.to_numpy();rows=[]
 for k in KS:
  src=[pd.read_parquet(OUT/'02_hotspots'/sec/f'{v[0]}_k{k}.parquet') for v in VARIABLES[:2]];cells={j:pd.read_parquet(OUT/'02_hotspots'/sec/f'CELL_{j:02d}_k{k}.parquet') for j in NON}
  for track in TRACKS:
   targets={('NICHE',f'niche{i}'):np.flatnonzero(t.accepted_label.to_numpy()==i) for i in range(9)}
   targets.update({('CELL',TYPES[j]):d.node_index.to_numpy()[d[f'{track}_hot'].to_numpy()] for j,d in cells.items()})
   for layer,sd in zip(LAYERS,src):
    valid=np.isfinite(sd.Z)&np.isfinite(sd.p_sim_native);E=sd.node_index.to_numpy()[valid];Hmask=np.zeros(n,bool);Hmask[sd.node_index.to_numpy()[sd[f'{track}_hot']]]=True;comps=np.full(n,-1,int);comps[sd.node_index]=sd[f'{track}_hot_component'];local=np.zeros(n,bool);local[sd.node_index]=sd.geometry_local
    for target_kind,target_id,j in TARGETS:
     tid=targets[(target_kind,target_id)];tm=np.zeros(n,bool);tm[tid]=True;d=nearest(xy,physical,E,tid);reason='NO_LOCAL_SOURCE_SUPPORT' if not len(E) else ('TARGET_NICHE_ABSENT' if target_kind=='NICHE' else 'NO_TARGET_HOTSPOT') if not len(tid) else 'NO_SOURCE_HOTSPOT' if not Hmask.any() else 'READY';targetvals=a['q05'][tid,j] if target_kind=='CELL' else np.array([]);tvmean=targetvals.mean() if len(targetvals) else np.nan
     for domain in (['MAIN','GEOMETRY_LOCAL_SUPPORT'] if k==10 else ['MAIN']):
      dom=np.ones(len(E),bool) if domain=='MAIN' else local[E];ee=E[dom];dd=d[dom];hh=Hmask.copy();hh[~local]=False if domain!='MAIN' else hh[~local];unreach=~np.isfinite(dd);reach=~unreach
      common=dict(dataset_id=meta['dataset_id'],section_id=sec,analysis_unit_uid=meta['analysis_unit_uid'],analysis_unit_type=meta['analysis_unit_type'],known_patient=meta['known_patient'],layer=layer,statistic='SOURCE_G',k=k,track=track,target_kind=target_kind,target_id=target_id,analysis_domain=domain,n_source_supported=len(ee),n_source_hot=int(hh[ee].sum()),source_hot_fraction=hh[ee].mean() if len(ee) else np.nan,source_formal_components=int(sd.loc[sd[f'{track}_hot'],f'{track}_hot_component'].nunique()),n_E_unreachable=int(unreach.sum()),n_H_unreachable=int(hh[ee[unreach]].sum()),reachable_E_fraction=reach.mean() if len(ee) else np.nan,target_relative_composition_mean=tvmean,LOW_ABSOLUTE_RELATIVE_COMPOSITION=bool(np.isfinite(tvmean) and tvmean<.001),n_target_physical_components=len(np.unique(physical[tid])),n_source_physical_components=len(np.unique(physical[ee])))
      for metric,domainmask in [('O',np.ones(len(ee),bool)),('D',reach),('D_out',reach&~tm[ee])]:
       ids=ee[domainmask];values=tm[ids].astype(float) if metric=='O' else dd[domainmask];why='TARGET_UNREACHABLE' if metric!='O' and reason=='READY' and not reach.any() else reason;rr=metric_record(metric,ids,hh,values,tm,t,comps,len(tid),why);rows.append(dict(**common,**rr))
 d=pd.DataFrame(rows);d.to_parquet(dest,index=False,compression='zstd');checkpoint('relations_'+sec);log('relations complete '+sec)

def family(kind,metric):
 if kind=='NICHE':return 'NICHE_OUTSIDE_SECONDARY' if metric=='D_out' else 'NICHE_PRIMARY'
 return {'D':'CELL_PRIMARY_DISTANCE','O':'CELL_OVERLAP_SECONDARY','D_out':'CELL_OUTSIDE_SECONDARY'}[metric]

KEY=['scope','layer','k','track','target_kind','target_id','metric','analysis_domain']
def units_and_tests(sections):
 allunits=[];alltests=[]
 for scope in SCOPES:
  sub=sections if scope=='ALL_DATASETS' else sections[sections.dataset_id.eq(scope)]
  for key,gg in sub.groupby(KEY[1:],sort=False):
   base=dict(zip(KEY,[scope]+list(key)));eligible=gg[gg.eligible & gg.effect.notna() & gg.analysis_unit_uid.ne('')];unitrows=[]
   for uid,uu in eligible.groupby('analysis_unit_uid',sort=False):
    r=dict(**base,analysis_unit_uid=uid,effect=uu.effect.mean(),n_sections=len(uu),known_patient=bool(uu.known_patient.all()),analysis_unit_type='|'.join(sorted(uu.analysis_unit_type.unique())),sections='|'.join(sorted(uu.section_id)),datasets='|'.join(sorted(uu.dataset_id.unique())),stratified_effect=uu.stratified_effect.mean(),component_equal_effect=uu.component_equal_effect.mean());unitrows.append(r);allunits.append(r)
   for subset in ['ALL_POLICY_UNITS','VERIFIED_PATIENTS_ONLY']:
    u=[x for x in unitrows if subset=='ALL_POLICY_UNITS' or x['known_patient']];stats=signed([x['effect'] for x in u]);stats['n_verified_units']=sum(x['known_patient'] for x in u);stats['n_other_units']=len(u)-stats['n_verified_units'];stats['n_sections']=sum(x['n_sections'] for x in u)
    if base['k']!=10 or base['analysis_domain']!='MAIN':stats.update(p=np.nan,W=np.nan,actual_method='DESCRIPTIVE_SENSITIVITY_NO_TEST',status='DESCRIPTIVE_SENSITIVITY' if len(u) else 'INSUFFICIENT_UNITS')
    alltests.append(dict(**base,subset=subset,family=family(base['target_kind'],base['metric']),planned_n=FAMILY_N[family(base['target_kind'],base['metric'])],n_sections_descriptive=int(gg.effect.notna().sum()),n_sections_source_hot=int(gg.n_source_hot.gt(0).sum()),n_sections_target_present=int(gg.n_target_full.gt(0).sum()),mean_source_hot_fraction=gg.source_hot_fraction.mean(),source_large_fraction_any=bool(gg.source_hot_fraction.gt(.8).any()),source_single_component_any=bool(gg.source_formal_components.eq(1).any()),n_sections_source_missing=int(gg.n_source_hot.eq(0).sum()),n_sections_any_unreachable=int(gg.n_E_unreachable.gt(0).sum()),mean_reachable_fraction=gg.reachable_E_fraction.mean(),mean_target_relative_composition=gg.target_relative_composition_mean.mean(),rare_target_any=bool(gg.LOW_ABSOLUTE_RELATIVE_COMPOSITION.any()),**stats))
 u=pd.DataFrame(allunits);r=pd.DataFrame(alltests);r['q_scope']=np.nan;r['q_seven_datasets']=np.nan;r['q_all_eight_scopes']=np.nan
 main=(r.k==10)&r.analysis_domain.eq('MAIN')
 for key,ii in r[main].groupby(['scope','track','subset','family']).groups.items():assert len(ii)==FAMILY_N[key[-1]];r.loc[ii,'q_scope']=bh(r.loc[ii,'p'],FAMILY_N[key[-1]])
 for key,ii in r[main].groupby(['track','subset','family']).groups.items():
  fam=key[-1];r.loc[ii,'q_all_eight_scopes']=bh(r.loc[ii,'p'],8*FAMILY_N[fam]);jj=[i for i in ii if r.loc[i,'scope']!='ALL_DATASETS'];r.loc[jj,'q_seven_datasets']=bh(r.loc[jj,'p'],7*FAMILY_N[fam])
 return u,r

def matched_layers(s):
 s=s[(s.k==10)&s.analysis_domain.eq('MAIN')];out=[];ur=[]
 for scope in SCOPES:
  sub=s if scope=='ALL_DATASETS' else s[s.dataset_id.eq(scope)]
  for key,gg in sub.groupby(['track','target_kind','target_id','metric'],sort=False):
   base=dict(zip(['track','target_kind','target_id','metric'],key));a=gg[gg.layer.eq(LAYERS[0])&gg.eligible];b=gg[gg.layer.eq(LAYERS[1])&gg.eligible];m=a.merge(b,on=['section_id','analysis_unit_uid'],suffixes=('_spot','_epi'),validate='one_to_one');g=m.groupby('analysis_unit_uid')[['effect_spot','effect_epi']].mean();n=len(g)
   for uid,row in g.iterrows():ur.append(dict(scope=scope,**base,analysis_unit_uid=uid,**row.to_dict()))
   x=signed(g.effect_spot);y=signed(g.effect_epi);out.append(dict(scope=scope,**base,family=family(base['target_kind'],base['metric']),n_matched_sections=len(m),n_matched_units=n,mean_spot=x['mean_effect'],mean_epi=y['mean_effect'],median_spot=x['median_effect'],median_epi=y['median_effect'],IQR_spot=x['IQR'],IQR_epi=y['IQR'],p_spot=x['p'],p_epi=y['p'],same_direction_units=int(((g.effect_spot*g.effect_epi)>1e-24).sum()),opposite_direction_units=int(((g.effect_spot*g.effect_epi)<-1e-24).sum()),nearzero_units=int(((abs(g.effect_spot)<=1e-12)|(abs(g.effect_epi)<=1e-12)).sum()),paired_mean_direction_agreement=bool(n and x['mean_effect']*y['mean_effect']>1e-24)))
 r=pd.DataFrame(out);r['q_spot']=np.nan;r['q_epi']=np.nan
 for key,ii in r.groupby(['scope','track','family']).groups.items():
  pp=np.concatenate([r.loc[ii,'p_spot'],r.loc[ii,'p_epi']]);qq=bh(pp,FAMILY_N[key[-1]]);r.loc[ii,'q_spot']=qq[:len(ii)];r.loc[ii,'q_epi']=qq[len(ii):]
 r['both_matched_layers_q05']=r.q_spot.lt(.05)&r.q_epi.lt(.05);tab('05_sensitivity/two_score_layers_matched_units.tsv',r);tab('05_sensitivity/two_score_layers_matched_unit_rows.tsv',ur);return r

def summarize():
 sec=pd.concat([pd.read_parquet(OUT/'03_niche_geometry'/f'{s}__all_effects.parquet') for s in COV.section_id],ignore_index=True);u,r=units_and_tests(sec);matched=matched_layers(sec)
 for kind,folder,prefix in [('NICHE','03_niche_geometry','niche'),('CELL','04_cell_hotspot_geometry','cell')]:
  tab(f'{folder}/{prefix}_section_effects.tsv',sec[sec.target_kind.eq(kind)]);tab(f'{folder}/{prefix}_unit_effects.tsv',u[u.target_kind.eq(kind)]);tab(f'{folder}/{prefix}_scope_tests.tsv',r[r.target_kind.eq(kind)])
  for track in TRACKS:tab(f'{folder}/{prefix}_scope_tests_{track}.tsv',r[r.target_kind.eq(kind)&r.track.eq(track)])
 tab('05_sensitivity/k_sensitivity.tsv',r[r.analysis_domain.eq('MAIN')]);tab('05_sensitivity/geometry_local_support.tsv',r[r.analysis_domain.eq('GEOMETRY_LOCAL_SUPPORT')]);tab('05_sensitivity/hotspot_threshold_sensitivity.tsv',r[r.k.eq(10)&r.analysis_domain.eq('MAIN')]);tab('05_sensitivity/verified_patients_only.tsv',r[r.subset.eq('VERIFIED_PATIENTS_ONLY')]);tab('05_sensitivity/stratified_background.tsv',sec[sec.k.eq(10)&sec.analysis_domain.eq('MAIN')]);tab('05_sensitivity/component_weighting.tsv',sec[sec.k.eq(10)&sec.analysis_domain.eq('MAIN')][['dataset_id','section_id','layer','track','target_kind','target_id','metric','effect','component_equal_effect','n_hot_components']]);tab('07_synthesis/unavailable_results.tsv',sec[~sec.eligible]);sec.to_parquet(OUT/'07_synthesis/all_section_effects.parquet',index=False);u.to_parquet(OUT/'07_synthesis/all_unit_effects.parquet',index=False);r.to_parquet(OUT/'07_synthesis/all_scope_tests.parquet',index=False)
 ev=r[r.k.eq(10)&r.analysis_domain.eq('MAIN')&r.subset.eq('ALL_POLICY_UNITS')].merge(matched,on=['scope','track','target_kind','target_id','metric','family'],suffixes=('','_paired'),validate='many_to_one');ev['statistic']='SOURCE_G';ev['direction']=np.where(ev.mean_effect.abs()<=1e-12,'NEAR_ZERO',np.where(ev.mean_effect>0,'POSITIVE','NEGATIVE'));ev.loc[ev.mean_effect.isna(),'direction']='NA';ev['interpretation']='NO_SUPPORTED_ASSOCIATION';sig=ev.q_scope.lt(.05)
 for metric,pos,neg in [('O','EXPLORATORY_OVERLAP_ENRICHMENT','EXPLORATORY_OVERLAP_DEPLETION'),('D','EXPLORATORY_SPATIAL_SEPARATION','EXPLORATORY_PROXIMITY'),('D_out','EXPLORATORY_EXTERNAL_SEPARATION','EXPLORATORY_EXTERNAL_PROXIMITY')]:
  ev.loc[sig&ev.metric.eq(metric),'interpretation']=np.where(ev.loc[sig&ev.metric.eq(metric),'mean_effect']>0,pos,neg)
 ev.loc[ev.p.isna(),'interpretation']=ev.loc[ev.p.isna(),'status'];ev['caveat']='Exploratory after v7; does not negate previous nonsignificant continuous-ring associations; no cell-contact/causal claim.'
 for kind,name in [('NICHE','niche_geometry_evidence_all9'),('CELL','cell_geometry_evidence_all29')]:tab(f'07_synthesis/{name}.tsv',ev[ev.target_kind.eq(kind)])
 tab('07_synthesis/dataset_evidence_all.tsv',ev);checkpoint('M06_RELATIONS_STATISTICS_COMPLETE')

def full():
 assert done('M03_ALL_HOTSPOTS_COMPLETE')
 from concurrent.futures import ProcessPoolExecutor
 with ProcessPoolExecutor(max_workers=2) as pool:list(pool.map(section_effects,COV.section_id))
 summarize()
if __name__=='__main__':full()
