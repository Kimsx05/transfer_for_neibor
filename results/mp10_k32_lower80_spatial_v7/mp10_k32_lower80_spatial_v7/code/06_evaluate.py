from common import *
from scipy.stats import spearmanr
O=R/'05_evaluation';labels=json.loads((R/'00_manifest/labels.json').read_text());bags=pd.read_csv(R/'02_pseudospots/bag_manifest.tsv',sep='\t');pan=pd.read_csv(R/'02_pseudospots/panel_references.tsv',sep='\t');truth=pd.read_csv(R/'02_pseudospots/truth_31.tsv',sep='\t').set_index('pseudospot_id')
def weight(d):
 w=1/d.groupby(['source_group','sample'])['sample'].transform('size').to_numpy(float);ng=d[['source_group','sample']].drop_duplicates().groupby('source_group').size();w/=d.source_group.map(ng).to_numpy(float);return w/w.sum()
def metrics(d,t,p,lo=None,hi=None):
 mask=np.isfinite(d[t])&np.isfinite(d[p]);q=d[mask];out={'n_spots_total':len(d),'n_evaluable':len(q),'n_samples':q['sample'].nunique(),'n_source_groups':q.source_group.nunique(),'n_datasets':q.dataset.nunique()}
 if not len(q):return out
 w=weight(q);x=q[t].to_numpy(float);y=q[p].to_numpy(float);xm=np.sum(w*x);ym=np.sum(w*y);var=np.sum(w*(x-xm)**2);slope=np.sum(w*(x-xm)*(y-ym))/var if np.ptp(x)>0 and var>1e-16 else np.nan;out.update(MAE=float(np.sum(w*abs(y-x))),RMSE=float(np.sqrt(np.sum(w*(y-x)**2))),bias=float(np.sum(w*(y-x))),true_mean=float(xm),predicted_mean=float(ym),calibration_slope=slope,calibration_intercept=float(ym-slope*xm) if np.isfinite(slope) else np.nan,spearman_descriptive=float(spearmanr(x,y).statistic) if np.ptp(x)>0 and np.ptp(y)>0 else np.nan)
 if lo:out['coverage90']=float(np.sum(w*((x>=q[lo])&(x<=q[hi]))));out['mean_interval_width']=float(np.sum(w*(q[hi]-q[lo])))
 return out
models=['DEV']+(['DEV_OLD'] if (outmap('DEV_OLD')/'RESULT_STATUS.json').exists() else [])
for model in models:
 path=outmap(model)
 if not (path/'RESULT_STATUS.json').exists():continue
 status=json.loads((path/'RESULT_STATUS.json').read_text());kl='Epi_K32' if model=='DEV_OLD' else K;ll='Epi_other' if model=='DEV_OLD' else L;d=bags.copy();mean=pd.read_csv(path/'abundance_mean.tsv.gz',sep='\t',index_col=0).reindex(d.pseudospot_id)
 for nm in ['mean','median','q05','q95']:
  ab=pd.read_csv(path/f'abundance_{nm}.tsv.gz',sep='\t',index_col=0).reindex(d.pseudospot_id);d['K32_'+nm]=ab[kl].to_numpy();d['background_reference_'+nm]=ab[ll].to_numpy()
 der=pd.read_csv(path/'posterior_derived.tsv.gz',sep='\t').rename(columns={'spot_id':'pseudospot_id'});d=d.merge(der,on='pseudospot_id',validate='one_to_one');d['vCAF_mean']=mean.vCAF.to_numpy();d['vCAF_true']=truth.loc[d.pseudospot_id,'vCAF'].to_numpy();d['numerical_status']=status['numerical_status'];d['model']=model;d['background_reference_definition']='historical all frozen nonK32 Epi_other' if model=='DEV_OLD' else 'nonK32 at or below TRAIN q80';wt(d,O/(model+'_spot_results.tsv.gz'));pr=pan[pan.mapped].merge(d,on=['pseudospot_id','sample','source_group','dataset','repeat','n_epi_design','n_dose_design'],validate='many_to_one');wt(pr,O/(model+'_panel_results.tsv.gz'))
 rows=[];targets=[('K32_count','n_K32_true','K32_mean','K32_q05','K32_q95'),('epithelial_count','n_epi_true','Epi_total_mean','Epi_total_q05','Epi_total_q95'),('K32_fraction_epi','K32_fraction_epi_true','I_state_mean','I_state_q05','I_state_q95'),('total_cells','n_total_true','total_cells_mean','total_cells_q05','total_cells_q95'),('vCAF_count','vCAF_true','vCAF_mean',None,None)]
 sets={'A_C_unique':d[d.pseudospot_id.isin(pr.loc[pr.panel.isin(['A','C']),'pseudospot_id'])],'ALL_unique':d,**{p:pr[pr.panel.eq(p)].drop_duplicates('pseudospot_id') for p in ['A','B','C']}}
 ac0=sets['A_C_unique'];sets['A_C_K32_positive']=ac0[ac0.n_K32_true.gt(0)];sets['A_C_K32_zero']=ac0[ac0.n_K32_true.eq(0)];sets['A_C_ratio_interior']=ac0[ac0.K32_fraction_epi_true.gt(0)&ac0.K32_fraction_epi_true.lt(1)]
 for scope,q in sets.items():
  for level in ['ALL','dataset','source_group','sample','n_epi_true','n_K32_true']:
   groups=[('ALL',q)] if level=='ALL' else q.groupby(level,observed=True)
   for value,g in groups:
    for target,t,p,lo,hi in targets:rows.append({'model':model,'scope':scope,'level':level,'unit':value,'target':target,'numerical_status':status['numerical_status'],**metrics(g,t,p,lo,hi)})
 wt(pd.DataFrame(rows),O/(model+'_recovery_metrics.tsv'))
 # Full 29-type accounting, by sample and source/dataset-aware weighting.
 rows=[]
 for lab in [x for x in labels if x not in [K,L]]:
  q=d.copy();q['true']=truth.loc[q.pseudospot_id,lab].to_numpy();q['pred']=mean[lab].to_numpy()
  for scope,ids in [('ALL_unique',d.pseudospot_id),('C',pr.loc[pr.panel.eq('C'),'pseudospot_id'])]:
   z=q[q.pseudospot_id.isin(ids)];rows.append({'model':model,'scope':scope,'label':lab,**metrics(z,'true','pred')})
 wt(pd.DataFrame(rows),O/(model+'_29type_allocation.tsv'))
 # Low-background distributions: distinct C bags with actual source results.
 c=sets['C'];neg=[]
 for fields in [['n_epi_true'],['dataset','n_epi_true'],['sample','source_group','dataset','n_epi_true']]:
  for vals,g in c.groupby(fields,observed=True):
   vals=vals if isinstance(vals,tuple) else (vals,);w=weight(g);order=np.argsort(g.K32_mean.to_numpy());cw=np.cumsum(w[order]);vsort=g.K32_mean.to_numpy()[order];wmed=vsort[min(np.searchsorted(cw,.5),len(vsort)-1)];wq95=vsort[min(np.searchsorted(cw,.95),len(vsort)-1)];neg.append({**dict(zip(fields,vals)),'level':'+'.join(fields),'n_bags':len(g),'n_groups':g.source_group.nunique(),'K32_mean_groupweighted':float(np.sum(w*g.K32_mean)),'K32_median_groupweighted':float(wmed),'K32_q95_groupweighted':float(wq95),'K32_median_technical':g.K32_mean.median(),'K32_q95_technical':g.K32_mean.quantile(.95),'K32_max':g.K32_mean.max(),'epi_pred_mean':float(np.sum(w*g.Epi_total_mean)),'epi_true':g.n_epi_true.iloc[0],'vCAF_pred_mean':float(np.sum(w*g.vCAF_mean)),'vCAF_true_mean':float(np.sum(w*g.vCAF_true))})
 wt(pd.DataFrame(neg),O/(model+'_low_background_distribution.tsv'))
 # A/B common backgrounds, no duplicated weighting of common zero-dose bag.
 ab=pr[pr.panel.isin(['A','B'])];keys=['sample','source_group','dataset','repeat','n_epi_design','n_dose_design'];pairs=ab[ab.panel.eq('A')].merge(ab[ab.panel.eq('B')],on=keys,suffixes=('_A','_B'),validate='one_to_one');assert (pairs.background_hash_A==pairs.background_hash_B).all();pairs['B_minus_A_K32_mean']=pairs.K32_mean_B-pairs.K32_mean_A;pairs['B_minus_A_I_state_mean']=pairs.I_state_mean_B-pairs.I_state_mean_A;wt(pairs,O/(model+'_A_B_matched_backgrounds.tsv.gz'))
 responses=[]
 for fields in [['panel','n_epi_design','sample','source_group','dataset'],['panel','n_epi_design','dataset']]:
  for vals,g in pr[pr.n_epi_design.isin([4,12])].groupby(fields,observed=True):
   vals=vals if isinstance(vals,tuple) else (vals,)
   for x in ['high_fraction_epi_true','mean_ucell_epi_true']:
    for y in ['K32_mean','I_state_mean']:
     mm=metrics(g,x,y);rr={k:v for k,v in mm.items() if k not in ['MAE','RMSE','bias']};rr['response_slope']=rr.pop('calibration_slope',np.nan);rr['response_intercept']=rr.pop('calibration_intercept',np.nan);responses.append({**dict(zip(fields,vals)),'level':'+'.join(fields),'x':x,'y':y,**rr})
 wt(pd.DataFrame(responses),O/(model+'_state_response.tsv'))
 # Matched end-minus-zero contrasts, source counts are explicit; no per-spot tests.
 endpoint=[]
 for (panel,sample,group,ds,rep,ne),g in pr[pr.panel.isin(['A','B'])&pr.n_epi_design.isin([4,12])].groupby(['panel','sample','source_group','dataset','repeat','n_epi_design']):
  lo=g[g.n_dose_design.eq(0)];hi=g[g.n_dose_design.eq(ne)]
  if len(lo)==1 and len(hi)==1:
   lo=lo.iloc[0];hi=hi.iloc[0];endpoint.append({'panel':panel,'sample':sample,'source_group':group,'dataset':ds,'repeat':rep,'n_epi_design':ne,'K32_mean_high_minus_zero':hi.K32_mean-lo.K32_mean,'I_state_high_minus_zero':hi.I_state_mean-lo.I_state_mean,'true_high_fraction_change':hi.high_fraction_epi_true-lo.high_fraction_epi_true,'mean_ucell_change':hi.mean_ucell_epi_true-lo.mean_ucell_epi_true,'background_identical':hi.background_hash==lo.background_hash})
 ep=pd.DataFrame(endpoint);wt(ep,O/(model+'_endpoint_contrasts.tsv'));ek=['sample','source_group','dataset','repeat','n_epi_design'];epair=ep[ep.panel.eq('A')].merge(ep[ep.panel.eq('B')],on=ek,suffixes=('_A','_B'),validate='one_to_one');wt(epair,O/(model+'_endpoint_contrasts_A_B_matched.tsv'))
 erows=[]
 for ne,q in epair.groupby('n_epi_design'):
  ww=weight(q);row={'n_epi':ne,'n_samples':q['sample'].nunique(),'n_groups':q.source_group.nunique(),'n_paired_repeats':len(q)}
  for col in ['K32_mean_high_minus_zero_A','K32_mean_high_minus_zero_B','I_state_high_minus_zero_A','I_state_high_minus_zero_B']:
   v=q[np.isfinite(q[col])];row[col]=float(np.sum(weight(v)*v[col])) if len(v) else np.nan;row[col+'_n_evaluable_pairs']=len(v)
  erows.append(row)
 wt(pd.DataFrame(erows),O/(model+'_matched_endpoint_groupweighted.tsv'))
js(O/'EVALUATION_STATUS.json',{'time_epoch':time.time(),'models':models,'binary_threshold':'NONE_BY_DESIGN','weights':'equal samples within conservative source group, then equal source groups; spots within sample equal','targets':'exact frozen K32 and broad MP10high state kept distinct; B withheld not assigned K32 truth','source_independence':'not verified patients; no per-spot significance tests','old_v6_comparison':'historical only unless paired DEV_OLD generated on identical v7 counts'})
