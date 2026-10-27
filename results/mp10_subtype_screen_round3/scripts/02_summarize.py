from common import *
edge=read('total_node_edges.tsv.gz');cov=read('node_coverage.tsv');order=read('node_order.tsv');obj=read('object_index.tsv');rows=[];contrasts=[];urows=[];partners=[]
for o in obj.itertuples():
 z=edge[edge.object_id.eq(o.object_id)];cc=cov[cov.object_id.eq(o.object_id)].set_index('node');univ=order[order.network_level.eq(o.network_level)].sort_values('node_order').node.tolist();actual=cc.index[cc.computed].tolist();tme=[x for x in actual if x not in EPI];mat={k:z.pivot(index='source',columns='target',values=k).reindex(index=univ,columns=univ) for k in ['W','C','R']}
 for node in univ:
  cv=cc.loc[node];rec=dict(object_id=o.object_id,scope=o.scope,network_level=o.network_level,node=node,parent=cv.parent,computed=cv.computed,status=cv.status,n_TME_ranked=len(tme),n_computable_partners=len(actual) if cv.computed else 0,n_matched_specimens=cv.n_matched_specimens,n_matched_recorded_IDs=cv.n_matched_recorded_IDs,n_raw_cells=cv.n_raw_cells,n_sampled_cells=cv.n_sampled_cells,n_contributing_specimens=cv.n_contributing_specimens,n_contributing_recorded_IDs=cv.n_contributing_recorded_IDs)
  for measure,a in mat.items():
   inc=a.loc[actual,node].sum(min_count=1) if cv.computed else np.nan;out=a.loc[node,actual].sum(min_count=1) if cv.computed else np.nan;self_=a.loc[node,node] if cv.computed else np.nan
   for k,v in [('incoming_all',inc),('outgoing_all',out),('self',self_),('incoming_other',inc-self_),('outgoing_other',out-self_)]:rec[k+'_'+measure]=v
  rec['incoming_self_fraction']=rec['self_W']/rec['incoming_all_W'] if rec['incoming_all_W']>0 else np.nan;rec['outgoing_self_fraction']=rec['self_W']/rec['outgoing_all_W'] if rec['outgoing_all_W']>0 else np.nan;rows.append(rec)
 for dr in ['epi_to_TME','TME_to_epi']:
  U={ep:(mat['W'].loc[ep,tme].sum(min_count=1) if dr=='epi_to_TME' else mat['W'].loc[tme,ep].sum(min_count=1)) for ep in EPI}
  urows.append(dict(object_id=o.object_id,scope=o.scope,network_level=o.network_level,direction=dr,n_common_TME=len(tme),U_low=U[EPI[0]],U_high=U[EPI[1]],U_difference=U[EPI[1]]-U[EPI[0]]))
  for node in [x for x in univ if x not in EPI]:
   cv=cc.loc[node];rec=dict(object_id=o.object_id,scope=o.scope,network_level=o.network_level,node=node,parent=cv.parent,direction=dr,computed=cv.computed,status=cv.status,n_TME_ranked=len(tme),n_matched_specimens=cv.n_matched_specimens,n_matched_recorded_IDs=cv.n_matched_recorded_IDs,n_target_specimens=cv.n_contributing_specimens,n_target_recorded_IDs=cv.n_contributing_recorded_IDs,U_low=U[EPI[0]],U_high=U[EPI[1]])
   for measure,a in mat.items():
    low=a.loc[EPI[0],node] if dr=='epi_to_TME' else a.loc[node,EPI[0]];high=a.loc[EPI[1],node] if dr=='epi_to_TME' else a.loc[node,EPI[1]];rec['low_'+measure]=low;rec['high_'+measure]=high;rec['delta_'+measure]=high-low
   rec['low_share']=rec['low_W']/U[EPI[0]] if U[EPI[0]]>0 else np.nan;rec['high_share']=rec['high_W']/U[EPI[1]] if U[EPI[1]]>0 else np.nan;rec['delta_share']=rec['high_share']-rec['low_share'];rec['W_R_direction_agree']=np.sign(rec['delta_W'])==np.sign(rec['delta_R']) if cv.computed else np.nan;rec['transition']='NE' if not cv.computed else 'zero_to_supported' if rec['low_W']==0 and rec['high_W']>0 else 'supported_to_zero' if rec['low_W']>0 and rec['high_W']==0 else 'both_zero' if rec['low_W']==0 and rec['high_W']==0 else 'both_supported';contrasts.append(rec)
 if o.network_level=='fine':
  for node in univ:
   for dr in ['incoming','outgoing']:
    src=node if dr=='incoming' else 'Treg';tar='Treg' if dr=='incoming' else node;rec=dict(object_id=o.object_id,scope=o.scope,direction=dr,partner=node,category=node if node in EPI else 'Treg_self' if node=='Treg' else 'other_TME',computed=bool(cc.loc['Treg','computed'] and cc.loc[node,'computed']))
    for measure,a in mat.items():rec[measure]=a.loc[src,tar]
    partners.append(rec)
g=pd.DataFrame(rows);c=pd.DataFrame(contrasts)
for field in ['incoming_all_W','outgoing_all_W','incoming_other_W','outgoing_other_W','incoming_all_C','outgoing_all_C','incoming_other_C','outgoing_other_C']:
 g['rank_'+field]=g[~g.node.isin(EPI)&g.computed].groupby('object_id')[field].rank(method='min',ascending=False)
for field in ['low_W','high_W','delta_W','low_C','high_C','delta_C','low_share','high_share','delta_share']:
 c['rank_'+field]=c.groupby(['object_id','direction'])[field].rank(method='min',ascending=False)
save(g,'global_incoming_outgoing.tsv');save(c,'epi_allTME_total_contrasts.tsv');save(urows,'epi_state_totals.tsv');pp=pd.DataFrame(partners);save(pp,'treg_partner_breakdown.tsv');pc=pp.groupby(['object_id','scope','direction','category'])[['W','C','R']].sum(min_count=1).reset_index();save(pc,'treg_partner_categories.tsv')
tr=g[g.network_level.eq('fine')&g.node.eq('Treg')].copy()
for dr in ['epi_to_TME','TME_to_epi']:
 a=c[c.network_level.eq('fine')&c.node.eq('Treg')&c.direction.eq(dr)].drop(columns=['scope','network_level','node','parent','direction','computed','status','n_TME_ranked','n_matched_specimens','n_matched_recorded_IDs']);a=a.rename(columns={k:dr+'_'+k for k in a if k!='object_id'});tr=tr.merge(a,on='object_id',validate='one_to_one')
save(tr,'treg_total_summary.tsv')
# Separate the five dataset directions from the pooled ALL overview, never average probabilities.
c['strength_direction']=np.where(~c.computed,'NE',np.where(c.delta_W.gt(0),'high_greater',np.where(c.delta_W.lt(0),'low_greater','zero_difference')));c['share_direction']=np.where(c.delta_share.isna(),'NE',np.where(c.delta_share.gt(0),'high_greater',np.where(c.delta_share.lt(0),'low_greater','zero_difference')));save(c,'cross_dataset_direction_long.tsv')
summary=[]
for keys,a in c.groupby(['network_level','node','parent','direction']):
 ds=a[a.scope.ne('ALL')];pool=a[a.scope.eq('ALL')].iloc[0];summary.append(dict(zip(['network_level','node','parent','direction'],keys))|dict(n_computable_datasets=int(ds.computed.sum()),n_high_greater=int(ds.strength_direction.eq('high_greater').sum()),n_low_greater=int(ds.strength_direction.eq('low_greater').sum()),n_zero_difference=int(ds.strength_direction.eq('zero_difference').sum()),n_NE=int((~ds.computed).sum()),n_datasets_at_least_two_matched_IDs=int(ds.n_matched_recorded_IDs.ge(2).sum()),n_high_share=int(ds.share_direction.eq('high_greater').sum()),n_low_share=int(ds.share_direction.eq('low_greater').sum()),n_W_R_direction_discordant=int((ds.computed & ds.W_R_direction_agree.eq(False)).sum()),high_datasets=';'.join(ds.loc[ds.strength_direction.eq('high_greater'),'scope']),low_datasets=';'.join(ds.loc[ds.strength_direction.eq('low_greater'),'scope']),pooled_ALL_delta_W=pool.delta_W,pooled_ALL_delta_share=pool.delta_share,pooled_ALL_is_independent_replication=False))
save(summary,'cross_dataset_direction_summary.tsv')
# Check LR contribution sums against supported W differences, not raw R differences.
lr=read('treg_LR_contribution_decomposition.tsv.gz');decomp=lr.groupby(['object_id','scope','direction']).agg(sum_contribution=('contribution','sum'),sum_raw_difference=('raw_difference','sum')).reset_index();decomp=decomp.merge(c[c.network_level.eq('fine')&c.node.eq('Treg')][['object_id','direction','delta_W','delta_R']],on=['object_id','direction']);decomp['W_absolute_error']=abs(decomp.sum_contribution-decomp.delta_W);decomp['R_absolute_error']=abs(decomp.sum_raw_difference-decomp.delta_R);assert decomp.W_absolute_error.max()<1e-10 and decomp.R_absolute_error.max()<1e-10;save(decomp,'treg_contribution_parity.tsv')
# Full partner categories must sum to the single shared Treg global totals.
for row in tr.itertuples():
 for dr in ['incoming','outgoing']:
  val=pc.loc[pc.object_id.eq(row.object_id)&pc.direction.eq(dr),'W'].sum(min_count=1);assert np.isclose(val,getattr(row,dr+'_all_W'),equal_nan=True)
# Supplement: a single LR never validates the aggregate network at patient level.
sup=read('supplementary_TIGIT_receivers_network.tsv');save(sup[sup.TME_type.eq('Treg')],'supplementary_NECTIN2_TIGIT.tsv');ex=read('supplementary_NECTIN2_TIGIT_expression.tsv');ev=ex[ex.coverage_eligible&ex.genes_observed_in_dataset&ex.expression_delta.notna()];patients=ev.groupby(['dataset','patient_uid_recorded']).agg(n_specimens=('sample_uid','size'),expression_delta=('expression_delta','mean'),aux_product_delta=('aux_product_delta','mean'),low_support=('low_support','max'),high_support=('high_support','max')).reset_index();save(patients,'supplementary_NECTIN2_TIGIT_patient_differences.tsv');pr=[]
for ds,z in patients.groupby('dataset'):
 vals=z.expression_delta.to_numpy();loo=[np.delete(vals,i).mean() for i in range(len(vals))] if len(vals)>1 else []
 pr.append(dict(dataset=ds,n_specimens=int(z.n_specimens.sum()),n_recorded_IDs=len(z),n_positive=int((vals>0).sum()),n_negative=int((vals<0).sum()),n_zero=int((vals==0).sum()),mean_expression_difference=vals.mean(),median_expression_difference=np.median(vals),LOO_min=min(loo) if loo else np.nan,LOO_max=max(loo) if loo else np.nan,LOO_direction_stable=bool(loo and (min(loo)>0 or max(loo)<0))))
save(pr,'supplementary_NECTIN2_TIGIT_dataset_expression.tsv')
tig=read('supplementary_TIGIT_receivers_expression.tsv');tg=tig[tig.n_cells.gt(0)].groupby(['dataset','node']).agg(n_specimens=('sample_uid','nunique'),n_recorded_IDs=('patient_uid_recorded','nunique'),n_cells=('n_cells','sum'),specimen_equal_TIGIT_mean=('TIGIT_mean','mean'),specimen_equal_TIGIT_detection=('TIGIT_detection','mean'),n_specimens_target30=('n_cells',lambda x:x.ge(30).sum())).reset_index();save(tg,'supplementary_TIGIT_receiver_summary.tsv')
par=read('aggregation_parity.tsv');assert par[['max_W_vs_v2','max_C_vs_v2','max_R_vs_v2']].to_numpy().max()<1e-10
assert edge.loc[edge.status.ne('computed'),['W','C','R']].isna().all().all()
print('ALL Treg summary');print(tr[['scope','incoming_all_W','outgoing_all_W','rank_incoming_all_W','rank_outgoing_all_W','n_TME_ranked','n_matched_recorded_IDs']].to_string(index=False));print('Aggregation parity maxima');print(par.select_dtypes('number').max().to_string())
