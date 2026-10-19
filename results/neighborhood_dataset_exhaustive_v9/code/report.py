from core import *

def summarize(final=False):
 c=pd.read_csv(OUT/'combinations.tsv',sep='\t',dtype={'combo_id':str});summary=[];primary=[];aux=[];rings=[]
 for row in c.itertuples():
  f=OUT/'combinations'/row.combo_id;state=json.loads((f/'state.json').read_text()) if (f/'state.json').exists() else {'stage':'REGISTERED'};sel=json.loads((f/'selection.json').read_text()) if (f/'selection.json').exists() else {};complete=(f/'COMPLETE.json').exists();record=dict(combo_id=row.combo_id,datasets=row.datasets,n_datasets=row.n_datasets,contains_inhouse=row.contains_inhouse,execution_order=row.execution_order,status='COMPLETE' if complete else state['stage'],K=sel.get('K'),n_accepted=sel.get('accepted_n'),stability_min=sel.get('stability_min'),seconds_refit=sel.get('seconds'),bytes=bytes_used(f) if f.exists() else 0);summary.append(record)
  for branch,name in [('REFIT','spotted_tests.parquet'),('FROZEN_LABEL','frozen_tests.parquet')]:
   if not (f/name).exists():continue
   t=pd.read_parquet(f/name);overall=t[t.scope.eq('ALL_DATASETS')].copy();overall['combo_id']=row.combo_id;overall['K']=sel.get('K') if branch=='REFIT' else 9;overall['contains_inhouse']=row.contains_inhouse;overall['datasets']=row.datasets
   ih=t[t.scope.eq(INHOUSE)].set_index(['layer','track','target_kind','target_id','metric','subset'])
   for col in ['n_units','n_sections','mean_effect','median_effect','n_positive','n_negative']:
    overall['inhouse_'+col]=[ih.loc[(r.layer,r.track,r.target_kind,r.target_id,r.metric,r.subset),col] if (r.layer,r.track,r.target_kind,r.target_id,r.metric,r.subset) in ih.index else (0 if col in ['n_units','n_sections','n_positive','n_negative'] else np.nan) for r in overall.itertuples()]
   overall['inhouse_estimable']=overall.inhouse_n_units.gt(0);overall['inhouse_same_direction']=overall.inhouse_mean_effect*overall.mean_effect>1e-24;overall['inhouse_actual_support']=overall.inhouse_estimable&overall.inhouse_same_direction
   mask=(overall.branch=='REFIT')&(overall.track=='SOURCE_P05')&(overall.target_kind=='NICHE')&overall.metric.isin(['O','D'])&(overall.subset=='ALL_POLICY_UNITS')
   primary.append(overall[mask]);aux.append(overall[~mask])
  if (f/'ring_tests.parquet').exists():
   t=pd.read_parquet(f/'ring_tests.parquet');t=t[t.scope=='ALL_DATASETS'].copy();t['combo_id']=row.combo_id;rings.append(t)
 ss=pd.DataFrame(summary);finished=int(ss.status.eq('COMPLETE').sum());all_done=finished==99;
 if all_done:
  for cid in ss.combo_id:
   for stage in ['NICHE_FROZEN','RING_DONE','SPOTTED_MAIN_DONE','FROZEN_DONE']:assert valid(OUT/'combinations'/cid,stage)
 correction_status='FINAL' if all_done else 'PROVISIONAL';table(OUT/'combination_summary.tsv',ss)
 if primary:
  p=pd.concat(primary,ignore_index=True);p['p_raw']=p.p;p['q_global_BH']=correction(p.p);p['q_global_BY']=correction(p.p,True);p['global_status']=correction_status;p['M_primary_current']=len(p);p['endpoint']=p.metric;p['target']=p.target_id;p['effect']=p.mean_effect;p['N']=p.n_units;p['GLOBAL_BY_SUPPORTED']=p.q_global_BY.lt(.05);p['result_class']=np.where(p.q_global_BH.lt(.05),'GLOBAL_BH_HIT',np.where(p.q_within_BH.lt(.05),'WITHIN_COMBO_EXPLORATORY','NO_GLOBAL_HIT'));p.loc[p.p.isna(),'result_class']='NOT_ESTIMABLE'
  if all_done:assert len(p)==int((ss.K*4).sum())
  table(OUT/'primary_relations.parquet',p);hits=p[p.q_global_BH.lt(.05)].sort_values('q_global_BH');within=p[p.q_within_BH.lt(.05)&~p.q_global_BH.lt(.05)].sort_values('q_within_BH');table(OUT/'global_primary_hits.tsv',hits);table(OUT/'within_combo_exploratory_hits.tsv',within);table(OUT/'inhouse_priority_candidates.tsv',p[p.contains_inhouse].sort_values(['q_global_BH','inhouse_actual_support'],ascending=[True,False]));table(OUT/'without_inhouse_controls.tsv',p[~p.contains_inhouse].sort_values('q_global_BH'))
  qc=p.groupby('combo_id').agg(primary_relations=('p','size'),testable=('p','count'),within_hits=('q_within_BH',lambda a:int(a.lt(.05).sum())),global_hits=('q_global_BH',lambda a:int(a.lt(.05).sum())),BY_hits=('q_global_BY',lambda a:int(a.lt(.05).sum())))
  ss=ss.merge(qc,on='combo_id',how='left');table(OUT/'combination_summary.tsv',ss)
 else:p=pd.DataFrame();hits=within=pd.DataFrame()
 if aux:
  a=pd.concat(aux,ignore_index=True);a['q_global_BH']=np.nan;a['q_global_BY']=np.nan
  for fam,ix in a.groupby('family_id').groups.items():a.loc[ix,'q_global_BH']=correction(a.loc[ix,'p']);a.loc[ix,'q_global_BY']=correction(a.loc[ix,'p'],True)
  a['global_status']=correction_status;table(OUT/'auxiliary_relations.parquet',a)
  if all_done:assert len(a[(a.branch=='FROZEN_LABEL')&(a.track=='SOURCE_P05')&(a.target_kind=='NICHE')&a.metric.isin(['O','D'])&(a.subset=='ALL_POLICY_UNITS')])==3564
 if rings:
  r=pd.concat(rings,ignore_index=True);r['q_global_BH']=np.nan;r['q_global_BY']=np.nan
  for family,ix in r.groupby('family_id').groups.items():r.loc[ix,'q_global_BH']=correction(r.loc[ix,'p']);r.loc[ix,'q_global_BY']=correction(r.loc[ix,'p'],True)
  r['global_status']=correction_status;table(OUT/'ring_auxiliary_relations.parquet',r)
 used=bytes_used();cl=pd.read_csv(OUT/'cleanup_manifest.tsv',sep='\t');cleaned=int(cl.bytes.sum()) if len(cl) else 0
 review_done=(OUT/'review/COMPLETE.json').exists();state='EXECUTION_COMPLETE' if all_done and review_done else 'EXECUTION_MAIN_COMPLETE' if all_done else 'SEARCH_INCOMPLETE';cmd=f"cd '{ROOT}' && PYTHONDONTWRITEBYTECODE=1 '{V8}/env/bin/python' '{OUT}/code/run.py' --resume"
 js(OUT/'RESUME.json',dict(status=state,completed=finished,total=99,global_status=correction_status,unsettled=ss[ss.status.ne('COMPLETE')][['combo_id','status']].to_dict('records'),resume_command=cmd,config_hash=config_hash(),input_hash=sha(OUT/'audit/input_manifest.tsv'),used_bytes=used,cleaned_bytes=cleaned,time=time.time()))
 lines=['# v9 穷举运行汇总',f'执行状态：{state}。已完成 {finished}/99；含 inhouse {int(ss[ss.contains_inhouse].status.eq("COMPLETE").sum())}/57，不含 {int(ss[~ss.contains_inhouse].status.eq("COMPLETE").sum())}/42。',f'全局校正状态：{correction_status}。当前 REFIT 主关系 {len(p)} 条；全局 BH 阳性 {len(hits)} 条；仅组合内 BH 阳性 {len(within)} 条。未完成时这些数量不能当最终搜索结论。',f'新增文件逻辑字节 {used:,}；累计清理本轮可再生临时文件 {cleaned:,} 字节。预算 3 GiB。','主分析 family 为 REFIT × SOURCE_P05 × k10 × ALL_POLICY_UNITS × 两层 × 全部冻结 niche × O/D。NA 仅表示已成功计算后的支持不足；技术失败保留未结算状态。','inhouse 纳入与实际支持分列；实际支持表示有可估计单位且其单位平均效应与总体同向，不等于独立显著或已核实独立患者。','99 个组合共享样本，属于穷举探索，不能称为独立重复验证。BH/BY 不消除赢家选择偏差，未核实身份和空间零模型限制仍然存在。','O>0 重叠增加，O<0 重叠减少；D<0 更近，D>0 更远。距离单位 spot pitch。不能解释为细胞接触、招募或主动免疫排斥。']
 if len(hits):lines+=['','当前主命中（完整表见 global_primary_hits.tsv）：',hits[['combo_id','layer','target_id','metric','mean_effect','n_units','q_within_BH','q_global_BH','q_global_BY','inhouse_n_units','inhouse_actual_support']].head(20).to_csv(sep='\t',index=False)]
 if (OUT/'review/selected_relations.tsv').exists():
  selected=pd.read_csv(OUT/'review/selected_relations.tsv',sep='\t',dtype={'combo_id':str});lines+=['','重点关系（未经赢家选择偏差校正；完整单位与敏感性表在 review）：',selected[['combo_id','target_id','layer','metric','n_units','mean_effect','q_within_BH','q_global_BH','q_global_BY','inhouse_n_units','inhouse_mean_effect','inhouse_actual_support']].to_csv(sep='\t',index=False)]
  for cid in selected.combo_id.unique():
   dest=OUT/'review'/cid
   if (dest/'leave_one_dataset_out.tsv').exists():
    loo=pd.read_csv(dest/'leave_one_dataset_out.tsv',sep='\t');rev=loo[loo.removed_dataset.ne('NONE')&loo.same_direction.eq(False)];lines.append(cid+'：删除单dataset后方向改变或不可估计的诊断 '+str(len(rev))+' 条；这是固定标签影响诊断，不是独立验证。')
    if len(rev):lines.append(rev[['target_id','layer','metric','removed_dataset','n_units','mean_effect']].to_csv(sep='\t',index=False))
   if (dest/'leave_one_unit_out.tsv').exists():
    uu=pd.read_csv(dest/'leave_one_unit_out.tsv',sep='\t');lines.append(cid+'：删除单单位后方向改变或不可估计 '+str(int(uu.same_direction.eq(False).sum()))+' 条。')
    top=uu.sort_values('absolute_mean_change',ascending=False).drop_duplicates(['target_id','layer','metric']);lines.append('最大单位均值影响（不把均值影响等同于 signed-rank P 的单一驱动）：\n'+top[['target_id','layer','metric','removed_unit','mean_remaining','absolute_mean_change','same_direction']].to_csv(sep='\t',index=False))
   if (dest/'refit_frozen_comparison.tsv').exists():
    z=pd.read_csv(dest/'refit_frozen_comparison.tsv',sep='\t');lines.append(cid+'：组成匹配后的新旧标签方向一致 '+str(int(z.same_direction.sum()))+'/'+str(len(z))+'；两分支源域与有效N分别保留。')
 lines.append('详细图表状态：'+('COMPLETE' if review_done else '等待全部主分析完成后自动生成；或生成中，见 review。'))
 (OUT/'RUN_SUMMARY.md').write_text('\n\n'.join(lines)+'\n');return all_done
if __name__=='__main__':summarize()
