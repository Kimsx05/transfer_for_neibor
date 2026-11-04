from aggregate import *
def briefrow(r):return f'{r.candidate_type} / {r.program}: r={r.pearson_r:.3f}, n={r.n_spots}, {r.distance_min_um:.0f}–{r.distance_max_um:.0f} μm, spatial CI [{r.spatial_ci_low:.3f}, {r.spatial_ci_high:.3f}], P={r.p_pearson:.3g}, across-section q={r.q_pearson_across_sections:.3g}'
def main():
 cov=pd.read_csv(O/'full_sections.tsv',sep='\t').set_index('section_id');context=pd.read_parquet(O/'tables/context85_effects.parquet');fr=pd.read_csv(O/'tables/front_only_niche_all.tsv.gz',sep='\t');targets=pd.read_csv(O/'tables/candidate_targets.tsv.gz',sep='\t');out=[]
 for sec in M.section_id:
  r=cov.loc[sec];a=pd.read_parquet(O/'tables'/f'{sec}_pearson_q.parquet');de=pd.read_parquet(O/'tables'/f'{sec}_delta_q.parquet') if (O/'tables'/f'{sec}_delta_q.parquet').exists() else pd.DataFrame();z=context[context.section_id.eq(sec)];normal=r.source_context=='adjacent_normal_source';pop='normal_Epi' if normal else 'P95';v=a[a.population.eq(pop)&a.range_type.eq('full_range')&a.score_method.eq('AddModuleScore_SCT_data')&a.program.isin(CORE)&a.pearson_r.notna()];v=v[v.n_spots.ge(10)];text=f'# {sec}\n\n来源: {r.source_context}；analysis unit: {r.analysis_unit_uid}；patient: {r.patient_uid}。\n\n组织spots {r.n_spots}，有效Epi {r.n_Epi}，内部Epi {r.n_internal}，冻结P95 {r.n_P95}。精确前沿Epi {r.n_front_exact_Epi}，前沿内部Epi {r.n_front_internal}，R12标注Epi {r.n_niche_labelled_Epi}，未标注 {r.n_niche_unlabelled_Epi}。\n\n## 候选距离原始Pearson\n\n'
  neg=pos='不可估计'
  if len(v):
   lo=v.sort_values('pearson_r').iloc[0];hi=v.sort_values('pearson_r').iloc[-1];neg=briefrow(lo);pos=briefrow(hi);text+=f'- 核心面板较低r：{neg}。\n- 核心面板较高r：{pos}。\n\n'
  else:text+='此片没有达到n≥10的核心面板相关展示；全部可计算n≥2结果仍保留在完整长表，不因覆盖不足删除。\n\n'
  dd=de[de.score_method.eq('AddModuleScore_SCT_data')&de.range_type.eq('common_support')&de.program.isin(CORE)&de.delta_r.notna()].copy() if len(de) else pd.DataFrame()
  if len(dd):
   dd['abs']=dd.delta_r.abs();row=dd.sort_values('abs',ascending=False).iloc[0];text+=f'核心面板最大绝对共同范围差异：{row.candidate_type} / {row.program}，r_P95={row.r_P95:.3f} (n={row.n_P95})，r_nonP95={row.r_nonP95:.3f} (n={row.n_nonP95})，delta={row.delta_r:.3f}，空间CI [{row.spatial_ci_low:.3f}, {row.spatial_ci_high:.3f}]，空间P={row.p_spatial_fisher_z:.3g}，跨片空间q={row.q_spatial_across_sections:.3g}，共同范围{row.common_min_um:.0f}–{row.common_max_um:.0f} μm。效应选择是探索性展示，不按显著性筛片。\n\n'
  if normal:text+='邻近正常不定义肿瘤P95或其组间差；全Epi候选距离与程序作为参照保留。\n\n'
  text+='## MP10中心及邻域\n\n'
  cc=z[z.analysis.isin(['B1_R0','B1_R1','B1_R2','B1_R3_R5','B1_R6_R10'])&z.adjustment.eq('adjusted')&z.program.isin(CORE)&z.effect.notna()]
  for ring in ['R0','R1','R2','R3_R5','R6_R10']:
   zz=cc[cc.analysis.eq('B1_'+ring)]
   if len(zz):
    lo=zz.loc[zz.effect.idxmin()];hi=zz.loc[zz.effect.idxmax()];text+=f'- {ring}: 核心面板较高 {hi.program}, β={hi.effect:.3f}；较低 {lo.program}, β={lo.effect:.3f}。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。\n'
  text+='\n## 前沿及niche\n\n'
  ff=z[z.analysis.isin(['C1_all_Epi_continuous','C2_P95_enrichment','C_active_front_MP10','A_front_vs_internal'])&z.adjustment.eq('adjusted')&z.program.isin([RAW,'P95_probability','PURAM_PEMT_TOP30','FLORES_CRISPR_PRO_INVASION_10'])]
  for rr in ff.itertuples():text+=f'- {rr.analysis} / {rr.program}: β={getattr(rr,"effect",np.nan):.3f}, CI [{getattr(rr,"low",np.nan):.3f}, {getattr(rr,"high",np.nan):.3f}], n={rr.n}, 新85上下文族q={getattr(rr,"q_context85_within_section",np.nan):.3g}。\n'
  nf=fr[fr.section_id.eq(sec)&fr.difference_percentage_points.notna()].copy()
  if len(nf):
   nf['abs']=nf.difference_percentage_points.abs();rr=nf.sort_values('abs',ascending=False).iloc[0];text+=f'\n仅前沿带[-200,0] μm，最大绝对R12比例差为niche {int(rr.niche)}：{rr.difference_percentage_points:.2f}个百分点，CI [{rr.ci_low_percentage_points:.2f}, {rr.ci_high_percentage_points:.2f}]，标注P95 n={rr.n_P95}、非P95 n={rr.n_nonP95}，跨片q={rr.q_spatial_across_sections:.3g}。所有九类和未标注覆盖同时保留；这里没有重新定义Top5%。\n'
  text+='\n## 限制与入口\n\n'
  tt=targets[targets.section_id.eq(sec)];weak=tt[tt.absolute_abundance_median.lt(1)];text+=f'29候选中 {len(weak)} 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。'
  if r.n_P95<30 and not normal:text+=' 本片P95少于30，相关取值仍报告，但其稳定性有限，不能据较大|r|认定样本特异表型。'
  text+=' 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。\n\n三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。\n'
  (O/'reports'/f'{sec}.md').write_text(text);out.append(dict(section_id=sec,core_low_r_example=neg,core_high_r_example=pos,n_core_negative_in_report_population=int(v.pearson_r.lt(-.1).sum()) if len(v) else 0,n_core_positive_in_report_population=int(v.pearson_r.gt(.1).sum()) if len(v) else 0,n_candidate_low_estimated_abundance=len(weak),candidate_distance_status='complete_29_with_reasons',neighborhood_status='previous80_plus_new5',front_only_niche_status='computed_or_explicit_coverage_NA',score_status='85_SCT_AMS_and85_UCell'))
 c=cov.reset_index().drop(columns=[c for c in pd.DataFrame(out).columns if c!='section_id' and c in cov],errors='ignore').merge(pd.DataFrame(out),on='section_id',validate='1:1');save(c,O/'full_sections.tsv');print('63 narrative notes complete',flush=True)
if __name__=='__main__':main()
