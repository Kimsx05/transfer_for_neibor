from run_v3_3 import *
S=pd.read_csv(O/'tables/section_summary.tsv',sep='\t');B=pd.read_csv(O/'tables/tissue_component_summary.tsv',sep='\t');sm=pd.read_csv(O/'tables/small_roi_decisions.tsv',sep='\t');a=pd.read_csv(O/'tables/author_registration_summary.tsv',sep='\t');weak=sm[sm.v3_2_weak_roi.eq(True)];W=weak.groupby(['decision','final_locations']).agg(rois=('old_roi_id','size'),spots=('n_spots','sum')).reset_index();W.to_csv(O/'tables/weak_roi_disposition_summary.tsv',sep='\t',index=False)
dom=[];support=[]
for sec in M.index:
 p=O/'by_dataset'/M.loc[sec,'dataset_id']/sec;d=pd.read_parquet(p/'spot_regions_v3_3.parquet');e=pd.read_csv(p/'lesion_boundary_edges.tsv.gz',sep='\t');q=e[e.boundary_type.eq('candidate_observed_outer')].groupby('boundary_support').size()
 for k,v in q.items():support.append({'section_id':sec,'boundary_support':k,'n_edges':v})
 if sec.startswith('GSE319536_'):
  for did,z in d.groupby('domain_component_id'):
   dom.append({'section_id':sec,'domain_component_id':did,'n_spots':len(z),'v3_1_inside':z.v3_1_lesion_location.eq('inside').sum(),'v3_3_inside':z.lesion_location.eq('inside').sum(),'coarse_inside_n':z.author_coarse_reference.eq('author_coarse_inside').sum(),'coarse_outside_n':z.author_coarse_reference.eq('author_coarse_outside').sum(),'reference_boundary_conflict_n':z.reference_status.eq('possible_author_conflict_in_registration_boundary_zone').sum(),'changed_n':z.lesion_location.ne(z.v3_1_lesion_location).sum()})
 if sec=='GSE238145_B16':
  z=d[d.barcode.eq('AGATAACACGGTATAA-1_4')];z[['section_id','barcode','array_row','array_col','v3_1_lesion_location','v3_2_lesion_location','lesion_location','domain_local_evidence_score','domain_shared_evidence_score','assignment_basis','lesion_roi_id']].to_csv(O/'tables/B16_singleton_correction.tsv',sep='\t',index=False)
pd.DataFrame(dom).to_csv(O/'tables/author_domain_review.tsv',sep='\t',index=False);pd.DataFrame(support).to_csv(O/'tables/candidate_boundary_support.tsv',sep='\t',index=False)
# Focused table so local support is not collapsed into an apparent validation score.
F=S.set_index('section_id');rr=[]
for sec in ['GSE319536_BC01','GSE319536_BC13','GSE319536_BC18','GSE285715_T3','GSE238145_B16']:
 x=F.loc[sec];rr.append(f'| {sec} | {x.inside:,} | {x.outside:,} | {x.author_inside_to_outside:,} | {x.author_outside_to_inside:,} | {x.candidate_edges:,} | {x.exact_estimable:,} | {x.band_200_determined:,} |')
size=sum(p.stat().st_size for p in O.rglob('*') if p.is_file())/1024**3
report=f'''# v3.3 实际结果

已完成 **{len(S)}/63张、7个dataset、{S.n_spots.sum():,}个原spot** 的候选位置修订、几何重算和PNG/PDF图。[全队列入口](index.html)。本轮没有重跑上游模型，没有新增MP10/TME生物学分析，v3.1/v3.2保留。

**恢复的是可检查的计算候选分区，不是把几何可计算升级为病理外缘确证。** 红蓝具有相同的计算估计地位，支持不足保存在独立字段/图中。最终{S.candidate_edges.gt(0).sum()}张有可观测候选界面；N1、N2、N3没有被强造margin。

## 灰化如何修正

| 项目 | 最终数量 |
|---|---:|
| inside候选 | {S.inside.sum():,} |
| outside候选 | {S.outside.sum():,} |
| 原有unlocalized保留 | {S.unlocalized.sum():,}（{100*S.unlocalized.sum()/S.n_spots.sum():.3f}%） |
| 从v3.2撤销状态恢复到v3.1候选outside、主位置未改变v3.1 | {S.restored_from_v3_2_only.sum():,} |
| 相比v3.1真正改变的位置 | {S.changed_from_v3_1.sum():,} |
| 其中局部小ROI修正 | {S.local_small_roi_changed.sum():,} |
| 其中作者粗参考校准 | {S.author_changed.sum():,} |
| 可估计精确候选距离 | {S.exact_estimable.sum():,} |
| 可判断200μm分带 | {S.band_200_determined.sum():,} |
| 其中仅分带可判断、精确距离仍NA | {S.band_only_200.sum():,} |

上述“恢复”是标签历史归因，不是覆盖率优化目标；恢复点中可另外具有作者支持。剩余局部参考冲突通过支持图表达，不再把普通基质/免疫组织统一变灰。原有unlocalized/QC状态没有机械恢复。[逐类变更摘要](tables/position_changes.tsv)

## 弱借域小ROI与B16

87个点逐对象审查后，**27个ROI/43点**的完整观测外邻域全为候选outside，自身局部证据不足，撤销独立ROI并归outside；**30个ROI/44点**因邻域缺测/混合等不足以可靠改位，保留v3.1 inside，标为不稳健小候选。本轮没有把它们统一改灰，也不将小候选称为卫星灶。其余有局部依据的小岛保留。局部决定、最终位置及附近margin计数见[逐ROI结果](tables/small_roi_decisions.tsv)。

B16 `AGATAACACGGTATAA-1_4` 已从v3.2新增inside恢复为outside，独立ROI取消。坐标(89,175)，局部分数0.463318、旧共享分数0.415650，六邻在v3.1均为outside。**上一轮“完整inside邻域支持”的解释错误；实际是改用局部分数过线造成的新孤点。** 原分数保留，不把取消远处借分变成新增单点病灶准入规则。[该点记录](tables/B16_singleton_correction.tsv) · [B16图形](by_dataset/GSE238145/GSE238145_B16/index.html)

## 作者参考实际参与双侧检查

16张GSE319536完成H&E对应尝试；14张达到固定粗配准标准。BC05、BC17未达到完整标准，不转移标签，候选位置仍保留。配准使用原H&E结构，未用位置标签拟合；图像边缘、细碎参考区域、Others及非二元区室不自动归内/外。校准只采用完整误差半径内仍保持粗区域共识的核心，**最终数量以本报告为准，较早试跑尚未完整排除配准误差带的计数不属于最终结果**。

最终作者校准共改变{S.author_changed.sum():,}点，包括 **{S.author_inside_to_outside.sum():,}点原inside→outside** 和 **{S.author_outside_to_inside.sum():,}点原outside→inside**。这些是实际修改数量，不是作者barcode真值一致率，未计算Dice/IoU。

| 切片 | 最终inside | 最终outside | 作者支持内→外 | 作者支持外→内 | 候选边段 | 精确距离点 | 200μm分带点 |
|---|---:|---:|---:|---:|---:|---:|---:|
'''+'\n'.join(rr)+f'''

- **BC01**：缓冲后的39个正常粗参考点均与原候选outside一致，已实质接入位置支持字段；不能为了制造红→蓝数量而越过配准误差带。46个旧outside按作者肿瘤相关核心改为inside。另1点为局部弱ROI修正。图边缘的红区/正常参考潜在冲突仍明确标记，未伪精确转移。
- **BC13**：53个旧红点改为outside，131个旧蓝点改为inside。作者的Normal Stroma核心与肿瘤相关核心同时参与，纠正了双向遗漏/覆盖；尚在误差带或非二元类别中的红区仍是候选，并非默认正确。
- **BC18**：179个旧outside改为inside，分别来自D8的121点、D3的53点和D15的5点，均由新配准后的粗区域核心支持；不再沿用503点D8特例，也不把剩余区域不经参考检查全部当作已证实正常。每个未转移点保留候选估计及作者边缘/冲突状态，Others不赋予作者二元真值。大尺度图仍支持部分原分界位于肿瘤相关区室内，当前候选边缘不能全部称宏观病灶外缘。
- **T3**：主体范围沿用v3.1，仅1个弱小ROI点按完整邻域修正。最终1895 inside/1446 outside。恢复可计算候选边界，但标注它可能反映病灶内部组成/状态差异；广泛EPI及bladder cancer来源不能证明整片inside。没有把Tumor-419a伪映射为T3。

[作者配准表](tables/author_registration_summary.tsv) · [域级双侧复核](tables/author_domain_review.tsv) · [来源/限制简表](tables/reference_comparison.tsv)。原作者图是本轮校准信息，不是独立验证；EPI一致性也是共享表达证据检查。

## 哪些边界可用、哪些仍需限制解释

60张的候选几何可以计算。**使用对象应称“当前计算候选ROI外缘”，不是已确认组织学侵袭前沿。** BC01/BC13/BC14/BC15具有作者肿瘤和正常粗区域信息，对区域身份提供更强背景；BC18主要有肿瘤相关区域背景。严格局部350μm双侧粗参考邻域条件下，本轮没有边段达到双侧粗参考近邻条件，因而未把任何边段称为作者精确轮廓验证。这不阻断候选距离。

其他队列以既有分子/空间候选为主；T3需特别限制宏观外缘解释，NPJ样本别名和Inhouse来源未被自行补全。逐片候选边段数、距离可估数与分带数见[切片汇总](tables/section_summary.tsv)，支持程度见[边段支持简表](tables/candidate_boundary_support.tsv)。

单一区域块单独保留：全inside例子包括CanDis Sample7:T2（202点）、Sample6:T3（39点）、Sample12:T3（37点），以及BC11:T4/T8/T12（51/86/36点）；均无观察到的外缘，距离NA，不称deep core。N1整片1417点outside；N2/N3保留既有未定位部分。其他全outside组织块也是候选位置，不等于新证明其全部正常。[完整单一区域清单](tables/single_region_components.tsv)

## 交付与验证

每片5套核心PNG/PDF，GSE319536另有配准和作者对照；7个dataset总览，逐片HTML入口。表达使用原共同色标，KRT18缺失仍unavailable。完整逐点Parquet保留原始分子值、三版本位置、支持字段、新ROI/边段/距离；另有[全队列标签](tables/all_spot_labels_v3_3.parquet)与TSV.gz。

检查原barcode/坐标/数值、原unlocalized、同块邻接、重建边段与最近距离、无候选界面不产生margin、B16指定点及所有HTML图片链接；身份映射沿用v3.1修正，B5/B18保持同P1，T1/T2/T3保持未核验独立切片单位。验证日志及`tables/checks.tsv`记录实际结果。

复现命令见[METHODS](METHODS.md)，实际脚本位于`scripts/`。本目录约{size:.2f}GiB；两片并发、每进程内部线程为1，不加载全队列表达大矩阵。
'''
(O/'RESULTS.md').write_text(report)
print('Report written; totals',S.select_dtypes('number').sum().to_dict())
