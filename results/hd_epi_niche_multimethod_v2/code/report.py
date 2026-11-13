"""Build a lightweight delivery from completed results; never fits a model."""
from common import *
import zipfile
from datetime import datetime
from zoneinfo import ZoneInfo

def tab(d):
 d=d.copy()
 for c in d.select_dtypes(include='floating'):d[c]=d[c].round(4)
 d=d.fillna('NA')
 return '| '+' | '.join(d.columns)+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join(map(str,row))+' |' for row in d.itertuples(index=False,name=None))
def read(name):
 f=R/'tables'/name
 return pd.read_csv(f,sep='\t') if f.exists() else pd.DataFrame()

m=read('selected_models.tsv');cov=read('coverage.tsv');matrix=pd.read_csv(R/'run_matrix.tsv',sep='\t')
for i,v in matrix.iterrows():
 ref='both' if v.method=='CellCharter' else v.reference;tag=f'{v.method}_{ref}_r{v.radius}_{v.representation}'
 matrix.loc[i,'status']='full_assignments_complete' if len(m) and tag in set(m.tag) else 'not_complete'
matrix.to_csv(R/'run_matrix.tsv',sep='\t',index=False)
valid_table=read('loader_validation.tsv')
checks=['IDs_aligned','original_counts_unchanged','original_QC_unchanged','nonanchors_not_labeled','loader_executed']
valid=bool(len(valid_table)==5 and valid_table[checks].astype(str).apply(lambda x:x.str.lower().isin(['true','1'])).to_numpy().all())
complete=bool((matrix.status=='full_assignments_complete').all() and valid and (R/'INTERPRETATION.md').exists())
wide=cov.pivot_table(index=['tag','section','anchors','evaluable'],columns='source',values='n',aggfunc='sum',fill_value=0).reset_index()
wide['assigned']=wide.anchors-wide.get('unresolved',0)
wide.to_csv(R/'tables/coverage_wide.tsv',sep='\t',index=False)
section_rows=[]
for s in SECTIONS:
 d=wide[wide.section==s];aa=d[d.tag.str.startswith('hellinger_kmeans')&d.tag.str.endswith('_NONE')]
 row={'section':s,'anchors':int(d.anchors.iloc[0]),'niche_train':int(d.native_train.iloc[0]),'evaluable_min':int(d.evaluable.min()),'evaluable_max':int(d.evaluable.max()),'A_NONE_assigned_min':int(aa.assigned.min()),'A_NONE_assigned_max':int(aa.assigned.max())}
 for r in [50,200]:
  cc=d[d.tag.eq(f'CellCharter_both_r{r}_scVI')].iloc[0]
  row[f'CC{r}_assigned']=int(cc.assigned);row[f'CC{r}_unresolved']=int(cc.unresolved)
 section_rows.append(row)
ss=pd.DataFrame(section_rows);ss.to_csv(R/'tables/section_completion.tsv',sep='\t',index=False)
recs=[]
for v in m.itertuples():
 if v.tag.startswith('hellinger_kmeans'):
  status='working_coarse_with_reference_sensitivity' if v.tag.endswith('NONE') else 'sensitivity_only'
  why='K2 composition axis; preserve continuous P and paired references/scales' if v.tag.endswith('NONE') else 'Transductive; lower source association changes partition; iteration cap reached'
 elif v.tag.startswith('adapted_leiden'):
  status='diagnostic_not_endorsed_fine_niche';why='Exact duplicate split fixed; fragmentation/source dependence or degenerate selected partition remains'
 elif v.tag.startswith('MENDER'):
  status='radial_support_sensitivity' if 'per_ring' in v.tag else 'exploratory_radial_representation_not_endorsed_partition'
  why='Support-driven pure-Epi splitting; selected 200um model is near-singleton partition'
 else:
  status='feature_sensitivity_only' if 'without' in v.tag else 'exploratory_molecular_environment'
  why='Raw neighbor molecular support; feature/source sensitivity prevents validated fine identity claim'
 recs.append(dict(tag=v.tag,n_clusters=int(v.n_clusters),recommendation=status,reason=why))
recommendations=pd.DataFrame(recs);recommendations.to_csv(R/'tables/model_recommendations.tsv',sep='\t',index=False)
summary=dict(updated=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),all_requested_branches_complete=complete,predefined_branches=24,conditional_branches=1,selected_result_branches=len(m),branch_count_note='Selected full-assignment results; auxiliary K2/candidate fits are not additional biological branches',counts_and_QC_reloaded=valid,section_coverage=section_rows,coverage_detail='tables/coverage_wide.tsv',selected_models='tables/selected_models.tsv',recommendations='tables/model_recommendations.tsv',working_coarse='A NONE K2 with continuous P and both references',exploratory_molecular_K={'50um':5,'200um':9},identical_full_inputs_split=0,MPBC2='never niche fitting; WEAK transductive; CellCharter official query adaptation',unfinished=matrix[matrix.status!='full_assignments_complete'].to_dict('records'),remaining_limitations=['reference dependence','MPBC tissue masks unknown','fine graph partitions not biologically endorsed','MENDER200 selected partitions degenerate','Harmony iteration cap','CellCharter feature and source sensitivity'],output_root=str(R),MP10_association_executed=False)
intro='已完成四方法的 24 个预设结果分支及 1 个条件敏感性，均输出五片全部冻结 anchors 的赋值或明确未定状态。主要工作层仍建议 A NONE K=2 + 连续全组成；CellCharter 作为分子辅助，B/MENDER 的细分保留诊断意义。'
lines=['# Visium HD 上皮相关 niche 四方法比较','',summary['updated'],'',intro,'','## 完成范围与固定分母','',tab(ss),'','evaluable_min/max 是全部25分支的可评价范围；A_NONE_assigned_min/max 是双参考×双尺度四组的已赋值范围。训练外数量=已赋值数−该片训练数；不同方法的具体 source 与未定原因见 coverage_wide.tsv 和逐中心压缩表。niche 描述 anchors/bin 环境，不是真实细胞数。','',f'五片原始对象重读和挂接验证：{valid}。非 anchors 保持 niche 缺失，原 counts/QC/坐标未修改。原525,552个anchor ID不变。','', '## 结果解释与工作建议','',(R/'INTERPRETATION.md').read_text() if (R/'INTERPRETATION.md').exists() else '空间复核尚未完成。','','## 实际运行矩阵与参数','',tab(matrix),'','selected_parameter 是有限网格中的计算候选；其生物学使用资格由 model_recommendations.tsv 单独记录，不因入选就认可为工作 niche。','',tab(m[[c for c in ['tag','selected_parameter','n_clusters','n_train','evaluable','unresolved','seed_ARI','ARI_to_full_on_buffer_test'] if c in m]]),'','## 输入、实现与验证入口','', 'METHODS.md 描述每种原生步骤及适配边界；config.yaml、versions.json、input_manifest.tsv 记录实际环境、来源和参数。旧METHODS的unpruned文字与代码不一致，本轮继承实际pruned字段。T_v2统一80000个，旧BLCA50仅1个不在共同可评价集合内的中心同片替换；training_overlap.tsv记录与四旧名单的差异。','', 'implementation_checks.json、real_ROI_MENDER_equivalence.json、adapted_graph_pilot.json、final_validation.tsv、loader_validation.tsv 保存针对半径、中心排除、Unknown分母、等价图/环/latent聚合及最终ID挂接的核验。initial scVI 400-epoch warm-up不适配有限epochs的问题已修正，旧结果隔离在diagnostics/pre_warmup_correction，不属于正式结果。详见logs/warmup_correction.json。','', '重计算在23:00–07:00完成；白天仅做读对象验证、图例修正、解释和打包。原counts/Seurat/参考均未复制。绘图使用实际Arial；PNG300dpi，PDF文字/图例矢量、密集点层栅格。原始单基因/H&E图链接复用上游，不声称是新生成数据。','', '## 文件导航','', '| 文件/目录 | 内容 |\n| --- | --- |\n| tables/frozen_anchors.parquet, T_v2.tsv, T_expr.parquet | 中心、共同niche训练、表达训练名单 |\n| all_anchor_assignments/ + tables/assignment_index.tsv | 25个逐中心结果分区；source明确训练/预测/投影/未定 |\n| tables/coverage_wide.tsv, section_completion.tsv | 全部组合逐片覆盖与精简汇总 |\n| models/ | 工作模型、必要预处理、训练表示及K2辅助层 |\n| tables/candidate_metrics.tsv, common_K2_comparisons.tsv | 参数稳定性及共同K2对照 |\n| tables/method_comparisons.tsv, partition_cross_tables.tsv | 同T与全量对应，附各自分母 |\n| tables/niche_composition.tsv, MENDER_ring_composition.tsv | 全部群原始累计/径向组成，未删小群 |\n| tables/quality_marker_support.tsv, technical_association.tsv | 原始邻居marker与质量/来源诊断 |\n| tables/identical_input*.tsv, native_graph_diagnostics.tsv | 重复完整输入和图诊断 |\n| tables/integration_scope.tsv, expression_fit_diagnostics.tsv | MPBC2参与范围、表达训练状态 |\n| tables/fixed_ROI_index.tsv, raw_ROI_figure_sources.tsv | 固定ROI坐标及上游原始证据入口 |\n| figures/ | 五片双参考四方法×双尺度、NONE/WEAK、ROI和诊断 |', '', '图中每个模型的数字ID独立。全片图绘制全部中心和全部标签；群数>12时PNG键列出该片每面板最大4群，并明确提示完整legend-only PDF。组成热图最多显示30个最大群，覆盖率见figure_display_scope.tsv；全量数值无截断。ROI使用相应全片图色键，空ROI未换选。','', '## 读取与复跑','', '```r',f'root <- "{R}"', 'ext <- file.path(dirname(root), "hd_bin8_banksy_extension_v1")','source(file.path(ext, "code/load_result.R"))','x <- load_extension(ext, "UC_SCC")','source(file.path(root, "code/load_result.R"))','x <- attach_hd_niches_v2(x, root, "UC_SCC") # x: 既有RNA对象，不覆盖任何旧列','```','', '```python','from importlib.util import spec_from_file_location, module_from_spec',f'root = {str(R)!r}', 'spec = spec_from_file_location("niche_loader", root + "/code/load_result.py")','loader = module_from_spec(spec); spec.loader.exec_module(loader)','metadata = loader.load_niche_metadata(root, "UC_SCC")','adata = loader.attach_niches(adata, root, "UC_SCC")','```','', 'loader按唯一bin ID挂接niche_v2_前缀字段；全片原对象的非anchors保持label缺失。全量压缩表本身的label=-1为未定，candidate仍保留。依赖上游项目原对象与本目录.venv，未另存25份表达对象。','', '`./rerun.sh` 核对缓存后按夜间窗口续跑；`./rerun.sh --report-only` 仅更新已有结果的报告和review，不启动拟合。review.zip不含完整表达、参考、全图邻接或训练缓存，未上传外部仓库。','']
(R/'REPORT.md').write_text('\n'.join(lines))
(R/'completion_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False))
(R/'COMPLETION_CN.md').write_text('# 中文完成摘要\n\n'+intro+'\n\n'+tab(ss)+'\n\n- 四方法共用80,000个niche训练中心；MPBC2不参加niche拟合。\n- 完整输入重复拆分为0，但图方法细碎分群仍不推荐直接命名生物niche。\n- MENDER200选定结果近乎单群；MENDER50存在环支持量敏感性。\n- CellCharter50/200最终K=5/9；去除发现特征后ARI=0.543/0.575，只作为探索性分子证据。\n- 弱Harmony降低来源关联但改变原始分区，且达到迭代上限；保留NONE主线。\n- 五片原对象挂接与counts/QC验证通过。未运行MP10关联或通讯。\n\n完整覆盖、预测/投影/未定见tables/coverage_wide.tsv；方法解释和局限见REPORT.md。\n')
files=[p for p in R.iterdir() if p.is_file() and p.suffix in ['.md','.json','.yaml','.tsv','.sh']]
files += [p for p in (R/'code').iterdir() if p.suffix in ['.py','.R','.sh']]
files += [p for p in (R/'tables').iterdir() if p.suffix in ['.tsv','.json'] and p.stat().st_size<20e6]
key_png=[f'{s}_blca_four_methods.png' for s in SECTIONS]+['D1_Johnson_hpca_four_methods.png','native_duplicate_diagnostic.png','Harmony_Top20_partition_sensitivity.png','technical_association.png','candidate_stability.png','coverage.png','D1_Johnson_ROI_conflict_focus.png','Cervilla_ROI_far_UMI_hotspot.png','MPBC2_ROI_high_information.png','D1_Johnson_raw_evidence_conflict_focus.png','D1_Johnson_raw_evidence_HE_conflict_focus.png','Cervilla_raw_evidence_far_UMI_hotspot.png','UC_SCC_raw_evidence_typical_information.png','MPBC2_raw_evidence_high_information.png','CellCharter_both_r50_scVI_raw_neighbor_markers.png','CellCharter_both_r200_scVI_raw_neighbor_markers.png','hellinger_kmeans_blca_r50_NONE_blca_composition.png','MENDER_adapted_blca_r50_ring_counts_radial_composition.png']
files += [R/'figures'/n for n in key_png if (R/'figures'/n).exists()]
files=list(dict.fromkeys(files))
pd.DataFrame([dict(path=str(p.relative_to(R)),bytes=p.stat().st_size) for p in files]).to_csv(R/'tables/review_contents.tsv',sep='\t',index=False)
files.append(R/'tables/review_contents.tsv')
files=list(dict.fromkeys(files))
def pack():
 with zipfile.ZipFile(R/'review.zip','w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for p in files:z.write(p,p.relative_to(R))
pack()
size=sum(p.stat().st_blocks*512 for p in R.rglob('*') if p.is_file() and not p.is_symlink())
summary['allocated_GiB']=round(size/2**30,3);summary['review_MiB']=round((R/'review.zip').stat().st_size/2**20,1)
(R/'completion_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False))
with (R/'REPORT.md').open('a') as f:f.write(f'\n新增占用（含隔离CUDA依赖）约 {summary["allocated_GiB"]} GiB；review.zip 约 {summary["review_MiB"]} MiB。\n')
pack()
if complete:
 runtime=R/'logs/runtime.json';old=json.loads(runtime.read_text()) if runtime.exists() else {}
 if old.get('state')!='complete':(R/'logs/runtime_before_delivery.json').write_text(json.dumps(old,indent=2))
 runtime.write_text(json.dumps(dict(state='complete',time=summary['updated'],validation='five original objects reloaded; Pandas loader boolean compatibility fixed',model_branches=25,previous_execution_state='logs/runtime_before_delivery.json'),indent=2))
print('REPORT_PACKAGED',complete,summary['review_MiB'],flush=True)
