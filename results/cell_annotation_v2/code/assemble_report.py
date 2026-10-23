from common import *
import pandas as pd,numpy as np,anndata as ad
names=['Cervilla','D1_Johnson'];summaries={};ann={};pred={};obs={}
for n in names:
 r=ROOT/n;summaries[n]=json.loads((r/'metadata/reconstruction_qc_summary.json').read_text());ann[n]=json.loads((r/'metadata/annotation_summary.json').read_text());pred[n]=json.loads((r/'metadata/prediction_summary.json').read_text());obs[n]=pd.read_parquet(r/'cell_labels.parquet')
def table(df):
 def fmt(v):
  return (f'{v:.6g}' if isinstance(v,(float,np.floating)) else str(v)).replace('|','\\|').replace('\n',' ')
 rows=[list(df.columns),['---']*len(df.columns)]+[list(x) for x in df.itertuples(index=False,name=None)]
 return '\n'.join('| '+' | '.join(fmt(v) for v in row)+' |' for row in rows)
def link(path,label):return f'[{label}](<{path}>)'
rows=[]
for n in names:
 s=summaries[n];o=obs[n];rows.append(dict(section=n,图像核=s['segmentation_nuclei'],有非零表达的重建细胞=s['reconstructed'],QC保留=s['retained'],低信息=int(o.low_information.sum()),低信息占比=f'{o.low_information.mean():.2%}',原始counts分配率=f"{s['assigned_counts']/s['source_counts']:.2%}",训练与UMAP核心=int(o.scvi_training_core.sum()),未映射但保留=int((~o.umap_mapped).sum())))
count_table=table(pd.DataFrame(rows));pd.DataFrame(rows).to_csv(ROOT/'metadata/two_section_summary.tsv',sep='\t',index=False)
maj=[];status=[];fine=[]
for n in names:
 o=obs[n]
 for label,k in o.major_label.value_counts().items():maj.append(dict(section=n,major=label,n=int(k),fraction=float(k/len(o))))
 for label,k in o.annotation_status.value_counts().items():status.append(dict(section=n,status=label,n=int(k),fraction=float(k/len(o))))
 for label,k in o.loc[o.annotation_status.eq('subtype_supported'),'final_label'].astype(str).value_counts().items():fine.append(dict(section=n,subtype=label,n=int(k)))
major=pd.DataFrame(maj);states=pd.DataFrame(status);sub=pd.DataFrame(fine);major.to_csv(ROOT/'metadata/major_composition.tsv',sep='\t',index=False);states.to_csv(ROOT/'metadata/uncertainty_composition.tsv',sep='\t',index=False);sub.to_csv(ROOT/'metadata/subtype_support.tsv',sep='\t',index=False)
majorwide=major.pivot(index='major',columns='section',values='n').fillna(0).astype(int).reset_index();statewide=states.pivot(index='status',columns='section',values='n').fillna(0).astype(int).reset_index();finewide=sub.pivot(index='subtype',columns='section',values='n').fillna(0).astype(int).reset_index() if len(sub) else pd.DataFrame({'result':['No fine subtypes supported']})
scvi=json.loads((ROOT/'metadata/scvi_config.json').read_text());comp=pd.read_csv(ROOT/'metadata/Cervilla_v1_v2_comparison.tsv',sep='\t');validation=json.loads((ROOT/'metadata/delivery_validation.json').read_text());visual=json.loads((ROOT/'metadata/final_visual_review.json').read_text())
coverage=[]
for n in names:
 s=summaries[n];coverage.append(dict(section=n,输入bins=s['source_bins'],图像未覆盖bins=s['image_uncovered_bins'],分配bins=s['assigned_bins'],未分配bins=s['unassigned_bins'],输入counts=s['source_counts'],分配counts=s['assigned_counts'],未分配counts=s['unassigned_counts'],QC保留counts=s['retained_counts']))
report=f'''# 两张 Visium HD 细胞重建与亚群注释 v2 审阅版

2026-10-23。已分别完成 Cervilla、D1_Johnson 的新图像核分割、真实2µm原始计数聚合、宽松QC、冻结参考候选预测与本轮深度校准，再完成两张共享scVI、表达聚类和marker审查。旧v1及原始输入保留。本目录是可修改的审阅版，没有新设冻结文件；本轮计算止于注释。

## 实际细胞数量与处理范围

{count_table}

“重建细胞”仅统计分配后counts>0的对象。另有Cervilla {summaries['Cervilla']['zero_expression_nuclei_with_bins']}、D1 {summaries['D1_Johnson']['zero_expression_nuclei_with_bins']}个核虽被分配到提供的bins，但合计counts=0，未称为有表达细胞；它们在排除元数据中。未分配任何bin的核只保存几何和原因。QC排除对象分别为 {summaries['Cervilla']['excluded']}、{summaries['D1_Johnson']['excluded']}，包括上述零表达对象。

{table(pd.DataFrame(coverage))}

Cervilla使用恢复的11418×8964同片RGB，原生0.8636122748 µm/px；D1使用27318×26836捕获区RGB，原生0.2478437293 µm/px，原fullres裁剪 `[23764,35851,51082,62687)`。D1局部坐标仅减一次裁剪原点。模型mpp约0.5，重采样不会创造原图没有的细节，也不扩大表达测量视野。模型像素中心变换为 `(fullres−cropOrigin+0.5)×scale−0.5`；物理坐标为原fullres×原生mpp。四套坐标分别保存于 `spatial`、`spatial_image`、`spatial_model`、`spatial_um`；详见每张 `metadata/transform.json`。没有counts的图像区域不能解释为零表达组织。

## 分割、扩张与QC

StarDist `2D_versatile_he`，prob=0.05、nms=0.3；先查看致密、间质/混合与组织边缘ROI，再分块2048推断。D1初选edge实际仍位于组织内部，目视检查后改为跨组织边界的窗口。各张 `pilot_review.json`、`segmentation_pilot.png` 保留实际审查。Bin2cell最大扩张2个2µm array bins（4µm），种子123、PCA30用于冲突处理；没有GEX-only补细胞。H&E核及扩张只是细胞近似，并非实测细胞膜。

正式X来自原始Space Ranger整数counts逐bin求和。destripe只影响临时冲突处理视图，汇总前恢复原始矩阵；没有将归一化值取整为counts。完整18,085个输入features均保留，HVG只用于模型视图。逐基因分配计数守恒；未分配counts保持未分配。

保留条件为n_counts≥10且n_genes≥10；10–30 genes单独标记但保留。未套200/500 genes或统一线粒体删除阈值。5MAD counts/genes/面积、触图像边缘等仅标记并影响表达核心选择，不单项删除。面积与counts未作为自动删除条件；实际代表图中确认的4个D1非单核血样聚集区假分割被排除（均原本不在scVI核心），另8个大而不规则/多叶轮廓保留但标为形态歧义。Cervilla其中1个在审阅前已进入表达核心，历史training标记保留，其最终身份不作为确定标签。代表形态审查：{visual['morphology_summary']}

## 参考、证据与建议标签

继续使用 `cell2location_v2_frozen_20260916__scRNA32`，字段 `reference_celltype_v2` 的32类词表；参考文件名含30type不改变实际词表。三套已训练CellTypist模型复用，没有重训全参考。两张输入features顺序与映射已验证：14,329个映射features合并为14,327参考基因。按这些实测counts做10k归一化、log1p与训练时相同scaler及截断；没有把scVI latent或估计表达传给分类器。

冻结留出参考4745细胞复用于深度校准，按每张新细胞共享基因的真实中位counts稀释：Cervilla普通/低信息={pred['Cervilla']['targets_actual_shared_counts']}；D1={pred['D1_Johnson']['targets_actual_shared_counts']}。类别至少20个接受预测、经验精度≥0.8才获得可用margin门槛，结合full-depth是否可校准；这不是独立HD真值准确率。保存候选、sigmoid模型支持分数、margin及校准接受状态；支持分数不等同于经HD真值校准的身份概率。

全局表达簇检查主要谱系；在预测谱系内重建表达图和分群，再按切片分别审阅实测marker频率。模型候选仍逐细胞保留。细胞自己的marker组合或相对一致表达群的成组marker均可提供证据，不再以每个细胞检出两个低表达特异marker作为统一否决。群体证据只支持模型自身候选一致的细胞，不将标签传播给簇中所有细胞；详见 `metadata/annotation_rules.json` 与每张 `cluster_annotation_review.tsv`、`cluster_marker_evidence.tsv`。参考和cluster证据不是独立真值。

主要谱系计数（包括保留的mixed对象原有主要谱系建议；mixed身份应结合状态排除于确定类型统计）：

{table(majorwide)}

最终注释状态，分母为该张所有QC保留细胞：

{table(statewide)}

获得亚群级支持的建议标签：

{table(finewide)}

`predicted_subtype`为模型候选，`final_label`为建议标签，`annotation_status/evidence`保留支持路径及不确定原因。主要谱系明确但亚群不足时保留对应父类unresolved。Epithelial、Plasma、Mast在本参考粒度下保留broad identity；上皮不推断恶性状态。GC B需B谱系及成组GC证据，不能仅凭BCL6/cycling。`low_information`与身份不确定分开；混合嫌疑来自实际跨谱系marker共存，不等同于已证实doublet。

实际重点谱系/亚群审阅结论：{visual['annotation_summary']}

低置信结果按表达簇与谱系集中审阅，16行谱系概览见 `metadata/uncertain_lineage_review.tsv`，每张一行一个表达群的概览见 `metadata/uncertain_group_review.tsv`；详细证据见每张 `cluster_annotation_review.tsv`、`lineage_cluster_top_markers.tsv` 与marker dotplot；没有创建逐个Unknown人工任务。描述性Welch排序只用于浏览实测log表达差异，不是独立确认性检验。

## 两张共同表达空间

scVI实际输入 **2 sections**：{scvi['structure_core_by_section']}，合计{scvi['total_training_cells']:,}核心细胞、{scvi['training_genes']}个HVG；section_id为batch，NB，latent20、hidden128、2层，seed123，batch1024，最多50epochs、KL预热10epochs、early stopping patience10，实际{scvi['trained_epochs']}epochs。GPU为{scvi['device']}，PyTorch {scvi['torch']} / CUDA {scvi['cuda']}。训练约{scvi['training_seconds']/60:.1f}分钟。HVG按两张各自变异排名等权选择并优先共同HVG；每张保留全基因主对象。

有HVG计数的保留细胞保存latent；零HVG计数的Cervilla {int((~obs['Cervilla'].scvi_latent_mapped).sum()):,} / D1 {int((~obs['D1_Johnson'].scvi_latent_mapped).sum()):,}个细胞latent明确NA，均保留在全基因主对象及空间分母。初始核心中676个零HVG对象在输入复核后转为仅预测，初次短训练被替代；最终模型按修正核心训练。表达近邻、Leiden、UMAP仅用核心，非核心UMAP明确NA但空间坐标和总数分母保留。另保存未整合log-count PCA基线。最终模型验证ELBO从112.328降至93.319（最小值约93.312），最后5轮范围约0.030，KL权重为1；达到50轮上限，未再次因未成熟的预热早停。这里的表达近邻不是物理邻域，未建立跨切片物理邻接。每来源仅一张，section、队列与真实组织差异混杂，不能把充分混合作为验收标准。没有为UMAP外观反复训练。一次默认400epochs的KL预热与50epochs上限不匹配，导致早停时KL权重仅0.025；查看实际ELBO与重建损失后，将预热设为10epochs重新训练，相关轨迹与原因保留于scvi_warmup_correction.json和日志。

## Cervilla v1/v2共同覆盖比较

{table(comp[['version','reconstructed_full','retained_full','low_information_full','common_reconstructed_fully_covered','common_retained_fully_covered','common_assigned_fraction','common_median_counts','common_median_genes']])}

共同覆盖按同一原始bin表两版image_covered交集定义；counts只统计交集bins，细胞比较保留所有已分配bins均在交集的对象，不按cell_id或行号配对。v1旧报告重建数106,934含212个零表达对象；本表统一按非零表达计算106,722。完整组成/状态对照见 `metadata/Cervilla_v1_v2_composition.tsv`。D1旧结果为NA，不记作0。

本轮图像、模型图像归一化视图和注释证据使用方式均有变化，不能把全部标签变化归因于图像恢复。分割/计数变化与注释决策分别列出；不以细胞数增加、Unknown下降或genes中位数提高单项判断成功。相同物理ROI的分割对照见 `Cervilla_v1_v2_segmentation.png`。实际比较审查：{visual['comparison_summary']}

## 交付与验证

主对象：
- {link(ROOT/'Cervilla/cells.h5ad','Cervilla/cells.h5ad')}
- {link(ROOT/'D1_Johnson/cells.h5ad','D1_Johnson/cells.h5ad')}

每张同时保存 `cell_labels.parquet`、`excluded_cell_metadata.parquet`、压缩核mask、2µm array cell mask、bin_to_cell映射、未获bins核几何及实际变换。mask与bin-to-cell保存重建时全部数值标签（包含QC排除对象），通过source_cell_id连接最终保留表，不把所有mask标签都当作保留细胞。没有将原图或重复全基因矩阵嵌入对象。最终共享模型为 `models/scvi_two_sections/`；参考模型只读复用v1。GPU同版本软件覆盖层复用原环境其余依赖，没有修改v1环境。

每张 `figures/` 包含原图/核/扩张ROI、QC、major/subtype/不确定状态UMAP、全捕获空间图、重点谱系图、marker dotplot和代表大核。跨切片同一配色、透明背景；绘图分母见标题。共享表达基线和新旧比较位于 `metadata/`。

验证通过{validation['checks_passed']}项，重点为原始整数X、完整features、行counts/genes、映射/标签对齐、原始/局部/模型/µm坐标及低信息分母；每张随机25个细胞加最高counts细胞独立从原2µm H5重汇总比对全部基因，非重复全量复算。结果见 `metadata/delivery_validation.tsv`。图像与marker是计算审查，不替代病理专家或独立细胞真值验证。

## 补充数据支线

复用既有65行候选审计，更新于 `metadata/candidate_datasets_v2.tsv`。本轮没有新增就绪HD切片、没有新增平台全流程：
- GSE326226另外四张Xenium5K仅各读131,072字节H5 features，实际均5001个gene targets、固定Top20覆盖9/20；与此前A2R2一起五张现均有实际面板核对。只登记覆盖，不替补、不计算MP10分数。相同研究/跨平台不重复计作独立患者。
- UNC LCCC_553_Xenium读取嵌套ZIP目录与少量Zarr元数据（cells约0.78MB、transcripts约3.17MB）：作者细胞数258,581，cell/nucleus centroid、area及边界schema、µm单位已确认；转录本gene/quality/location字段已确认。未找到独立原始cell count矩阵，转录本到细胞分配尚未建立，不能称为已就绪细胞counts。未下载9.1GB整包。
- MPBC1/2仍需作者细胞原始counts+空间映射，或同片图像/positions/scales；UC_SCC仍缺现成2µm及形态输入，未下载22.69GB FASTQ。STT0000103与Oh/Jung仍属受限访问。其他已核候选沿用前轮证据，没有把未复查项写成新近复验。

人工操作集中于 {link(ROOT/'metadata/MANUAL_INPUT_ACCESS.md','MANUAL_INPUT_ACCESS.md')}，列明完整链接、所需文件、阻碍、操作、落盘位置与可解锁切片。未联系作者/发邮件/索取凭据，未让这些缺项阻塞本次两张。

## 后续是否具备物理邻域分析基础

已有两张独立cell_id、真实原始counts、经检查物理坐标、分割映射、建议谱系/亚群和明确不确定性字段，具备从本审阅版设计后续分析的技术基础。血液/碎片区和广泛基质背景仍需针对性审查，不能把本版等同于已证明准确的单细胞真值。应按切片分别构建物理图，预先约定低信息、未定类及mixed的处理，对分割/扩张和标签不确定性做敏感性检查；稀少或不支持的亚群不可强行解释，两个来源各一张不足以分离患者/批次/组织效应。

本轮没有运行MP10/胆固醇评分、neighborhood/niche、spottedpy、空间阳性检验、通讯或CNV。交付后停止，不自动启动下游。
'''
candidates=pd.read_csv(ROOT/'metadata/candidate_datasets_v2.tsv',sep='\t').fillna('')
for n,cohort in [('Cervilla','Cervilla'),('D1_Johnson','GSE300692')]:
 hit=candidates.cohort.eq(cohort);candidates.loc[hit,'status']='Two-section v2 reconstruction and annotation review complete';candidates.loc[hit,'cell_outputs']=str(ROOT/n/'cells.h5ad');candidates.loc[hit,'next_step']='Review suggested identities; downstream not started';candidates.loc[hit,'v2_incremental_update']='New original image segmentation, raw2um aggregation, permissive QC, reference+cluster annotation and shared2section scVI completed'
candidates.to_csv(ROOT/'metadata/candidate_datasets_v2.tsv',sep='\t',index=False)
(ROOT/'ANNOTATION_REVIEW_REPORT.md').write_text(report);print(count_table,flush=True)
