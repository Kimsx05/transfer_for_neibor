from common import *
import numpy as np,pandas as pd,anndata as ad,importlib.metadata
from datetime import datetime

o=pd.read_parquet(ROOT/'cell_annotations.parquet');review=pd.read_csv(ROOT/'ALL_CLUSTER_REVIEW.tsv',sep='\t');methods=json.loads((ROOT/'metadata/lineage_structure_methods.json').read_text());mapping=json.loads((ROOT/'metadata/reference.json').read_text())['major_mapping'];checks=[]
def check(name,value,details=''):
 checks.append(dict(check=name,passed=bool(value),details=details));assert value,(name,details)
check('cell_id_unique',o.cell_id.is_unique and o.index.is_unique);check('row_index_matches_cell_id',np.array_equal(o.index,o.cell_id));check('all_retained_cells',len(o)==308549)
for inp in json.loads((ROOT/'metadata/inputs.json').read_text()):
 sec=inp['section'];source=Path(inp['source']);st=source.stat();check(sec+'_source_unchanged',st.st_size==inp['source_stat']['size'] and st.st_mtime_ns==inp['source_stat']['mtime_ns']);a=ad.read_h5ad(source,backed='r');z=o[o.section_id.eq(sec)].set_index('cell_id').reindex(a.obs_names);check(sec+'_coverage',len(z)==a.n_obs==inp['cells'] and z.proposed_major.notna().all());check(sec+'_coordinates_unchanged',np.allclose(z[['x_um','y_um']].to_numpy(),a.obsm['spatial_um'],atol=1e-6,rtol=0));check(sec+'_low_info_preserved',np.array_equal(z.low_information.to_numpy(),a.obs.low_information.to_numpy()));check(sec+'_core_preserved',np.array_equal(z.scvi_training_core.to_numpy(),a.obs.scvi_training_core.to_numpy()));a.file.close()
check('hierarchy_unknown',not ((o.proposed_major=='Unknown')&(o.proposed_subtype!='Unknown')).any());fine=o[o.proposed_subtype.ne('Unknown')];check('hierarchy_fine_parents',np.array_equal(fine.proposed_subtype.map(mapping),fine.proposed_major));check('unknown_has_specific_reason',o.loc[o.proposed_subtype.eq('Unknown'),'unknown_reason'].ne('not_applicable').all());check('review_is_not_freeze',o.review_status.eq('computational_proposal_pending_user_review').all())
for head in ['major','subtype']:check('SingleR_'+head+'_all_cells_actually_ran',o['SingleR_'+head+'_run_status'].eq('ran_cell_level').all() and o['SingleR_'+head+'_candidate'].notna().all())
check('every_global_group_reviewed',set(o.global_cluster_v21)==set(review.loc[review.analysis_level.eq('global'),'cluster'].astype(str)));check('every_lineage_group_reviewed',set(o.loc[o.proposed_major.ne('Unknown'),'lineage_cluster_v21'])==set(review.loc[review.analysis_level.eq('lineage'),'cluster']));check('all_actual_lineages_reviewed',set(o.proposed_major)-{'Unknown'}==set(methods));pd.DataFrame(checks).to_csv(ROOT/'metadata/delivery_checks.tsv',sep='\t',index=False)
compact=['section_id','cell_id','n_counts','n_genes','low_information','scvi_training_core','old_major_label','old_final_label','predicted_major','predicted_subtype','SingleR_major_candidate','SingleR_major_pruned_label','SingleR_major_delta_median','SingleR_subtype_candidate','SingleR_subtype_pruned_label','SingleR_subtype_delta_median','global_cluster_v21','lineage_cluster_v21','proposed_major','proposed_subtype','state','annotation_support','evidence_source','unknown_reason','quality_flag','artifact_review_flag','major_changed_from_v2','subtype_changed_from_v2','review_entry','review_status'];o[compact].to_csv(ROOT/'cell_annotations.compact.tsv.gz',sep='\t',index=False,compression='gzip')
def table(df):
 df=df.reset_index();return '| '+' | '.join(map(str,df.columns))+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(str(x).replace('|','/') for x in row)+' |' for row in df.itertuples(index=False,name=None))
def ct(col):return pd.crosstab(o[col],o.section_id)
levels=[]
for l in methods:
 for sec in ['Cervilla','D1_Johnson']:
  z=o[o.proposed_major.eq(l)&o.section_id.eq(sec)];labs=z.loc[z.proposed_subtype.ne('Unknown'),'proposed_subtype'].value_counts();states={c[6:]:int(z[c].sum()) for c in z if c.startswith('state_') and z[c].dtype==bool and z[c].any()};levels.append(dict(lineage=l,section=sec,n=len(z),fine_supported=int(z.proposed_subtype.ne('Unknown').sum()),subtype_unknown=int(z.proposed_subtype.eq('Unknown').sum()),subtypes='; '.join(f'{k} ({v})' for k,v in labs.items()) or 'none',states='; '.join(f'{k} ({v})' for k,v in states.items()) or 'none'))
level=pd.DataFrame(levels);level.to_csv(ROOT/'metadata/all_lineage_resolution.tsv',sep='\t',index=False)
changes=[]
for sec in ['Cervilla','D1_Johnson']:
 z=o[o.section_id.eq(sec)];changes.append(dict(section=sec,n=len(z),old_major_unknown=int(z.old_major_label.eq('Unknown').sum()),new_major_unknown=int(z.proposed_major.eq('Unknown').sum()),old_display_unknown=int(z.old_final_label.eq('Unknown').sum()),old_fine_supported=int(z.annotation_status.astype(str).eq('subtype_supported').sum()),new_fine_supported=int(z.proposed_subtype.ne('Unknown').sum()),major_changed=int(z.major_changed_from_v2.sum()),new_mixed_flag=int(z.mixed_observed_programs.sum()),artifact_review=int(z.artifact_review_flag.sum())))
reg=[]
for sec in ['Cervilla','D1_Johnson']:
 z=o[o.section_id.eq(sec)];reg.append(dict(section=sec,old_predicted_Myeloid_to_T=int(((z.predicted_major=='Myeloid')&(z.proposed_major=='T/NK')).sum()),old_predicted_Myeloid_to_B=int(((z.predicted_major=='Myeloid')&(z.proposed_major=='B')).sum()),old_predicted_Fibroblast_to_Mural=int(((z.predicted_major=='Fibroblast')&(z.proposed_major=='Mural')).sum()),C1QC=int(z.proposed_subtype.eq('C1QC TAM').sum()),SPP1=int(z.proposed_subtype.eq('SPP1 TAM').sum()),C1Q_and_SPP1_measured_programs=int(z.C1Q_SPP1_coexpression.sum()),fine_supported_old_calibration_failed=int((z.proposed_subtype.ne('Unknown')&~z.subtype_calibrated).sum())))
reg=pd.DataFrame(reg);reg.to_csv(ROOT/'metadata/regression_evidence_checks.tsv',sep='\t',index=False)
oldnew=o.groupby(['section_id','old_major_label','proposed_major'],observed=True).size().rename('n').reset_index();oldnew.to_csv(ROOT/'metadata/major_transition_v2_to_v21.tsv',sep='\t',index=False)
priority=review[(review.scope=='joint')&(review.analysis_level=='lineage')].sort_values(['review_priority','secondary_lineage_signal_fraction'],ascending=[True,False]).groupby('lineage',sort=False).head(2)[['lineage','cluster','n','suggested_dominant_label','review_question']]
versions={k:importlib.metadata.version(k) for k in ['anndata','scanpy','numpy','pandas','scipy','scikit-learn','pynndescent','leidenalg']};js(ROOT/'metadata/python_method_versions.json',versions)
text=f'''# 两张 Visium HD 全谱系细胞注释 v2.1 审阅材料

生成时间：{datetime.now().isoformat(timespec='seconds')}（服务器本地时间）。本轮计算止于细胞注释。所有标签均为 **computational_proposal_pending_user_review**，不是人工确认或最终冻结。

## 交付与范围

- `cell_annotations.parquet`：全体 308,549 个 v2 保留细胞，一行一个 ID，含原 QC/候选/旧校准、SingleR 大类和细类完整精简诊断、新建议、状态、具体 Unknown 原因与审阅索引。
- `cell_annotations.compact.tsv.gz`：便于一般表格工具读取的精简副本。
- `ALL_CLUSTER_REVIEW.tsv`：所有全局簇和全部谱系局部簇的联合及逐切片证据/命名审阅总表。包含多数标签同时也包含完整组成，不应将多数标签赋给所有成员。
- `metadata/*_markers.parquet`：实测 marker 表，含平均 log1p10k、检出率、与同分析集合其余细胞的差、描述性 z/FDR；`nearest_cluster_marker_comparisons.parquet` 提供局部 PCA 最近竞争簇的表达和检出率差。
- `figures/`：两切片共享标签配色、透明背景；所有大类与各谱系已生成 known/data-driven marker 图及空间/表达图。全局 UMAP 只画既有 core；空间图包含全部保留细胞。局部图标明 informative 数量。
- `code/attach_annotations.py`：按 cell_id 对齐，默认只验证；指定新的输出路径才生成附加 v21 字段的对象，拒绝覆盖原文件。

原始 counts、v1/v2 对象、分割和坐标均未更改。本轮没有重新训练 scVI、分割、恢复图像、下载新的生物数据，也未做 MP10/胆固醇、空间邻域、通讯、CNV/恶性推断。

## 各切片大类建议

{table(ct('proposed_major'))}

## 全谱系分辨率和状态

{table(level.set_index(['lineage','section']))}

Epithelial、Plasma、Mast 也重新检查了谱系内表达结构及状态；参考对这些群体本来只有 broad 标签。没有将上皮 basal/luminal/squamous、浆细胞 IgA/IgG 或 mast 活化程序自动升格为新身份。保留大类而细类 Unknown 是可用的父类结论，不代表未分析。D1 的2个 Erythroid 建议属于参考之外的实测血红蛋白程序对象，必须结合 artifact flag 审查，不据此宣称发现成熟红细胞核或新群体。小群体不强制聚类。

## 实际表达发现与需要区分的解释

- 上皮：Cervilla 的 luminal 程序较多，D1 的 basal/squamous 程序较多；程序可以重叠。Epithelial::0 富集 KRT5/KRT14/KRT6A，::5 有 OLFM4/HPGD/MUC1/PIGR，::2 主要由低深度/MT 差异区分。保留 Epithelial 身份并另列状态，不推断恶性或机制。
- Fibroblast/CAF：Matrix_CAF 有 COL11A1、MMP11、POSTN、CTHRC1/FAP 等支持，而不是仅用 COL1A1。局部 ::1/::3/::5 的数据驱动差异混有 KRT/上皮程序，::6 带 VWF/RGS5/血管信号，不能把这些局部簇当成纯 CAF 新类型。逐细胞竞争转录本另列 secondary_lineage_signal。
- 髓系：Myeloid::5 是 D1 的43个细胞，CHIT1/CHI3L1/APOE/SPP1/GPNMB 同向支持 SPP1 相关巨噬程序；更多 C1Q 细胞分布于其他簇。::2/::3 有 S100A8/A9 和上皮差异，S100A8/A9 本身不决定 Classical Mono，需 FCN1/CCR2/S100A12 等更具体支持；炎症程序另存 state。C1Q/SPP1 并存保留，不把缺乏某个 marker 当互斥证据。
- T/NK：::0 有 CXCL13、SELL/GIMAP7，::1 有 CD3D/E，但并非整簇均可定 Tfh 或同一 T 细类；::2 受 Ig 信号、::3 受炎症/上皮信号影响。多数保留 T/NK，具体 Treg/Tfh/naive-memory/NK 等建议仅对应有证据的成员。
- B：两个局部群主要是 BANK1 等 B 程序与 Ig 丰富背景差异，少部分支持 Naive/Memory；本轮没有足够证据提出 GC B。极少 Cervilla B 和个位数细类不构成稳定群体的证明。
- 内皮：局部两簇并未自然分成四种血管身份，背景差异明显；只对有区分 marker 的成员提出 Arterial/Venous/Capillary/Lymphatic EC，其余保留 Endothelial。
- Mural：RGS5/CSPG4/PDGFRB 与 MYH11/CNN1/MYOCD/LMOD1 用于区分 Pericyte/SMC；局部两簇仍有上皮/Ig 背景差异，不能凭簇号等同于这两个身份。原 Fibroblast 到 Mural 的纠正与剩余 CAF 均分别审阅。
- Plasma：两簇有 IgG/IgA 相关表达差异，但只命名 Plasma 并单列免疫球蛋白状态；游离 Ig 表达不能单独建立 Plasma 身份。
- Mast：两个小表达簇以一般表达差异为主，没有可靠细身份分离；保留 Mast。活化程序支持极少，不能由 marker 缺失推断其他细胞静息。
- 参考外红系信号：仅2个 D1 对象满足协同血红蛋白程序，均带 artifact_review_flag；这是待核查对象，不是已确认红系群体。没有支持独立 Neural_glial 群体；不能把稀疏信号缺失写成生物学不存在。

这些命名是可修订的计算建议。主要身份有支持但另有竞争转录本的细胞，可能受环境 RNA 或重建边界影响；本轮不能确定其来源。局部 HVG/PCA 揭示了这些差异，但表达聚类本身并不证明新的细胞身份。

## 与 v2 的变化

{table(pd.DataFrame(changes).set_index('section'))}

{table(reg.set_index('section'))}

完整转移见 `metadata/major_transition_v2_to_v21.tsv` 和 `label_transition_v2_to_v21.tsv`。old_predicted 的回归检查与 old_major_label（v2 已接受标签）是不同基线，勿混用。Mural 从旧 Fibroblast 父类拆出是层级修正的一部分。C1Q 与 SPP1 程序允许在同一对象共现，两个状态不互斥；提出的 TAM 标签不能被解释成互斥稳定机制。旧校准失败可被实际 marker/表达支持覆盖，但这并不消除不确定性。

## Unknown 与质量的具体原因

{table(ct('unknown_reason'))}

`major_unknown_reason` 对大类未定单列原因；大类清楚但细类不清保留父类。low_information、v2 mixed、当前竞争程序和 artifact flag 均独立保存，彼此可重叠，不可相加当互斥细胞数。当前 mixed 仅表示同对象可比的跨谱系表达，不等于已证明 doublet。血红蛋白背景、宽泛 ECM 和可能分割问题需要区分；没有为减少 Unknown 改变分母。低信息总数仍为 Cervilla 12,573、D1 30,342。

## 全簇审查与需集中判断的问题

共 {int((review.scope=='joint').sum())} 个联合 review 条目：全局 {int(((review.scope=='joint')&(review.analysis_level=='global')).sum())} 个（33 个表达簇 + 一个明确标出的 noncore 集合），谱系 {int(((review.scope=='joint')&(review.analysis_level=='lineage')).sum())} 个。low_information_unclustered 是审阅集合而非表达簇。所有条目均有逐切片支持；细胞很少或某切片没有成员时不制造空群。

{table(priority.set_index('cluster'))}

全局 1/11/12/14/18/21/25 等含不同谱系的成员，不能凭一个多数 marker 整簇覆盖。多个小全局簇的 top genes 主要为 MT 基因且深度低，已明确标出，不命名为新身份。`ALL_CLUSTER_REVIEW.tsv` 的 competing_lineage_markers 必须结合细胞级 mixed 标记读取：整簇出现不同谱系不等于每个成员混合。上表仅为集中查看入口，全部谱系和有把握簇也在总表中。

## 方法及解释边界

1. 使用原 v2 重构细胞 raw integer X。marker 使用全部输入基因，映射符号并合并重复后 18,082 个；参考不包含的实测基因仍可用于 marker 检查。实测检出率和 log1p(10,000/细胞总 counts) 的差为描述性证据，不以 scVI 去噪/插补表达发现 marker。
2. 原共享 scVI 涉及 **2 sections**，core 252,601（81,789 + 170,812）；本轮只复用其 latent20、全局簇和 core UMAP。大类首先来自身份 marker/竞争证据，允许跨旧标签纠正；弱成员可通过 marker 支持种子的表达近邻投射，具体阈值、距离范围和证据路径见 `major_rules.json`、`major_projection.json`。没有物理邻居标签传播。
3. 每个有足够细胞的纠正后谱系重新选 section-aware HVG（至多 2,000，去 MT/RPL/RPS，补充实测身份 marker），稀疏归一化和方差缩放，PCA30、表达 kNN15、Leiden0.65，seed123；小群体仅做 marker 审阅。未强制来源混合。仅每来源一张切片，来源/技术与生物差异无法完全分离。局部 marker 以同谱系为对照，并记录最近局部 PCA 竞争簇。
4. SingleR **确实运行**：R4.5.3、SingleR2.12.0，原冻结 `reference_celltype_v2` 32 类，原参考300,033细胞按标签及 dataset/specimen 轮转选至多400/类，12775参考细胞（vCAF全375）。specimen不假定独立患者。共享14,327去重基因、实测 counts→log1p10k；wilcox/de.n20，aggr.ref=TRUE，seed123，8线程；大类9类（SMC/Pericyte归Mural）和细类32类分别训练一次，所有 query 按5000分块逐细胞预测，按整切片统一prune。保存 candidate/pruned/selected score/delta/最初第一第二名；score/delta 是相似性诊断，不是正确概率。within-major initial 候选来自同一细类初始分数，不称为额外分层微调运行。安装缺失 scrapper1.4.0 及其依赖；服务器 R 启动路径优先级使最终依赖位于现有用户库，未重装项目环境，日志和 sessionInfo 保存。
5. 复用 CellTypist，不重训全参考。两工具同参考一致只作为算法间支持；任何单工具或旧校准均不是真值。直接 marker、特异 anchor、参考兼容和表达簇支持是不同路线，没有全细胞通用的双工具一致/旧校准通过/固定 marker 数/旧候选纯度 AND 门槛。实现仍有显式可审查的描述性阈值，未声称通过独立 HD 真值校准。
6. 亚群证据要求相近身份的区别；近似竞争时保留父类。state 是实测程序，不是机制或恶性结论。单 marker 未检出不是不存在证据。GC B 不预设存在；单独 BCL6/cycling 不足以定 GC。
7. 大量细胞不能替代生物学重复，z/FDR 仅用于排序和描述；nearest 比较主要看效应与检出率。分割核+局部扩张不是精确膜边界；本轮复用 v2 已有形态审阅，8 个原形态疑问仍留表，没有新人工确认。原 H&E/分割例图位于 `../cell_annotation_v2/{{Cervilla,D1_Johnson}}/figures/representative_large_nuclei.png`。

官方实现依据：[SingleR 单细胞参考](https://www.bioconductor.org/books/release/SingleRBook/sc-mode.html)、[参考聚合与复用训练](https://www.bioconductor.org/books/release/SingleRBook/advanced-options.html)、[诊断和剪枝](https://bioconductor.org/books/release/SingleRBook/annotation-diagnostics.html)。

## 复现和必要检查

Python：`../cell_annotation_freeze_v1/.venv/bin/python`，版本见 `metadata/python_method_versions.json`；R版本/依赖见 `metadata/SingleR_session.txt`。依次运行 code 下 prepare_inputs.py → run_singleR.py；expression_evidence.py → major_proposals.py → lineage_structure.py → subtype_proposals.py → review_tables.py → make_review_figures.py → validate_and_report.py。SingleR 使用现有模型时固定训练不变；全部路径从脚本目录解析。

必要一致性检查 {len(checks)} 项通过，见 `metadata/delivery_checks.tsv`：ID 唯一与完整、按 ID 对齐、两张分母正确、low_info/core 保留、空间坐标未改、源文件大小/修改时间未变、标签层级、Unknown 原因、SingleR 实际全覆盖、全簇/全谱系审阅覆盖。没有额外全文件哈希或多份完整表达矩阵副本。

本轮停在集中审阅材料；建议标签尚未经用户确认，不进入下游分析。
'''
(ROOT/'ALL_LINEAGE_ANNOTATION_REVIEW.md').write_text(text);print('DELIVERY_CHECKS_PASSED',len(checks),flush=True);print(ct('proposed_major').to_string());print('Report written',flush=True)
