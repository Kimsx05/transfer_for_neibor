import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import pathlib,json,pandas as pd,numpy as np
P=pathlib.Path.cwd();B=P/'results/mp10_hypothesis_followup_v2';T=B/'tables'
def rd(n):return pd.read_csv(T/n,sep='\t')
def put(d,n):d.to_csv(T/n,sep='\t',index=False,na_rep='NA')
def doc(n,s):(B/n).write_text(s.strip()+'\n')
r=pd.concat([rd('score_model_results.tsv.gz'),rd('gsea_model_results.tsv.gz'),rd('recovered_SCENIC_model_results.tsv')],ignore_index=True)
e=rd('EC_M20_GSEA_registry.tsv').merge(rd('EC_M20_Hallmark_diagnostic.tsv.gz'),on=['source_model','pathway','mask'],how='left');e['FDR_family']=e['family'];e['P']=e['pval'];e['y']=e['pathway'];e['x']='Epi MP10 Q50';e['kind']='gsea';e['setting']='EC_M20_diagnostic';e['analysis_role']='v2_diagnostic';e['primary']=False;e['status']=np.where(e.P.notna(),'completed','not_evaluable');e['reason']=np.where(e.P.notna(),'','outside_GSEA_size_15_500');e['n_specimens']=12;e['n_recorded_patient_units']=12;e['n_datasets']=4;e['residual_df']=7;e['method']='original edgeR rank / fgseaMultilevel';e['statistical_unit']='specimen / recorded patient unit';e['min_cells']=50;e['patient']='main';e['testable_n']=e.P.notna().sum();e['source_DE']=str(P/'results/mp10_immune_endothelial_mechanisms_v1/tables/de/M0073.tsv.gz')
# Explicitly restore registered but too-small pathways as NE, never zero.
h=rd('hallmark_gene_sets.tsv.gz');de=pd.read_csv(e.source_DE.iloc[0],sep='\t');top=set((P/'neighborhood_mp10_calibration_v5/01_score_fidelity/top20.txt').read_text().split());universe=set(de.gene)-top
for i,z in e.iterrows():e.loc[i,'n_mapped']=len(set(h[h.signature==z.pathway].gene)&universe)
e['n_expected']=e.pathway.map(h.groupby('signature').gene.nunique());put(e,'EC_M20_Hallmark_complete.tsv.gz')
r=pd.concat([r,e],ignore_index=True);r['source_result_table']=np.where(r.model_id.str.startswith('V2E_'),'EC_M20_Hallmark_complete.tsv.gz',np.where(r.model_id.str.startswith('V2S_'),'recovered_SCENIC_model_results.tsv',np.where(r.kind=='gsea','gsea_model_results.tsv.gz','score_model_results.tsv.gz')));put(r,'model_results.tsv.gz')
reg=pd.concat([rd('model_registry.tsv'),e[['model_id','kind','setting','FDR_family','x','y','mask','min_cells','patient','primary','planned_n','statistical_unit','analysis_role']]],ignore_index=True);put(reg,'model_registry_complete.tsv');put(r.groupby(['FDR_family','planned_n'],dropna=False).agg(registered=('model_id','size'),testable=('P','count'),completed=('status',lambda x:(x=='completed').sum())).reset_index(),'FDR_family_summary.tsv')
assert len(r)==1328 and not r.model_id.duplicated().any() and set(reg.model_id)==set(r.model_id)
assert r.loc[r.status!='completed',['P','q']].isna().all().all()
sel=r[(r.state=='H05')&(r.setting=='main')];put(sel,'H05_primary_and_masked_results.tsv.gz')
notes={
'CHOL_SYNTH':'胆固醇总体合成；含共享甲羟戊酸上游和调控成员；转录程序不等于胆固醇通量。',
'BLOCH':'下游 Bloch 胆固醇合成；小集合 GSEA NE 与样本评分可检验分开。',
'FA_ACYL_COA':'脂肪酰 CoA 生物合成涵盖合成、延长、活化；不能简称纯新生合成。',
'FA_ELONG':'脂肪酸延长；不是总脂肪酸合成或释放。',
'FA_UNSAT':'不饱和脂肪酸生物合成；不是纯 SCD/FADS 去饱和酶面板。',
'LCFA_IMPORT':'长链脂肪酸进入细胞；含调控因子，不能由全成员正向表达推断净通量。',
'CHOL_IMPORT':'胆固醇摄入；接收端主评分覆盖不足，不能缩分母。',
'CHOL_EFFLUX':'胆固醇外排；包含非编码和正负调控成员，必须保留完整分母；不代表游离脂肪酸释放。',
'TG_SYNTH':'甘油三酯合成/储存相关；不代表外排。',
'TG_CATAB':'甘油三酯分解；不直接证明胞外游离脂肪酸释放。',
'FAO':'线粒体脂肪酸 beta 氧化，与 OXPHOS 不同。',
'PRENYLATION':'冻结版本实际为 GO:0097354 GOBP_PRENYLATION；不冒充 GO:0018342 下载集合，转录不等于修饰量。',
'LXR':'NR1H2/NR1H3 相关响应；不能确定具体氧固醇、来源或转移方向。',
'LCFA_TRANSPORT':'长链脂肪酸运输；不是跨内皮递送的特异量度。'}
s=rd('signature_registry.tsv');s['interpretation_boundary']=s.signature.map(notes);s['applicable_cells']='Epi / Treg / Capillary EC subject to fixed lineage coverage';put(s,'signature_interpretation.tsv')
dec=[('SREBF2 final regulon','有支持','原始六模块均被映射门槛排除；唯一符号修正恢复 top50 模块，36 targets；H05 Full/M20 均升高。','技术修正；单 ranking、同批细胞，非独立验证。','核内成熟 SREBF2 与靶位点结合/扰动'),('Epi mTORC1 beyond Top20','有支持','H05 repaired NES 3.217 / masked 2.744；masked leading-edge 55，其中43非正式脂质集合成员。','signature 不等于 mTORC1 磷酸化；状态和 RNA 含量仍可能参与。','p-S6/p-4EBP1 与 SREBF 成熟蛋白并行测量'),('Cholesterol relative transcriptional bias','有支持','CHOL_SYNTH−FA_ACYL_COA D=.679 Full、.808 M20；固定基因匹配支持正方向。','仅所指定标准化程序；Full Bloch−FA 未获支持，不能推广所有甾醇或绝对通量。','同位素胆固醇/脂肪酸合成通量对照'),('Epi cholesterol efflux upregulation','未获支持','H05 Full q=.423，M20 q=.215；正式集合未出现明确升高。','未获支持不等于不存在特定转运蛋白或条件性输出。','测外排而非只测合成'),('Epi output to Treg uptake','无法评价','两项主关联均因 Treg 摄取签名覆盖不足 NE；没有通用 FFA 释放读数。','真实 profile 交集18，不能补零或扩大背景救显著性。','稳定同位素供体上皮→接收 Treg 的质量守恒/转移测量'),('H05 lipid programs to Treg cycling','未获支持','两项可估计 quasi-binomial 主模型 q=1。','18记录患者单位、5数据集；区间仍宽。','同标本共培养下 cycling 与功能分开读出'),('H05 lipid programs to Treg OXPHOS/G2M','支持有限','主 HC3 singleton leverage=1 为 NE；预设17单位子集可估计但未获支持。','不能以子集替换主分析，也不能据转录分数推断实测呼吸。','直接呼吸测量及 cycling 匹配'),('Capillary EC mTORC1','支持有限','独立 EC RNA 的 Q50 DE/GSEA 正向，M20 NES1.742 q=.000499。','12单位4数据集，7来自 PRJNA662018；评分98/200 NE、HC3 singleton；未做环境 RNA 校正。','同患者内皮纯度核验与蛋白/代谢读数'),('EC nutrient delivery','无法评价','五项预登记 EC 路线因固定评分覆盖不足 NE。','运输与摄取不是跨内皮流向，不能声明营养递送。','跨内皮示踪与屏障方向性实验'),('Spatial mTORC1','支持有限','混合 spot 的距离负相关；组成/深度调整后跨数据集 q=.223/.163。','不是上皮专属 RNA；条件区间与数据集检验不同。','HD 同细胞类型/同切片邻域验证，pending')]
put(pd.DataFrame(dec,columns=['relationship','decision','evidence','limitation','next_discriminating_test']),'decision_table.tsv')
status=[('M0','completed','冻结 H05/L50 主分析、1206 初始模型及条件性/技术诊断追加登记；保留原冻结文件和修订。'),('M1_original_SCENIC','completed','8577 候选模块重建，12个 SREBF 模块逐项诊断；未重跑 GRNBoost2。'),('M1_symbol_correction_AUCell','completed','唯一历史别名修正；SREBF2 36 targets；116448 Epi 重评分；48条件性配对模型。'),('M1_alternative_ranking','not_evaluable','没有预先选定且配套的本地第二 ranking；未按救出 SREBF2 的结果选择新库。'),('M1_current_HGNC_download','technically_blocked','当前 HGNC 外网下载受网络限制；已用安装的 HGNC/NCBI 资源完成核心核验，不阻断主模块。'),('M2_sources_filter_QC','completed','三类来源、EC 2x2过滤、原始3905 EC计数和标记、设计与混合spot审计；旧结果只读。'),('M2_EC_score_HC3','not_evaluable','固定背景 mTORC1 98/200<50%；另有 singleton dataset h=1，分别报告。'),('M3_official_sets','completed','13正式过程+1探索运输集合，2026.1.Hs，完整成员/非编码/重叠/覆盖。'),('M4_H05_lipid_processes','completed','24配对、5数据集；Full/M20、H10/H50、min30、a/b、LODO、质量/组成、共同样本。'),('M4_small_set_GSEA','not_evaluable','小集合未达到15基因门槛；样本级程序可检验结果单独保留。'),('M4_direct_D_matching','completed','72个直接D模型和每mask200次固定种子基因匹配；非生物学重复。'),('M5_primary_Treg_uptake','not_evaluable','主两项：Treg CHOL_IMPORT4/12，LCFA_IMPORT7/19；均不满足固定覆盖。'),('M5_EC_exploratory','not_evaluable','五项主探索路线因sender或receiver固定集合覆盖不足；不是零效应。'),('M5_secondary_states','completed','适用的quasi模型及预登记去singleton子集完成；主连续HC3 NE另列，阳性/未支持/NE均保留。'),('delivery_validation','completed','完整模型表1328行；独立核验、原counts校验、18图62文件渲染；只读结果表重画入口。')]
pd.DataFrame(status,columns=['module','status','reason']).to_csv(B/'module_status.tsv',sep='\t',index=False)
# Coverage summary and exact findings for report.
put(r.groupby(['kind','status'],dropna=False).size().rename('n').reset_index(),'model_status_summary.tsv')
D=r[(r.kind=='direct_D')&(r.state=='H05')&(r.setting=='main')];put(D[['model_id','mask','x','beta','CI_low','CI_high','P','q','FDR_family','n_specimens','n_datasets']],'sterol_FA_primary_summary.tsv')
doc('README.md',f'''# MP10 hypothesis follow-up v2

实际结果目录：`{B}`。代码目录：`{P}/mp10_hypothesis_followup_v2/scripts`；结果目录内 `scripts/` 为同版本轻量交付副本。

本轮以冻结 **H05 vs L50** 为上皮主比较，H10/H50 为敏感性。M0 后连续完成所有可运行模块；旧结果未修改，未发布远程仓库。

## 主要发现

1. 原 SREBF2 六个候选模块均在 ranking 成员映射门槛被排除，并非没有 TF 表达、GRN 边或 motif 支持。PRAG1 的可核实历史符号 SGK223 修正使 top50 从40/51到41/51，通过原门槛，恢复36成员 regulon。H05 活性增加，去 Top20 后仍有支持。
2. 正式集合支持 H05 胆固醇合成；直接与脂肪酰 CoA 程序比较，D=0.679 (95% CI 0.442–0.917)，M20 D=0.808 (0.537–1.079)。这是固定尺度的转录变化，不是通量比；Full Bloch 比较不支持更大变化。
3. 上皮 mTORC1、EC 的 Q50 关联和空间混合 spot 结果来自不同表达单位。EC 固定评分覆盖98/200不足，空间调整后的跨数据集支持减弱。
4. 上皮胆固醇外排程序未获明确支持；两项上皮→Treg 摄取主关联及五项 EC 探索关联均因覆盖不足为 NE。可估计的 cycling 主关联未获支持，不能据此证明脂质转移。

## 阅读顺序

- `RESULTS.md`：逐问题结果、效应与局限。
- `figures/F01_SREBF_module_mapping.pdf`、`F01_recovered_SREBF_activity.pdf`：技术原因与修正后活性。
- `figures/F02_mTORC1_singlecell_sources.pdf`、`F02_EC_filter_difference.pdf`、`F02_mTORC1_spatial_source.pdf`：三类来源及过滤。
- `figures/F03_H05_formal_lipid_effects.pdf`、`F04_sterol_FA_specimen_deltas.pdf`：正式过程与直接比较。
- `figures/F05_prespecified_routes_status.pdf`：预登记路线的 NE；次要状态见 F05 其余两图。
- `tables/model_results.tsv.gz`：1328条完整登记结果；`model_registry_complete.tsv`、`FDR_family_summary.tsv`、`decision_table.tsv` 配套。
- `METHODS.md`、`module_status.tsv`、`FIGURE_INDEX.md`、`qc/`：实际方法、模块状态与核验。

## 复现

从项目根目录运行：

```bash
cd '{P}'
mp10_immune_endothelial_mechanisms_v1/venv/bin/python mp10_hypothesis_followup_v2/scripts/run.py --resume
mp10_immune_endothelial_mechanisms_v1/venv/bin/python mp10_hypothesis_followup_v2/scripts/run.py --replot
```

`--resume` 核对输入、代码、冻结配置及输出哈希；兼容完成结果复用。缺失产物可按冻结管线恢复，发生不兼容变化时停止覆盖并要求新版本。`--recompute` 在同一输入/配置下重算受管理阶段；不会重跑 GRNBoost2、CellChat 或 NicheNet。`--replot` 仅依赖结果表和绘图配置，不读 Seurat、不重新统计。无原始数据的接收者可直接执行 `Rscript --vanilla scripts/11_replot.R`，设置 `MP10_RESULTS_DIR` 指向本结果目录，并安装 `logs/sessionInfo_scores.txt` 对应绘图依赖。

所有运行限制单进程计算线程1，R使用 `--vanilla`；顺序复现不超过总4线程。输入绝对路径、大小、时间、SHA256见 `tables/input_manifest.tsv`；原大对象仅引用，不复制。交付体积见 `tables/resource_usage.tsv`。未运行湿实验、HD、环境RNA校正、全量GRN重建、第二套完整seed456 motif网络或新排名库搜索。
''')
doc('METHODS.md',f'''# METHODS — actual execution

## Design and units

冻结时间、版本和技术修订见 `config/analysis_plan.yaml`、`M0_FROZEN.json`、`technical_revision_*.json`。这是已看过发现性结果后的定向内部检验，不是前瞻性独立验证。H05主，H10/H50敏感性；全部沿用原 exact c2l reference membership、未舍入分组和 L50，不重新划组，不修改评分。主profile≥50，敏感性≥30；≥4记录患者单位、满秩、残差df≥3、有效暴露变化。主T4463209a、替换b敏感性；其余映射保留项目 provisional 状态，不能将记录单位当作已确认独立患者。

原始计数配对DE继续引用 edgeR TMM robust QL `~ specimen + state` 的全基因结果；不为本轮候选重做两基因FDR。新增样本程序模型对 High−Low 差值用 OLS HC3，双侧t检验和同一HC3残差df的95% CI。跨细胞模型为 `receiver ~ dataset + standardized H05 sender`；H05−L50 sender差值及其他设置为登记敏感性。连续模型 h≥1−1e−10 时 HC3 NE；cycling为quasi-binomial，使用各模型离散度和残差df t区间。NE保留NA，真实零不替代缺失。

## Scores, gene sets, masking

复用 v1 `objects/scoring_matrices.qs` 的真实全谱系 TMM-logCPM 背景与reference profiles；没有截取候选并集重归一化。GSVA 2.4.9 ssGSEA，alpha=.25，normalize=FALSE。每集合在main-a Epi_all≥50或接收端原参考集固定mean/SD；所有High、LODO及患者替换均沿用该尺度，不中心化每个delta。

正式13过程+1EC运输补充均来自冻结MSigDB2026.1.Hs/Reactome资源；完整成员和GO非编码项不静默移除。程序评分≥5可用成员且≥50%预期集合；GSEA minSize15/maxSize500，signed全基因DE排序，seed123，nPermSimple10000，eps1e−10，scoreType std。EC整套Hallmark因数值诊断统一提高至100000，保留10000结果和警告；49可检验+1小集合NE全部登记，BH计划50。极小P为计算下限，不解释精确尾概率。

Full为正式过程描述主结果；M20仅从下游集合及排序/评分宇宙移除Top20，不改变RNA原counts、library、TMM或分组。接收端Full主，M20仅明确敏感性。旧自编面板来源单列，不冒称正式集合。基因标识ACAT2未替换成SOAT2，NECTIN2/PVRL2不重复计数。

## Direct sterol comparison

同一H05/L50配对中，D=(zCHOL_H05−zCHOL_L50)−(zFA_H05−zFA_L50)。主CHOL_SYNTH对FA_ACYL_COA；Bloch补充。每个mask两项一个家族。相同评分算法和预设参考；两个集合各自参考SD而非通量单位。200次固定seed123匹配：可用基因数相同，在全Epi参考背景标准化的平均logCPM和检出比例上最近邻无放回匹配，随机匹配顺序，所有结果保留。它检验集合构成敏感性，不是200独立重复。

## SCENIC stage audit and correction

原pySCENIC0.12.1 seed123 adjacency及精确训练cache复用；8577模块重建而不再跑GRNBoost2。rho=.03、非dropout-mask、仅activating、min_genes20；importance .75/.90、top50targets、top5/10/50TF-per-target；ranking5000、AUC .05、NES≥3，原hg38/refseq-r80/mc9nr 10kb ranking及匹配v9 annotation，motifsimilarityFDR .001。逐模块ID含TF/context/targets哈希，原始分母保留。

实际软件 `missing_fraction >= .20` 排除，即正好80%也排除。所有motif先在同输入下补算诊断再逐门核验。唯一无歧义HGNC/NCBI历史别名与ranking链接可修正；保留未映射、歧义和原名，不删基因降低分母。安装HGNC20170418与NCBI2025-09-24参考来源和hash均登记；PRAG1→SGK223获得同一NCBI157285支持。没有为SREBF2挑新ranking或放宽阈值。

修正后regulon逆映射回表达原名，保留原TF-self成员语义。AUCell在原20216训练基因顺序、完整116448 Epi normalizedRNA上按2000cell chunks计算，seed123/aucThreshold .05/1worker。Full与M20；原SREBF1前3000cells逐数值校验maxabs0。活性差异是样本配对证据，不是ChIP、成熟蛋白或独立数据复现。

## mTOR sources and QC

Epi读取v1修复后的M0031完整Hallmark家族；EC M0073为EC自身rawRNA计数对Epi Q50的dataset调整关联，保留历史primary标签；空间B3为混合spotRNA UCell，Epi位置仅决定query，不能将RNA归属上皮。空间条件区间和dataset-level q独立列出。

EC过滤交叉诊断固定两组样本与两种filterByExpr设计，逐基因输出原检出/CPM/保留/leading-edge/Top20。对原3905 EC用原RNA counts核对内皮标记、EPCAM/keratin共表达、原QC字段及单细胞计数贡献。该标记不是doublet判定；没有执行ambient RNA校正。单个代谢基因表达不被定义为污染。

## FDR, reproducibility and scope

完整登记1328行（1278主登记+50EC Hallmark诊断），所有family计划数见 `FDR_family_summary.tsv`。每High/mask/setting程序家族13；直接D家族2；跨细胞primary2、EC探索5、Treg次要6；恢复SCENIC每setting/mask/state2。敏感性另成family，NE保持NA而BH使用计划n；不选最小P。旧RNA/TF/GSEA的原q及来源另存。

`qc/validation.tsv`独立Python重算HC3、quasi、BH和约束；`tables/counts_membership_validation.tsv`校验原整数计数及精确成员；`qc/render_audit.tsv`核对18图62文件，PDF矢量Arial、约12pt、300dpiPNG，带图例分离PDF。未运行新的HD、CNV、CellChat、NicheNet或全GRN，未改变原始或历史目录。

实际项目：`{P}`。版本见 `logs/sessionInfo_registry.txt`、`sessionInfo_scores.txt`、`environment_versions.tsv`；代码、参数、成员与参考哈希见checkpoint及各缓存metadata。原始大对象只读；工作峰值RSS来自恢复AUCell约44.61GiB，内部线程1，总并发计算不超过4。
''')
doc('RESULTS.md','''# 结果解释

## 1. SREBF2 为什么没有 regulon，修正了什么？

原 SREBF2 在 Epi 训练基因中，自身可映射、属于 TF predictor，也有 GRNBoost2 边。对8577个原始候选模块精确重建后，SREBF1/2各6个唯一模块；不能用日志警告条数代替模块数。SREBF2六模块均存在 NES≥3 且带匹配TF注释的motif，但原成员映射比例为67.5%–78.4%，全部在 `missing_fraction>=0.20` 门被跳过，因而没有进入最终regulon合并/AUCell。不是转换阶段莫名丢失，也不是已经计算出零活性。

最接近门槛的是top50：原始51成员中40可映射。PRAG1在ranking使用历史名称SGK223，唯一HGNC/NCBI映射校正后为41/51 (80.39%)，原NES/annotation/覆盖门槛均不变，保留14个motif并恢复36成员SREBF2 regulon。其余5模块仍不通过。SREBF1原top50为42/51并保留21成员，修正后仍21。双方其他成员的命名/参考覆盖原因有逐成员台账，不能称所有未映射成员都已解决；当前参考库覆盖局限仍在。

修正后的SREBF2含9个Top20，M20后27。H05 24配对、5数据集中，Full AUCell差0.1123 (95% CI 0.0758–0.1487)，q约3.3×10^-6；M20差0.02764 (0.01573–0.03955)，q=0.000153。SREBF1 Full差0.04647，M20差0.00482、q=0.0211。与原SREBF1逐细胞实现核对的Full/M20误差均0。

结论由“原参考下SCENIC不可评价”更新为“版本化命名修正后有regulon及配对活性支持”。旧NE仍保留，不能追溯改写旧结果，也不能称新计算为独立验证。单个历史符号使模块跨门槛说明参考敏感性；seed456只做过关键边复核，没有完整第二套motif网络。仍缺成熟核内SREBF2、结合或扰动证据。

## 2. 三类 mTORC1 证据及边界

**上皮自身RNA**：H05 M0031修复后 Full NES=3.217、q=1.26×10^-9；M20 NES=2.744、q=1.67×10^-9，极小P受eps下限限制。Full leading-edge62中10为Top20、12为其他正式脂质集合成员、40为其他响应；M20 leading-edge55中12脂质、43其他。支持超出Top20直接重叠。新增配对评分Full效应1.111 SD，M20 .669 SD；分别加入RNA质量或原cluster0组成后方向保留。它是状态相关转录程序，不是激酶磷酸化证据。

**毛细内皮自身RNA**：旧M0073使用Epi Q50连续暴露，不能叫H05比较。12记录患者单位/4数据集，CNP0000460=3、GSE326225=1、Inhouse=1、PRJNA662018=7；旧Full NES=1.864、q=1.58×10^-5。新增整个Hallmark家族M20诊断NES=1.742、q=.000499，说明Top20不是全部驱动。它仍是小样本且主要由一个数据集贡献；两个singleton导致HC3 h=1，不能换标准误强行得到P。

旧DE保留176个mTOR基因，新固定评分背景仅98/200。2×2复核：旧样本+实际设计176；旧样本+单组过滤102；新reference样本+单组过滤98；新样本+实际设计157。主要差别来自过滤设计，同时有样本集合贡献；不是名称丢失或把同一个背景错算。98/200未达50%，原评分NE保留，不以旧宽背景救活。

原3905个Capillary EC中3771至少检出一个内皮标记，105检出EPCAM，74 (1.90%) 同时EPCAM与keratin阳性。每个关键基因逐标本的原始counts、检出比例、最高1%/5%细胞及标记细胞贡献均输出。共表达提示需要核验，不等于确定污染；没有环境RNA校正，也不能用共享代谢基因本身诊断污染。

**空间混合表达**：B3在Epi富集位置查询的是混合spotRNA UCell，距局部MP10 P95参考点0–1000µm，6数据集53分析单位。Full/noMP10原r约−.276/−.250，dataset-search q=.0456/.0328；加入上皮比例和RNA深度后partial r约−.146/−.110，dataset q=.223/.163。条件区间仍可不跨零，但跨数据集检验支持减弱。不能把位置中心归属当成上皮专属RNA，也没有新启动HD。

## 3. 正式脂质过程支持什么？

13个正式过程均来自冻结2026.1.Hs，另有1个EC运输探索集合。保留旧面板名称和来源性质；正式PRENYLATION实际为GO:0097354。H05主分析24对/5数据集；H10为33对、H50为39对；≥30的H05为30对，属于嵌套敏感性。

H05正式评分显示胆固醇合成、Bloch、脂肪酰CoA、延长、不饱和脂肪酸合成、摄入和甘油三酯合成等多个程序上升，并非只有甾醇。胆固醇合成Full效应2.749参考SD、q=3.49×10^-14，M20效应1.989、q=4.52×10^-10；对应GSEA NES=3.789/3.469。脂肪酰CoA Full效应2.070、M20 1.181，均获支持。

对同24对标本直接检验D=Δcholesterol−Δfatty-acyl-CoA：Full D=.679 (95% CI .442–.917)，q≈9.9×10^-6；M20 D=.808 (.537–1.079)，q≈5.5×10^-6。H10/H50、≥30、a/b、共同样本和质量/组成设置均在表中；固定种子200次基因数/表达/检出匹配的D全部正向，Full中位约.666，M20约.719。这支持指定标准化转录程序相对偏胆固醇，且不完全由MP10定义重叠造成；匹配重抽不是200生物学重复。

限制同样明确：Full Bloch−FA D=−.113 (−.341,.114)，q=.312，未支持下游Bloch变化更大；M20后该比较正向。不能由此宣称全部甾醇反应高于全部脂肪酸通量。脂肪酰CoA集合包含活化和延长；所有分数均不是绝对通量单位。

**没有得到支持的过程**：胆固醇外排Full .0817、q=.423，M20−.122、q=.215；FAO与prenylation也未获主比较支持。LXR Full有变化而M20减弱。小集合GSEA因<15而NE，不等于阴性；可评分的小集合单独报告。没有可靠通用签名证明上皮向外释放游离脂肪酸，脂解/囊泡不替代这一缺口。

## 4. 同标本发送—接收关联到哪一步？

发送端主指标确实为每标本H05程序分数；接收端是真实单份Treg/Capillary EC profile，不复制成环境High/Low。两项预登记Treg主问题均NE：胆固醇摄入只有4/12固定背景基因；长链脂肪酸摄取7/19，未达到5基因和50%规则。真实profile交集18记录患者单位/5数据集，但不能因有样本就忽略基因覆盖。此处是无法评价，不能称关联不存在。

五项EC探索路线也因覆盖不足NE：EC mTORC1为98/200、LCFA摄入2/19、运输11/65、胆固醇摄入4/12。没有建立内皮营养递送证据；运输程序本身也不能确定跨内皮流向。

有限Treg次要状态中，18单位的连续OXPHOS/G2M主模型因singleton dataset出现h=1而HC3 NE；两项可估计cycling quasi-binomial的efflux和FA系数分别−.0529、.1075，95% CI分别[−.4595,.3538]和[−.0669,.2820]，BH q均1，未获支持。预设去singleton的数据集子集17单位/4数据集可估计，全部次要q约.873，未支持清晰关联；不能把子集当主分析替代。单基因检出/表达表只作描述，不拿其替代NE程序模型。

## 5. 更新假说与剩余未知

现在更支持“MP10高分尾部伴随广泛脂质重塑，指定比较中有超出Top20的胆固醇合成偏向，SREBF2参考命名修正后也有调控活动支持”。这仍不建立mTORC1→SREBF2的直接因果，也没有证明上皮脂质输出或向Treg转移。

上皮生产/合成这一端证据较强；外排→摄取→接收端功能的中间链条缺乏可估计或支持性证据。EC自身mTOR相关状态可作为限定背景的支线保留，但受覆盖、singleton和数据集构成限制。NECTIN2旧RNA证据保留，不因脂质分支强而宣布接触轴已被排除。

最有区分力的下一步是同位素示踪的上皮—Treg共培养：分别扰动上皮胆固醇合成与ABCA1/ABCG1外排，并测供体、上清及Treg中的标记胆固醇；以接触/隔膜条件区分接触与可扩散机制，再独立测cycling、呼吸和免疫调节功能。它能区分“共同脂质状态”与真正生产—输出—摄取；仅测更多转录marker不足以填补该缺口。该建议尚未实施，HD局部接触和细胞类型特异读数也仍pending。

## 完成、NE和技术限制

全部可运行审计/补充分析已完成；1328条完整结果保留阴性及NE，18组图及直接绘图数据。独立数值核验4654项、原counts/成员核验170项、渲染62文件通过。当前HGNC外网下载受限，以安装参考和NCBI身份核验完成修正；未下载第二套ranking，也未完成完整seed456 ctx复现，不能称参考/种子稳定性已验证。原大对象仅按需读取，峰值RSS约44.61GiB，未复制完整表达或数据库。
''')
doc('HYPOTHESIS_UPDATE.md','''# 假说更新

项目相关证据支持：H05上皮的mTORC1相关转录、SREBF1/2活动及多条脂质程序升高；胆固醇合成相对指定脂肪酰CoA程序的变化，在去Top20后仍得到直接比较支持。SREBF2从旧SCENIC NE变为历史名称修正后可评价，属于技术修复，不是独立验证。

尚未成立的连接：mTORC1激酶活化→SREBF成熟/结合；脂质合成→外排；上皮外排→Treg摄取；内皮运输→定向营养递送；接收端转录→实测呼吸或免疫抑制。Full Bloch直接比较、正式胆固醇外排及可估计cycling模型未获支持；摄取关联多数因覆盖或HC3条件NE，不能等同反证。

保留共同环境、RNA含量/细胞组成、数据集差异、反向作用和参考库敏感性。记录患者单位包含provisional映射，不外推为确认独立患者。

优先验证：供体上皮稳定同位素胆固醇合成与转移，合成和ABCA1/ABCG1外排分开扰动；检测供体/上清/接收Treg质量守恒。加入接触/隔膜条件和NECTIN2分支扰动，可区分接触候选与代谢环境候选。独立测Treg cycling、呼吸、免疫调节功能，不将RNA分数当功能。内皮先补同患者覆盖及纯度、再做跨内皮示踪。

HD pending：同切片、同细胞类型内的MP10/NECTIN2、Treg cycling及邻接距离；以患者/切片为统计结构，控制局部组成和深度。这里只登记问题，没有启动、修改或等待HD分析。

逐条判断与下一项实验见 `tables/decision_table.tsv`。所有结论来自同批数据定向内部检验，多算法/阈值/LODO不计为独立生物学重复。
''')
print('DELIVERY TABLES',len(r),r.status.value_counts().to_dict())
