from common import *
import zipfile,shutil

def md(d):
 d=d.astype(object).where(d.notna(),'NA');lines=['| '+' | '.join(map(str,d.columns))+' |','| '+' | '.join(['---']*len(d.columns))+' |']
 for vals in d.itertuples(index=False,name=None):lines.append('| '+' | '.join(f'{v:.4g}'if isinstance(v,float)else str(v).replace('|','/')for v in vals)+' |')
 return '\n'.join(lines)
summ=pd.read_csv(R/'tables/section_summary.tsv',sep='\t');main=pd.read_csv(R/'tables/association_primary.tsv',sep='\t');cover=pd.read_csv(R/'tables/measurement_coverage.tsv',sep='\t');allres=pd.read_csv(R/'tables/association_all_results.tsv',sep='\t');readiness=pd.read_csv(R/'tables/communication_readiness.tsv',sep='\t')
methods='''# 方法与执行记录

## 输入、分母与研究单位

本轮纳入 Cervilla、D1_Johnson、UC_SCC、MPBC1、MPBC2。定向核查 data、已有 asset manifest 和加载入口，没有发现可另行纳入的独立 HD 切片；不同分辨率文件不是新增切片。普通 Visium63 未运行。输入表为 `tables/input_manifest.tsv`，每片原始路径与文件身份保存在 `prepared.json`。5 张切片不等于5个已确认独立患者，患者映射仍未知；MPBC为两个额外MIBC微乳头/尿路上皮混合标本，不能按论文标题归为NMIBC。

复用完整 supplied 8 µm raw integer counts，共18,085基因。Cervilla/D1/UC有原组织范围；MPBC1/2来源为上轮按固定原点从实际2 µm counts聚合的8 µm，原组织掩膜保持未知。没有将 filtered 文件名或表达阈值当作组织掩膜。聚合前未删除低UMI bins。确认组织外才排除；未知组织归属保留。所有结果是区域测量，非真实Treg细胞数或密度。

坐标、物理尺度和Cervilla/D1图像变换复用已验证的节点表与 v2 transform。UC使用已有hires图及scalefactors；MPBC仅相对阵列坐标。x向右、y向下；未提供图像/表达位置保持缺测。支持面积是实际成员bins×64 µm²，MPBC不能称为组织面积。输入连通组件不是病理分区。

## 区域构建与来源分离

T主测量在原阵列0点的32 µm非重叠格内按8连通、同输入组件汇总；16 µm为唯一尺度敏感性。保存全部区域，不只旧T/NK或有标签区域。每个独立分支核对逐基因和全局 `included+excluded=input`。未知、实测零、测量不合格独立记录。

上皮支持使用未经平滑的10基因面板(EPCAM,KRT7,KRT19,UPK2,KRT5,KRT14,KRT6A,KRT6B,KRT17,COL17A1)：至少2个检出且任一冻结参考预测Epi，或至少3个检出。无Treg/MP10参与此选择。仅将支持bins在固定64 µm格内连通汇总，不跨缺口。名称为“上皮富集区域MP10”，不是已分离纯上皮RNA。所有成员映射以唯一bin ID保留。

位置为成员bin等权质心。50/200 µm为精确累计欧氏半径；不跨输入组件，一个Epi区域只贡献一次。主终点排除任何共享原始8 µm bins的T/Epi区域对；包含共享来源单独作为共同检测诊断。不能把评分复制到成员后当作独立观测。邻域需至少1个独立可评分Epi区域，并且圆内实际合格输入网格支持面积/圆面积≥0.5(截顶1)；这是可观察输入支持，不是已证实的组织覆盖率。50 µm缺少目标时不扩大半径补足。

## Treg signature与参考校准

使用最终30类型精细参考 `subtype/v2/reference_training_v2/11_final_spatial_reference/reference_v2_final_30type_regression_input.h5ad` 的raw counts；复用旧冻结参考成员20,446，18,289训练/2,157留出，91个留出Treg。留出按orig.ident样本，不称患者独立验证。类型背景保留Tfh/CXCL13-like、普通/活化CD4、CD8、NK、上皮及其他非T谱系。

只在训练部分选基因：候选中Treg检出率≥0.08，平均每万UMI表达>1.2倍最高竞争T/NK类，且>1.5倍合并非T背景。最终11基因：FOXP3,IL2RA,CTLA4,IKZF2,RTKN2,CCR8,TNFRSF18,TNFRSF4,ICOS,LAYN,MAGEH1。LRRC32未通过；排除MP10及预定通讯读数TIGIT,NECTIN2,CD96,NECTIN1,PVR,CD226。未按HD关联选基因。

文献复核为De Simone et al., Immunity2016 Table S4的完整309基因，来自实际补充Excel；保存原始表、来源及逐基因修改。剔除预定读数/MP10并要求参考和HD共同可测后，评分版本283基因。原文为CRC/NSCLC肿瘤Treg，不是BLCA专属分类器；含活化和增殖相关基因，单独解释。不是自行拼接marker冒充原文signature。

参考验证含样本留出、预期500/2000/5000 UMI thinning、leave-one-gene-out、无Treg的70%上皮+30%竞争T RNA阴性混合以及10%/30% Treg RNA混合。每组100次。混合比例指RNA贡献而非细胞比例。训练期特征冻结后才用于留出；不将验证集重新选择特征。

## UCell与测量资格

实际R4.5.3、UCell2.14.0；软件清单及sessionInfo附后。按五片T支持32 µm区域的最低中位nFeature向下取整100，限制500–3000，得到共同maxRank=1300。保留平均并列秩，包括零值；使用完整基因背景。没有逐片调到最有利参数。

原生公开ScoreSignatures_UCell在80个完整表达示例核验。全量使用UCell2.14 `u_stat_signature_list`原生评分core、明确maxRank=1300，对在完整18,085背景计算的精确秩分块执行。原因是公共precalc.ranks会从缓存最大值推断cap，精选基因/半整数并列秩缓存不适于该推断；不是换用自定义评分器。数值比较修复了浮点32位约3e-8差异，double实现与原生core/历史MP10达到<1e-12。正式区域字段使用原生UCell结果。未保存全基因标准化或稠密排名矩阵。

T可测量：nFeature≥325、≥2个T-lineage基因、T-lineage UMI≥3。这不是单细胞QC删除规则，而是本轮区域评分资格；所有原始区域保留。Treg相关连续分数为主，合格T区中的计算零值也进入主分析。描述性Treg-like层级另需主signature≥3基因、核心FOXP3/IL2RA/CTLA4/IKZF2中≥1、分数超过训练竞争T的Q95=0.199553892。这个层级不等于可信单细胞身份，不是每片强取最高百分比。主分析不要求达到这个Treg-like层级。

## MP10与冻结顺序

Top20严格复用项目冻结文件：HMGCS1,MSMO1,INSIG1,FDFT1,SQLE,FDPS,IDI1,ACSS2,ELOVL5,LDLR,MVD,SCD,HMGCR,ACAT2,SLC20A1,ACSL3,ACSL1,ADIPOR2,MXD1,NPC1。历史公式：全18,085基因平均并列秩，sum(max(905-rank,0))/sum(905-1,...,905-20)。无默认AUCell替代、无平滑或HVG裁剪。每个评分区域≥905检出基因、≥5个Top20成员检出；原始Top20未更换。仅为区域测量资格，不把旧“5细胞”改成“5 bins”。合格计算零保留，不合格不是Low。

单次半深度测试和分数-检测量关系用于测量限制记录。全片合格上皮区域type7分位数生成L50≤Q50、H50>Q50、H10>Q90、H05>Q95，不拆并列。连续值仍为主。

配置于2026-11-17 00:28 CST冻结，早于正式空间关联；参考与几何先行。ROI使用T证据、空间分散及T强度匹配对照，无MP10关系参与。24个最终ROI清单在关联计算已启动但结果尚未查看时保存，明确记录在models/roi_freeze.json；不是在关联图上追选重叠热点。少量候选(Cervilla/MPBC各2)全部展示；D1/UC按预定分数范围及空间分布选4并配2对照。

## 逐片统计

因变量为独立Epi评分区域MP10等权均值，自变量为连续主Treg UCell。效应是标准化平均秩OLS斜率，不是表达倍数或细胞比例。未调整和固定协变量调整均保留。固定顺序：总T score、T UMI、Epi平均UMI、Epi支持面积、输入圆覆盖、原始Epi-Fib表达背景。常量或VIF≥10的后续变量删除，具体实际协变量写每行，不按方向或显著性选择。Tfh/Activated_T附加条件作为单独敏感性。

空间400 µm区块根据最大邻域及边际相关诊断事先固定，800 µm为唯一更大区块敏感性。至少50个T区域、12个区块及足够不同分数才做推断。CI为499次空间区块pairs bootstrap；P为限制零假设模型的HC2残差、区块Rademacher wild bootstrap(499次，seed123)。不是普通bootstrap系数尾部比例，也不全片逐bin打乱。区块模型仍是片内探索性推断，边界重叠和未测混杂不能完全消除。

主家族为5片×2半径，10行全部保留、7行可检验，BH对非NA P执行。文献、16 µm、包含共享来源、总T、Tfh、活化T、反向分析、更大区块各自单列家族；全结果保留。反向分析以Epi区域为单位，不当独立重复。无法确认患者独立性且测量支持不同，不输出合并总体P或meta效应。

## 旧细胞、SpottedPy与通讯前检查

Cervilla复用上轮1,209个细胞限定Epi区域MP10，不重跑分割或细胞评分。通过原2 µm bin-to-cell映射及原阵列floor/4恢复到8 µm来源，保守排除任何共享8 µm格。除原覆盖分析外，使用完全相同247个T中心和相同协变量比较两种MP10来源；50 µm共同仅7个，描述性。其他片没有可直接复用的同类完整细胞限定MP10区域结果，不补造。

SpottedPy0.1.7源代码锁定SHA，原样提取calculateDistancesHelper执行，不安装/替换其旧依赖；这是原生距离helper适配调用，不声称运行Gi*。预定Treg-like位置与Epi H50/H10/L50位置作物理µm最近距离，每次一个源位置、同输入组件、排除共享原始bins，未找到或距离超出可观察边缘记NA/截尾，不填最大值。仅辅助可视化，不代替主统计。

冻结scRNA候选为TIGIT_NECTIN2、CD96_NECTIN1，来源mp10_cellchat_high_definition_sensitivity_v1 H05表，CellChatDB.human1.6.1。保留数据库编码的ligand/receptor字段，但TIGIT/CD96受体型分子与上皮结合伙伴的生物学方向不能由该字段顺序证明。只核对raw表达、来源分离、50/200 µm邻近和MP10组描述，不算配体受体得分。没有运行CellChat/LIANA/NicheNet/TF等。

## 绘图、资源、验证

采用已安装jins-routine-plotting规范。其常规模板不足以覆盖本次空间图册/森林图布局，使用matplotlib自定义布局，Arial约12 pt、白底300dpi PNG、PDF文字和图例矢量、密集点及H&E栅格。有图例输出full/panel_only/legend_only PDF。黑比例尺64 µm，红圆50 µm；ROI256 µm。橙色空心方块标中心T区域实际成员bins，青色空心方块标可评分Epi区域成员bins；不画跨缺口的凸包或插值边界。空白是无输入，灰底为supplied，不可评分在覆盖图单列；分数零不等于缺测。全片底图可抽稀显示，所有分析和评分使用全量合格区域。

默认总CPU亲和性4，BLAS/OpenMP/Numba1、R --vanilla及运行时线程控制，稀疏计数/空间索引。用户明确“没有必要等，直接算就可以了”后本轮允许日间继续，记录runtime_authorization.json；没有升级共享环境或重复上游BANKSY/SingleR。最大记录单进程RSS约5.1 GiB，原生UCell全量14分09秒。资源表列实际计时步骤，不把多进程峰值相加冒充实测总峰值。

实际重读所有区域/成员/邻域、核对唯一ID、坐标质心、counts来源、逐基因守恒、评分来源、缺测、共享bin排除、50/200几何与同组件、统计区域不复制；Python/R loader实际读回五片。图像检查见ROI_VISUAL_REVIEW.md。上游原counts/QC不修改。

## 方法来源

- [UCell参数](https://bioconductor.org/packages/release/bioc/vignettes/UCell/inst/doc/UCell_parameters.html)：以本机2.14实现核对，不套用默认maxRank。
- [De Simone2016](https://pubmed.ncbi.nlm.nih.gov/27851914/)；[完整XML与附表入口](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5119953/fullTextXML)。本地sources/DeSimone_tables.xlsx Table S4为309基因来源。
- [Oh2020 BLCA T细胞](https://pmc.ncbi.nlm.nih.gov/articles/PMC7321885/)：背景参考，不将自组panel称其完整signature。
- [SpottedPy作者源码](https://github.com/secrierlab/SpottedPy)：本地源码快照与实际helper执行记录。
'''
(R/'METHODS.md').write_text(methods)
results='''# 全部HD切片：Treg相关区域与邻近上皮MP10

**五片全量可执行分析已完成；共同bin区域主线没有检出可跨片重复的正向Treg–MP10关联。** 这不是证明不存在联系：Treg特征稀疏、竞争状态和上皮评分来源均限制解释。主分析涵盖有T-lineage支持的连续Treg相关分数，而不是已确认Treg细胞群。

## 覆盖与可评估性

'''+md(summ[['section','T_measurement_regions','Treg_like_regions','Epi_MP10_regions','comparable_T_regions_50','comparable_T_regions_200','mask_status']])+'''

Treg-like为描述性多基因证据层，不是Treg细胞数。五片均有可测T区域和Epi区域，50 µm的Cervilla/MPBC1/MPBC2支持不足以推断；200 µm全部可做片内探索性分析。MPBC的统计分母仅是输入/表达支持范围，未伪造组织掩膜或真实Treg/mm²。所有实际零分和未测量状态分别保留。

## 主效应：32 µm T区，排除共享原始bins，调整总T/深度/覆盖/背景

'''+md(main[['section','radius_um','n','n_blocks','effect','low','high','P','BH_q']])+'''

效应为标准化秩斜率，low/high为空表示未进行有效区间/P推断。7个可检验主终点均未达到BH q<0.05。D1 200 µm为弱负向(beta=-0.041，P=.014，q=.098)，不是预期正向；Cervilla接近零，UC/MPBC为小正值但区间包含零。样本身份独立性未知，未将五片合成“5患者重复验证”。

## Signature与测量质量

11基因参考主signature在样本留出中区分Treg/其他T类AUROC约0.893；与Tfh单独比较0.819，不能宣称完全特异。500 UMI稀疏化后约0.830/0.776；参考结果支持组合有信息，但不把单细胞参考性能当HD区域身份准确率。单基因去除、阴性上皮+活化T混合和少量Treg RNA混合结果完整保留。低比例Treg在混合背景中的识别能力有限。

De Simone完整309基因可追溯，本次283基因修改版单独保存；竞争T AUROC约0.882、Tfh约0.777。文献版含活化/增殖信号，其空间结果不能视为独立金标准。

T32半深度排序rho约0.55–0.66；Epi MP10约0.57–0.81，属于有限稳定性。Top20在五片均完整可测，历史K905/平均并列秩数值核验通过。合格Epi仍可能混有其他谱系RNA，所以主结果称“上皮富集区域MP10”。

## 调整、来源分离与敏感性

D1未调整主分数在50/200 µm为负向约-0.091/-0.136，控制总T、检测量等后减弱。200 µm使用800 µm区块时仍负向(P=.006，敏感性家族q=.042)，加入Tfh/活化T后约-0.034(P=.036，q=.252)。文献signature在同片200 µm方向转为小正值但无明确支持，提示状态/定义依赖，不能把任何一个更显著版本选作统一结论。

允许共享表达来源会增加可比较位置，UC 50 µm从418增至851，效应从约0.045增至0.080(P=.026，但该敏感性家族q=.117)。独立来源主结果不支持这一名义正向信号。其他全部敏感性、16 µm、反向、一般T和竞争状态结果见association_all_results.tsv，未只保留显著行。

Cervilla已有细胞限定MP10给出了值得复核的测量来源差异：在完全相同247个T中心、相同协变量、200 µm下，cell-limited beta=0.210，CI[0.065,0.335]，P=.018；bin-enriched beta=0.080，CI[-0.070,0.221]，P=.310。800 µm区块cell-limited P=.034，单列敏感性家族q=.068。此辅助结果提示上皮表达分离/支持范围会影响关联，不能据此替换五片预定主终点或宣称跨片阳性。50 µm共同仅7个位置，不做推断。旧细胞结果也未提供本轮可确认的Treg单细胞位置，因此目前没有更强的Treg细胞身份验证。

## ROI与原始表达

24个候选/配套对照ROI全部输出。Cervilla代表候选处有密集核和T表达，Treg核心UMI零散、部分伴CXCL13/PDCD1；不能由H&E或单个FOXP3证明Treg。D1展示的候选中有较多IL2RA/CTLA4/IKZF2/CCR8组合支持，但存在Tfh或活化竞争。UC下采样图只支持粗组织核验，Treg核心信号稀疏；MPBC仅相对网格，不能补造形态结论。ROI中Epi MP10面板空白常表示该256 µm窗口没有合格Epi评分区域，不是MP10生物学零。

SpottedPy原生距离helper已对预定证据位置与Epi H50/H10/L50运行，真实µm、同组件、共享来源排除及观察边缘截尾。目标少且截尾多，作为辅助定位，不产生额外阳性检验。未声称运行完整Gi*热点模型。

## 通讯前候选与下一步

预定TIGIT–NECTIN2和CD96–NECTIN1原始两侧表达/邻近支持已完成。200 µm下，D1分别有1787/1626个一般T测量区具有两侧检测，其中仅4/7个达到本轮Treg-like层级；Cervilla两轴各2个Treg-like位置，其他三片对应支持更弱。数据库的编码方向不能直接解释为已证实信号方向。没有计算LR配对分数。

D1适合继续做区域候选轴空间支持与目标身份核验；Cervilla应优先复核细胞限定MP10与bin区域测量差异；UC/MPBC当前主要支持一般T细胞层面。所有片的“直接进入Treg特异正式通讯网络”均为False。signature-high混合区不能全部重标为Treg输入CellChat。

## 完成与未进行事项

- 已完成：5片所有supplied输入的32/16µm区域、上皮64µm连通区域、原生UCell、历史MP10、参考验证、50/200几何、全部合格区域的主/敏感性统计、24ROI、SpottedPy距离辅助、两条候选轴、Cervilla旧cell对照、加载和守恒验证。
- 已计算但不可有效推断：Cervilla/MPBC1/MPBC2主50µm，及表中明确标记的部分敏感性终点。不是漏跑；P保留NA。
- 输入不支持：MPBC组织形态/真实组织面积与患者级合并推断；除Cervilla外无可直接复用的同定义完整cell MP10区域结果。未伪造。
- 按范围未进行：新分割/全量SingleR/BANKSY/四方法niche、普通Visium63、正式通讯网络/NicheNet/TF。不存在排队冒充完成的分支。

所有上游文件保持只读。详细参数、数据字典、路径和资源见METHODS/README及tables。没有为获得正向结果更改半径、signature、ROI或选择切片。
'''
(R/'RESULTS.md').write_text(results);(R/'REPORT.md').write_text(results)
readme=f'''# hd_treg_mp10_spatial_all_v1

实际完成五片全量可执行分析。先读 [RESULTS.md](RESULTS.md) 和 [METHODS.md](METHODS.md)。主结果未检出跨片可重复正向关联；部分终点可评估性不足，详见逐行状态。

## 目录与依赖

结果根目录：`{R}`。
原8µm矩阵、空间节点、旧细胞结果与冻结参考均通过 `tables/input_manifest.tsv`、每片 `prepared.json` 和代码引用，没有复制整套counts/RDS/H5AD。
Python环境：`{A}/hd_epi_niche_multimethod_v2/.venv/bin/python`；R库路径见 `code/04_native_UCell.R` 引用的项目common.R，不修改共享库。

- 每片 `T32_regions.parquet`：所有32µm区域，分数、原始支持基因UMI/名单、资格和证据层。
- `T16_regions.parquet`：唯一测量尺度敏感性。
- `Epi64_regions.parquet`：上皮支持连通区域、MP10及合格分组。
- `*_members.parquet`：唯一原始bin→区域映射；不将区域分数复制成独立bin观察。
- `T32_neighborhoods.parquet`/`T16_neighborhoods.parquet`：50/200µm，两种来源规则、覆盖/资格/邻域均值。
- `region_pairs.parquet`：精确距离、共享来源标记；不包含全空间N×N图。
- `bin_evidence.parquet`：定位及原始目标panel表达，不是完整表达副本。
- `*_measure.parquet`：用于精确复跑的选定基因raw counts/完整背景秩缓存，非全基因矩阵。
- `tables/association_primary.tsv`、`association_all_results.tsv`：全部终点、效应、CI、P、BH及实际协变量。
- `tables/measurement_coverage.tsv`、`section_summary.tsv`：可评估性分母。
- `figures/`、`FIGURE_INDEX.md`：全图及索引；review.zip不带原矩阵和训练参考。

## 读取

```python
import sys
sys.path.insert(0, r"{R}/code")
from load_result import load_section
r = load_section("Cervilla", include_members=True, include_pairs=True)
t = r["T32"]
# 统计使用区域表；按bin ID挂回只附成员ID。
```

```r
root <- "{R}"
source(file.path(root, "code/load_result.R"))
t <- load_treg_mp10(root, "Cervilla", "T32")
e <- load_treg_mp10(root, "Cervilla", "Epi64")
# existing <- attach_treg_mp10_membership(existing, root, "Cervilla")
```

## 复跑与重绘

`bash rerun.sh` 断点续跑，完成步骤跳过；输入身份与冻结配置保留检查。`bash redraw.sh` 仅从已保存表/原图像重绘，不重跑统计。完整强制复跑 `FORCE=1 bash rerun.sh`，按已冻结参数执行，不重新挑signature/ROI/最佳关联。上游依赖必须存在；不能用review包代替原输入。
本轮用户已授权立即计算，runner设置HD_ALLOW_DAYTIME=1；默认CPU总亲和性4、数值线程1。其他任务的夜间约束不由此改变。

## 字段与限制

`measurement_eligible`是本轮用途资格，原QC不改；`evidence`区分Treg-like、多基因一般T、信息不足。`MP10`/各UCell值对不合格区域也可保留计算诊断，但仅在资格True时用于正式终点；不可测量不等于Low。`tissue_status=unknown`仍为未知。
`effect`是标准化秩回归效应，`status`标明片内探索性或仅描述性；不是细胞比例、因果效应或多患者验证。
`support_area_um2`为成员网格面积，MPBC不能解释为病理组织面积。

检出基因名单和UMI直接来自raw区域计数。原始symbol/Ensembl完整映射在tables。未评价正式通讯，所有formal_communication_ready=False。
'''
(R/'README.md').write_text(readme)
idx=['# Figure index','', 'PNG 300dpi；PDF文字/坐标/图例矢量，H&E和密集点为位图。带图例图另有panel_only/legend_only PDF。','', 'ROI黑线=64µm，红圈=50µm，坐标y向下。橙色空心方块为中心T测量的真实成员bin位置；青色空心方块为可评分Epi区域的真实成员位置，不是插值轮廓。只有显示底图抽稀，分析未抽样。灰色supplied背景不等于生物零；MPBC为相对坐标、组织未核验。','', '| File stem | Data / interpretation |','|---|---|']
for f in sorted((R/'figures').glob('*.png')):
 stem=f.stem
 desc='每片所有合格区域；仅显示质心分数'if 'whole_scores'in stem else 'supplied/合格/不可评分的不同状态'if 'measurement_coverage'in stem else '固定ROI；H&E或网格、raw T、主/文献score、raw核心、Tfh、MP10、深度；中心和选择见roi_manifest'if 'atlas'in stem else '对应tables同主题数据；统计单位为区域，非bin重复'
 idx.append(f'| {stem} | {desc} |')
(R/'FIGURE_INDEX.md').write_text('\n'.join(idx)+'\n')
print('REPORTS WRITTEN')
