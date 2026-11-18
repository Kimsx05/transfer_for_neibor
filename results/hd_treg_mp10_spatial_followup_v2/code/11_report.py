from fu_common import *
import re
read=lambda f:pd.read_csv(R/'tables'/f,sep='\t')
def md(d):
 def f(x):
  if isinstance(x,(float,np.floating)):return 'NA'if pd.isna(x)else f'{x:.4g}'
  return str(x).replace('|',' / ').replace('\n',' ')
 return '| '+' | '.join(d.columns)+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join(f(x)for x in row)+' |'for row in d.itertuples(index=False,name=None))
pair=read('paired_effects.tsv');step=read('stepwise_decomposition.tsv');base=read('baseline_reproduction.tsv');flow=read('qualified_HD_flow.tsv');ev=read('five_section_evidence_changes.tsv');orig=read('v1_five_section_main.tsv');decision=yaml.safe_load((R/'evidence_decision.yaml').read_text());support=read('common_support_paired_effects.tsv');lobo=read('leave_one_block_out.tsv');cov=read('source_coverage_summary.tsv');parts=read('mixture_partitions.tsv');rules=read('training_negative_thresholds.tsv');worst=rules.loc[rules.candidate_threshold.idxmax()]
resource=[]
for p in (R/'logs').glob('*.log'):
 s=p.read_text(errors='replace');rss=re.findall(r'Maximum resident set size \(kbytes\): (\d+)',s);elapsed=re.findall(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([^\n]+)',s);exitcode=re.findall(r'Exit status: (\d+)',s)
 if rss:resource.append(dict(log=p.name,peak_RSS_GiB=int(rss[-1])/1024**2,elapsed=elapsed[-1]if elapsed else'NA',exit_code=int(exitcode[-1])if exitcode else-1))
put(pd.DataFrame(resource),'resource_usage.tsv')
roi=read('ROI_evidence_notes.tsv');rn=read('ROI_neighborhood_coverage.tsv');rtext=['# 固定 ROI 逐区复核','保留原24个ID及中心，14候选、10对照。512 µm只扩大显示视野；所有统计半径仍为50/200 µm累计圆，热点排除仍使用原256 µm候选矩形。原始单基因计数以8 µm supplied bins显示，评分以T32/上皮区域质心及成员显示。配色截顶见色条三角，灰色不是MP10-low。']
for rr in roi.itertuples():
 q=rn[rn.section.eq(rr.section)&rn.roi.eq(rr.roi)];q50=q[q.radius_um.eq(50)].iloc[0];q200=q[q.radius_um.eq(200)].iloc[0]
 rtext+=['',f'## {rr.section} / {rr.roi}',f'中心 `{rr.center_id}`。主signature实测基因：{rr.raw_main_genes}；核心：{rr.raw_core_genes}。原/新证据：{rr.old_evidence} / {rr.new_evidence}。',f'Treg UCell={rr.Treg_main:.4f}，Tfh={rr.Tfh:.4f}，Activated T={rr.Activated_T:.4f}；分数不是身份概率，不能直接按不同signature分数大小裁决类型。',f'排除共享表达来源后，50/200 µm可评分bin上皮区域分别为 {q50.bin_score_regions}/{q200.bin_score_regions}；200 µm平均区域UMI={q200.bin_mean_depth:.1f}，观测输入比例={q200.input_observed_fraction:.3f}。'+(f'Cell来源邻居分别为 {int(q50.cell_score_regions)}/{int(q200.cell_score_regions)}，200 µm平均区域UMI={q200.cell_mean_depth:.1f}。'if rr.section=='Cervilla'else'本分支无可复用的cell限定MP10结果，图中记未测量，未重新分割。'),('多基因证据允许区域Treg-related解释；仍不能把混合区域全部bins视为Treg细胞。'if rr.old_evidence=='Treg_like_multigene'else'T-lineage/竞争状态和稀疏性仍是解释；未达到联合规则不能等同于无Treg。'),('200 µm仍无可评分上皮，邻域MP10缺失，未移动中心或扩大统计半径。'if q200.bin_score_regions==0 else'完整200 µm中的上皮评分来自图示蓝色成员与邻居质心；详见逐邻居来源表。'),('Cervilla两候选相距约32 µm，属于同一热点邻接测量，不能当独立重复。'if rr.section=='Cervilla'and rr.roi.startswith('candidate')else'')]
 if rr.section.startswith('MPBC'):rtext+=['无H&E/原组织掩膜；这里只展示相对坐标及表达覆盖，不声称组织内面积或细胞边界。']
 elif rr.section=='UC_SCC':rtext+=['H&E仅下采样hires，用于粗结构，不作单细胞边界判读。']
 rtext+=[f'[512 µm 图册](figures/{rr.section}_{rr.roi}_512atlas.png)']
(R/'ROI_REVIEW.md').write_text('\n\n'.join(rtext))
report=f'''# Treg–MP10 v2 探索性来源复核

本轮三分支已执行：Cervilla共同中心/共同设计矩阵及计数级共同空间支持敏感性；五片混合测量校准；全部24个固定ROI的512 µm图册。它是看过v1后的探索性复核，**不是独立验证或确认性研究**。未运行通讯、重新分割、BANKSY、Harmony或新niche。

## 1. 基线、中心集合与来源

本地冻结配置、旧共同中心代码及五片主效应表，与指定提交 `e16e986fc776329cbb6f77b686131d1aecc98c42` 的相应文件逐字一致。真实原始矩阵仍在服务器，未将review包当作原始表达。基线效应复现，详见 `baseline_reproduction.tsv`。

旧共同中心200 µm结果是bin β=0.07959、cell β=0.21027，但两模型名为Epi_nCount的协变量**数值随来源改变**。本轮在247个相同T中心，用同一数值设计矩阵同时纳入bin与cell来源深度；共同面积不可定义，均不填0或强行加入。

{md(step[['step','unit','n','effect','low','high']])}

限制中心集合而保留原完整bin模型，仅使β从−0.0067变为0.0046；改为共同协变量后的bin β为0.0709。最后同设计矩阵的cell−bin差约0.1407。不能把整段变化都解释为分割或纯度改善。A/B的秩标准化按各自集合重算，因此包含集合改变的测量尺度影响。

## 2. 配对来源效应与稳健性

{md(pair[['variant','block_um','n','n_blocks','beta_bin','beta_cell','delta','delta_low','delta_high','P_delta','BH_q_auxiliary']])}

唯一主要新增检验是共同模型200 µm的Δβ=0：**Δβ=0.1407，配对400 µm区块95%区间[−0.0211,0.2641]，受限零假设空间wild bootstrap P=0.105**。这是标准化秩关联差，不是MP10活性的增量。cell自身区间不跨0，不能替代来源差值检验。

800 µm重抽样的百分位区间略高于0，但对应零假设wild P=0.088；二者使用不同的抽样分布，不要求恰好反演，不选其中更有利的结论。竞争Tfh/活化T调整后Δ约0.1393；与Treg主signature无基因交叠，当前共同中心分数相关也不强，未发现严重共线性。原热点范围排除12个T中心后，Δ约0.1058，区间仍跨0。热点排除仅改中心，剩余邻域与结局不变。

400 µm LOBO差值范围 {lobo[lobo.block_um.eq(400)].delta.min():.4f}–{lobo[lobo.block_um.eq(400)].delta.max():.4f}；800 µm为 {lobo[lobo.block_um.eq(800)].delta.min():.4f}–{lobo[lobo.block_um.eq(800)].delta.max():.4f}。正向来源差不只由单个区块产生，但统计不确定性和覆盖偏差仍在。50 µm仅7个共同中心，保留覆盖/描述，不输出可靠推断。

新增辅助P统一登记并BH校正，含800 µm、竞争状态、热点和B模型诊断；不把旧v1 P/q混入新家族。见 `new_test_registry.tsv`。

## 3. 同中心并非同上皮范围

247中心的两来源8 µm支持位置Jaccard中位数为 {cov.median_Jaccard.iloc[0]:.3f}，无中心具有完全一致来源范围。raw MP10相互相关不等于测量来源可交换。

严格保留两边完整旧评分区域，在共同footprint内找不到同时可用的区域。因此另外执行了一次有限的**原始counts共同范围重建**：仅保留全部分配位置落在共同范围内的完整可信上皮细胞；固定64 µm原始网格内合并接触的来源footprints；bin端使用相同8 µm位置的全部raw counts。每个T中心内每个细胞及bin各计一次，不按面积切细胞、不复制旧分数。494条中心×来源账目均逐基因守恒，记录几何排除UMI与共同评分资格排除UMI；重叠邻域的账目不能跨中心相加。几何核查发现旧cell评分区域可能共用8 µm位置，已按几何合并以消除重复计数。

最终1864个几何评分单元中585个两侧均达到原MP10评分资格，覆盖170个T中心。两来源采用同一设计矩阵，cell深度因与已纳入协变量残差方差≤0.1在两模型中共同省略，见设计文件。结果：

{md(support[['block_um','n','n_blocks','beta_bin','beta_cell','delta','delta_low','delta_high']])}

共同范围差值约0.0393，区间跨0，明显小于原共同中心结果。但这同时改变评分区域、可测量中心和覆盖；不能据此证明哪种来源是真值，也不能将减少归因于某一个机制。相同8 µm footprint仍不是相同RNA来源：bin含该位置全部RNA，cell只含可信上皮对象RNA。本敏感性不新增正式P值。

## 4. Treg混合测量的漏检与单一校准

完整重放2700个v1 70%上皮背景模拟，选定基因的全背景秩及检测量与旧缓存一致。扩展模拟固定五种背景、三种T竞争背景、三档总UMI和四档Treg/T RNA比例，共21600条。Treg占总RNA与占T RNA分别保存，并额外保存固定模拟下实现的RNA来源量；均不是细胞比例。

{md(parts)}

参考按样本留出，非患者独立验证；同一cell及所有派生profile不能跨分区。激活CD4以固定训练分位数对既有CD4 profiles的Activated_T表达作分层诊断，未伪称参考有独立“活化CD4”最终标签。高活化留出背景数量偏少，不能给出强泛化保证。

实际HD中，漏检不只是最后的高分门槛。合格T32的顺序流失如下（最后一列是通过该阶段后仍保留数，失败原因不可直接相加）：

{md(flow[['section','qualified_n','stage','first_loss','survives']])}

五片分别有78、3443、606、49、341个**T测量合格且检出signature，但UCell=0**的区域；逐例复核发现所有入选marker平均并列秩均在1300窗口之外，属于秩截断，而非没有检测到RNA。多数区域还先受多基因/核心支持不足限制。降低高分门槛无法恢复这些连续测量损失。

训练负混合物分层误报目标生成唯一候选门槛 **{decision['candidate_threshold']:.5f}**（旧为0.19955）。最严格背景为 `{worst.background}` / `{worst.competitor}` / {int(worst.depth)} UMI。候选门槛反而更高。留出正混合物完整联合规则检出率从 **{decision['heldout_old_positive_sensitivity']:.2%} 降至 {decision['heldout_new_positive_sensitivity']:.2%}**；候选最大条件误报{decision['max_heldout_negative_FPR']:.2%}，但不满足预定检出改善/条件无损标准。**不采用新规则**，不再尝试第二套门槛。

这些性能是有限参考profiles条件下的模拟检出/误报，不是HD FDR或患者泛化性能。500–5000 UMI是候选校准范围；范围外不补造性能。分层完整性能、原70%重放、顺序与重叠失败均已保存。

## 5. 五片完整保留

{md(ev[['section','all_T32','T_eligible','old_Treg_like','new_Treg_like','eligible_detected_but_zero','eligible_outside_range']])}

旧/新标签不变；主11基因、文献对照、连续UCell、T32资格、五片主分析集合及MP10定义均未改。只有Cervilla公平来源比较增加了明确的新共同模型和敏感性测量；ROI视野从256扩大为512 µm，统计半径不变。

{md(orig[['section','radius_um','n','effect','low','high','P','BH_q']])}

表中全部P/q来自v1，本轮未重跑五片连续关联。方向不一致，没有跨五片一致的正向证据；D1 200 µm负向线索旧q=0.098。Cervilla/MPBC两片50 µm支持不足，仅描述。MPBC组织范围仍未知，患者跨来源独立性仍未知，不合并成五患者验证。

## 6. ROI与结论边界

24/24原固定ROI全部绘为512×512 µm；保留14候选、10对照，并显示50与200 µm累计圆、原始四核心基因计数、T-lineage、主UCell/Tfh/活化T、上皮评分范围、实际来源成员及共享来源排除。见 [ROI_REVIEW](ROI_REVIEW.md) 与 [图索引](FIGURE_INDEX.md)。

Cervilla候选1有FOXP3、IL2RA等6个主基因；候选2只有IKZF2、TNFRSF18、TNFRSF4三基因。两者相距约32 µm，是同一热点内邻接测量；两者都有Tfh相关信号，区域Treg-related解释不能升级为已确认Treg细胞身份。两者50 µm均无可评分上皮邻居；200 µm分别只有2个bin上皮区域与1个cell限定区域，主要在图中右侧边缘的上皮支持范围。更大显示视野解释了“信号在哪里”，没有产生新的独立观测。

UC_SCC下采样H&E只支持粗结构；MPBC无H&E/原组织掩膜，空白是缺少来源覆盖而不是MP10-low。未测量、未达到资格、合格计算零值与有分数低值均分开显示。

## 7. 下一步与实际完成状态

- 来源公平对照、400/800区块、竞争状态、LOBO、原热点排除、原始counts共同支持敏感性：已计算。
- 五片原规则流失、2700基线重放、21600有限混合、唯一门槛训练/留出决策：已计算；新规则未通过，不采用。
- 24ROI及全部五片证据/旧主效应总览：已交付。
- 50 µm共同来源推断：已核对，但7个中心不足，未作有效检验。
- 旧完整评分区域的严格共同支持：已计算，零可用；以不拆细胞的计数级重建完成可行敏感性。
- 其余四片cell限定MP10：本分支无可复用正式结果，不补做分割；MPBC病理边界不能恢复。

目前适合停留在**区域关联及局部身份核验**。可优先核对Cervilla同一热点的Treg与Tfh/活化T竞争身份，并评估小淋巴细胞表达归属，但没有证据支持为追求关联而全片重分割或把cell主线普遍判为更优。需要改善连续Treg测量时，应另立测量校准任务；单改高分标签不能“挽救”原连续关联。所有 `formal_communication_ready=False`，本轮不运行通讯。

## 8. 交付与资源

加载示例见 [README](README.md)。本轮仅保存表、映射、图和代码，counts通过原路径读取。配置、逐项验证和软件版本分别见 `frozen_config.yaml`、`tables/validation.tsv`、`tables/software_versions.tsv`。运行日志记录实测资源；4核CPU亲和性、数值线程1，无并行任务嵌套。沿用本会话用户“直接算、不必等待”的明确授权完成本轮有限复核；没有把排队视作完成。实测峰值RSS约{max(x['peak_RSS_GiB']for x in resource):.2f} GiB。最终磁盘占用与压缩包大小见 `DELIVERY.md`。
'''
(R/'REPORT.md').write_text(report)
methods='''# 实际方法与复跑边界

本轮为v1结果后的探索性复核。原counts、QC、评分、参考与结果只读；配置于2026-11-17 14:18 CST冻结，之后只进行了记录在common_support_reconstruction.yaml中的几何去重修复。没有基于MP10关联调signature、半径、ROI或高分门槛。

## 输入身份与版本

`tables/input_manifest.tsv`列出五片真实H5、坐标节点、18085输入基因、数据集和未知患者关系。Cervilla native8输入来自上一轮native8目录；cell计数来自segmentation_assignment_v2/Cervilla/cells.h5ad，旧可信上皮来源由hd_mp10_treg_endo_precommunication_v1/tables/epithelial_region_membership.tsv.gz限定。源2 µm归属映射用于恢复cell表达footprint，转换为row8=floor(row2/4)、col8=floor(col2/4)，不移动原点。

远程仅核对指定基线提交的三个小文件（旧共同中心脚本、配置、主表），通过GitHub只读接口；未下载新原始数据、未上传数据。精细参考使用subtype/v2/reference_training_v2/11_final_spatial_reference/reference_v2_final_30type_regression_input.h5ad及v1固定成员。混合参考为17987个唯一symbol，完整参考背景秩；HD为18085原基因背景。差异继承v1，不将参考/HD基因背景混称相同。

## A：配对设计与统计

共同中心T32必须在两个200 µm来源下都有可用独立上皮结局及共同协变量，无共享原始T32 bins。50 µm不放宽资格。协变量固定顺序：T-lineage、T侧nCount、输入圆内可观察比例、原T区域上皮/成纤维背景平衡、bin上皮深度、cell上皮深度。两侧深度同时进入两套模型。面积无法可比，均不使用。竞争状态加入Tfh及Activated_T。

平均并列秩→总体标准差标准化（ddof0），在固定共同集合上一次完成；bootstrap不重新计算秩。按固定顺序计算每列对已纳入设计残差方差，>0.1才保留，相当于预定VIF<10原则。共同模型使用完全相同的X、相同Treg向量、相同中心顺序；未按P或方向删变量。缺失采用两结局及共同协变量联合complete-case，缺面积不填0。

区块由原component及floor(x/400),floor(y/400)定义，另作800 µm。999次有放回抽取区块；两模型使用完全相同索引，直接保存每次βcell、βbin、差值。CI为配对百分位区间。唯一主要新增P针对Ydiff=zrank(cell)-zrank(bin)，在不含Treg的受限共同模型下拟合零假设，用HC2残差、每区块一个Rademacher符号进行999次wild bootstrap，P=(1+|null|>=|observed|次数)/1000。它是依赖区块内相关、区块间近似独立与对称误差的片内探索性零模型，非患者检验、非普通bootstrap尾部P。百分位CI与该零假设P不互为严格反演。

800 µm、竞争状态、原热点排除和新增B原完整模型P进入统一辅助BH家族；来源各自β仅CI，不用两个P之差判断来源差。LOBO不重新选择变量/重做秩。原热点定义是v1两个候选256 µm矩形并集；本轮512 µm显示框不参与排除。

共同支持重建将整个可信上皮cell限定在两来源实际8 µm成员交集及单个原始64 µm网格内，删除与T32共享的来源。触碰的cell footprints按8邻接合并，避免两个历史cell区域在同一8 µm位置重复使用bin counts。整cell counts从原h5ad求和；bin端从相同footprint的原H5求和；两端不按面积分摊。170中心及其585对评分单元均保存原始成员。各T中心内部来源成员唯一；不同重叠邻域可引用同一评分来源，不被解释为独立重复。与原共同中心相比测量单位/中心发生变化，因此结果仅作有边界的敏感性。评分继承18085背景、平均秩、K905及20基因归一化，nFeature>=905、MP10检出>=5、上皮marker>=2；两来源同时满足才进入均值。

## B：混合校准

11基因、UCell maxRank1300、平均并列秩、T32资格及旧联合规则均不改变。旧门槛0.1995538920。无MP10、TIGIT、NECTIN2、CD96、NECTIN1、PVR、CD226用于定义Treg。文献signature只引用v1，未新选基因。

先按v1 RNG顺序重放原70%上皮2700混合，核对选定原始计数和全背景秩。扩展21600混合来自固定18,289训练/2,157留出cell、46/12样本。每表达profile按类型取5cell合并原counts再归一化为RNA概率，参考抽样不代表组织比例。背景Epi40/70/90%，及40%Epi+30%Fib或Myeloid；其余T RNA。竞争普通CD4、Tfh/CXCL13、耗竭CD8；Treg占T RNA为0/0.1/1/3/1。总UMI500/2000/5000，60profile重复。训练/留出来源完全分离，派生混合不跨分区。保留预设概率比例和用独立随机流条件归因的实现RNA来源量；后者只补充真值记录，不改变抽样表达、评分或门槛。

资格、原UMI、检出宽度、核心基因、秩窗口、旧/候选门槛及完整联合规则分别保存，并同时给顺序首次失败与单规则失败。UCell0例子逐项核查检测marker秩是否>=1300；不能把秩窗口截断等同无RNA。

唯一候选门槛是在训练无Treg混合中，各背景×竞争×深度满足T资格/多基因/核心的子集95%上取整分位的最大值。如此固定全局门槛是保守操作规则，不等于真实HD FDR。训练分层及留出完整规则性能全部输出。接受标准事先冻结：所有有支持负条件FPR<=5%、条件Wilson上界<=15%、整体阳性检出至少提高5个百分点且每条件不降低超过2个百分点。候选未通过，保留旧规则。没有第二次阈值搜索。500–5000 UMI之外不外推校准。

活化CD4补充检查使用训练常规CD4 profiles的Activated_T均值q75固定分层，在留出重复评价，不新造细胞类型。Wilson区间反映有限模拟的条件变异；父参考profiles重复使用，不能解释为患者泛化CI。高激活留出样本少，明确有限性。

## C：绘图和缺失语义

原24中心、候选/对照身份及统计半径完全不变；显示边长512 µm，质心累计圆50/200 µm。原始RNA计数来自v1已验证bin_evidence；全基因矩阵不重跑。T32与Epi64/旧cell区域成员按ID对齐，来源footprint和共享排除明确显示。每ROI的来源数量、深度、覆盖及邻居距离另存表，可直接复核。

白底Arial，PNG300dpi，PDF文字/坐标/图例矢量，H&E与密集点层栅格。提供完整、主体、图例PDF；统计图无图例者仅完整版本。灰色supplied背景不是低分，紫色叉为资格不足，色条0是合格计算零，白色无supplied，旧cell缺项明确未测量。截顶用色条上端三角及cap文本；比例尺64 µm。UC_SCC仅hires，MPBC无H&E与原掩膜，不能读单细胞形态或组织面积。

## 复用与依赖

软件版本与实际日志见tables/software_versions.tsv和logs。数值线程1、最多4核亲和性；不修改共享环境。所有新区域结果和来源映射为Parquet/TSV，没有新增全量表达对象。R加载器经Python桥读取轻量区域表，不要求再保存Seurat。旧counts通过manifest路径访问。

方法接口沿用已核验v1实现；本轮不切换到AUCell原生默认值，也不重新跑UCell参数选择。参考官方定义可见[UCell参数](https://bioconductor.org/packages/release/bioc/vignettes/UCell/inst/doc/UCell_parameters.html)，具体实现以冻结代码和数值核对为准。可复跑入口支持检查点；--force仅重算本轮输出，旧结果保持只读。replot.sh从表及已有图像重绘，不运行模型。
'''
(R/'METHODS.md').write_text(methods)
readme=f'''# hd_treg_mp10_spatial_followup_v2

先读 [REPORT.md](REPORT.md)。三分支均已执行；未采用候选新证据门槛，五片v1主结果保留。

```r
source("{R}/code/load_result.R")
r <- load_followup("{R}", "Cervilla", "regions")
p <- load_followup("{R}", "Cervilla", "paired")
```

```python
import sys
sys.path.insert(0, "{R}/code")
from load_result import load_followup, attach_region_evidence
regions = load_followup("MPBC2")
members = load_followup("Cervilla", "membership")
```

区域分数不复制为独立bin观测。加载器仅给原T32区域表追加fu2_证据字段，或返回bin→T32成员ID；不会覆盖旧labels/counts/QC。通用AnnData/Seurat对象仍通过上游loader读取。

- `tables/Cervilla_paired_analysis.parquet` + `design_r200_common.parquet`：共同247中心、两来源结局与同一数值X。
- `paired_effects.tsv`、`new_test_registry.tsv`：主/辅助检验来源与家族。
- `common_support_*`：有限counts重建、整数总量与成员，非全片新细胞对象。
- `expanded_mixtures.parquet`/`mixture_rule_results.parquet`、`mixture_parent_membership.parquet`：混合设计、真值、单规则/联合结果与参考成员。
- `*_evidence.parquet`：五片全部T32旧/新证据及秩窗口解释。
- `ROI_*`、`ROI_REVIEW.md`：24区域原始支持、完整50/200邻居及图册。
- `v1_five_section_main.tsv`：仅转录v1全五片主结果，不是本轮重新估计。

依赖根为同级 `hd_treg_mp10_spatial_all_v1`、`hd_mp10_treg_endo_precommunication_v1`、`segmentation_assignment_v2` 与manifest的真实counts/图像/参考。移动根目录时保持analysis同级结构；原始绝对路径可通过manifest检查，包内没有假装自包含的表达矩阵。

执行 `bash rerun.sh` 续跑；`bash rerun.sh --force` 重算本轮各步；`bash replot.sh` 只重绘。默认沿用项目23:00–07:00窗口；明确授权即时计算时使用 `HD_ALLOW_DAYTIME=1`。配置变化时不静默使用旧结果，应另建版本并检查输入。review.zip不包含counts/参考/大型对象。
'''
(R/'README.md').write_text(readme)
files=sorted((R/'figures').glob('*.png'));index=['# 图件索引','所有PNG为300dpi；对应同名PDF。含图例图另有_panel_only.pdf与_legend_only.pdf。H&E及密集图层是嵌入栅格。','| 文件 | 内容 |','|---|---|']
for p in files:
 desc='固定512 µm ROI；50/200 µm圆与原始来源'if '512atlas'in p.name else'五片v1主结果未改'if 'forest'in p.name else'完整联合规则条件检出/误报'if 'mixture'in p.name else'来源/证据诊断（见REPORT）';index.append(f'| [{p.name}](figures/{p.name}) | {desc} |')
(R/'FIGURE_INDEX.md').write_text('\n'.join(index));print('REPORTS_UPDATED',stamp(),flush=True)
