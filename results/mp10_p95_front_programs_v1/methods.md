# MP10 中心、邻域与正式前沿 (实施前固定)

输入为 full_v1 的63片、冻结 query、Top20 AUCell、全部55项UCell，正式v3.4 signed_distance_um (内部负、外部正)。有效组织沿用full_v1 posterior mean/QC；不新增上皮纯度或CAF门槛。肿瘤Epi指冻结query与v3.4 inside交集，仅是分析区域身份，未声称每点病理证实恶性。P95在该片全部合格肿瘤Epi中取线性分位数，保留并列；正常来源无肿瘤参考集合，保留全Epi中心分析和独立描述。

新增Puram Top30、Flores正向10，以及公开完整Puram100敏感性。原符号与候选映射保存；UCell2.14.0，maxRank1500、ties average、chunk50、单核，完整RNA宇宙重复symbol按counts求和。官方55项不重算；所有与Top20重叠的22项另计算去重叠敏感性。文献集合采用UCell，非完整复现文献原始算法。

真实array六邻接只连接观测有效组织且同组织块，R0、R1、R2、R3–5、R6–10为最短图距离的互斥层。每中心邻域可重叠，均分保留有效邻居数和物理范围；每指标分别按可用性纳入。P95集合最小图距离形成互斥集合环。d_MP10为同块最近P95中心的欧氏距离，nominal Visium 100μm间距，缺同块参考则NA；不插值。d_front仅正式严格signed_distance_um；原观测外缘距离、距离带和缺失原因一并保留，不以带中点替代。局部Epi–间质接触距离不参与本轮主检验。

前沿带[-200,0]、近内部[-500,-200)、深内部<-500。pEMT_high_front由前沿带及同片全部肿瘤Epi Top30 P75定义，不用MP10。只作标签，不把定义的高pEMT本身算独立验证。外侧组织程序单独描述。

中心MP10—各环程序、P95距离—程序、前沿距离—程序/MP10均输出原始效应及组成/深度调整。连续距离另用固定二次项作整体关系检验，报告分箱实际距离曲线；不择优选环或多项式次数。前沿轴不调整掉前沿位置。P95程序和niche对照采用同片同块非P95肿瘤Epi，对log1p上皮mean丰度和log1p(nCount)标准化最近匹配、不放回、固定seed20261103。第二套再要求相同前沿带且精确距离差≤200μm；缺精确距离不伪造匹配。匹配使用真实丰度和深度，不匹配程序/niche；保存逐对差异及实际覆盖。控制前沿的连续MP10模型加入标准化d_front及平方项。

统计族预设为：A前沿分子证据(全部57程序)、B0中心程序、B邻域/P95距离程序(全部程序与预设层共同校正)、C MP10与前沿、D niche与前沿(九类)、D MP10–niche控制前沿(九类)。原始、调整、匹配和去重叠/Top100敏感性标记分开。所有结果保留。

逐片模型使用标准化连续响应/预测量，二分类比较保留0/1预测量；niche和P95比例效应保留概率差单位。组成调整为Epithelial fraction和log1p(nCount_RNA)，中心邻域模型另控制相应环平均Epithelial fraction；不把所有细胞比例入模。常数协变量移除，模型需≥20点。空间不确定性使用组织块内固定方格聚合回归影响量的cluster-robust sandwich、有限簇修正及t区间，前沿500μm、重叠邻域2000μm；至少8簇才给片内P/CI。块间近似独立是限制，尤其远环。二次模型联合Wald F检验整体距离关系。此法无置换最小P分辨率问题，不使用普通spot自由度。

主要生物学汇总先片内效应，再已有analysis_unit内等权，再dataset等权。本轮主汇总使用标准化回归系数/均值差及比例差直接等权；Pearson与残差相关保留在片内明细，不混作β。输出全部切片描述、有效单位、已核验患者三套层级，缺身份不制造ID。dataset效应零均值t检验(≥3datasets)，层级bootstrap1999次估计区间，BH覆盖完整预设族；患者方向及leave-one-dataset-out稳定性保留。来源正常/肿瘤分开，spot数不作权重。统计为探索性，不以组成调整或位置匹配宣称因果。

## 版本来源与核对
正式 `accepted_label` 直接来自full_v1，后者从v8 nodes连接SAMPLE_R12_LEIDEN/labels.parquet；0–8共126,075个标签全部落在冻结query中。标签缺失保留NA，绘图原色板对应旧显示niche1–9是标签+1，仅是名称映射。没有引入dataset-local/K12。本轮不重算q05资格；上游位置接口的analysis_eligible_legacy是既有posterior mean Epithelial/all30≥0.25，R12/common_anchor/E0用途不同，继续原样保存，不以其覆盖条件替换主query。两个代表片新图R1与原v8 R1在相同有效组织/组织块限制下完全相同。

文献依据：[Puram Cell 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5878932/)、[公开100基因原表](https://zenodo.org/records/3260758/files/pemt_signature.txt?download=1)、[Flores 2026](https://www.nature.com/articles/s41420-026-03202-y)。Flores原文已实际取得并核对十个pro-invasion hits，与本任务提供成员一致；其功能筛选证据不能自动转译为表达高低的侵袭强弱。

单细胞使用项目既有reference_v2原始counts和reference_celltype_v2，在每sample×subtype固定随机最多60细胞，共19,053细胞，Epi、CAF、SMC、Myeloid分别检查。完整基因宇宙评分，不训练新注释；按sample×subtype汇总，再在lineage图中按sample先汇总，不作空间样本配对。检测支持中top1_count_share为签名总counts中最高单基因的份额，不能声称是UCell秩贡献分解。

所有同类空间分数采用全队列98%分位显示上限，下限0；前沿距离显示±1000μm、P95距离0–1000μm，显示裁切不改变原值或统计。图中灰色为未显示/无覆盖，红圈分别表示P95或pEMT-high front，H&E橙线为原正式观测整体外缘。所有可分析切片均输出，不按显著性筛图。

### 实际覆盖与显示层级
全部63片完成，不按方向筛片。59片具有至少一个内侧肿瘤Epi P95参考点；GSE171351_D1仅一个候选，保留其P95(实际比例100%)，不把它作为有足够覆盖的推断结果。GSE238145_B13没有冻结query。3张邻近正常没有本定义下的肿瘤P95/前沿对象，但已保留全Epi连续中心分析和空间程序图。正式单位汇总不包含缺analysis_unit_uid的切片；这些仍完整保留在section及dataset描述中。C1同时输出全部Epi跨内外侧及肿瘤Epi内部版本。

曲线首先在片内对有效spots取均值，重复切片按已有单位等权，再dataset等权；同箱要求≥5目标点。统计表保留每箱spots/sections/units/verified patients/datasets。不同距离箱覆盖可变，因此曲线为观察分布，不因曲线视觉差异声称配对总体效应；正式匹配比较另以同一高分中心及其控制中心的对应环均分差作为片内效应。曲线阴影为dataset均值的标准误，不是层级bootstrap 95% CI；正式95% CI见效应表和森林图。

前沿对近内部、前沿对深内部两项预设距离带比较完整保留，与主前沿对全部内部比较属于同一A/C检验族。包含位置协变量的原始与调整比较使用相同有精确位置的可用集合，避免将覆盖变化当作调整效应。

主MP10是冻结Epi查询点上的mixed-spot EpiMP10_Top20_AUCell。已有inverse-c2l派生spot/epi模型层及epi_score_valid仅原样保留供溯源，未替换本轮主评分，未称为直接测得的上皮特异表达。
