exec(open(__file__.replace('07_finalize.py','01_prepare.py')).read().split('# Archive contents')[0])
import html,base64,datetime
D=B/'04_summary';now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(timespec='seconds')
read=lambda p:pd.read_csv(B/p,sep='\t')
models=read('03_clinical/model_results.tsv');lines=read('02_cell_lines/overall_effects.tsv');deltas=read('02_cell_lines/line_deltas.tsv');epi=read('01_epi/sample_effects.tsv');patients=read('03_clinical/patient_inclusion_endpoint_audit.tsv');coverage=pd.concat([read('00_audit/'+g+'_matrix_audit.tsv') for g in ['GSE171537','GSE169455','GSE212810']]);resources=read('00_audit/resource_usage.tsv');validation=read('00_audit/validation_checks.tsv');sc=read('03_clinical/score_standardization.tsv')
summary=[dict(module='01 Epi',status='PARTIAL / patient inference BLOCKED',unit='36 eligible tumor specimens; independent patient N unknown',finding='UP median within-specimen rho=0.2516; DN=0.3187, DN opposite expected direction',evidence='transcriptional specimen association only; deterministic average-tie AUCell-like',limit='No confirmed patient crosswalk; S.Score/G2M.Score unavailable'),dict(module='02 GSE171537',status='COMPLETED',unit='5 cell-line backgrounds,60 RNA measurements',finding='CisR2 positive/3 negative; mean delta0.0433,95% bootstrap CI[-0.2083,0.2949],sign-test P1',evidence='cell-model association; no consistent elevation',limit='3 passages per derivative,not independent inductions;5 backgrounds'),dict(module='03 GSE212810',status='COMPLETED with endpoint missingness',unit='39 matched;4 GCb excluded;35 cisplatin-coded;28 endpoint-evaluable',finding='nonresponse OR0.5756 [0.2172,1.5255],P0.2667;7 endpoint unknown',evidence='pretreatment pathological-response association',limit='small N; observational; author GP code interpreted in cisplatin context; no treatment-benefit prediction'),dict(module='03 GSE169455',status='PARTIAL / strict cisplatin BLOCKED',unit='149 expression;125 NAC +24 induction',finding='125 NAC regimen-unverified sensitivity OR0.7674 [0.5335,1.1040],P0.1537',evidence='NAC pathological-response sensitivity association only',limit='No individual regimen/cycle data; cannot designate strict cisplatin main'),dict(module='04 integration',status='COMPLETED for available evidence',unit='modules retained separately',finding='No convergent evidence for a general MP10 cisplatin-resistance direction',evidence='exploratory conclusion',limit='No causality,cholesterol flux or predictive-treatment-benefit claim')]
save(pd.DataFrame(summary),'04_summary/evidence_summary.tsv')
status=dict(status='EXECUTED_WITH_DOCUMENTED_INPUT_BLOCKERS',updated=now,configuration_sha256=(B/'config/frozen_v1.sha256').read_text().strip(),validation_checks=len(validation),validation_all_pass=bool(validation['pass'].all()),core_tables_replay_identical=bool(read('00_audit/reproducibility.tsv').identical.all()),peak_RSS_GiB=float(resources.peak_RSS_KiB.max()/1024**2),scheduler_job_id=None,report=str(D/'中文结果报告.html'),blockers=str(B/'00_audit/blockers.tsv'))
(D/'status.json').write_text(json.dumps(status,ensure_ascii=False,indent=2))
def table(x):return x.to_html(index=False,border=0,na_rep='NA',float_format=lambda x:f'{x:.4g}',escape=True)
def fig(name,caption):
 data=base64.b64encode((B/'figures'/f'{name}.png').read_bytes()).decode();return f'<figure><img src="data:image/png;base64,{data}"/><figcaption>{caption} <a href="../figures/{name}.pdf">矢量 PDF</a></figcaption></figure>'
def link(path,text=None):return f'<a href="../{path}">{html.escape(text or path)}</a>'
raw=epi[(epi['subset']=='tumor')&(epi.adjustment=='unadjusted')];es=raw.groupby('signature',sort=False).agg(n_specimens=('rho','count'),median_rho=('rho','median'),n_positive=('rho',lambda x:(x>0).sum()),n_negative=('rho',lambda x:(x<0).sum())).reset_index()
clinical_cols=['cohort','scope','endpoint','adjustment','n_enrolled','n','events','missing','OR','CI_low','CI_high','P','FDR','method']
body=f'''<h1>MP10 与含铂治疗反应：首轮计算验证</h1>
<p>生成时间：{now} (Asia/Shanghai)。分析配置于 <b>2026-11-09 21:44:39 +08:00</b> 冻结，早于任何评分关联分析。结果目录：<code>{B}</code>。</p>
<div class="key"><p><b>当前结论：可执行分析没有形成“MP10 普遍对应顺铂耐药”的一致证据。</b>五种细胞系的 CisR 效应为2升3降；GSE212810 病理主终点 OR 约0.58，95% CI跨1；GSE169455 方案未核实的 NAC 补充结果也未显示较差反应关联。Epi 的 UP 和 DN 多数同为正相关，其中 DN 与预期耐药方向相反。阴性或不一致结果不等于证明不存在作用。</p>
<p>这不是全部主分析均已完成的声明：<b>Epi 患者级独立推断、指定细胞周期调整、GSE169455 严格顺铂主分析仍有明确输入阻塞</b>。它们保留 NA/BLOCKED，没有以细胞数、推导患者编号或未知治疗方案补足。</p></div>
<h2>1. 实际完成范围与证据层级</h2>{table(pd.DataFrame(summary))}
<p>没有重新运行 NMF，没有修改 MP10 Top20、优化阈值、训练分类器、预测 IC50、增加耐药 signature 或用生存替代疗效。外部队列分别预处理、评分和建模；未合并表达矩阵或做跨队列 ComBat，也没有用两个不可比人群合并 OR。</p>
<h2>2. 输入、评分与预先冻结的规则</h2>
<p>固定 Top20 (按原顺序)：<code>{', '.join(cfg['genes']['MP10'])}</code>。MP10可用基因至少16/20 (80%) 是本轮操作阈值，不是公认生物学界值。完整列表与参数见 {link('config/frozen_v1.json')}；SHA-256：<code>{status['configuration_sha256']}</code>。</p>
{table(coverage[['cohort','n_original','n_exact_or_alias','n_background','n_present','missing','preprocessing']])}
<p>GSE171537和GSE212810原始counts按基因符号合并后，采用edgeR TMM、log2 CPM (prior.count=1)。GSE169455已经过作者RMA、分位数归一化、批次校正和探针中位数合并，原值范围0.662–13.542；没有再次log或重新批次校正。它缺少IDI1、ACAT2，但18/20达到阈值。三队列均用GSVA ssGSEA (alpha=0.25，normalize=TRUE)，使用完整提供的表达基因背景；不是HVG或仅20个目标基因。ssGSEA归一化并不保证数值在0–1，各队列绝对值也不直接比较。</p>
<p>映射先保留作者/注释库精确SYMBOL，再采用org.Hs.eg.db唯一ALIAS；重复counts取和，重复已归一化表达取中位数；未映射背景保留原ID。完整逐基因映射、矩阵覆盖、检测率与软件注释版本均已保存。细胞表达为0没有被误作“基因不在矩阵”，也没有逐患者改变signature。</p>
<p>Epi复用已核实的Top20 <b>平均秩AUCell-like</b>：原RNA counts、90011基因、K=4501、ties取平均、归一化分母为Σ(K−j)。新增集合采用同一实现。它不是官方AUCell随机ties实现；原项目比较表已显示二者不同，本报告明确这个方法偏离。全部123699细胞重算Top20的最大误差为3.33×10⁻¹⁶。Epi和bulk ssGSEA的绝对值、阈值不能互换。</p>
<p>Tsunoda采用官方MSigDB现时快照及文件hash锁定，UP15、DN50映射基因 (DN网页53个原始accession映射50个symbol)。原文两张表核实了集合方向，以及C1S在T24和KK47模型中相反；没有机械按网页倒置的Exact source交换集合。INSIG1是MP10与UP共有成员。没有把UP−DN包装成临床耐药指数。</p>
<h2>3. Epi：完成样本内转录关联，患者推断受阻</h2>
<p>原分析对象为123699个Epi/既有上皮谱系细胞、57个样本、5个dataset。保留Epi命名，不将全体Epi称为恶性细胞。肿瘤来源42个样本101263细胞，按每样本≥50细胞纳入36个样本101198细胞；6个样本因2–21个细胞排除。正常来源15个样本22436细胞，13个样本22418细胞可描述。既有CNV支持的malignant_like_main子集仅作敏感性：23个合格肿瘤样本66584细胞；没有新推断CNV。</p>
{table(es)}
<p>上表是<b>样本内Spearman相关的样本分布</b>，不是患者等权总体效应。原始UP中位数0.2516，36/36正相关；原始DN中位数0.3187，35/36正相关、1/36负相关，DN与预期负方向相反。去INSIG1后UP中位数降至0.1384；同时去INSIG1/C1S为0.1171 (33正、3负)。去C1S对DN几乎无影响。CNV支持子集中UP/DN中位数为0.2687/0.2972，仍未形成UP正、DN负的一致耐药方向。</p>
<p>样本映射中37项为derived，20项为provisional，没有confirmed/verified患者ID。公开ENA确认PRJNA662018有患者编号，但缺少原始编号到本地重命名样本的已验证交叉表；正常组织也不能强行配对。因此没有把orig.ident当患者，患者级Fisher-z汇总、整患者bootstrap、双侧符号检验及UP/DN主FDR均保留NA。计算入口已实现这些方法并在身份条件不满足时拒绝推断。</p>
<p>nCount_RNA/nFeature_RNA可用，但S.Score/G2M.Score缺失，所以指定rank-residual调整和替换深度指标敏感性未完成；没有以depth-only模型冒充。正相关可能同时涉及测序检测深度、表达可检出性、细胞状态和真实生物学，不能据未完成的调整宣称其为假象，也不能凭这些分数标记“已证实耐药细胞”。</p>
{fig('A_specimen_descriptive_original','图A的可完成部分：每点一个肿瘤Epi样本内效应；不是患者效应。患者级图本轮因身份阻塞无法生成。')}
{fig('A_specimen_sensitivity','全部预设基因去除版本与既有CNV支持子集，保留所有可估计样本。没有患者级P值。')}
<h2>4. 独立耐药细胞模型 GSE171537</h2>
<p>核实60个GSM/测定：253J、5637、KU1919、T24、TCCSUP各有parental、CisR、GemR、GemCisR及3个RNA采样。原文说明RNA来自三个不同传代；不能作为三次独立耐药诱导。GSM培养说明为取样前数代不加抗生素、顺铂或吉西他滨，采样为Untreated；未记录可独立拆分的测序批次。MP10发现数据的五个本地dataset不包含此队列，但无法证明未记录的历史筛选行为不存在。</p>
{table(deltas[['cell_line','comparison','parental','resistant_mean','delta']])}
{table(lines)}
<p>主要CisR−parental的五个背景效应分别为：253J −0.0345、5637 +0.3642、KU1919 −0.1357、T24 −0.3677、TCCSUP +0.3903。均值+0.0433，但中位数−0.0345，说明均值方向不能代表普遍升高。以完整细胞系背景等权bootstrap10000次，均值95%区间[−0.2083,0.2949]；双侧精确符号检验P=1。n=5下最小可能双侧P为0.0625，不能以P&lt;0.05作成功门槛。</p>
<p>GemCisR也是2升3降，GemR为1升4降；三项比较共享亲本，并不独立。三项符号检验一起BH后的FDR均为1。差值是ssGSEA分数变化，不是基因表达fold change；n=5的bootstrap区间仅作探索描述，不能精确刻画所有BLCA模型的异质性。</p>
{fig('B_cell_line_RNA_replicates','每个圆点为原始RNA传代采样，菱形为line×state均值；连接线只用于比较状态均值，不表示耐药演化轨迹。')}
{fig('B_cell_line_deltas','五个独立细胞系背景的全部预设差值，0线表示与亲本相同。')}
<h2>5. 临床：病理三态编码与独立建模</h2>
<p>共同主终点为≤ypT1且ypN0；nonresponse=1表示未达到。严格pCR=ypT0N0为次要终点，non_pCR=1表示未达pCR。本项目预先固定这一顺序以便跨队列比较，即使GSE169455原文强调pCR也不交换主次。OR均为MP10每增加一个固定队列SD的效应，OR&gt;1才是较差病理反应优势比升高。</p>
<h3>GSE212810</h3>
<p>39个表达病例全部通过GSM、Case_#和作者患者编号连接到Table S6A；原标签21R/18NR完整保留，未用更大临床总体替代样本数。作者编码30 GP、5 ddMVAC、4 GCb；按原研究顺铂方案语境将GP归为gemcitabine/cisplatin，GCb为卡铂并严格排除，原始编码保留以便审核；未推断具体给药剂量。主入组35例，其中良好病理反应13、未反应15、未知7，实际主模型n=28。</p>
<p>作者S6A中777表示“无残余肿瘤”，0表示“非浸润或CIS”，999表示淋巴结未评估。7例ypN未知的作者R病例没有被静默算成良好反应；其中1例非浸润/CIS仍可明确排除pCR，因此pCR可评估n=29 (11 pCR、18非pCR)，未知6。原文NR句子中的ypN0与“淋巴结转移”矛盾，采用可追溯逐例病理，不机械照抄。两个作者表中2周期病例与GEO“至少3周期”总体说明不一致；主分析保留，≥3周期单独敏感性。</p>
<p><b>共同主终点：OR=0.576，95% CI 0.217–1.526，P=0.267，预设临床家族保守BH(n=2) FDR=0.533。</b>未观察到较差病理反应方向；区间很宽，不能宣称保护作用。cT调整OR=0.508 (0.154–1.668)。含GCb全NAC可评估32例OR=0.689 (0.300–1.580)；≥3周期严格组n=26，OR=0.560 (0.207–1.512)。原作者R/NR全39例的独立补充分析OR=0.562 (0.248–1.274)，不替代病理三态终点。</p>
<p>pCR次要模型OR=0.530 (0.177–1.586)，P=0.256。cT调整模型触发稀疏/分离，按预设规则使用Firth的Jeffreys惩罚似然和profile区间：OR=0.687 (0.187–1.447)，P=0.371。估计方法没有按P值选择；本地实现已通过已知完全分离2×2数据的解析Firth OR=441数值核验，该单元测试从未混入科学分析数据。</p>
<h3>GSE169455</h3>
<p>149个表达病例中125 NAC、24 induction，与GEO逐例字段一致。表达为治疗前TURB肿瘤；术后pTN在治疗后情境中重建ypTN。125 NAC均可判定共同主终点，54良好反应、71未反应；42 pCR、83非pCR。induction独立描述：23例共同终点可判定，1例pTCisNx为未知；该例可以排除pCR，故induction的pCR n=24。病理解析覆盖IS/CIS、a/b后缀及Nx，相关修正有留痕。</p>
<p><b>严格顺铂主分析仍BLOCKED：</b>GEO不含逐例方案和周期；出版社/preprint补充材料多次403/401，后续原作者开放补充表为基因统计，未给出可与GSM连接的方案表。系列总体“cisplatin-based”描述没有替代逐例审计。作为明确降级的125例NAC方案未核实补充分析，OR=0.767 (0.533–1.104)，P=0.154；cT调整0.811 (0.553–1.191)，亚型调整0.801 (0.542–1.185)。这些结果不能标成严格顺铂独立验证。</p>
<p>pCR补充分析OR=0.776 (0.531–1.134)，P=0.191；induction共同终点OR=1.133 (0.426–3.011)，与NAC点估计方向不同但区间宽。没有按最显著亚型筛选结果，也没有进行足以断言亚型交互的分层推断。方案、周期、实际肿瘤含量调整受数据限制；微切割/宏切割不能视为100%肿瘤纯度，没有做“评分除以纯度”。</p>
<h3>固定标准化与全部模型</h3>{table(sc)}
<p>每队列只使用一次固定均值和样本SD。GSE212810在35例顺铂方案入组者中标准化，同一z用于28/29例不同终点及含GCb补充模型；GSE169455在125例NAC中标准化，induction沿用此参考尺度。没有按R/NR分别标准化。所有P双侧，主模型为未调整logistic；下表保留实际n、事件、缺失和估计方法。</p>
{table(models[clinical_cols])}
<p>临床主FDR家族预设两队列：仅GSE212810严格分析可估计，保留总家族大小2进行保守BH，GSE169455缺失项保持NA；不将方案未核实敏感性偷偷换成主检验。pCR另列次要家族；敏感性不作显著性菜单。没有合并两队列OR。</p>
{fig('C_GSE212810_nonresponse','严格顺铂方案入组35例，每患者一个点；未知7例单独显示。P来自28例可评估者的连续logistic模型，不是图中三组检验。')}
{fig('C_GSE169455_nonresponse','125例NAC方案未核实的补充分布，不能称为严格顺铂主验证。P来自连续logistic模型。')}
{fig('D_clinical_OR_nonresponse','共同主终点及预设调整/敏感性分别标注；严格GSE169455主结果保留BLOCKED。')}
{fig('D_clinical_OR_non_pCR','严格pCR为独立次要终点；三角点使用触发后的Firth profile推断。')}
<h2>6. 不一致、解释边界与下一步</h2>
<p>跨层级不一致是本轮最重要的发现：Epi中UP正相关，但DN同样正相关；去共享INSIG1后UP明显衰减。细胞系升降混合，5637/TCCSUP升高而T24/KU1919/253J下降。临床两个可分析NAC人群点估计偏向OR&lt;1，induction点估计偏向OR&gt;1，均无精确证据；不是“独立临床验证成功”。患者间异质性目前无法独立估计，样本差异不能伪称患者差异。</p>
<p>MP10涉及甲羟戊酸/胆固醇相关转录过程，但表达不直接等于胆固醇含量或通量，不能证明胆固醇本身导致耐药。单臂NAC病理反应关联不能等同获得性耐药；没有治疗比较及交互分析，不能宣称预测顺铂治疗获益。临床bulk还混有组成因素，Epi统计层级与bulk也不同。</p>
<ol><li><b>优先补齐低成本关键审计资料：</b>Epi样本到确认患者ID的交叉表及既有细胞周期分数；GSE169455逐例方案/混合换药/周期表。补齐后仍按当前冻结主终点与Top20重算，不改方向和阈值。</li><li><b>暂降低“普适顺铂耐药MP10”的机制投入优先级。</b>当前不能以大规模TME/空间关联扩展来补救不一致的耐药证据。</li><li>若继续功能验证，优先在升高的5637或TCCSUP和下降的T24/KU1919各选背景，测真实顺铂剂量反应并设计独立诱导/培养重复；结合通路干预、救援及胆固醇/甲羟戊酸通量测量，检验背景依赖而非预设普遍方向。此为下一步建议，本轮未实施。</li></ol>
<h2>7. 阻塞项、可复现性与交付入口</h2>{table(read('00_audit/blockers.tsv'))}
<p>GSE87304仅审计系列与样本：305例pre-NAC/TURBT说明可核实，但抽查GSM无逐例方案或术后ypTN，没有获得可追溯患者终点表，故未满足条件扩展纳入；未下载6GB CEL，更未用OS代替疗效。</p>
<p><b>{len(validation)}项风险相关检查全部通过；7张核心表完整重跑hash一致。</b>全流程各步骤退出码0，最高实测RSS为{status['peak_RSS_GiB']:.2f} GiB，最多2个data.table线程、BLAS/OMP1线程，重跑入口设12GiB虚拟内存上限，远低于48CPU/200GB上限。没有加载完整Seurat、复制大对象或构建全细胞排名矩阵。无需大任务或调度器，未提交任何job。用户23:00–07:00 Asia/Shanghai限制用于大任务；本轮完成的是低资源小任务。</p>
{table(resources)}
<p>图件为白底300dpi PNG及矢量PDF、真实Arial约12pt、线宽≥1pt。有图例的图同时给出panel_only/legend_only PDF。配对与森林图无现成技能模板，采用统一主题和导出器单独实现。图A患者级版未伪造，以准确标注的样本描述图提供可完成部分。报告内图片嵌入，可离线独立阅读。</p>
<ul><li>{link('scripts/run_all.sh','完整重跑入口')}；执行：<code>bash '{B}/scripts/run_all.sh'</code>。</li><li>{link('scripts/download_sources.py','冻结公共输入核验/缺失时恢复')}；上游hash改变将停止，不自动覆盖。</li><li>{link('00_audit/source_manifest.tsv','公共来源与hash')}；{link('00_audit/local_input_manifest.tsv','本地输入及原对象路径')}；{link('00_audit/implementation_amendments.json','解析/Firth修订记录')}。</li><li>{link('01_epi/cell_scores.tsv.gz','Epi细胞分数')}；{link('01_epi/sample_effects.tsv','全部样本原始/敏感性效应')}；{link('01_epi/patient_summary.tsv','患者级BLOCKED/NA')}。</li><li>{link('02_cell_lines/sample_mapping.tsv','细胞系GSM映射')}；{link('02_cell_lines/line_state_scores.tsv','line×state分数')}；{link('02_cell_lines/line_deltas.tsv','配对差值')}。</li><li>{link('03_clinical/patient_inclusion_endpoint_audit.tsv','临床逐患者入组/终点审计')}；{link('03_clinical/model_results.tsv','全部临床模型')}；{link('03_clinical/model_patient_usage.tsv','模型实际患者使用清单')}。</li><li>{link('04_summary/evidence_summary.tsv','跨模块证据汇总')}；{link('04_summary/figure_index.tsv','图件索引')}；{link('00_audit/validation_checks.tsv','验收检查')}；{link('00_audit/reproducibility.tsv','重跑hash比对')}；{link('logs/sessionInfo.txt','R sessionInfo')}。</li></ul>
<h2>8. 原始来源</h2>
<ul><li><a href="https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/TSUNODA_CISPLATIN_RESISTANCE_UP.html">MSigDB UP</a>、<a href="https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/TSUNODA_CISPLATIN_RESISTANCE_DN.html">DN</a>；<a href="https://www.nature.com/articles/1208313/tables/1">Tsunoda Table1</a>、<a href="https://www.nature.com/articles/1208313/tables/2">Table2</a>。</li><li><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171537">GSE171537</a>、<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC11094426/">原论文与RNA传代说明</a>。</li><li><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169455">GSE169455</a>、<a href="https://doi.org/10.1016/j.eururo.2021.10.035">Sjödahl等原研究</a>。</li><li><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212810">GSE212810</a>、<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC10591055/">原论文及Table S6A</a>，本地原始附件为00_audit/sources/mmc9.xlsx。</li></ul>
'''
style='body{max-width:1280px;margin:30px auto;padding:0 24px;font:16px/1.7 Arial,"Noto Sans CJK SC",sans-serif;color:#20242a;background:white}h1{font-size:28px}h2{margin-top:2em;border-bottom:1px solid #ddd}h3{margin-top:1.5em}.key{padding:14px 22px;background:#f4f6fa;border-left:5px solid #377eb8}table{display:block;overflow-x:auto;border-collapse:collapse;font-size:13px;line-height:1.4;margin:16px 0}th,td{border:1px solid #ddd;padding:7px;max-width:470px;overflow-wrap:anywhere}th{background:#eef1f5}code{overflow-wrap:anywhere;font-size:13px}a{color:#245b9a}figure{margin:22px 0}img{width:100%;height:auto}figcaption{font-size:14px;color:#444}'
report=D/'中文结果报告.html';report.write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>MP10首轮含铂反应验证</title><style>'+style+'</style><body>'+body+'</body></html>')
# Plain-text synopsis and entry point, no Page writes or markdown persistence required.
(B/'README.txt').write_text(f'''MP10 resistance validation v1
中文完整报告: {report}
状态: 已执行全部当前可行分析，患者映射/细胞周期/部分逐例方案仍BLOCKED。详见00_audit/blockers.tsv。
完整重跑: bash '{B}/scripts/run_all.sh'
核验公共输入: python '{B}/scripts/download_sources.py'
配置已冻结且有SHA256；不重新运行00_freeze.py，不修改原对象。
结果解释: Epi样本转录关联、细胞模型关联、临床病理反应关联分开；没有一致支持普适顺铂耐药。
所有统计/图表均来自真实数据；Firth数值单元测试单独存在scripts/firth.R，不用于科学结果。
''')
# Code identities and reviewable artifact index (exclude raw downloaded matrices from duplicate bundles).
files=[]
for p in sorted(B.rglob('*')):
 if p.is_file() and ('sources' not in p.parts) and p.name not in ['artifact_manifest.tsv']:
  files.append(dict(path=str(p.relative_to(B)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest() if p.stat().st_size<100_000_000 else 'not_rehashed_large'))
save(pd.DataFrame(files),'04_summary/artifact_manifest.tsv');print('REPORT',report,'status',status['status'])
