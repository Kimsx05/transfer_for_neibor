from common import *
import numpy as np,pandas as pd

o=pd.read_parquet(ROOT/'cellcharter_assignments.parquet');old,z,xy=load();s=pd.read_csv(ROOT/'metadata/niche_summary.tsv',sep='\t');c=pd.read_csv(ROOT/'metadata/major_composition.tsv',sep='\t');m=pd.read_csv(ROOT/'metadata/representative_markers.tsv',sep='\t');sens=pd.read_csv(ROOT/'metadata/sensitivity.tsv',sep='\t');per=pd.read_csv(ROOT/'metadata/per_niche_sensitivity.tsv',sep='\t');graph=pd.read_csv(ROOT/'metadata/graph_summary.tsv',sep='\t');hop=pd.read_csv(ROOT/'metadata/hop_scale_summary.tsv',sep='\t');auto=pd.read_csv(ROOT/'metadata/autok_stability.tsv',sep='\t');selection=json.loads((ROOT/'metadata/K_selection.json').read_text());K=selection['K'];checks=[]
def check(name,ok):checks.append(dict(check=name,passed=bool(ok)));assert ok,name
check('all_retained_unique',len(o)==308549 and o.cell_id.is_unique and o.cell_id.equals(old.cell_id));check('working_labels_unchanged',o.working_major.equals(old.working_major) and o.working_subtype.equals(old.working_subtype));finite=np.isfinite(z).all(1)&np.isfinite(xy).all(1);check('primary_actual_matrix_mask',np.array_equal(o.PRIMARY_included,finite&~o.low_information));check('expanded_actual_matrix_mask',np.array_equal(o.EXPANDED_LOWINFO_included,finite));check('conflict_complete_exclusion',np.array_equal(o.NO_CURRENT_CONFLICT_included,o.PRIMARY_included&~o.review_strong_competition));disputed=o.review_strong_competition|o.morphology_identity_uncertain;check('interpretation_disputed',o.loc[disputed,'interpretation_major'].eq('Disputed').all() and o.loc[disputed,'interpretation_subtype'].eq('Unknown').all());check('spatial_unchanged',np.array_equal(o[['x_um','y_um']],xy))
for version in ['PRIMARY','EXPANDED_LOWINFO','NO_CURRENT_CONFLICT']:check(version+'_labels_and_exclusion',o.loc[o[version+'_included'],version+'_niche'].notna().all() and o.loc[~o[version+'_included'],version+'_niche'].isna().all())
check('all_primary_niches_reviewed',len(s[s.version.eq('PRIMARY')])==len(o[o.PRIMARY_included].groupby(['section_id','PRIMARY_niche'],observed=True)));pd.DataFrame(checks).to_csv(ROOT/'metadata/delivery_checks.tsv',sep='\t',index=False)
def table(df):
 d=df.reset_index();return '| '+' | '.join(map(str,d.columns))+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join(str(v).replace('|','/') for v in r)+' |' for r in d.itertuples(index=False,name=None))
interpret=[]
for _,r in s[s.version.eq('PRIMARY')].iterrows():
 cc=c[(c.version=='PRIMARY')&(c.section_id==r.section_id)&(c.niche==r.niche)].sort_values('n',ascending=False);mm=m[(m.section_id==r.section_id)&(m.niche==r.niche)].sort_values('difference_log2',ascending=False);flags=[]
 if r.Disputed_fraction>.15:flags.append('争议身份占比较高')
 if r.Unknown_fraction>.4:flags.append('未知身份占比较高')
 if r.median_genes<100:flags.append('深度偏低')
 if r.n<20:flags.append('成员极少，不作独立生物群解释')
 if mm.empty:flags.append('没有满足当前检出率/正表达差规则的代表基因，不等于无表达')
 other=s[(s.version=='PRIMARY')&(s.niche==r.niche)&(s.section_id!=r.section_id)];source=r.n/(r.n+other.n.sum())
 if source>.95:flags.append('主要由本来源贡献，技术/生物差异混杂')
 interpret.append(dict(section=r.section_id,niche=r.niche,n=int(r.n),section_pct=round(100*r.section_fraction,2),interpretable_pct=round(100*r.interpretable_coverage,2),composition_all='; '.join(f'{v.major} {v.fraction_all:.1%}' for v in cc.head(5).itertuples()),representative_measured_genes='; '.join(f'{v.gene}({v.detection:.0%},Δ{v.difference_log2:.2f})' for v in mm.head(6).itertuples()),median_genes=r.median_genes,Unknown_pct=round(r.Unknown_fraction*100,2),Disputed_pct=round(r.Disputed_fraction*100,2),review='；'.join(flags) or '组成/实测marker/位置共同审阅，暂不赋机制性名称'))
interpret=pd.DataFrame(interpret);interpret.to_csv(ROOT/'metadata/niche_interpretation.tsv',sep='\t',index=False);cover=o.groupby('section_id')[[v+'_included' for v in ['PRIMARY','EXPANDED_LOWINFO','NO_CURRENT_CONFLICT']]+['review_strong_competition','morphology_identity_uncertain']].sum();worst=per.sort_values('Jaccard').groupby(['version','section_id']).head(3);primarys=s[s.version=='PRIMARY'];kshared=primarys.groupby('niche').section_id.nunique();shared=kshared[kshared==2].index.tolist();unique=kshared[kshared==1].index.tolist()
text=f'''# CellCharter 两片联合探索 v1：待用户验收

本轮已完成 PRIMARY、EXPANDED_LOWINFO 和 NO_CURRENT_CONFLICT 的实际计算。只到空间簇发现、全谱系解释和两个敏感性版本；未计算 MP10/胆固醇、spottedpy、通讯、CNV 或正式 TLS。原标签不是人工确认，空间簇也不是机制性 niche 的证明。

## 输入与身份解释

全部保留308,549，输入对象、counts、分割、scVI和工作注释不修改。按cell_id对齐，实测latent20所有维度及空间坐标有限才有建模资格。三集合与给定预期完全一致：

{table(cover)}

PRIMARY=264,957，EXPANDED=302,271，NO_CURRENT_CONFLICT=241,000；6,278缺有效latent者未补零。PRIMARY不是旧UMAP/core集合，noncore不一刀切排除。37,314个有latent的low_info仅进入expanded。三集合都保留Unknown；模型输入未使用细胞身份标签。

当前强竞争23,957（1,947/22,010）全部标为Disputed，不是只处理其中2,402个仍有已知大类的子集。此子集中的1,038个继承细类解释也回退Unknown。8个形态身份疑问单独标Disputed（7个PRIMARY），不再新增拟合版本。历史artifact/mixed标志原样保存，不断言污染或doublet。

## 构图、尺度、孤立对象

每片分别以真实spatial_um进行SciPy Delaunay（三角剖分等价于Squidpy generic Delaunay底层），显式按section_id分组；聚合sample_key=section_id。逐版本断言跨切片边为0，重复坐标检查通过。预先声明PRIMARY每片99百分位与40µm取小值，后续两版沿用同一物理上限。

{table(graph.set_index(['version','section_id']).round(3))}

没有为了连接碎片而补边。主分析的最长100条保留边图已审阅：大片空白与分离组织块未被长连接跨越；距离裁剪无法保证识别所有小裂隙，仍是局限。Cervilla有多个大组织分量，D1主体更连续。孤立点保留拟合并标 `isolated_self_expression_only`：其邻居层无邻居，不能把分类称为完整邻域证据。

CellCharter官方aggregate_neighbors：n_layers=3、mean，中心20维+第1/2/3跳各20维，实际80维。以下为固定seed123每片300个中心的抽样，3-hop不是固定微米半径，也不是普通Visium ring：

{table(hop.set_index(['version','section_id','hop']).round(2))}

## K选择与计算方法

CellCharter0.3.8、Squidpy1.8.3；实际接口已核对。GPU复用现有CUDA PyTorch，8 CPU线程、GPU进程上限40%，batch16384，full covariance GMM，协方差正则1e-5，最多100epochs。AutoK基础seed123；选定PRIMARY模型的实际seed1514383052，两项敏感性均固定为同一实际seed。参数完整见 `metadata/parameters.json`。AutoK初始K2–12、max_runs5、官方稳定性收敛阈值0.01。该版本为评估相邻K稳定性还实际拟合边界邻居K1/13，非额外生物学候选。

{table(auto.set_index(['stage','K']).round(4))}

主K=**{K}**。初始内部峰={selection['peaks_initial']}；是否仅扩展一次到16：{selection['extended']}。优先内部稳定性峰，并检查最小簇≥0.5%以排查退化小簇，再比较官方稳定性；没有用身份标签、预期阳性或MP10选K。来源/深度/组织分布在下表逐簇解释；不是以完全混合为验收标准。如多个K接近，不把本K称为唯一真实数量。实际重复次数见runs列，独立运行种子见K_selection.json。保存的是选定K的最佳似然模型 `models/primary_gmm.pkl`，不是所有候选模型。

[AutoK曲线](figures/autok_stability.png)；[主分析两片并列](figures/PRIMARY_two_sections.png)。两片均出现的簇：{shared}；只在一片有成员：{unique}。出现不等于同等比例或同一机制。尤其Cervilla的niche7仅3个，不能据此宣称两片共享的可靠群体；niche10仅见D1。K11稳定性约0.777，K4约0.764，差距有限，当前K是探索性分辨率而非唯一最优真值。

## 每片每簇解释（分母包括Unknown/Disputed）

{table(interpret.set_index(['section','niche']))}

便于审阅的生物学描述如下，所有身份都属于工作解释：

| niche | Cervilla | D1_Johnson |
| --- | --- | --- |
| 0 | 成纤维占46.3%，COL3A1/COL1A1/TAGLN；小群790 | 成纤维占46.7%，COL1A1/COL1A2/LUM/TAGLN；仍含21.3% Unknown |
| 1 | 上皮占89.0%，S100P/PSCA/HPGD；主要来自本片 | 仅385个，上皮占88.3%，SPINK1/UPK1A/KRT7；敏感性较弱 |
| 2 | 上皮占99.8%，RAB3IP/MDM2/CPSF6 | 上皮占89.7%，S100A8/S100A9/S100P；与C片同编号不代表完全相同表达状态 |
| 3 | Ig/JCHAIN高，Unknown53.9%，并有成纤维/髓系/浆细胞 | Ig/JCHAIN高，Unknown54.1%，成纤维26.2%；两片均不宜命名纯浆细胞区或TLS |
| 4 | COL1A1/COL3A1/COL1A2，基质相关混合区；Unknown38.3% | SPARC/COL4A1/VIM，基质与血管相关成分；Unknown42.5% |
| 5 | 仅596个，Unknown66.4%，中位87genes，免疫marker伴低深度 | 39,319个，Unknown57.8%，TRAC/TRBC2/IGHM支持免疫混合成分；来源偏倚明显 |
| 6 | 上皮66.0%、成纤维9.8%，MDM2/RAB3IP与胶原共同支持混合邻域 | 仅944个，Unknown46.9%，KRT5/KRT14与胶原；不稳定小群 |
| 7 | 仅3个，不作生物群解释 | 上皮40.7%，KRT14/KRT5/KRT6A/KRT17；Disputed32.1%，不能视为纯上皮区 |
| 8 | 上皮81.3%，KRT13/KRT19；中位78genes，空间外围广泛分布 | Unknown60.6%，中位52genes；没有满足筛选规则的正向代表基因，深度解释优先 |
| 9 | Ig及胶原共同表达，含上皮/成纤维/未知身份 | 成纤维33.5%、Disputed31.8%，Ig及胶原；争议身份影响较大 |
| 10 | 无成员 | 上皮57.4%，S100A8/S100A9/KRT16/KRT5；Disputed18.9%，D1特有且敏感 |

空间图显示Cervilla上皮相关niche1/2/6占据不同连续区域，niche8多分布外围；D1的niche2主要位于图中下部，niche5偏上部/右部，niche7/10在主体中部互相邻接。上述位置是聚类图观察，尚无病理分区证实，不能据此确定侵袭边界或机制。niche8具有空间连续性，也不能排除深度驱动。全体逐片逐簇的数值证据见上表和marker文件。

本轮统一称niche0…niche{K-1}，表示中心及周围表达共同定义的空间簇。依据组成、实测marker和空间位置作描述，未强行命名TLS/GC/机制。全谱系均纳入：major_composition.tsv含全体建模分母与可解释大类分母两列；同时有interpretable_coverage，避免隐藏Unknown/Disputed。subtype_parent_composition.tsv明确沿用细类或父类，父类不冒充细类。

marker来自未改动v2整数counts，合并同名基因、log1p10k，逐片相对于其他PRIMARY空间簇计算平均log表达差和检出率。未使用插补、P值筛选或独立验证措辞；仅两张切片，细胞不等于患者重复。`representative_markers.tsv`有限代表基因，`measured_marker_evidence.parquet`同时保留全谱系已知marker。原feature缺失（如TPSAB1/TPSB2、KRT18等）不同于检测为零，缺失明细沿用上轮 `../joint_umap_major_review/metadata/review_panels.json`。

每片质量图分别显示Unknown/Disputed、n_genes、n_counts与强竞争，同原物理坐标；深度采用对数颜色。以上review列指出深度、争议身份、未知覆盖或来源偏倚，不自动删掉这些簇。

## 固定K的两个敏感性版本

两个版本均真正重新Delaunay构图、按主分析物理阈值裁剪、3层聚合、固定K重新拟合，实际seed1514383052与选定PRIMARY模型一致。早期seed123敏感性试算已被此统一种子的结果替换；以下为最终数值。相同种子约束初始化随机性，但不同样本集合仍会改变初始化与优化路径，不能把差异全归于生物学。不重新AutoK，也不是只投射新增细胞。

{table(sens.set_index(['version','section_id']).round(4))}

ARI及交叉表只在每版与PRIMARY共有ID上计算，排除者不当不一致。显示编号在两片联合交叉表上用Hungarian匹配一次，映射文件保留，每片不独立改名；raw_cluster仍保留。NO_CURRENT_CONFLICT不是无污染真值。

两版比较中的具体结构变化（百分比均以共有ID为分母）：

- 较稳定：两片niche3的Jaccard分别为0.909/0.837（expanded）与0.865/0.860（去冲突）；D1 niche5为0.893/0.911。Cervilla niche8为0.892/0.948，但稳定性不消除其深度偏倚。D1 niche2两版分别有99.3%/99.4%成员保留在对应簇。
- 纳入low_info后，Cervilla niche1保留63.3%，另有21.8%进入显示niche10、14.2%进入niche2；niche6的65.8%进入niche1、23.4%进入niche0，仅10.1%保持匹配niche6，显示混合邻域边界重组。D1 niche6有92.6%保留，但该匹配簇共同成员从944增至12,880，Jaccard仅0.067，属于显著扩张而非稳定边界。
- 去冲突后，Cervilla niche2有97.0%共同成员进入显示为niche6的新簇；旧niche1的50.5%及旧niche6的61.0%也进入该簇。旧niche2内部大部分保持共同归属，但与其他上皮/混合邻域的区分消失，属于合并，不能仅凭对角Jaccard接近0就说全部结构消失。D1 niche6有96.3%进入niche0，同样是小群并入更大区域。
- D1 niche10两版中保留到同一匹配簇的比例为65.1%/59.4%，Jaccard0.648/0.462；D1 niche7为89.1%/78.8%，Jaccard0.801/0.683。niche7比niche10更稳定，但高Disputed比例仍限制其身份解释。
- Cervilla niche7原本只有3个，去冲突后共有2个；不对这些极小分母的稳定性作生物解释。Cervilla niche0/4也会合并或吸收其他成员，需结合组成变化审阅。

[逐簇转移图](figures/sensitivity_transitions.png)同时显示非对角的大块归属；最大去向比例只说明共同归属，不能代替Jaccard证明簇边界稳定。变化最明显的簇（完整逐簇表见per_niche_sensitivity.tsv）：

{table(worst.set_index(['version','section_id','niche']).round(4))}

[纳入low_info空间图](figures/EXPANDED_LOWINFO_two_sections.png)；[去当前强竞争空间图](figures/NO_CURRENT_CONFLICT_two_sections.png)。请结合逐簇retention/Jaccard、两版组成/规模和hop范围比较，整体ARI不能掩盖小群变化。低信息增点会改变邻接和聚合；去强竞争也会改变局部密度，因此改变不等于某版一定错误或正确。结果不稳定处保留待讨论，没有新增筛选版本。

## 是否进入下游

目前值得把这些结果作为下一步研究设计和核对材料，但**尚不宜直接作为确认性MP10/共定位结论的固定niche**。先验收空间结构及敏感性较稳定的簇，明确是否把Disputed/孤立点仅用于建模而限制解释；对高Unknown、强来源或逐簇Jaccard偏低的结构尤其保留不确定性。两片无独立患者重复，marker是同数据描述性支持。B/Plasma聚集不自动等于TLS，GC B未标出不证明组织无GC。任何后续评分/共定位均等待用户明确授权。

## 文件、复现和检查

主入口 `cellcharter_assignments.parquet` 覆盖全部308,549个ID，各版未纳入者niche为空且有排除原因；工作标签、解释标签、历史/current质量字段均保留。`code/prepare.py → code/fit.py → code/summarize.py → code/refine_marker_plots.py → code/review_sensitivity.py → code/report.py`；每步用 `code/run.sh code/脚本名.py`（以上路径相对于本结果目录）执行；运行科学流程用 `code/run.sh` 复用旧GPU环境和本轮隔离软件依赖。主聚合表示只保存一份 `metadata/primary_aggregate.npy`（约85MB），其余邻接和聚合只在内存。图件为透明PNG、真实µm尺度；没有新h5ad或全体两两距离。

必要交付检查{len(checks)}项通过（delivery_checks.tsv）。计数、ID、纳入掩码、原工作标签/空间坐标、完整强竞争排除、解释层级及每版标签空值核对通过；运行中验证所有版本跨片边为0和80维聚合有限。

官方依据：[教程](https://cellcharter.readthedocs.io/en/stable/notebooks/cosmx_human_nsclc.html)、[聚合API](https://cellcharter.readthedocs.io/en/stable/generated/cellcharter.gr.aggregate_neighbors.html)、[AutoK API](https://cellcharter.readthedocs.io/en/stable/generated/cellcharter.tl.ClusterAutoK.html)。安装环境版本可能较网页实现更新，实际参数以本地代码和日志为准。

本轮已停在CellCharter空间簇及敏感性结果验收，等待用户审阅。
'''
(ROOT/'CELLCHARTER_REVIEW.md').write_text(text);print('REVIEW_COMPLETE',K,len(checks));print(sens.to_string(index=False))
