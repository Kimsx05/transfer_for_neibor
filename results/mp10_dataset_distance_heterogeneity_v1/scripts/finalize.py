from common import *
import json,shutil

def tab(d,cols):
 d=d[cols].copy()
 def f(x):
  if isinstance(x,(float,np.floating)):return 'NA' if not np.isfinite(x) else f'{x:.4g}'
  return str(x).replace('|',' / ')
 return '| '+' | '.join(cols)+' |\n| '+' | '.join(['---']*len(cols))+' |\n'+'\n'.join('| '+' | '.join(f(x) for x in row)+' |' for row in d.itertuples(index=False,name=None))+'\n'
def main():
 reg=read('endpoints.tsv.gz');mapping=read('canonical_programs.tsv');disc=read('discussion_relations.tsv');sens=read('discussion_sensitivity.tsv.gz');h=read('heterogeneity_configurations.tsv.gz');hs=read('heterogeneity_input_scope.tsv');co=read('representative_coordinates.tsv');loads=read('representative_loadings.tsv');st=read('stability_summary.tsv');rs=read('representative_method_range_sensitivity.tsv');eco=read('ecology_associations.tsv.gz');loo=read('ecology_leave_dataset_out.tsv');es=read('ecology_sensitivity.tsv.gz');comb=read('combinations.tsv').set_index('combination');ds=read('dataset_endpoint_summary.tsv.gz');meta=read('unit_metadata.tsv');scope=read('search_scope.tsv');bg=read('unit_background.tsv.gz');ur='D119'
 coverage=read('input_coverage.tsv');coverage=coverage[~coverage.module.str.startswith('background::')]
 br=[]
 for name,z in bg.groupby('input_family'):
  v=z[z.value.notna()];known=v[~v.record_id.str.startswith('section_only::')];br.append(dict(module='background::'+name,statistic='observed_fraction_or_score_summary',score_method='existing_small_summary',range_type='as_named',signature='separate_full_or_noMP10_as_named',datasets='|'.join(sorted(v.dataset_id.unique())),n_datasets=v.dataset_id.nunique(),n_sections=np.nan,n_analysis_units=known.record_id.nunique(),n_programs=v.feature.nunique(),n_candidates=0,n_missing_rows=z.value.isna().sum(),missing_reason='zero labelled denominator or absent applicable region remains NA'))
 save(pd.concat([coverage,pd.DataFrame(br)],ignore_index=True),'input_coverage.tsv')
 # Preserve source q names rather than inventing a pooled old q.
 ix=read('input_index.tsv');ix.loc[ix.module.eq('C1')&~ix['columns'].str.contains('q_pearson_within_section',na=False),'columns']+='; q_pearson_within_section, q_pearson_across_sections retained at source';ix.loc[ix.module.eq('C2')&~ix['columns'].str.contains('q_spatial_within_section',na=False),'columns']+='; q_spatial_within_section, q_spatial_across_sections retained at source';save(ix,'input_index.tsv')
 refs=disc[disc.selection.eq('all_identified_units')].copy();selected=disc[disc.selection.eq('exploratory_selected_subset')].copy();selected['datasets']=selected.combination.map(comb.datasets);selected['includes_Inhouse']=selected.datasets.str.contains('Inhouse')
 save(selected[['endpoint_id','module','canonical_program','candidate_type','combination','datasets','includes_Inhouse','estimate','n_analysis_units','q_conditional_across_search','q_between_across_search']],'discussion_selected_datasets.tsv')
 text='# MP10 数据集组合、距离关系和样本异质性：本轮实际结果\n\n'
 text+='本轮已完成。新增的是**小表层面的全组合搜索、真实r/partial r的统一单位汇总、跨输入异质性比较以及niche生态轴对距离结果的解释模型**。没有重跑MP10/UCell/AddModuleScore、cell2location、点级距离或niche聚类，没有恢复点级parquet。既有空间niche标签仍为R12 0–8；本轮Ward分组是分析单元层的候选分组，不是新的spot niche。\n\n'
 text+=f'登记63片(60肿瘤来源、3正常参照)。全Epi中心端点实际可用58张肿瘤切片，对应54个已有analysis units加3个缺单位映射的切片记录，共57条合并记录；不能把这57条全部称为已核实患者或独立分析单元。单位推断主参照因此为6个dataset；七dataset描述保留GSE285715的三张T切片。不同端点有不同覆盖，完整记录见input_coverage和section_manifest。\n\n'
 text+=f'主名称61个，按完整成员合并Budden/Reactome胆固醇别名后为{len(mapping)}个canonical主程序：{int(mapping.n_shared.gt(0).sum())}个有现成noMP10版本，{int(mapping.n_shared.eq(0).sum())}个原集合本就无共享成员。相同成员的别名不在同一矩阵中重复加权；以Budden实际评分为统一代表，没有把两个AddModuleScore数值强行当成完全相同。Top100仅替换Top30作敏感性。\n\n'
 text+=f'共登记{len(reg):,}个统计端点配置，保存{int(scope.n_rows.sum()):,}条数据集组合汇总。七dataset所有127个非空组合均在支持它们的端点中遍历；仅6个dataset有单位身份的矩阵有63个非空组合。第一层异质性覆盖{len(hs)}类输入配置、{len(h):,}条组合/K或状态记录。包括单dataset、无需Inhouse的组合；没有在dataset内部按效应挑片。\n\n'
 text+='## 1. 新增和直接复用的部分\n\n'
 text+='A使用旧80与新增5各自真实保存的R0/R1/R2 r、partial r，按origin去重；B使用B3_Epi_distance/B3_tissue_distance的原逐spot r/partial r，0–1000 μm和零距离保持；C1/C2复用已算出的AMS/UCell候选距离相关、共同支持及联合空间delta；C3/C4复用距离摘要和全Epi MP10相关。新结果的条件P/q来自真实单位层，不借用旧二次模型、β或原始Pearson的P/CI。原q列仍在input_index所指源文件中，未改名冒充本轮校正。\n\n'
 text+='每个统计量先按冻结单位合并重复切片，r在Fisher-z上等权；配对delta直接平均各切片的配对delta，不由两组汇总r相减。对于缺单位映射的记录，只保留section-only描述；正常来源不进入肿瘤异质性建轴。未找到以dataset_reanalysis/mp10_dataset_reanalysis/mp10_inhouse_pca_sweep命名的当前本地缓存，已由原始效应小表重建，不以目录名称缺失阻塞。\n\n'
 text+='## 2. 两类距离轴及全数据参照\n\n'
 text+='下表均为全部六个有身份dataset的等权参照。B3是到**MP10-P95 Epi参考集合**的距离；C1是P95 Epi自身程序到**候选细胞相对富集核心**的距离。二者不是同一个问题。conditional区间/P以所列dataset固定、单位为重复；between检验以dataset均值为重复。所有q_across_search覆盖同端点族的全部组合扫描，不能将条件均值支持等同于可推广的跨队列证据。\n\n'
 text+=tab(refs,['module','canonical_program','candidate_type','estimate','conditional_low','conditional_high','n_analysis_units','unit_direction_fraction','q_conditional_across_search','q_between_across_search'])
 text+='\nB3非脂代谢方向值得关注：Epi中的UPR随离MP10高分参考点更远而总体降低(r=−0.2044，53单位、6dataset，48/53方向为负)；mTORC1 r=−0.2758，51/53为负。全组织MYC_TARGETS_V1 r=−0.3381，54单位，51/54为负。它们是混合spot程序的总体线性关联，不能据此推断所有距离段单调或细胞内代谢通量。\n\n'
 text+='C1中Classical Mono–TNFA/NFKB为负(r=−0.3082，47单位，43/47负向)，Neutrophil–OXPHOS为正(r=+0.3090，45单位，41/45正向)。前者是MP10-P95 Epi内部的炎症程序梯度，后者意味着远离该候选核心时OXPHOS程序总体较高；不能把常见细胞组成关系自动称为MP10特异机制。它们的条件q分别约6.8×10⁻⁵和0.0339，但dataset-replicate跨扫描q约0.170/0.171，跨队列外推证据弱于条件单位层证据。\n\n'
 text+='## 3. P95内部梯度与P95/nonP95差异必须分开\n\n'
 text+='相同Classical Mono–TNFA关系的六dataset配对delta约−0.0107，条件区间[−0.0796,0.0582]，跨搜索q=0.979；Neutrophil–OXPHOS的delta约+0.0996，区间[−0.0155,0.2148]，q=0.833。故其P95内部梯度不能升级为已证实的P95特有关系。\n\n'
 text+='C2存在选定dataset组合的探索性差异：CanDis、GSE171351、GSE238145(D007)的Plasma–OXPHOS配对delta=+0.1993，21单位，条件区间[0.1629,0.2356]，跨该扫描族q=0.0009256。但以3个dataset均值为重复的跨搜索q=0.9462，全部六dataset也没有C2主面板q<0.05(最低约0.1149)。这是**指定组合中的条件差异候选**，不是独立验证或全队列P95差异。所有单片旧联合空间区间及其字段在源表保留；本轮条件区间是单位层变异，不冒充原空间区间。\n\n'
 text+='C3并未给出普遍“P95更靠近Treg”的证据。全六dataset均值差为+63.49 μm(P95更远)，但跨搜索q≈0.352，不能宣称排斥。均值差、中位数、零比例差和匹配摘要分别保留；匹配源结果只有描述性时，没有借用未匹配距离的空间CI/P。\n\n'
 text+='## 4. 哪些完整dataset组合更支持？\n\n'
 text+=tab(selected,['module','canonical_program','candidate_type','combination','datasets','estimate','n_analysis_units','q_conditional_across_search','q_between_across_search'])
 text+='\n上表按效应、单位/数据集方向一致性及覆盖综合排序，至少两个dataset和15单位用于优先讨论；这只是展示规则，完整搜索仍包括所有单dataset和小组合。Classical Mono–TNFA的较强组合D085为CanDis、GSE238145、GSE319536、NPJ，**不包含Inhouse**；其r=−0.3578，但dataset层q=0.198。UPR较强组合D115只排除GSE238145；mTORC1较强组合D113排除GSE171351/GSE238145。这些排除并未改变全数据参照，也不能被解释成某dataset“错误”。七dataset含section-only记录的描述与六dataset单位推断分列，不以缺失填0。\n\n'
 text+='## 5. 哪种输入刻画异质性更清楚？\n\n'
 names=co.input_family.unique();best=[]
 for name in names:
  z=h[h.input_family.eq(name)&h.combination.eq(ur)].sort_values('silhouette',ascending=False)
  if len(z):best.append(z.iloc[0])
 best=pd.DataFrame(best);text+=tab(best,['input_family','n_units','n_features','PC1_variance','PC2_variance','PC1_dataset_eta2','PC1_log_coverage_r','K','silhouette','cluster_sizes'])
 text+='\n在深入检查的五类代表输入中，**全Epi九类niche组成最稳定且最容易解释**。54单位，PC1/PC2解释39.85%/24.17%，PC1的dataset解释比例约7.6%。PC1主要是标签0正向(载荷+0.736)与标签1负向(−0.656)的生态构成对比。原R1参考中，标签0的Epithelial比例约82.1%；标签1约59.1%，同时B比例约18.6%。这是原邻域profile背景，不是本轮中心spot纯度，更不能叫全组织细胞组成。\n\n'
 text+='Ward的描述性最优K为4，组大小20/14/4/16，silhouette≈0.349；80%单位留出平均ARI≈0.833、特征扰动≈0.705、K-means敏感性≈0.765。这支持连续生态轴和若干内部稳定候选组，尚不支持正式生物学亚型。组内niche组成图使用建组变量本身，只是解释，不是独立验证。\n\n'
 text+='R0 partial谱虽有较高silhouette≈0.567，却是53+1；单独单位是CanDis Sample12(n_Epi=405)，K-means一致性低，不能包装为稳定二分类。B3 Epi partial谱和CAF候选梯度更接近连续变化，候选聚类稳定性有限。全组织九大谱系的PC1同时与dataset(η²≈0.344)和Epi点数(r≈0.637)相关，解释时需警惕组成/覆盖驱动。其余核心及扩展60程序分块、P95 R0状态、匹配状态差异及internal/front/near/deep niche均完整进入第一层表，没有因结构弱而删除。\n\n'
 text+=tab(st[st.check.isin(['unit_subset_80percent','feature_subset_80percent','kmeans_same_K'])],['input_family','check','count','mean','min','max'])
 text+='\n## 6. niche生态背景能解释MP10距离关系吗？\n\n'
 text+='已实际拟合29,664个生态轴–端点关系，生态PC1/PC2只用all_Epi niche比例构建，未使用距离r/delta筛轴。模型使用真实分析单位，控制dataset项，另检查技术覆盖、全组织Epithelial比例、可用空间SE权重和dataset留出。**没有关系通过全部生态搜索BH q<0.05**。\n\n'
 ee=eco[eco.endpoint_id.isin([5847,13691])&eco.ecology_axis.eq('PC1')].merge(reg,on='endpoint_id');text+=tab(ee,['module','candidate_type','canonical_program','statistic','n_units','beta','low','high','p_HC3','q_within_endpoint_family','q_across_ecology_search','between_dataset_r','composition_adjusted_beta','technical_adjusted_beta'])
 text+='\n两个可保留的预验证线索：vCAF–UPR的**partial r**随生态PC1增加更负(每1 SD PC1的Fisher-z系数−0.0959，47单位，区间[−0.1365,−0.0553])；NK–DNA_REPAIR配对delta随PC1增加更正(β=+0.1139，48单位，区间[0.0586,0.1693])。六次dataset留出分别保持负向/正向，但全搜索q为0.549/0.974，不能声称生态背景已经解释了这些差异。前者在换回原始r、UCell或其他距离范围时减弱；partial r没有可借用的原始r空间SE。后者可用空间SE加权在同点子集约+0.101，仍属于内部敏感性。\n\n'
 text+='两例的between-dataset相关分别约−0.304和+0.038，不能将内部生态轴关系直接等同于dataset均值差异；逐dataset单位关联和留出结果已单独保存。已有两PC和dataset模型的增量R²是样本内拟合解释量，不是外部预测性能。\n\n'
 text+='## 7. signature、调整、评分、距离和单dataset敏感性\n\n'
 bz=sens[sens.module.eq('B3_Epi_distance')&sens.combination.eq(ur)];text+=tab(bz,['canonical_program','signature','statistic','estimate','conditional_low','conditional_high','q_conditional_across_search','q_between_across_search'])
 text+='\nmTORC1去共享后原始r仍为−0.2497；再控制组成/深度后的partial r减至−0.1101，条件q=0.000100，但dataset层q=0.163。UPR无共享成员，调整后由−0.2044减至−0.0475，dataset层q=0.263。故部分梯度与组成/深度共同变化，不能仅凭noMP10结果宣称细胞内独立机制。\n\n'
 text+=tab(rs,['input_family','variant','n_paired_units','n_common_features','K','adjusted_rand','absolute_PC1_correlation'])
 text+='\nCAF候选梯度的连续PC1在AMS/UCell之间仍有一定对应(|r|≈0.81)，但离散分组ARI接近0；范围和partial处理同样改变分组。n≥30检查在保留下来的单位上比较原完整拟合标签与新拟合，不能将两个完全相同的子矩阵ARI=1当稳定性。Top100替换单一Top30主题对整体轴影响很小，这只是内部冗余敏感性。\n\n'
 text+='C4另有明确的稀疏影响：原全部Epi汇总中，GSE171351_D1仅3个查询点，近±1的r在Fisher-z上对dataset均值影响很大。没有删除原主结果；单列总Epi n≥20敏感性后，Venous EC r由约0.394降至0.116，Matrix_CAF由约0.368降至0.083；进一步partial后约0.043/0.040，跨该敏感性搜索q均>0.1。这是全Epi覆盖问题，不是P95门槛，也不应据原大r给样本命名特有表型。\n\n'
 text+='## 8. 下一轮优先验证的五个方向\n\n'
 text+='1. **MP10-P95周围的mTORC1/UPR梯度。** 多dataset方向支持，mTORC1去共享仍存在；需要新队列预设验证和更明确的上皮来源表达/组成控制，以区分混合spot背景。\n2. **P95 Epi中的Classical Mono–TNFA/NFKB关系。** 目前支持P95内部梯度，不支持P95相对nonP95的特殊性；下一轮需同时测量两组、控制空间位置并验证候选核心定位。\n3. **Neutrophil–OXPHOS的正向距离关系。** 可作为与炎症梯度不同方向的对照；需核查稀有细胞定位、评分方法和低氧/组织结构，不能推断代谢转移。\n4. **Plasma–OXPHOS的指定dataset组合delta。** 条件单位层支持，但跨dataset推断和全数据不支持泛化；需在未参与搜索的队列中固定组合对应问题，不能重复挑阳性dataset。\n5. **R12连续生态轴及候选四组。** 先在独立空间队列验证生态组成轴/分组，再预先限定检验vCAF–UPR、NK–DNA_REPAIR等少数端点；当前FDR未通过，应明确保持探索性。\n\n'
 text+='## 9. 可选缺项和最小提取需求\n\n'
 text+=tab(read('optional_inputs.tsv'),['input','status','minimum_needed'])
 text+='\n已纳入的是真实内侧Epi UCell中位数，以及原P95参考R0均值、P95匹配R0程序差；不是14个Epi MPs的活性替代品。全组织组成来自分区内9大谱系均值按有效spot数精确加权，不用candidate core丰度或niche比例冒充全组织组成。若后续需要全30类型组成，最小提取为section_id、有效组织n、各类型mean fraction/abundance及覆盖；若需要空间组织，需额外输出正式R12邻接计数及边分母。没有验证的临床/分子分型或scRNA–空间对应关系不用于本轮解释。\n\n'
 text+='本轮没有声称重新评估所有样本的曲线形状，也没有生成或虚构逐spot散点。图中的点是已有分析单位效应，缺身份的section-only记录明确区别；曲线摘要不完整不妨碍已有真实r的组合遍历。\n\n'
 text+='## 10. 文件入口和体积\n\n'
 text+='- `tables/input_coverage.tsv`、`canonical_programs.tsv`、`section_manifest.tsv`: 输入、成员映射及全部切片状态。\n- `tables/unit_endpoints.tsv.gz`、`unit_background.tsv.gz`: 唯一统一单位/背景输入；含无法形成单位的section-only记录。\n- `tables/endpoints.tsv.gz`、`combinations.tsv`: 统计端点和完整dataset组合字典。\n- `tables/combination_associations.tsv.gz`: 全部组合效应及两层P/q；字段明确条件推断与dataset-replicate推断。\n- `tables/heterogeneity_configurations.tsv.gz`: 核心/扩展面板全组合PC1/PC2、带符号主载荷、K/规模及技术关系。\n- `tables/ecology_associations.tsv.gz`、`ecology_leave_dataset_out.tsv`、`ecology_sensitivity.tsv.gz`: 完整生态模型及代表敏感性。\n- `tables/discussion_relations.tsv`、`discussion_selected_datasets.tsv`、`discussion_sensitivity.tsv.gz`: 精简讨论入口。\n- `tables/figure_index.tsv`、`figures/`: 九组、16张独立代表图，Arial白底、矢量PDF及300dpi PNG；有图例时三版PDF。\n- `scripts/run_all.sh`、`config.json`、`METHODS.md`: 可重启运行入口和实际定义。\n\n'
 text+='未创建ZIP或重复上传目录。全部源结果只读。本轮不恢复点级缓存、不复制表达/SCT/反卷积对象，只保留必要矩阵和汇总。@@SIZE@@\n'
 p=O/'REPORT.md';p.write_text(text)
 for cache in (O/'scripts').rglob('__pycache__'):shutil.rmtree(cache)
 size=sum(x.stat().st_size for x in O.rglob('*') if x.is_file());p.write_text(text.replace('@@SIZE@@',f'新增文件实际合计约 **{size/1e6:.2f} MB ({size/1024**2:.2f} MiB)**，低于100 MB；未删除任何旧项目结果。'))
 print('Report complete',size,'bytes',flush=True)
if __name__=='__main__':main()
