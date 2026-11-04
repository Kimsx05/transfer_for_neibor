from aggregate import *
import zipfile,json,shutil

def table(df,cols,rename=None):
 z=df[cols].copy()
 if rename:z=z.rename(columns=rename)
 def fmt(x):
  if isinstance(x,(float,np.floating)):return 'NA' if not np.isfinite(x) else f'{x:.4g}'
  return str(x).replace('|','/')
 return '| '+' | '.join(z.columns)+' |\n| '+' | '.join(['---']*len(z.columns))+' |\n'+'\n'.join('| '+' | '.join(fmt(x) for x in row)+' |' for row in z.itertuples(index=False,name=None))+'\n'

def main():
 frozen=pd.concat([pd.read_csv(OLD/'tables'/f'{sec}_coverage.tsv',sep='\t') for sec in M.section_id]).set_index('section_id');sections=pd.read_csv(O/'full_sections.tsv',sep='\t');sections.loc[sections.n_P95.eq(1),'P95_status']='single_P95_spot_not_inferential';sections['MP10_P95_threshold']=sections.section_id.map(frozen.P95_threshold);sections['pEMT_P75_frozen']=sections.section_id.map(frozen.pEMT_P75);sections['n_matched_base_frozen']=sections.section_id.map(frozen.n_match_base);sections['n_matched_front_frozen']=sections.section_id.map(frozen.n_match_front);save(sections,O/'full_sections.tsv')
 from section_notes import main as section_notes
 section_notes()
 c=json.loads((O/'tables/result_counts.json').read_text());cov=pd.read_csv(O/'full_sections.tsv',sep='\t');e=pd.read_csv(O/'tables/illustrative_Pearson_examples.tsv',sep='\t');dd=pd.read_csv(O/'tables/illustrative_delta_examples.tsv',sep='\t');hier=pd.read_csv(O/'tables/context85_hierarchical_summary.tsv',sep='\t');a=pd.read_parquet(O/'tables/pearson_all.parquet');de=pd.read_parquet(O/'tables/delta_all.parquet');h=pd.read_csv(O/'tables/dataset_heterogeneity.tsv.gz',sep='\t');rec=pd.read_csv(O/'tables/dataset_priorities.tsv',sep='\t');niche=pd.read_csv(O/'tables/front_only_niche_all.tsv.gz',sep='\t');source='tumor_study_specimen_not_region_diagnosis';s=hier[hier.source.eq(source)&hier.population.eq('units')&hier.adjustment.eq('adjusted')];ctr=s[s.analysis.eq('B1_R0')];focus=['HALLMARK_UV_RESPONSE_UP','HALLMARK_MTORC1_SIGNALING','HALLMARK_MTORC1_SIGNALING__NO_MP10','HALLMARK_UNFOLDED_PROTEIN_RESPONSE','HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY','HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT','PURAM_PEMT_TOP30','FLORES_CRISPR_PRO_INVASION_10']+NEW
 covsummary=cov.groupby('source_context').agg(sections=('section_id','size'),tissue_spots=('n_spots','sum'),Epi=('n_Epi','sum'),internal_Epi=('n_internal','sum'),P95=('n_P95','sum'),exact_front_Epi=('n_front_exact_Epi','sum')).reset_index()
 text='# MP10逐切片Pearson与样本差异：实际结果\n\n'
 text+=f'本轮已登记并完成全部 {c["n_sections"]} 张切片，{c["n_tissue"]:,} 个组织spot，85个程序(61主、24敏感性)和29个非Epi候选。所有切片均有三页报告与模块状态；没有依据效应方向、显著性或前沿覆盖筛片。MP10仍为冻结Top20 AUCell，新AddModuleScore不改变P95或pEMT-high front。\n\n'+table(covsummary,list(covsummary.columns))
 text+='\n## 1. MP10高的Epi周围有哪些功能变化？\n\n以下为与上一轮兼容的UCell中心/邻域模型，不是本轮候选距离Pearson。旧80项点级结果复用，新增5项同法补齐，并在85程序面板内重新进行族校正。效应先在既有analysis unit内等权，再在dataset间等权；95%区间为既有层级bootstrap，P来自dataset均值检验，因此二者不是同一检验的对偶。\n\n'
 z=ctr.set_index('program').reindex(focus).reset_index();text+=table(z,['program','effect','low','high','q_context85','n_units','n_datasets','positive_units'])
 text+='\n非脂代谢中心信号中，UV_RESPONSE_UP的调整标准化β为0.04894 (95%区间0.01895–0.07787，85面板族q=0.01356；54单位、6个dataset中41单位为正)，效应较小但具有重复性。UPR、ROS和去MP10重叠mTORC1也偏正，但该层级q分别约0.070、0.076和0.093，暂作为支持不足的后续线索。Puram pEMT中心关系偏负(β=−0.05084，区间−0.08385至−0.02510，q=0.01032)，不能把MP10高直接等同于pEMT高。新增S1P和TG程序没有一致的总体中心升高。\n\n'
 text+='\n原始脂代谢集合与MP10的共享基因会增加一致性；Budden胆固醇27与Top20重叠9个，去重叠18成员单列。mTORC1去重叠版本与原版也分别保留；OXPHOS原集合没有Top20重叠，因此不制造额外“去重叠”集合。程序名称不是代谢通量、细胞活性或暴露本身的直接测量。\n\n'
 rz=s[s.analysis.isin(['B1_'+r for r in ['R0','R1','R2','R3_R5','R6_R10']])&s.program.isin(['HALLMARK_MTORC1_SIGNALING__NO_MP10','HALLMARK_UNFOLDED_PROTEIN_RESPONSE','PURAM_PEMT_TOP30','BUDDEN_S1P_RESPONSE_18','BUDDEN_TNF_EMT_RESPONSE_15'])];text+=table(rz,['program','analysis','effect','low','high','q_context85','n_units'])
 text+='\n这些是中心MP10与对应环平均程序的联系。R0不混入邻域平均；缺环不填0。旧重叠邻域仍采用原2000 μm空间块处理，未硬套新点级500 μm方案。P95/匹配对照、d_MP10曲线与前沿模型也在context85完整表中；二次项不解读成全程单调斜率。\n'
 text+='\n## 2. P95 Epi自身程序随候选距离如何变化？\n\n'
 text+=f'AddModuleScore主结果共有 {c["n_P95_computable"]:,} 个可计算的P95全范围相关(含敏感性程序)。在主61程序且n≥30、贡献空间块≥8的组合中，r<−0.1有 {c["n_negative"]:,} 个，r>0.1有 {c["n_positive"]:,} 个，|r|≤0.1有 {c["n_near_zero"]:,} 个。这是重叠候选/程序扫描的描述性计数，不能当作独立生物学发现数。\n\n'
 text+='下列覆盖充分的例子分别保留负向、正向和接近零，按效应和覆盖展示，并非只列显著结果。距离为候选相对丰度P95核心，而非到MP10-P95集合。\n\n'
 text+=table(e,['example_type','section_id','candidate_type','program','pearson_r','n_spots','distance_min_um','distance_max_um','spatial_ci_low','spatial_ci_high','p_pearson','q_pearson_across_sections'])
 text+='\n负r表示总体负线性关系，不证明每个距离段都下降；近零也可能含弯曲或平台。报告中的散点为原始spot自身AddModuleScore，负分保留，完整范围不因结果裁切。普通P/q基于独立spot假设；空间区间单独列，不能以q替代空间相关处理。\n'
 text+='\n## 3. P95与非P95是否不同？\n\n'
 text+=f'共同距离范围内、两组n≥30且各≥8空间块的组合中，{c["n_delta_CI_supported"]:,} 个delta_r的联合块bootstrap区间不跨0；其中空间差异P的跨切片同族BH q<0.05为 {c["n_delta_spatial_q05"]:,} 个(主61程序 {c["n_delta_spatial_q05_main"]:,} 个)。计数含29候选及相关程序的重叠扫描，区间与Fisher-z空间检验不是完全相同的检验，因此支持数不必相同。\n\n'
 text+=table(dd,['section_id','candidate_type','program','program_family','n_P95','n_nonP95','common_min_um','common_max_um','r_P95','r_nonP95','delta_r','spatial_ci_low','spatial_ci_high','p_spatial_fisher_z','q_spatial_across_sections'])
 text+='\n同一空间块抽样同时用于两组，重复抽中块保留重复权重。delta_r是相关强度差，不是单位距离下降速度；没有以“一组显著、另一组不显著”替代比较。full_range delta仅作描述，正式结果使用common_support，跨样本另提供common_support∩[0,1000] μm。两组实际距离方差、断档、零距离比例及n均保留。\n'
 text+='\n## 4. 同dataset内的切片/样本异质性\n\n'
 hz=h[h.program.isin(CORE)&h.population.eq('P95')&h.n_units.ge(3)&h.n_positive.gt(0)&h.n_negative.gt(0)].copy();hz['spread']=hz.r_max-hz.r_min;hz=hz.sort_values('spread',ascending=False).drop_duplicates('dataset_id').head(7);text+=table(hz,['dataset_id','candidate_type','program','n_sections','n_units','r_min','r_max','n_negative','n_positive','n_near_zero','p_heterogeneity','q_heterogeneity'])
 text+='\n上表使用统一0–1000 μm内各片两组共同支持范围。异质性Q使用空间bootstrap Fisher-z标准误，未使用n_spot−3精度权重；描述性单位均值仍按切片等权Fisher-z合并。身份缺失者保留切片描述，不补造patient。表中出现双向效应说明值得逐片核查，只有直接区间/异质性检验有支持时才讨论统计差异。\n\n'
 # Direct real section examples for the same relationship, preserving min/max from the common window.
 opposite=[]
 for rr in hz.itertuples():
  z=a[a.dataset_id.eq(rr.dataset_id)&a.candidate_type.eq(rr.candidate_type)&a.program.eq(rr.program)&a.score_method.eq('AddModuleScore_SCT_data')&a.population.eq('P95')&a.range_type.eq('common_support_0_1000')&a.pearson_r.notna()];opposite.extend([z.sort_values('pearson_r').iloc[0],z.sort_values('pearson_r').iloc[-1]])
 op=pd.DataFrame(opposite);save(op,O/'tables/within_dataset_opposite_examples.tsv')
 if len(op):text+=table(op,['dataset_id','section_id','candidate_type','program','pearson_r','n_spots','distance_min_um','distance_max_um','spatial_ci_low','spatial_ci_high','q_pearson_across_sections'])
 text+='\n样本单位与已核实patient分别汇总，B5/B18按已有P1关系合并。重复切片有限，以上优先称“切片/样本间异质性”，不能仅凭一张切片命名稳定患者表型。GSE285715缺analysis_unit映射的切片保留全部效应及图册，不混进患者推断。\n'
 # A fixed discussion example retains negative, positive and near-zero sections.
 discussion=a[a.section_id.isin(['GSE319536_BC11','GSE319536_BC12','GSE319536_BC13','GSE319536_BC18'])&a.candidate_type.eq('vCAF')&a.program.eq('BUDDEN_TNF_EMT_RESPONSE_15')&a.score_method.eq('AddModuleScore_SCT_data')&a.population.eq('P95')&a.range_type.isin(['full_range','common_support_0_1000'])].copy()
 save(discussion,O/'tables/discussion_vCAF_TNF_Pearson.tsv')
 text+='\n### 可用于机制对照的具体例子：GSE319536 / vCAF / Budden TNF–EMT response\n\n'+table(discussion,['section_id','range_type','pearson_r','n_spots','distance_min_um','distance_max_um','spatial_ci_low','spatial_ci_high','p_pearson','q_pearson_across_sections'])
 text+='\nBC13完整范围r=−0.469 (n=45，100–900 μm)，非P95共同范围r=+0.286 (n=776)，delta_r=−0.755，联合空间区间[−1.043,−0.380]，空间P=3.5×10⁻⁵、跨片同族q=0.00416。BC12也负相关，但P95与非P95差异区间跨0，不能称为P95特有。BC18完整范围为弱正向且空间区间跨0；限定预设共同0–1000 μm窗口时r=+0.233 (n=90)，区间[0.040,0.422]，提示相反方向但对距离范围敏感。BC11接近零。这四片共同说明局部关系存在异质性，不能泛化为所有MP10高Epi都呈同样梯度。\n\n'
 text+='同dataset统一窗口的vCAF–TNF程序P95相关异质性Q=42.116，空间标准误检验P=6.3×10⁻⁵、族q=0.000627 (15单位，14个有可用空间区间)。BC12与BC18的直接Fisher-z差异为−0.922，区间[−1.275,−0.570]，固定核心面板配对校正q=0.000289。这个对照优先于仅比较两片各自是否显著。BC13候选核心posterior mean丰度中位数1.76，但仅6个P95点有精确前沿距离、其中3个在前沿带，因此不能判断其关系是否独立于前沿位置。\n'
 text+='\n## 5. 前沿和niche能否解释？\n\n'
 f=s[s.analysis.isin(['A_front_vs_internal','C1_all_Epi_continuous','C2_P95_enrichment','C_active_front_MP10','C3_P95_state'])&s.program.isin(['PURAM_PEMT_TOP30','FLORES_CRISPR_PRO_INVASION_10',RAW,'P95_probability']+NEW)];text+=table(f,['analysis','program','effect','low','high','q_context85','n_units','n_datasets'])
 text+='\n跨单位/数据集没有支持所有前沿都具有更高pEMT或Flores分数：前沿对内部的调整β分别−0.0184 (q=0.928)和−0.0442 (q=0.557)。MP10连续前沿距离、P95前沿富集及分子活跃前沿的MP10总体比较亦无一致支持(q分别0.884、0.665、0.665)。这不排除局部边段或单片现象，但不能把全队列边界统一命名为已验证侵袭前沿。\n\n'
 text+='\n原有精确v3.4前沿、pEMT-high定义与R12标签保持冻结。缺精确前沿距离只限制位置分析，不限制候选距离相关。补充partial r先控制Epi比例与深度，再在同一精确前沿子集比较加入d_front前后；不能把换了spot集合导致的变化叫作前沿调整效应。\n\n'
 text+='\n仅前沿带内仍存在片别的niche差异，例如CanDis_2026_Sample4的R12标签2在P95中低28.47个百分点 (标注P95/非P95 n=22/288，区间−42.55至−14.81个百分点，空间跨片q=0.0185)；GSE238145_B14的标签7低22.04个百分点 (n=12/313，区间−35.49至−9.84，q=0.0378)。这是同一前沿带内的比例差，不代表已匹配带内每一点的几何位置，也不是所有切片的共同niche偏好。\n\n'
 n=niche[niche.n_P95.ge(10)&niche.n_nonP95.ge(10)&niche.blocks_P95.ge(8)&niche.blocks_nonP95.ge(8)].copy();n['abs_difference']=n.difference_percentage_points.abs();n=n.sort_values('abs_difference',ascending=False).drop_duplicates('dataset_id').head(7);text+='本轮额外完成严格[-200,0] μm带内的原P95 vs非P95九类niche比较，分母为各组已有正式标签的Epi，不在前沿内重新取Top5%。以下展示实际比例差，缺少空间区间的行仅为描述：\n\n'+table(n,['section_id','niche','n_P95','n_nonP95','fraction_P95','fraction_nonP95','difference_percentage_points','ci_low_percentage_points','ci_high_percentage_points','q_spatial_across_sections'])
 ss=pd.read_csv(O/'tables/sc_lineage_sample_source.tsv.gz',sep='\t');sm=ss.groupby(['program','lineage'])['median'].median().unstack().reset_index();text+='\n新增文献程序的非配对单细胞来源支持(样本层UCell中位数的中位数)：\n\n'+table(sm,list(sm.columns))
 text+='\nS1P响应和TNF–EMT响应集合可在CAF/SMC等来源显著表达，空间相关不能直接称为Epi专属信号。胆固醇去重叠版本单列用于区别共享基因一致性；这些UCell来源分数与空间SCT AddModuleScore不能直接比较绝对数值。\n'
 if len(dd):
  kk=['section_id','candidate_type','program'];adj=a[a.population.eq('P95')&a.range_type.eq('common_support')&a.score_method.eq('AddModuleScore_SCT_data')].merge(dd[kk],on=kk,validate='1:1');text+='\n上文共同范围示例中的组成/深度及前沿敏感性(前沿加入前后严格同点)：\n\n'+table(adj,['section_id','candidate_type','program','pearson_r','partial_r_epi_depth','n_partial','partial_r_front_subset_before','partial_r_front_subset_after','n_front_subset'])
 text+='\n协变量完整点数小于10的partial r单独标记为稀疏描述；特别是5–6个点在加入多个协变量后可因剩余自由度极少接近±1，不能用来证明调整稳健性。原始Pearson仍按可计算点数保留，不因无法稳定调整而删除。\n\n'
 text+='\n全部九类及未标注覆盖保留在niche_counts_all；原前沿位置匹配和一般前沿分布保留原名称，未改名为仅前沿带比较。R12解释使用原R1/R2 profile，不用中心spot组成冒充邻域组成。\n'
 text+='\n## 6. 优先审阅与证据限制\n\n'+table(rec,['dataset_id','n_sections','n_P95_ge30','n_units','n_verified_patients','n_exact_front_sections','n_core_pairs_with_opposite_section_directions'])
 text+='\n优先顺序按P95数量、共同距离覆盖、同dataset重复单位及可用于位置解释的前沿覆盖安排；具体候选–程序和切片见section_relationship_priorities。方向相反的切片是机制对照线索，未因“难看”删除。罕见类型的P95只是相对排名，候选核心绝对丰度和实际目标比例见candidate_targets，不能仅凭相对P95称大量存在。\n\n'
 text+=f'主要覆盖限制：{int(cov.n_Epi.eq(0).sum())} 张无冻结Epi；{int(cov.n_P95.eq(1).sum())} 张仅1个内部P95 Epi；{int(cov.n_P95.ge(30).sum())} 张达到P95≥30。GSE238145_B13保留为无Epi模块状态；GSE171351_D1仅1个内部Epi；GSE319536_BC15有4个内部Epi、仅1个P95，两片均不能计算P95相关。三张邻近正常只作合格Epi参照，没有伪造肿瘤P95或前沿。\n\n'
 text+='所有85项AddModuleScore以完整可用SCT data基因池评分，SCT拟合要求至少5个细胞检出，故少数低检测成员可能不进入该切片SCT池；各片实际成员/缺失和RNA检测数完整记录。新5项UCell按完整RNA排名宇宙补齐。单细胞来源支持复用上一轮19,053个分层抽样细胞，按真实样本和原注释汇总；与空间不配对者仅支持可能来源。混合spot信号不能直接写成Epi特异表达、脂质转移、S1P浓度或代谢通量。空间500 μm块近似、共同范围缺口和少量高端P95仍限制精度。\n\n'
 text+='63片评分、四项候选比较、旧邻域/前沿/niche整理、新前沿带内niche比较、单细胞来源及逐片/数据集图册均已实际运行。没有重跑cell2location、niche、CNV或通讯。旧结果保持不变。\n'
 (O/'RESULTS.md').write_text(text);(O/'summary.md').write_text(text)
 # Member table with actual RNA/SCT alias mapping rather than pretending an alias is missing.
 cover=pd.read_csv(O/'tables/program_gene_coverage.tsv.gz',sep='\t');alias={'CCN2':'CTGF','CAVIN3':'PRKCDBP','ACKR3':'CXCR7','GSDME':'DFNA5','P3H2':'LEPREL1'};maps=[]
 for rr in M.itertuples():
  genes=set(Path(rr.genes_path).read_text().splitlines());v=cover[cover.section_id.eq(rr.section_id)].set_index('program')
  for p,gs in G.groupby('gs_name'):
   used=set(str(v.loc[p,'used_symbols']).split('|'))
   for row in gs.itertuples():
    g=row.gene_symbol;actual=g if g in genes else alias.get(g,g);maps.append(dict(section_id=rr.section_id,program=p,original_symbol=row.original_symbol,accepted_symbol=g,RNA_symbol=actual if actual in genes else None,RNA_available=actual in genes,SCT_available=actual in used,MP10_overlap=row.MP10_overlap,panel=row.panel))
 save(maps,O/'tables/program_symbol_mapping.tsv.gz');save([dict(program=p,panel=gs.panel.iloc[0],original_size=gs.original_symbol.nunique(),mapped_size=gs.gene_symbol.nunique(),MP10_overlap_n=gs.loc[gs.MP10_overlap,'gene_symbol'].nunique(),overlap='|'.join(gs.loc[gs.MP10_overlap,'gene_symbol']),source=gs.source.iloc[0]) for p,gs in G.groupby('gs_name')],O/'tables/program_index.tsv')
 # Original matching and detection diagnostics retain old names/meaning.
 for filename in ['matching.tsv.gz','signature_detection_summary.tsv','gene_detection.tsv.gz','summary_effects.tsv','bands.tsv.gz','niche.tsv.gz']:
  if (OLD/'tables'/filename).exists():shutil.copy2(OLD/'tables'/filename,O/'tables'/('reused_previous80_'+filename))
 profile=R/'neighborhood_ring_harmony_v3/03_main/SAMPLE_R12_LEIDEN/profiles_major.tsv';shutil.copy2(profile,O/'tables/R12_original_major_profiles.tsv')
 readme='''# MP10逐切片Pearson分析\n\n本目录包含全部63片实际计算，主分析为SCT data AddModuleScore + 原始逐spot Pearson；主MP10仍为冻结AUCell。\n\n- `full_overview.pdf`: 全队列总览。\n- `RESULTS.md` / `summary.md`: 中文解释、实际数值与限制。\n- `full_sections.tsv`: 每片身份、输入、覆盖及模块状态。\n- `reports/<section>.pdf`: 每片3页主报告；同名MD为覆盖及具体关系。\n- `atlas/<dataset>.pdf`: 全部切片与同dataset样本比较。\n- `figure_index.tsv`: 主图、逐片图位置；各图有300dpi PNG和有图例时的三个PDF。\n- `tables/pearson_all.parquet`: 完整85×29×两评分×各人群/范围长表。\n- `tables/delta_all.parquet`: full_range描述、common_support正式比较、0–1000共同支持。\n- `tables/MP10_candidate_all.tsv.gz` / `closeness_all.tsv.gz`: 比较四/一。\n- `tables/context85_effects.parquet` / `context85_hierarchical_summary.tsv`: 旧80+新5的UCell上下文模型，β不等于Pearson。\n- `tables/unit_correlations.tsv.gz` / `dataset_heterogeneity.tsv.gz`: 真实身份层级与空间不确定性异质性。\n- `tables/front_only_niche_all.tsv.gz`: 仅前沿带新增九类niche比较。\n- `tables/dataset_priorities.tsv` / `section_relationship_priorities.tsv.gz`: 覆盖和效应驱动的审阅建议。\n- `tables/program_gene_coverage.tsv.gz` / `program_symbol_mapping.tsv.gz`: 实际成员、别名、检测与SCT覆盖。\n- `spots/<section>.parquet`: 冻结身份、三距离轴、29候选、85 UCell及85 AMS绘图源数据，无重复表达矩阵。\n- `scores/`: 紧凑评分与参数；不保存大型SCT对象。\n- `scripts/run_all.sh`: 可重启入口；`scripts/redraw.py`: 任意候选/程序散点。\n\n所有普通Pearson P/q都有独立spot假设限制；空间区间/P另列。灰色NA不等于0。相邻正常与肿瘤来源分开。完整Parquet可用pandas/pyarrow读取，主程序与敏感性在program_family列明确区分。review.zip只放说明、总览、关键小表和数据集森林图；大型全组合表、逐点缓存、表达数据不放入包。\n'''
 (O/'README.md').write_text(readme)
 ix=pd.read_csv(O/'figure_index.tsv',sep='\t');known=set(ix.pdf);added=[]
 for path in (O/'figures').rglob('*.pdf'):
  if path.stem.endswith(('_panel_only','_legend_only')) or str(path) in known:continue
  added.append(dict(section_id='',dataset_id=path.parent.name if path.parent.name not in ['overview','datasets'] else '',page='supplement',pdf=str(path),png=str(path.with_suffix('.png')) if path.with_suffix('.png').exists() else None))
 ix=pd.concat([ix,pd.DataFrame(added)],ignore_index=True)
 for suffix,col in [('_panel_only.pdf','panel_pdf'),('_legend_only.pdf','legend_pdf')]:ix[col]=[str(Path(p).with_name(Path(p).stem+suffix)) if Path(p).with_name(Path(p).stem+suffix).exists() else None for p in ix.pdf]
 save(ix,O/'figure_index.tsv')
 files=[O/'README.md',O/'METHODS.md',O/'RESULTS.md',O/'full_sections.tsv',O/'figure_index.tsv',O/'full_overview.pdf']
 small=['candidate_mapping.tsv','candidate_targets.tsv.gz','program_index.tsv','result_counts.json','dataset_priorities.tsv','dataset_display_selection.tsv','illustrative_Pearson_examples.tsv','illustrative_delta_examples.tsv','within_dataset_opposite_examples.tsv','front_only_niche_all.tsv.gz','context85_hierarchical_summary.tsv','dataset_equal_summary.tsv.gz','section_relationship_priorities.tsv.gz']
 files += [O/'tables'/n for n in small if (O/'tables'/n).exists()];files += list((O/'figures/datasets').glob('*.pdf'));files += [O/'figures/discussion_GSE319536_vCAF_TNF.pdf',O/'tables/discussion_vCAF_TNF_Pearson.tsv']
 with zipfile.ZipFile(O/'review.zip','w',zipfile.ZIP_DEFLATED) as z:
  for p in files:z.write(p,str(p.relative_to(O)))
 # Disposable sparse scRNA extract is not a scientific output; source and sample IDs remain indexed.
 for n in ['sc_temporary.mtx','sc_genes.txt']:
  p=O/'scores'/n
  if p.exists():p.unlink()
 print('Final reports and review package complete',flush=True)
if __name__=='__main__':main()
