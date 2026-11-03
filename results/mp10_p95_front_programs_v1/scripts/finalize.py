from plots import *
import zipfile,shutil,textwrap
S=pd.read_csv(O/'tables/summary_effects.tsv',sep='\t');T=S[S.source.eq('tumor_study_specimen_not_region_diagnosis')&S.population.eq('units')];C=pd.read_csv(O/'full_sections.tsv',sep='\t')
def row(an,p,adj='adjusted'):
 a=T[T.analysis.eq(an)&T.program.eq(p)&T.adjustment.eq(adj)];return a.iloc[0] if len(a) else None
def fmt(v,d=3):return 'NA' if pd.isna(v) else f'{v:.{d}f}'
def eff(an,p,adj='adjusted'):
 a=row(an,p,adj)
 return '不可估计' if a is None else f"β={a.effect:.3f} (95% CI {a.low:.3f}–{a.high:.3f}), FDR={a.q:.3g}; {a.n_units}单位/{a.n_datasets}datasets"
def mdtable(frame):
 f=frame.copy();return '| '+' | '.join(map(str,f.columns))+' |\n| '+' | '.join(['---']*len(f.columns))+' |\n'+'\n'.join('| '+' | '.join(map(str,r))+' |' for r in f.itertuples(index=False,name=None))
# Four complete module files; no positive-result selection.
e=pd.read_csv(O/'tables/effects.tsv.gz',sep='\t')
for module,prefix in [('A_front','A_'),('B_neighborhood','B'),('C_MP10_front','C'),('D_niche','D')]:save(e[e.analysis.str.startswith(prefix)],O/'tables'/f'{module}_statistics.tsv.gz')
# Discussion shortlist contains cautions, not a new statistical selection rule.
rows=[]
for p,note in [('HALLMARK_UV_RESPONSE_UP','零MP10重叠；小效应。集合名称不代表UV暴露，需独立应激表型。'),('HALLMARK_MTORC1_SIGNALING','与MP10重叠10基因；去重叠效应约减半，尚未通过敏感性FDR。'),('PURAM_PEMT_TOP30','中心负向且随环减弱；混合spot与CAF/SMC来源阻止解释为MP10抑制上皮侵袭。'),('HALLMARK_UNFOLDED_PROTEIN_RESPONSE','无MP10重叠；小的正向趋势，尚未通过主族FDR。'),('HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY','无MP10重叠；小的正向趋势，需重复验证。'),('HALLMARK_ANDROGEN_RESPONSE','去除7个共享基因后关联基本消失，不作独立雄激素机制证据。'),('HALLMARK_TGF_BETA_SIGNALING','去除SLC20A1后关联基本消失，不作独立TGF-β机制证据。')]:
 a=row('B1_R0',p);b=row('B1_R0',p+'__NO_MP10');rows.append(dict(program=p,effect=a.effect,low=a.low,high=a.high,FDR=a.q,positive_units=a.positive_units,negative_units=a.negative_units,n_units=a.n_units,no_overlap_effect=b.effect if b is not None else np.nan,no_overlap_FDR=b.q if b is not None else np.nan,interpretation=note))
save(rows,O/'tables/discussion_candidates.tsv')
# Descriptive predefined band effects are shown in full, including near/deep comparisons.
z=T[T.analysis.isin(['A_front_vs_near','A_front_vs_deep'])&T.adjustment.eq('adjusted')&T.program.isin(MAIN)];heat(z.pivot(index='program',columns='analysis',values='effect'),'01_front_near_deep','Adjusted front-minus-interior effect (SD)',.4)
# Section-level heterogeneity for main molecular evidence and MP10.
x=e[e.analysis.isin(['A_front_vs_internal','C1_front_vs_internal'])&e.adjustment.eq('adjusted')&e.program.isin([PE,'FLORES_CRISPR_PRO_INVASION_10',RAW])].merge(M[['section_id','dataset_id']],on='section_id');save(x,O/'tables/front_section_effects.tsv')
idx=pd.read_csv(O/'figure_index.tsv',sep='\t');idx=idx[~idx.name.isin([i['name'] for i in index])];idx=pd.concat([idx,pd.DataFrame(index)],ignore_index=True);save(idx,O/'figure_index.tsv')
# Update software-independent coverage checks that affect scientific correctness.
checks=[]
for sec in M.section_id:
 d=pd.read_parquet(O/'spots'/f'{sec}.parquet');f=pd.read_parquet(F/'spots'/f'{sec}.parquet',columns=['barcode','query',RAW]);c=pd.read_parquet(O/'centers'/f'{sec}.parquet');assert len(d)==len(f) and d.barcode.is_unique;assert d.barcode.equals(f.barcode) and np.array_equal(d['query'],f['query']) and np.allclose(d[RAW],f[RAW],equal_nan=True)
 assert d.loc[d.MP10_P95,'d_MP10'].eq(0).all();assert (d.loc[d.lesion_location.eq('inside'),'d_front'].dropna()<=0).all();assert (d.loc[d.lesion_location.eq('outside'),'d_front'].dropna()>=0).all();assert c.R0_n_neighbors.eq(1).all();assert d.accepted_label.dropna().isin(range(9)).all();assert np.allclose(c['R0__'+PE],d.set_index('barcode').loc[c.barcode,PE],equal_nan=True)
 checks.append(dict(section_id=sec,n_spots=len(d),n_query=d['query'].sum(),frozen_query_and_MP10_preserved=True,P95_self_distance_zero=True,R0_is_center=True,front_sign_preserved=True))
assert sum(x['n_spots'] for x in checks)==265002;save(checks,O/'tables/essential_checks.tsv')
# Concise research report, numeric values drawn from the final tables.
r0=[]
for p in ['HALLMARK_MTORC1_SIGNALING','HALLMARK_ANDROGEN_RESPONSE','HALLMARK_TGF_BETA_SIGNALING','HALLMARK_UV_RESPONSE_UP','HALLMARK_UNFOLDED_PROTEIN_RESPONSE','HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY','HALLMARK_ANGIOGENESIS',PE,'FLORES_CRISPR_PRO_INVASION_10']:
 a=row('B1_R0',p);b=row('B1_R0',p+'__NO_MP10');r0.append([p,fmt(a.effect),f'{a.low:.3f}–{a.high:.3f}',fmt(a.q),fmt(b.effect) if b is not None else '无重叠',fmt(b.q) if b is not None else '不适用',f'{a.positive_units}/{a.n_units}'])
r0=pd.DataFrame(r0,columns=['程序','R0调整β','95% CI','FDR','去重叠β','去重叠FDR','正向单位'])
front=[]
for p in [PE,'PURAM_PEMT_TOP100_SENSITIVITY','FLORES_CRISPR_PRO_INVASION_10','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION']:
 a=row('A_front_vs_internal',p);front.append([p,fmt(a.effect),f'{a.low:.3f}–{a.high:.3f}',fmt(a.q),f'{a.positive_units}/{a.n_units}'])
front=pd.DataFrame(front,columns=['程序','前沿减内部调整β','95% CI','FDR','正向单位'])
summary=f'''# MP10 邻域功能、正式前沿及 R12 niche

完成63张切片、265,002个组织spots；60张肿瘤研究来源和3张邻近正常分开。冻结Epi queries共143,461个，其中肿瘤研究来源141,152个；v3.4内侧合格肿瘤Epi 112,547个，P95参考点5,664个(59片)，前沿带Epi 12,517个(57片)，pEMT-high front 3,428个。基础匹配5,623对，精确前沿位置匹配1,979对(56片；达到≥5对的单位汇总45单位)。所有55原程序、新增Top30/Flores10及23项敏感性均计算并保留。

下列主要统计为肿瘤来源、已有analysis_unit合并后dataset等权。R0通常覆盖55片/54单位/6datasets；另外保留全部7datasets的section描述和27名核实患者/3datasets的身份子集。不把未知身份改称患者。β是标准化回归系数或标准化均值差，不是普通spot P值。95% CI为层级bootstrap，FDR来自dataset层检验的预设族校正；二者不是同一种检验的对偶关系。完整数据与负结果在tables中。

## 1. MP10 周围有哪些非脂代谢功能变化？

最清楚的结果位于中心R0。原分数中mTORC1、androgen、TGF-β等与MP10同向，但去共享基因改变了解释：androgen去7个成员、TGF-β去SLC20A1后效应接近0。mTORC1去10个重叠成员后仍为正向趋势，强度约减半，敏感性FDR未达0.05。因此不能把这些原关联全部写成独立于胆固醇程序的新机制。

{mdtable(r0)}

无Top20重叠的UV_RESPONSE_UP为小而较一致的中心正关联(41/54单位、6/6datasets同向；各dataset β约0.022–0.078)，集合名称不能解释为实际紫外线暴露。UPR、ROS和PI3K–AKT–mTOR是小效应候选，主族FDR尚不足。增殖E2F/G2M、Flores10、IFN等未获得一致中心上调证据。p-EMT、angiogenesis、广义EMT/ECM与中心MP10呈小的负关联；p-EMT Top100方向相同，不能据此推断MP10抑制侵袭。

中心到远端呈衰减而非稳定的全邻域增强：mTORC1调整β按R0/R1/R2/R3–5/R6–10为0.159/0.109/0.084/0.065/0.045；p-EMT为−0.051/−0.038/−0.028/−0.015/−0.002。实际平均物理距离的切片中位数约为0/100/187/381/746μm，不能把图环号当作精确距离。各环缺覆盖按自身可用点计算。

同前沿位置的P95匹配对照中，R0 mTORC1约+0.300 SD，androgen约+0.342 SD，但整个匹配程序×环检验族FDR约0.093；近邻和远环没有通过全族FDR的稳定非脂代谢增高。完整P95集合/对照的距离曲线、互斥图环、组织与Epi目标分别保留；视觉曲线不能替代配对效应。原始全部P95场和匹配后参考子集分别命名，未把匹配掉的参考点悄悄改为主P95。

核实患者子集中，p-EMT中心负向21/27、mTORC1正向26/27，但仅3个dataset，相关族FDR较宽；不能把匿名单位的较强汇总写成已经独立患者验证。

## 2. 现有前沿获得了什么分子支持？

**没有得到整个v3.4前沿一致升高的p-EMT或Flores证据。** 前沿带对内部的总体调整结果如下；前沿对近内部、对深内部另有完整预设比较，而不是看到阴性后更换边界。

{mdtable(front)}

存在局部差异：例如GSE319536_BC01的前沿p-EMT效应为正(调整β约+0.766)，CanDis_2026_Sample6为负(约−0.406)，全片空间图及全部逐片效应均保留。这些例子用于展示异质性，不用于挑选纳入切片。按预设P75定义的3,428个pEMT-high front点标明了分子活跃候选边段，不是独立验证。其EMT/ECM、血管生成和缺氧等背景较高；其中部分程序与p-EMT共享基因，重叠表单列，不能当成多条完全独立证据。

Flores在pEMT-high versus其他前沿的独立支持：{eff('A_active_front_support','FLORES_CRISPR_PRO_INVASION_10')}。未获得一致支持，不设置“通过两个集合”的门槛，也不因此停止MP10模块。

单细胞来源检查在19,053个分层抽样细胞、原注释Epi/CAF/SMC/Myeloid内进行。样本层p-EMT UCell中位数约CAF0.211、SMC0.191、Myeloid0.084、Epi0.054；这说明混合spot p-EMT具有明显的非上皮来源。Flores10在空间肿瘤Epi中的检测基因数切片中位数为4，约19%的点仅检出≤2个成员，最高单基因计数份额中位数约0.40。其阴性不能等同于缺乏侵袭功能；该集合是功能筛选来源的支持集合，并非已验证的表达量侵袭标尺。本轮功能程序评分均为混合spot UCell；MP10始终为冻结Top20 AUCell。

因此本轮保留“正式解剖前沿”“pEMT-high候选边段”的区别，不统一升级为病理或分子确认的侵袭前沿。没有精确距离也不意味着切片没有侵袭前沿。

## 3. MP10 与前沿有什么关系？

- 全部Epi的连续MP10：{eff('C1_all_Epi_continuous',RAW)}。未调整曲线有从内向外降低的总体趋势，组成/深度调整后明显减弱；标准化距离二次项及分箱曲线也已输出，不能只据一个线性斜率作结论。
- 前沿P95比例：{eff('C2_P95_enrichment','P95_probability')}。该效应单位为概率差，即约+0.87个百分点，区间约−0.67至+2.26个百分点，未支持总体前沿富集。分母为各带全部合格肿瘤Epi，不是P95点的距离直方图。
- 仅P95自身的程序–前沿关系：可估计主汇总约23单位/5datasets；包括p-EMT、Flores、增殖、缺氧和应激的完整面板，没有通过预设族FDR的可重复状态梯度。稀疏与覆盖限制保留，不能称为确定无生物学差异。
- pEMT-high front的MP10：{eff('C_active_front_MP10',RAW)}。原始均值偏低，控制上皮比例、深度及前沿内部精确位置后证据减弱，未支持MP10偏向分子活跃边段。

三个问题均独立完成：总体位置、P95前沿富集、P95内部状态不能相互替代。没有把局部CAF阈值接触线当作本轮前沿。

## 4. 九类 niche 中，哪些关系与 MP10 直接相关？

正式R12原始标签0–8全部分析；126,075个已标注点均是冻结Epi queries。各带组成分母为已标注合格中心，同时保存未标注比例；没有混入dataset-local或K12。

R12_2的总体前沿倾向为约+7.2个百分点(FDR约{row('D_front_enrichment','R12_2').q:.3f})，属于位置背景趋势，不是可靠的MP10特异结果。原R1 profile的主要非上皮成分包括B/Plasma和vCAF；这与本轮中心posterior mean来源不同。R12_3在相同前沿位置的P95对照中偏低：匹配概率差约−2.34个百分点(FDR约{row('D_matched_front','R12_3','matched').q:.3f})，其原R1非上皮组成包括Matrix_CAF/vCAF。但九类共同校正后，没有任何niche取得稳健的MP10特异偏好证据。

连续MP10、P95回归和逐对前沿匹配三套结果并列。R12_3的负向趋势可保留为描述性线索，不能据某些区间不跨0忽略全族FDR，亦不能将一般前沿组成变化写成MP10驱动的生态位重塑。全部九类组成、位置效应、活跃前沿分布和缺标签覆盖见完整表/图。

## 5. 哪些结果最值得后续推进？

1. **MP10中心的共享基因独立应激背景**：UV_RESPONSE_UP的小效应及mTORC1去重叠后的剩余趋势最值得核查；UPR/ROS为次级候选。验证应以真实上皮内部、样本重复和独立蛋白/功能读出判断，不能把Hallmark名字直接当作机制。
2. **MP10与p-EMT呈空间状态差别**：中心负关联较一致，向远环减弱；先区分混合组成与上皮细胞内部状态，不提出MP10因果抑制侵袭的结论。
3. **暂不优先提出普遍MP10-high侵袭前沿或固定MP10 niche**：前沿总体富集、P95状态梯度和九类niche特异性均缺少稳定全族支持。局部边段值得结合图册讨论，不替代总体结果。
4. **降低androgen/TGF-β独立机制解释优先级**：去共享基因后信号基本消失，这是本轮较有信息量的敏感性结果，而不是隐藏原阳性结果。

## 实际限制与交付状态

全部63片计算及五组图完成，原输入未修改。GSE238145_B13没有冻结Epi queries；GSE171351_D1只有3个query、1个内侧肿瘤Epi(保留P95并列规则，实际入选100%)，相关推断不足。3张正常来源没有本定义下的肿瘤P95/前沿对象，保留全Epi中心程序和独立空间图。

精确前沿距离覆盖58,556个肿瘤来源Epi queries，不能外推到所有未定位点。前沿匹配仅覆盖1,979/5,664个P95点，片内前沿绝对差中位数的切片中位数约65μm；3/56片匹配后上皮丰度或深度的绝对SMD>0.2，完整匹配质量保留。没有删除这些切片以改善结果。

层级汇总按单位、dataset等权，核实患者仅3datasets；空间方格簇间近似独立和不完整前沿覆盖都是限制。2000μm块用于重叠邻域，空间块<8时只输出片内效应，仍保留可用单位层级结果。没有使用海量spot独立P值，也没有采用最小置换P值导致的FDR分辨率截断。

入口为overview.pdf、figure_index.tsv、atlas/、tables/summary_effects.tsv、methods.md和scripts/run_all.sh。review.zip仅包含审阅所需报告、主图和核心表；原始RNA与完整Seurat/h5ad不复制入包。
'''
(O/'summary.md').write_text(summary)
# Preserve pEMT overlap context for interpretation of active-front background programs.
pegenes=set(G[G.gs_name.eq(PE)].gene_symbol)
save([dict(program=p,pEMT_overlap_n=len(pegenes&set(g.gene_symbol)),overlap='|'.join(sorted(pegenes&set(g.gene_symbol)))) for p,g in G[G.panel.eq('main')].groupby('gs_name')],O/'tables/pEMT_program_overlap.tsv')
# Source figure: one lineage summary per sample, after within-sample subtype medians.
ss=pd.read_csv(O/'tables/sc_sample_scores.tsv',sep='\t').groupby(['sample','lineage','dataset'])[[PE,'FLORES_CRISPR_PRO_INVASION_10']].median().reset_index()
f,axs=plt.subplots(1,2,figsize=(11,4));levels=sorted(ss.lineage.unique())
for ax,p in zip(axs,[PE,'FLORES_CRISPR_PRO_INVASION_10']):
 ax.boxplot([ss.loc[ss.lineage.eq(k),p].dropna() for k in levels],labels=levels,showfliers=False);ax.set_ylabel('Sample-level UCell (subtype medians)');ax.set_title(short(p));ax.tick_params(axis='x',rotation=25)
f.tight_layout();export(f,'01_single_cell_source')
# Preserve machine-readable protocol notes and original 55-gene-set coverage.
shutil.copyfile(F/'tables/gene_set_coverage.tsv',O/'tables/frozen55_gene_set_coverage.tsv')
# Cover page for the vector PDF overview (report prose remains Chinese Markdown).
f=plt.figure(figsize=(11.7,8.3));f.text(.07,.92,'MP10 centers, neighborhoods and formal fronts',fontsize=21,weight='bold');f.text(.07,.85,'63 sections | 265,002 tissue spots | 143,461 frozen Epi queries',fontsize=14)
lines=['57 primary programs + 23 sensitivity programs; original Top20 AUCell unchanged.','5,664 P95 tumor-Epi reference points; 1,979 exact-front-position matched pairs.','R0: modest stress/growth associations; shared-gene removal attenuates several signals.','Androgen and TGF-beta associations largely disappear without MP10 overlap genes.','p-EMT is modestly lower at MP10-high centers; spatial effects weaken with distance.','No consistent global p-EMT / Flores increase at the formal v3.4 front.','No robust P95 front enrichment or MP10-specific R12 preference after family FDR.','Main repeated-unit summary: 54 units / 6 data sets; verified subset: 27 patients / 3 sets.','Normal-source samples and unverified identities remain separate.','All effects, uncertainty, coverage, negative findings and missing modules are retained.']
for i,line in enumerate(lines):f.text(.07,.75-i*.054,line,fontsize=12)
f.text(.07,.10,'Read summary.md and methods.md before biological interpretation.\nPDF graphics are vector; H&E retains its original raster format.',fontsize=11);export(f,'00_reading_guide')
idx=pd.read_csv(O/'figure_index.tsv',sep='\t');new=pd.DataFrame(index);idx=idx[~idx.name.isin(new.name)];idx=pd.concat([idx,new],ignore_index=True).drop_duplicates('name');save(idx,O/'figure_index.tsv')
# Merge vector PDFs without rasterizing their pages.
def combine(files,dest):
 p=O/'logs/merge_paths.txt';p.write_text('\n'.join(str(x) for x in files));subprocess.run(['Rscript','-e','a<-commandArgs(TRUE);qpdf::pdf_combine(readLines(a[1]),output=a[2])',str(p),str(dest)],check=True,stdout=subprocess.DEVNULL)
for ds,g in M.groupby('dataset_id'):combine([O/'figures/sections'/ds/(s+'.pdf') for s in g.section_id],O/'atlas'/(ds+'_atlas.pdf'))
combine([O/'figures/sections'/r.dataset_id/(r.section_id+'.pdf') for r in M[M.source_context.eq('adjacent_normal_source')].itertuples()],O/'atlas/adjacent_normal_reference.pdf')
main_names=['00_reading_guide','00_module_coverage','02_center_rings_adjusted','02_R0_overlap_sensitivity','02_matched_front','02_P95_exclusive_rings_Epi','02_P95_exclusive_rings_tissue','03_P95_distance_effects','03_distance_coverage']+['03_P95_control_curves_'+str(i) for i in range(1,6)]+['01_front_programs_adjusted','01_front_near_deep','01_single_cell_source','01_active_front_support']+['04_front_curves_'+str(i) for i in range(1,6)]+['04_P95_front_proportions','04_MP10_front_effects','05_niche_front_composition','05_niche_effects','05_original_R12_profiles']
reps=M[M.pilot_representative.eq(True)];files=[O/'figures'/(n+'.pdf') for n in main_names]+[O/'figures/sections'/r.dataset_id/(r.section_id+'.pdf') for r in reps.itertuples()];combine(files,O/'overview.pdf')
save([dict(section_id=r.section_id,dataset_id=r.dataset_id,reason='Original pilot representative, fixed before current results') for r in reps.itertuples()],O/'tables/overview_sections.tsv')
# Only the temporary matrix created in this task is removed; original inputs remain untouched.
(O/'scores/sc_counts_temporary.mtx').unlink(missing_ok=True)
small=['summary.md','methods.md','README.md','input_index.tsv','parameters.json','software.txt','full_sections.tsv','figure_index.tsv','overview.pdf']
tables=['summary_effects.tsv','dataset_effects.tsv.gz','unit_effects.tsv.gz','effects.tsv.gz','bands.tsv.gz','niche.tsv.gz','matching_balance.tsv','signature_index.tsv','gene_sets.tsv','new_gene_set_coverage.tsv','literature_symbol_coverage.tsv','signature_detection_summary.tsv','pEMT_program_overlap.tsv','module_coverage.tsv','ring_coverage.tsv','curve_summary.tsv.gz','ring_summary.tsv.gz','front_program_bands.tsv.gz','discussion_candidates.tsv','front_section_effects.tsv','sc_sample_scores.tsv','sc_gene_source.tsv','full_sections.tsv','dataset_module_coverage.tsv','active_front_segments.tsv','niche_palette.tsv','overview_sections.tsv']
with zipfile.ZipFile(O/'review.zip','w',zipfile.ZIP_DEFLATED,compresslevel=6) as zf:
 for n in small:
  p=O/n
  if p.exists():zf.write(p,n)
 for n in tables:
  p=O/'tables'/n
  if p.exists():zf.write(p,'tables/'+n)
 for p in (O/'sources').glob('*.txt'):zf.write(p,'sources/'+p.name)
 for p in (O/'scripts').glob('*'):
  if p.is_file():zf.write(p,'scripts/'+p.name)
 assert zf.testzip() is None
print('Complete',len(idx),'figures;',len(files),'overview pages; review MiB',round((O/'review.zip').stat().st_size/2**20,1),flush=True)
