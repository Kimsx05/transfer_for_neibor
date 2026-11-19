from base import *
from source import kinds
from zipfile import ZipFile,ZIP_DEFLATED
import sklearn,scipy,anndata,platform,joblib
from PIL import Image
active=rd(R/'unit_comparison.tsv');summary=active[active.roi=='all_unique'];validation=rd(R/'reference_validation.tsv');hold=validation[validation.truth=='ALL'];nuc=rd(R/'cell_assignment_audit.tsv');sources=rd(R/'FOXP3_capture_localization.tsv');cov=pd.read_parquet(R/'source_overlap_and_coverage.parquet');manifest=rd(R/'input_manifest.tsv')
qc=[];grows=[];depth=[];history=[];part=[]
for s in sections:
 root=R/s;src=pd.read_parquet(root/'sources.parquet');genes=rd(root/'genes.tsv');assert set(T+epi+core+fib+mye).issubset(set(genes.symbol))
 pc=rd(root/'count_accounting.tsv');assert pc.max_per_gene_residual.eq(0).all()
 assert (pc.object_UMI+pc.excluded_UMI+pc.outside_selected_objects_or_unassigned_UMI==pc.all_UMI).all()
 frozen=rd(R/'mp10_frozen_check.tsv').query('section==@s')
 for col in [c for c in frozen if c.endswith('_exact') or c=='eligible_matches_v3_presence']:
  assert frozen[col].dropna().all()
 rctdfile=A/'hd_bin8_comparison_v1'/s/'rctd_ROI_results.tsv.gz';rc=rd(rctdfile);hist=pd.read_parquet(root/'historical_native8_context.parquet')
 rc=rc.drop_duplicates('bin_id').set_index('bin_id')
 for c in ['rctd_status','spot_class','first_type','second_type','weight_1','weight_2']:hist['old_BLCA_RCTD_'+c]=hist.bin_id.map(rc[c])
 hist.to_parquet(root/'historical_native8_context.parquet',index=False)
 for f,kind in [(rctdfile,'historical bin8 BLCA RCTD'),(A/'hd_annotation_diagnosis_singler_v1'/s/'diagnosis.tsv.gz','historical cell SingleR'),(Z/s/'T_regions.parquet','historical strict and FOXP3 T32'),(Z/s/'Epi_regions.parquet','v3 scored Epi64'),(Z/'tables/roi_manifest.tsv','ROI source'),(F/'tables/mixture_parent_membership.parquet','frozen mixture parent IDs'),(F/'tables/expanded_mixtures.parquet','frozen simulation conditions')]:
  history.append(dict(section=s,path=str(f),readable=f.is_file(),size=f.stat().st_size,counts=kind))
 prior=sum(hist.old_BLCA_RCTD_rctd_status.notna())
 for k in kinds:
  e=pd.read_parquet(root/(k+'_evidence.parquet'));m=pd.read_parquet(root/(k+'_members.parquet'));x=sparse.load_npz(root/(k+'.npz'));cc=cov[(cov.section==s)&(cov.branch==k)]
  assert e.id.is_unique and m.barcode.is_unique and x.shape==(len(e),len(genes)) and np.array_equal(np.asarray(x.sum(1)).ravel(),e.nCount)
  assert len(cc)==2*len(e)and not cc.duplicated(['id','radius_um']).any()and set(cc.id)==set(e.id)
  if k.startswith('T'):
   f=4 if k=='T8_from2'else 8;grp=(src.array_row//f).astype(str)+'_'+(src.array_col//f).astype(str);ag=src.groupby(grp).agg(supplied=('barcode','size'),confirmed_outside=('in_tissue',lambda u:(u==0).sum()),supplied_zero=('raw_counts',lambda u:(u==0).sum()))
   tab=e[['section','branch','id','n_supplied_usable','nCount']].join(ag,on='id');tab['nominal_positions']=f*f;tab['unprovided_positions']=f*f-tab.supplied;assert tab.unprovided_positions.min()>=0;part.append(tab)
  qc.append(dict(section=s,branch=k,n_units=len(e),unique_source_bins=len(m),features=len(genes),UMI=int(x.sum()),paired_Epi_ids_unique=True,verified_count_residual=0,legacy_BLCA_RCTD_context_bins=prior))
  covered=e[e.model_candidate|e.FOXP3_raw_screen].copy();gg=covered.merge(cc[['id','radius_um','eligible_epi_before','shared_excluded','independent_epi','nearest_independent_um']],on='id')
  grows.append(gg)
 for q in rd(R/'roi_manifest.tsv').query('section==@s').itertuples():
  im=R/'figures'/f'{s}_{q.roi}.png'
  with Image.open(im) as img:img.verify()
put(pd.concat([manifest,pd.DataFrame(history)]).drop_duplicates(['section','path'],keep='last'),'input_manifest.tsv')
put(pd.DataFrame(qc),'validation.tsv');pd.concat(part).to_parquet(R/'grid_input_support.parquet',index=False);pd.concat(grows).to_parquet(R/'candidate_evidence_and_opportunity.parquet',index=False)
model=joblib.load(R/'model.joblib');assert not set(model['features'])&set(exclude)
put(pd.DataFrame([dict(python=platform.python_version(),sklearn=sklearn.__version__,scipy=scipy.__version__,anndata=anndata.__version__,numpy=np.__version__,pandas=pd.__version__,threads_per_compute=1,allowed_cores=4)]),'versions.tsv')
def md(df):
 return '\n'.join(['| '+' | '.join(df.columns)+' |','|'+'---|'*df.shape[1]]+['| '+' | '.join(map(str,row))+' |' for row in df.itertuples(index=False,name=None)])
table=summary[['section','branch','units','median_nCount','FOXP3_and_T','model_candidates','supported']].copy();table.columns=['切片','单元','去重数','中位 UMI','FOXP3 与 T 共检出','原始模型候选','达到参考支持等级']
report='''# Treg 身份与空间单元比较

Cervilla 的 4 个和 D1_Johnson 的 6 个预定 ROI 已分别完成全部五路分析。D1 值得做针对性的细胞位置验证；目前两片均宜保留区域表达主线。核 ID 和来源映射澄清了 RNA 被分给哪个对象，证据仍不足以据此获得可靠的全片 Treg 细胞数。

'''+md(table)+'''

五个分支使用同一冻结分类器。每行独立去重，横向不能相加。FOXP3 与 T 共检出指同一单元内两类原始计数均非零。原始模型候选为多类别第一预测。唯一达到参考支持等级的是 D1 的 8 μm 格 304_622，它同时落在 candidate_4 与 control_1，表中只计一次。该格 233 UMI、3 个 T 面板基因共 3 UMI、IL2RA 3 UMI，FOXP3、CTLA4 和 IKZF2 均未检出。其模型支持 0.984，次选为上皮。该处支持 T 相关区域进一步复核，不能定位为一个已确认的 Treg。

## 检测量与亚型混淆

30 类模型复用 18,289 个训练对象、2,157 个留出对象。训练内选 C=0.1、候选阈值 0.95。原样本留出含 12 个样本、91 个 Treg 和 2,066 个其他类型。原深度得到 45/91 真阳性，500 UMI 为 23/91，100 和 25 UMI 均为 0/91；四种条件下均为 0/2,066 假阳性。前两种精度分别为 45/45、23/23，后两种无阳性调用，精度不可估。100 UMI 若只看首选类别，91 个真实 Treg 中有 36 个被选中，91 个 Tfh/CXCL13-like 中有 12 个也被选中，普通 CD4 为 8/96，耗竭 CD8 为 6/90。低检测量和竞争状态均限制解释。

70% 上皮 RNA、500 UMI 的留出模拟中，总 RNA 的 Treg 比例为 0%、3%、10%、30% 时，各 30 个混合物达到支持等级的数量依次为 0、0、0、1。模拟重复和复用的留出样本不能计作独立 HD 真值。混合格分类支持低时，无法区分真实缺乏与漏检。

## 同格历史基线

固定评价范围内，native8 与 supplied2 重聚合的每个 8 μm 父格均匹配。Cervilla 3,200 格共 508,127 UMI，D1 5,760 格共 1,257,112 UMI，两片在这些父格上逐格总 UMI 均相同。该核对范围仅限本轮 ROI。旧 T32 在同一评价范围内分别有 44 和 285 个满足原 T 资格，FOXP3-supported 为 9 和 63，strict 为 2 和 11；旧 T16 与 native8 对照另见 historical_unit_baselines.tsv。这里的 T32 资格与新小单元的模型支持等级有不同分母。

## 核内和核周表达来源

Cervilla 固定视野内有 20 个 FOXP3 阳性 2 μm 捕获格，0 个在同格检出 T 面板。D1 为 81 个，其中 5 个同格检出。以每个 FOXP3 格为中心查看 12 μm 内的 T 阳性捕获格，4 μm 规则下能找到同核归属 T 信号的格数分别为 5/20、50/81，而附近还存在其他已分配核的 T 信号的格数分别为 14/20、76/81。这两个条件可以重叠。核中心在 ROI 外仍可贡献 ROI 内的捕获格，反之亦然，因此此处分母与逐核表不同。

两片同核扩张均增加 FOXP3 与 T 共检出。D1 的 57 个 8 μm 外扩共检出对象中，47 个也有上皮或胶原面板 RNA，14 个从 4 扩到 8 μm 时胶原计数进一步增加。Cervilla 的 7 个中有 5 个共检出上皮或胶原。增加的转录本使区域信号更易观测，也增加跨核归属歧义。核内检出少、参考验证显示明显的深度依赖；现有证据不能把漏核列为低识别率的主要原因。没有独立膜、蛋白或细胞级真值，本轮没有依据更换分割算法。

## 空间与下一步

两片原始阵列父格中心与 native8 坐标差的第 99 百分位分别为 0.0044、0.0187 μm。两片 Epi64 成员、资格、坐标、MP10 及 H05/H10 直接沿用，逐 ID 核对见 mp10_frozen_check.tsv。Cervilla 保守 4 μm 的 1,958 个核中，在同组件排除共享父格后有 175 个在 50 μm 内有合格 Epi64，200 μm 内为 1,421 个；D1 分别为 165/4,166、1,224/4,166。此处按完整 Epi64 质心计机会，并保留来源成员。缺少独立邻居表示此处无可测机会。

Cervilla 两个邻接候选属于同一热点。candidate_1 中 T 信号有局部聚集，FOXP3 捕获格稀少，核级证据不足；继续保留 T32 的区域表述。D1 的 candidate_4 中 T 信号覆盖较广，FOXP3、IL2RA、CTLA4 分布在多个捕获格。该处优先做局部多重标记或人工核边界核验，以邻近的 control_1 重叠区作为同一位置处理。核内/4 μm 配对适于查来源；现有结果不足以支持扩展全片 cell 分类。

相比已有 T16、4 μm 和 SingleR，新增内容为同一 supplied 2 μm 输入上的五路配对、30 类 L2 模型的样本分组验证与外推、逐核/逐捕获格审计、来源排除后的 Epi64 机会及固定图册。冻结的 MP10 没有重算。研究停止于局部单位比较，未运行正式富集或通讯。

## 读取

主图册为 figures/ROI_atlas.pdf；各 ROI 的完整 PNG、PDF、主体与图例存于 figures。unit_comparison.tsv 给出去重合计和逐 ROI 结果，cell_assignment_audit.tsv 以核为行；candidate_evidence_and_opportunity.parquet 保留所有 FOXP3 筛查及模型候选，按半径一行。全体单位的分数、次选、原始计数和来源可用 code/load_result.py 读取。网格实际覆盖见 grid_input_support.parquet。单个捕获格只能定位到格内。历史 RCTD 的覆盖差异见 validation.tsv。

'''
(R/'REPORT.md').write_text(report)
