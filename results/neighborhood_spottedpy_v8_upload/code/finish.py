from common import *
import zipfile,tempfile,re,html,shutil

def mdtable(d):
 cols=list(d.columns);out=['| '+' | '.join(cols)+' |','| '+' | '.join(['---']*len(cols))+' |']
 for row in d.itertuples(index=False,name=None):out.append('| '+' | '.join(str(x).replace('|','/') for x in row)+' |')
 return '\n'.join(out)

def write_reports():
 ev=pd.read_csv(OUT/'07_synthesis/dataset_evidence_all.tsv',sep='\t');cov=pd.read_csv(OUT/'02_hotspots/hotspot_coverage_all.tsv',sep='\t');geo=pd.read_csv(OUT/'01_geometry/coordinate_and_graph_audit.tsv',sep='\t');execution=pd.read_csv(OUT/'02_hotspots/variable_execution_audit.tsv',sep='\t');ss=pd.read_parquet(OUT/'07_synthesis/all_section_effects.parquet');units=pd.read_parquet(OUT/'07_synthesis/all_unit_effects.parquet');validation=json.loads((OUT/'reports/validation_summary.json').read_text());pilot=pd.read_csv(OUT/'00_audit/pilot_checks.tsv',sep='\t')
 src=cov[cov.domain.eq('SOURCE')&cov.k.eq(10)];hot_summary=src.groupby(['dataset_id','variable','track','direction']).agg(n_formal_spots=('n_formal','sum'),n_sections_with_regions=('n_formal',lambda x:int((x>0).sum())),n_regions=('n_components','sum'),n_supported_total=('n_supported','sum')).reset_index();tab('07_synthesis/source_hot_cold_by_dataset.tsv',hot_summary)
 support=ev.groupby(['scope','track','target_kind','metric']).agg(n_planned=('p','size'),n_calculable=('p',lambda x:int(x.notna().sum())),n_q05=('q_scope',lambda x:int((x<.05).sum()))).reset_index();tab('07_synthesis/scope_support_summary.tsv',support)
 scope_units=units[units.k.eq(10)&units.analysis_domain.eq('MAIN')].groupby('scope').analysis_unit_uid.nunique();totalE=int(geo.n_E0.sum());nE=int(geo.n_E0.gt(0).sum());text=f'''# v8 实际执行报告

状态：`V8_COMPLETE_AWAITING_REVIEW`。旧输入保持只读；不自动启动任何后续分析。

## A. 执行身份与完成量

- 复用 v7 score_eligible 共同评分中心：{totalE:,} spots，实际有中心的切片 {nE}/63。63 张切片全部保留，包括无中心的状态记录。
- 冻结两层 MP10（spot_MP10_model_universe / epi_MP10_model_universe）、SAMPLE_R12_LEIDEN niche0–8；细胞目标为完整组织原始 q05/all30 的29种非Epi相对组成。
- 完整组织 {int(geo.n_tissue.sum()):,} spots，其中组成有效 {int(geo.n_valid_composition.sum()):,}；缺组成位置保留在物理组织图中。
- 已执行 63×31×3={len(execution):,} 个变量/切片/k 任务。统计可计算 {int(execution.status.eq('CALCULABLE').sum()):,}，不满足位置数/离散值门槛 {int(execution.status.ne('CALCULABLE').sum()):,}，全部保留状态。
- 主法 SOURCE_G：star=False、R 行标准化、k10；k6/k18全部完成；9999 permutations，seed123，两个切片worker，每个 n_jobs=1。
- 物理像素坐标以完整组织 R1 边中位长度换算 spot pitches；未声称 μm。{int(geo.pixel_jitter_present.sum())} 张切片有可见像素坐标抖动，六方向中位间距审计保留。
- 七dataset pilot全部通过：锁定源码统计逻辑一致、重复种子一致、keep_simulations一致、原生较小尾P及双侧转换、G/G*区别、空目标/跨组件/零距离核验。实际验证 {validation['passed']} 项通过。
- 点级统计 {validation['point_rows']:,} 行，以 section/variable/k parquet 分区存储；完整索引见 point_hotspot_statistics_manifest.tsv。

## B. 热点存在性

下表统计正式 R1≥5 成员区域；零表示没有检出的正式区域，并非计算失败。SOURCE_P05未做局部FDR；LOCAL_BH05使用双侧局部BH，两者不混接。

{mdtable(hot_summary)}

## C. 九niche几何结果

所有9个目标都保留。O为百分比点（正=重叠富集）；D和D_out为pitch（负=更近）。下列总体单元格为 `单位均值效应 / N / q_scope`，N是符合切片覆盖要求后归并的既有分析单位，不全部称为患者。
'''
 def formatval(row,metric):
  if pd.isna(row.mean_effect):return 'NA / '+str(int(row.n_units))+' / NA'
  val=row.mean_effect*(100 if metric=='O' else 1);q=f'{row.q_scope:.3g}' if pd.notna(row.q_scope) else 'NA';return f'{val:+.3g} / {int(row.n_units)} / {q}'
 for kind in ['NICHE','CELL']:
  if kind=='CELL':text+='\n## D. 全29类细胞目标\n\nD为主终点，O和D_out为辅助终点；其各自BH家族分开。细胞目标是相对组成热点，低绝对比例警告未删除任何类型。下表同样保留无可用结果。\n'
  targets=[f'niche{i}' for i in range(9)] if kind=='NICHE' else [TYPES[j] for j in NON]
  for track in TRACKS:
   for layer in LAYERS:
    z=ev[ev.scope.eq('ALL_DATASETS')&ev.target_kind.eq(kind)&ev.track.eq(track)&ev.layer.eq(layer)];rows=[]
    for target in targets:
     row={'target':target}
     for metric in ['O','D','D_out']:rr=z[z.target_id.eq(target)&z.metric.eq(metric)].iloc[0];row[metric]=formatval(rr,metric)
     rows.append(row)
    text+=f'\n### {layer} · {track}\n\n'+mdtable(pd.DataFrame(rows))+'\n'
 text+='''
## E. 稳健性与限制

- 两层公平比较使用同一实际可用切片，随后按原单位等权合并。单层q、配对两层q及方向分别存储，不能把不同N称为复现。
- k6/k18、局部BH、几何局部支持、技术条件分层、组件等权与已核实患者子集均完整输出。尺度敏感性不提供可择优的P。
- 3×3技术分层按固定规则退至2×2；不足时明确不可用，不以单层冒充调整。匹配敏感性比较还提供相同切片/单位集的方向和实际效应。
- Source G不包含中心自身；热点表示周围高值，不要求每个成员本身高分。自身分数分布在空间页展示。
- 单位层Wilcoxon检验是同片背景差的跨单位一致性，不是保留空间自相关的随机位置零模型。方向性最近距离不是细胞接触，也不能推出招募、排斥或机制。
- 模型Epi评分沿用v5来源恢复限制，不能据此证明上皮特异或恶性细胞特异。
- 本次为v7之后的热点几何探索。不得以热点阳性否定v7既有连续环关联的非显著结果，不删除阴性dataset后声称跨队列稳定。

## F. 完整文件

见 [REVIEW_INDEX.md](REVIEW_INDEX.md)、[全部dataset索引](INDEX_ALL_DATASETS.html)、[续跑说明](RESUME.md)。七dataset图包包含全切片、全目标及对应数值表。所有图为透明300dpi PNG及PDF；无距离曲线。主压缩包包含总体图和每dataset一张预先技术选择的overview。
'''
 primary_distances=ev[ev.metric.isin(['D','D_out'])];overall=ev[ev.scope.eq('ALL_DATASETS')];paired=overall[overall.track.eq('LOCAL_BH05')&overall.two_layers_paired_q05_same_direction&overall.metric.eq('O')].drop_duplicates(['target_kind','target_id']);keytext='\n## 主要观察\n\n'
 keytext+=f"主分析（全部既定政策单位）在总体及七dataset、两条轨道中的D/D_out，q_scope<0.05共{int(primary_distances.q_scope.lt(.05).sum())}项。不能据此确认更近、更远或区域外接近。\n\n"
 for tr in TRACKS:
  oo=overall[overall.track.eq(tr)&overall.metric.eq('O')&overall.q_scope.lt(.05)];keytext+=f"{tr}总体重叠终点：niche×评分层{int(oo.target_kind.eq('NICHE').sum())}项、celltype×评分层{int(oo.target_kind.eq('CELL').sum())}项过各自q_scope；这些是重叠差异，不是距离结论。\n\n"
 keytext+='严格LOCAL_BH05下，两层在匹配切片/单位上同向且各自配对家族q<0.05的重叠结果为：'+', '.join(paired.target_id)+'。完整效应、N与q见下表及匹配表。\n\n'
 verified=pd.read_parquet(OUT/'07_synthesis/all_scope_tests.parquet');verified=verified[verified.subset.eq('VERIFIED_PATIENTS_ONLY')&verified.metric.isin(['D','D_out'])&verified.q_scope.lt(.05)]
 for rr in verified.itertuples():keytext+=f"预定已核实患者子集另有{rr.target_id}的{rr.metric}结果：{rr.scope}、{rr.track}、{rr.layer}，效应{rr.mean_effect:+.3f} pitch，N={rr.n_units}，q_scope={rr.q_scope:.5f}。这是子集敏感性，不能替换全部单位的主结果。\n\n"
 text=text.replace('\n## A. 执行身份与完成量',keytext+'## A. 执行身份与完成量',1)
 (OUT/'reports/RUN_SUMMARY.md').write_text(text)
 methods='''# Methods and implementation adaptations

Locked SpottedPy commit: 9dc364f203b7617b6701aa8fb579cbbb35357b26. Original helper/main Git blobs and installed esda/libpysal source/version audit are under 00_audit. This is an explicit project adapter, not a claim of unmodified package reproduction.

Frozen domains: E0 is exactly v7 score_eligible; both source score layers share component-supported KNN axes. Niche targets include all accepted frozen common anchors with that label, not only E0. Cell targets use full-tissue q05 divided by all30 row sums, excluding Epithelial only from target list. Invalid-composition tissue nodes still determine physical connectivity. No score recalculation, expression reassignment, training or clustering occurred.

Geometry: canonical spot_uid ordering; fullres imagecol/imagerow divided by median full-tissue R1 edge length. Verified Visium offsets (0,±2),(±1,±1); ideal hex Euclidean mapping checked independently. Real image coordinates can have point jitter: directional medians and edge distributions are audited, with no smoothing of observed XY. KNN tie order is rounded distance (12 decimals), then canonical UID. Each domain component needs k+1 members. Directed binary KNN converted to R weights; no edges between independent physical tissue components. Randomization background remains the section's supported domain, not separately renormalized per physical component.

SOURCE_G: esda.G_Local(star=False, transform='R', permutations=9999, seed=123, n_jobs=1, keep_simulations=False). Local support needs ≥30 positions and ≥3 distinct values. Native p_sim in esda2.9 is directed smaller-tail conditional randomization with +1 correction; retained simulation audit confirms this. SOURCE_P05 requires native P<.05 and Z direction. LOCAL_BH05 uses min(1,2*p_native), then BH over all calculable local positions per section/variable/k, then q<.05 and Z direction. No separately positive-Z family. pvalue_convention is constant for all point partitions and stored both in every point table and the partition manifest. Original vec_permutations output is memoized by (max_card,n,permutations,seed), without changing draws or statistics; seven pilot source comparisons and an independent cell-domain comparison were bitwise identical. Pilot G* is audit only. Locked source parity was tested by AST-loading its hotspot function with the exact same W and injected G_Local; component detection was bypassed only in that statistics parity test because production intentionally uses induced full-tissue R1 components. Production retains all raw hot/cold points, including 1–4-member components; formal components have ≥5 members.

For each section/layer/k/track, E is the locally calculable source domain and H its formal hot members. O=mean target indicator in H minus E. D=mean nearest target distance in H minus E, restricted to components containing a target; D_out excludes target members from both groups. KDTree retains self-overlap zero and returns NA for unavailable components/targets, with no max-distance fill. Direction is source→target. Point-weighted means are primary, medians and component-equal hot means minus E are descriptive. Nonhot differences and the shrinkage identity Delta=(1-H/E)*(meanH-meanNonhot) are audited on each actual metric domain. Absent niche gives NA even for O, rather than treating obligatory zero as evidence.

Inference: all descriptive effects retained, but eligible O needs E≥30,H≥20,nonhot≥20,target-in-E≥20,rest≥20. D/D_out need actual domain E≥30,H≥20,nonhot≥20,target-full≥5. Equal section means within original analysis_unit_uid; ALL_DATASETS does not split shared units by dataset. Two-sided Wilcoxon signed-rank vs0 needs ≥5 units and ≥5 nonzero effects; allzero with≥5 returns P1. Values within1e-12 of zero are set to zero for testing only and counted. Exact only for n<50 without zeros/absolute ties; otherwise approximation with continuity correction and Wilcox zero handling. The mean is reported as an effect summary, not as the exclusive signed-rank estimand.

Fixed k10 BH families per scope and track: NICHE_PRIMARY36 (9×2 layers×O,D), CELL_PRIMARY_DISTANCE58, NICHE_OUTSIDE_SECONDARY18, CELL_OVERLAP_SECONDARY58, CELL_OUTSIDE_SECONDARY58. Planned n includes missing hypotheses. Seven-dataset and eight-scope joint q are additional sensitivity fields; datasets are not independent replications by assumption. Verified-patient tests have separate same-size families. k6/k18/geometry sensitivity does not produce competing P. Stars on relation figures are q_scope, never local P/q. Full 63/31/3 computation ledger includes unavailable tasks; zero regions are not errors.

Technical sensitivity: Epi fraction and log1p(center_nCount) quantile bins computed on actual E of each metric, including outside domain. Try3×3 with duplicate cuts dropped, then2×2 if any nonempty bin<20; require ≥2 actual nonempty strata; a constant covariate collapses to one level under the duplicate-boundary rule; otherwise unavailable. Hot stratum weights standardize background means. Covariate medians and strata status retained. This is descriptive, not a guarantee of deconfounding. Geometry sensitivity requires all k neighbors within full-tissue R≤3 and kth distance≤3 pitches, restricting H and E without redefining target hotspots. Scale and score-layer comparisons explicitly match actual eligible sections before unit averaging.

Source code: https://github.com/secrierlab/SpottedPy/tree/9dc364f203b7617b6701aa8fb579cbbb35357b26
PySAL reference: https://pysal.org/esda/v2.9.0/generated/esda.G_Local.html
'''
 (OUT/'reports/METHODS.md').write_text(methods)
 (OUT/'reports/WHAT_IS_STILL_UNCERTAIN.md').write_text('''# 尚未解决的问题

所有结果为探索性热点几何关系。尚未实施保留空间自相关的随机位置/空间平移/Moran谱零模型，也未完成额外病理审定、模型表达来源恢复或严格空间排斥验证。原始H&E注册沿用既有坐标，技术可视核验不等于病理确认。匿名政策单位不重命名为患者；共享来源不视为独立队列复制。

G_Local hot并不要求中心分数本身高；相对组成热点不等于绝对细胞数增加。低比例目标的警告、长KNN支持、单个大组件、局部BH无热点和各dataset无可估计结果都保留。两层不一致、技术背景敏感和尺度敏感是证据限制，不能选择其中显著的一套替换主结果。

不执行通讯、LIANA+、模型训练、重聚类、重算MP10或后续机制推断。本轮停止等待人工审阅。
''')
 resume=f'''# 续跑与核验

环境和代码均位于本项目，无需修改旧成果。任意目录执行：

```bash
"{OUT}/code/run_spottedpy_v8.sh" --resume
```

无参数同样按指纹复用 checkpoint，补齐缺失阶段。输入、spec或锁定源码/依赖身份变化会拒绝复用，不自动覆盖旧输入。

仅核验（不会重新进行置换）：

```bash
"{OUT}/code/run_spottedpy_v8.sh" --verify-only
```

运行入口为code/run.py；环境env/bin/python；每次日志保留于logs。--resume需要本服务器完整数据，review包用于离线审阅，不包含表达矩阵或依赖环境。

点级表分布在02_hotspots/<section>/<variable>_k<k>.parquet，包含UID、原始值、G、Z、native P、双侧P、local q、原始/正式hot/cold和组件ID；逻辑清单为02_hotspots/point_hotspot_statistics_manifest.tsv。*_components.parquet保留小区域。计算状态JSON包含输入/spec/锁定源与依赖指纹。

不以修改阈值、减少9999置换、填最大距离、跨组织块或自动挑选显著组合来补结果。完成状态为V8_COMPLETE_AWAITING_REVIEW。
''';(OUT/'reports/RESUME.md').write_text(resume)
 return ev

def build_index():
 pilots=pd.read_csv(OUT/'00_audit/pilot_sections.tsv',sep='\t');lines=['# v8 离线审阅入口','','状态：`V8_COMPLETE_AWAITING_REVIEW`。先读完整报告，再核查全量证据和覆盖。','','- [实际结果与全部目标](RUN_SUMMARY.md)','- [方法](METHODS.md)','- [限制](WHAT_IS_STILL_UNCERTAIN.md)','- [续跑](RESUME.md)','- [自动验证](validation_results.tsv)','- [九niche证据](../07_synthesis/niche_geometry_evidence_all9.tsv)','- [29类细胞证据](../07_synthesis/cell_geometry_evidence_all29.tsv)','- [所有dataset证据](../07_synthesis/dataset_evidence_all.tsv)','- [全部不可推断结果及原因](../07_synthesis/unavailable_results.tsv)','- [点级分区索引](../02_hotspots/point_hotspot_statistics_manifest.tsv)','- [全部热点覆盖](../02_hotspots/hotspot_coverage_all.tsv)','- [两层匹配结果](../05_sensitivity/two_score_layers_matched_units.tsv)','- [匹配敏感性](../05_sensitivity/matched_sensitivity_effects.tsv)','','点级原始parquet留在服务器，不放入轻量包；索引所列相对path是服务器项目根目录相对路径。','','## 总体图','','图上星号为关系检验q_scope。局部热点P/q只用于定义区域。']
 for p in sorted((OUT/'06_figures/ALL_DATASETS').rglob('*.png')):lines.append(f'- [{p.stem}](../{p.relative_to(OUT).as_posix()}) · [PDF](../{p.with_suffix(".pdf").relative_to(OUT).as_posix()})')
 lines+=['','## 七dataset技术选择的空间overview','']
 for r in pilots.itertuples():p=Path('06_figures')/r.dataset_id/'sections'/r.section_id/'MP10_overview.png';lines.append(f'- [{r.dataset_id}：{r.section_id}](../{p.as_posix()}) · [PDF](../{p.with_suffix(".pdf").as_posix()})')
 lines+=['','## 各dataset完整图包','']
 for ds in SCOPES[1:]:lines.append(f'- `figures_{ds}_spottedpy_v8.zip`：解压后打开其 REVIEW_INDEX.md 或 INDEX.html；含该dataset全部切片/29目标及总体关系页。')
 (OUT/'reports/REVIEW_INDEX.md').write_text('\n'.join(lines)+'\n');htmltext=['<!doctype html><meta charset="utf-8"><title>v8 datasets</title><h1>v8 dataset overview</h1><p><a href="REVIEW_INDEX.md">Review index</a> | <a href="RUN_SUMMARY.md">Full report</a></p>']
 for r in pilots.itertuples():p=Path('06_figures')/r.dataset_id/'sections'/r.section_id/'MP10_overview.png';htmltext.append(f'<h2>{html.escape(r.dataset_id)}</h2><p>{html.escape(r.section_id)}; complete package: figures_{r.dataset_id}_spottedpy_v8.zip</p><a href="../{p.as_posix()}"><img loading="lazy" width="900" src="../{p.as_posix()}"></a>')
 (OUT/'reports/INDEX_ALL_DATASETS.html').write_text('\n'.join(htmltext));return pilots

def packages():
 pilots=build_index();manifest=[]
 # Full per-dataset packages preserve paths and include independent offline index.
 for ds in SCOPES[1:]:
  file=OUT/f'figures_{ds}_spottedpy_v8.zip';figs=sorted((OUT/'06_figures'/ds).rglob('*'));figs=[p for p in figs if p.is_file()];md=[f'# {ds} · complete v8 figures','','All sections and targets retained; NA means unavailable, not zero.',''];hx=[f'<!doctype html><meta charset="utf-8"><title>{ds}</title><h1>{ds}</h1><p>All figures: relation q_scope stars; local hot/cold thresholds labeled separately.</p>']
  with zipfile.ZipFile(file,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
   for p in figs:
    rel=p.relative_to(OUT).as_posix();z.write(p,rel)
    if p.suffix=='.png':md.append(f'- [{p.relative_to(OUT/"06_figures"/ds).as_posix()}]({rel}) · [PDF]({p.with_suffix(".pdf").relative_to(OUT).as_posix()})');hx.append(f'<p><a href="{rel}">{html.escape(p.relative_to(OUT/"06_figures"/ds).as_posix())}</a> | <a href="{p.with_suffix(".pdf").relative_to(OUT).as_posix()}">PDF</a></p>')
   for fn in ['dataset_evidence_all.tsv','unavailable_results.tsv','source_hot_cold_by_dataset.tsv']:
    d=pd.read_csv(OUT/'07_synthesis'/fn,sep='\t');col='scope' if 'scope' in d else 'dataset_id';d=d[d[col].eq(ds)];z.writestr('tables/'+fn,d.to_csv(sep='\t',index=False));md.append(f'- [{fn}](tables/{fn})')
   for fn in ['METHODS.md','WHAT_IS_STILL_UNCERTAIN.md']:z.write(OUT/'reports'/fn,'reports/'+fn);md.append(f'- [{fn}](reports/{fn})')
   z.writestr('REVIEW_INDEX.md','\n'.join(md));z.writestr('INDEX.html','\n'.join(hx))
  manifest.append(dict(package=file.name,size=file.stat().st_size,sha256=sha(file),dataset_id=ds));log('packaged '+ds)
 mainfiles=set()
 for folder in ['00_audit','01_geometry','03_niche_geometry','04_cell_hotspot_geometry','05_sensitivity','07_synthesis','reports','code']:
  for p in (OUT/folder).rglob('*'):
   if p.is_file() and p.suffix in ['.tsv','.json','.md','.html','.py','.sh','.txt'] and '__pycache__' not in str(p):mainfiles.add(p)
 mainfiles.add(OUT/'00_audit/pilot_G_vs_Gstar_points.parquet')
 for name in ['hotspot_coverage_all.tsv','hotspot_components.tsv.gz','point_hotspot_statistics_manifest.tsv','variable_execution_audit.tsv']:mainfiles.add(OUT/'02_hotspots'/name)
 for p in (OUT/'06_figures/ALL_DATASETS').rglob('*'):
  if p.is_file():mainfiles.add(p)
 for r in pilots.itertuples():
  for ext in ['png','pdf']:mainfiles.add(OUT/'06_figures'/r.dataset_id/'sections'/r.section_id/('MP10_overview.'+ext))
 mainfiles={p for p in mainfiles if p.relative_to(OUT).as_posix() not in ['reports/package_manifest.tsv','reports/package_validation.tsv']}
 file=OUT/'review_spottedpy_v8.zip'
 with zipfile.ZipFile(file,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for p in sorted(mainfiles):z.write(p,p.relative_to(OUT).as_posix())
  z.writestr('START_HERE.md','# v8 review\n\nOpen [review index](reports/REVIEW_INDEX.md). Full per-dataset figures are separate zip packages.\n')
 manifest.append(dict(package=file.name,size=file.stat().st_size,sha256=sha(file),dataset_id='ALL_DATASETS'));tab('reports/package_manifest.tsv',manifest)

def verify_packages():
 records=[]
 for package in [OUT/'review_spottedpy_v8.zip']+[OUT/f'figures_{ds}_spottedpy_v8.zip' for ds in SCOPES[1:]]:
  with zipfile.ZipFile(package) as z:
   assert z.testzip() is None
   with tempfile.TemporaryDirectory(dir=OUT/'logs') as td:
    z.extractall(td);root=Path(td);nlinks=0
    for p in root.rglob('*'):
     if p.suffix not in ['.md','.html']:continue
     txt=p.read_text();links=re.findall(r'(?:href|src)="([^"]+)"',txt) if p.suffix=='.html' else re.findall(r'\]\(([^)]+)\)',txt)
     for link in links:
      if '://' in link or link.startswith('#'):continue
      target=(p.parent/link.split('#')[0]).resolve();assert target.exists(),(package.name,p.relative_to(root),link);nlinks+=1
    records.append(dict(package=package.name,crc=True,extract=True,relative_links_checked=nlinks,passed=True))
 tab('reports/package_validation.tsv',records);return records

def finish():
 (OUT/'07_synthesis/source_hot_k10_interim.tsv').unlink(missing_ok=True);write_reports();tab('00_audit/implementation_code_manifest.tsv',[input_record(p,'v8 execution/review implementation') for p in sorted((OUT/'code').glob('*')) if p.is_file()]);js('reports/STATUS.json',dict(status='V8_COMPLETE_AWAITING_REVIEW',completed_variables=5859,completed_sections=63,completed_k=KS,tracks=TRACKS,finished_UTC=pd.Timestamp.now(tz='UTC').isoformat(),fingerprint=fingerprint()));packages();verify_packages();
 with zipfile.ZipFile(OUT/'review_spottedpy_v8.zip','a',zipfile.ZIP_DEFLATED) as z:z.write(OUT/'reports/package_validation.tsv','reports/package_validation.tsv')
 pm=pd.read_csv(OUT/'reports/package_manifest.tsv',sep='\t');pm.loc[pm.package.eq('review_spottedpy_v8.zip'),'sha256']=sha(OUT/'review_spottedpy_v8.zip');pm.loc[pm.package.eq('review_spottedpy_v8.zip'),'size']=(OUT/'review_spottedpy_v8.zip').stat().st_size;tab('reports/package_manifest.tsv',pm);checkpoint('V8_COMPLETE_AWAITING_REVIEW');log('V8_COMPLETE_AWAITING_REVIEW')
if __name__=='__main__':finish()
