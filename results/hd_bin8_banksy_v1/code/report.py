from pathlib import Path
import pandas as pd,json,zipfile,os,re
r=Path(__file__).resolve().parents[1]
read=lambda p:pd.read_csv(p,sep='\t')
fmt=lambda x:f'{x:.2%}'
def tab(d):
 cols=list(d.columns);return '| '+' | '.join(cols)+' |\n| '+' | '.join(['---']*len(cols))+' |\n'+'\n'.join('| '+' | '.join(str(v) for v in row)+' |' for row in d.itertuples(index=False,name=None))
diagnostics=[]
for section in ['Cervilla','D1_Johnson']:
 txt=(r/'logs'/f'models_{section}.log').read_text();parts=re.split(r'MODEL_START '+section+r' (l\d+) [^\n]*\n',txt)
 for i in range(1,len(parts),2):
  body=parts[i+1];rec=re.findall(r'Annoy recall = ([\d.]+)%',body);warn=re.findall(r'(\d+) smooth knn distance failures',body)
  if rec:diagnostics.append(dict(section=section,branch=parts[i],UMAP_Annoy_recall_percent=float(rec[-1]),UMAP_knn_calibration_warnings=int(warn[-1]) if warn else 0))
pd.DataFrame(diagnostics).drop_duplicates(['section','branch'],keep='last').to_csv(r/'tables/UMAP_runtime_diagnostics.tsv',sep='\t',index=False)
rows=[];coarse=[];structure=[]
for s in ['Cervilla','D1_Johnson']:
 o=r/s;v=json.loads((o/'validation.json').read_text());c=read(o/'coarsening_comparison.tsv');z=c[c.denominator=='both_candidates'].set_index('level');st=read(o/'spatial_structure_metrics.tsv');ne=json.loads((o/'neighbor_summary.json').read_text());row={'切片':s,'全量 QC bins':v['n'],'原始 UMI':v['UMI'],'原大类分歧 (双方可评价)':fmt(z.loc['original_major','disagree_fraction']),'四粗类分歧 (同分母)':fmt(z.loc['four_coarse','disagree_fraction']),'下降 (百分点)':round(100*(z.loc['original_major','disagree_fraction']-z.loc['four_coarse','disagree_fraction']),2),'λ0 / 0.2 / 0.8 簇数':' / '.join(str(st.set_index('branch').loc[t,'clusters']) for t in ['l00','l02','l08']),'推荐空间层':('λ0.2 主用，λ0.8 粗域对照' if s=='Cervilla' else 'λ0.8 组织域，λ0.2 局部对照'),'状态':'完成并实际加载验证'};rows.append(row)
 c['section']=s;coarse.append(c);st['section']=s;structure.append(st)
pd.DataFrame(rows).to_csv(r/'summary.tsv',sep='\t',index=False);pd.concat(coarse).to_csv(r/'tables/coarsening_all_denominators.tsv',sep='\t',index=False);pd.concat(structure).to_csv(r/'tables/spatial_structure_summary.tsv',sep='\t',index=False)
lines=['# 8 μm bins BANKSY 匹配对照','',tab(pd.DataFrame(rows)),'','两片均保留原有全部 QC bins、18,085 个输入基因及原 RNA counts/data。三分支为全量直接 PCA/聚类，非 sketch 投影。旧 cluster_projected、逐 bin SingleR、分割、归属和 RCTD 均未重跑或覆盖。区域解释与候选类型分层保存，未据此正式修订细胞身份。','', '## 主要结论','']
notes=r/'reviewed_conclusions.md';lines.append(notes.read_text() if notes.exists() else '详细空间复核结论尚未写入；不能将本状态视为最终交付。')
lines+=['','## 分母与粗分类','', '主表分母是同片两套参考均有候选的 bins。全 QC 和原 HPCA 覆盖类别限定的分母另存 `tables/coarsening_all_denominators.tsv`，不可混用。原生大类→四粗类只是解释粒度改变；Mural+Fibroblast 为 Mesenchymal，Mast/Plasma 等归 Immune 但原字段保留。HPCA 缺少 Mural/Mast/Plasma 的独立覆盖，不构成对这些类型不存在的证据。Unknown 是缺失/不可评价的候选；低信息量区域、未测试 RCTD 和没有输入覆盖另行区分。SingleR pruning 保留不等于高置信度。','']
for s in ['Cervilla','D1_Johnson']:
 o=r/s;z=read(o/'coarse_composition.tsv');lines+=['### '+s,''];# exact schemas are retained as linked tables; avoid silently selecting different denominators.
 lines+=['独立原始/剪枝后粗类的 bin 数、比例及 UMI 占比：`'+s+'/coarse_composition.tsv`；交叉表：`'+s+'/coarse_crosstabs.tsv`。','']
 st=read(o/'spatial_structure_metrics.tsv');st['adjacent_agreement']=st.adjacent_agreement.map(fmt);lines+=[tab(st),'','相邻网格同簇率和碎片数是描述性的空间连续性指标；BANKSY 本来就鼓励连续，不能把这些值解释为准确率。旧 sketch 聚类与新匹配 λ0 的实现/规模/resolution 均不同，旧→新差异不能全部归因于空间信息。','']
 ne=json.loads((o/'neighbor_summary.json').read_text());lines +=[f"空间邻域：{ne['components']} 个 supplied in-tissue 连通块；跨块边为 0；{ne['n_small_neighbor']} 个 bins 因所在块较小使用少于 50 个邻居。第 50/实际最远邻居的中位距离 {ne['actual_max_distance_quantiles'][1]:.2f} μm，99% 分位 {ne['actual_max_distance_quantiles'][4]:.2f} μm，最大 {ne['actual_max_distance_quantiles'][-1]:.2f} μm。块内孔洞/稀疏边缘仍可能被跨过，逐 bin 距离已交付；该尾部不是新 QC 删除标准。分组只限制邻接，不改变全片 scaling。",'']
lines+=['## 证据与文件','', '- `cluster_summary.tsv`、`cluster_reference_composition.tsv`：所有历史/三路簇的 counts、信息量、marker 检出、参考组成及边缘关系。','- `cluster_markers.tsv.gz`：全原始基因表达上的描述性效应/检出差，非只用 HVG；`panel_marker_statistics.tsv.gz` 保留成组 marker 均值和比例。','- `cluster_SingleR.tsv`、`cluster_SingleR_scores.tsv`：λ0.2/0.8 的原始表达聚合对照，保留两参考分歧；没有给整簇强行补同一种 bin 身份。','- `cluster_domain_interpretation.tsv`、`reviewed_region_overrides.json`：原 cluster ID、证据、暂定区域名及是否经过重点图像复核。未重点复核的小簇仍标为摘要级解释。','- `conflict_information_summary.tsv`、`conflict_cluster_distribution.tsv`、`conflict_ROI_summary.tsv`、`ROI_cluster_summary.tsv`：D1 正反方向分歧及一致区的固定分母对照。','- `fixed_rois.tsv`：沿用 ROI；D1 新增冲突热点的选择依据和坐标在 `additional_roi_selection.json`。新增热点未重跑 RCTD；原 RCTD 与 BLCA 共享参考，不是独立真值。','- `metadata.tsv.gz` 与 `l00/l02/l08_embeddings.rds`：增量结果，三路 ID 与坐标严格对齐；不含复制的表达矩阵。原详细候选、scores、marker 支持继续从源 Seurat 获取。','', '主要图：`figures/spatial_clusters_historical`、`spatial_coarse_references`、`spatial_clusters_l00/l02/l08`、`spatial_interpreted_three_lambda`、`spatial_recommended_regions_HE`、`ROI_types_*`、`ROI_clusters_*`、`ROI_raw_markers_*`；D1 另有 `D1_reference_conflict_location`、`D1_conflict_region_composition` 与 `D1_conflict_raw_marker_evidence`。图内无大标题；H&E/高密度点层为栅格，文字图例为矢量；白底 PNG 300 dpi，含图例图均有全图/主体/图例 PDF。','', '## 加载与复跑','', '依赖原对象：`../hd_bin8_comparison_v1/Cervilla/seurat.rds`、`../hd_bin8_comparison_v1/D1_Johnson/seurat.rds`；外部图像继续由原加载器解析。加载不会保存重复的完整 Seurat。','', '```r',"root <- '/home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/analysis/hd_bin8_banksy_v1'","source(file.path(root, 'code/load_result.R'))","x <- load_banksy(file.path(root, 'D1_Johnson'))","table(x$banksy_l08_region)","table(x$coarse_blca, x$coarse_hpca)","SeuratObject::Embeddings(x, 'banksy_l02_umap')","SeuratObject::LayerData(x, assay='RNA', layer='counts')",'```','', '`bash rerun.sh` 顺序复用已完成分支，检查固定参数，更新证据汇总和最后报告。软件版本在各片 `model_versions.json`、`sessionInfo.txt`；完整设置及与官方示例的区别见 `METHODS.md`、`parameters.json`。只进行一次最终加载验证，结果在各片 `validation.json`，并从实际重读对象绘制 `figures/reloaded_attachment`。','', '本轮未扩展 MP10、精细亚群、正式邻域生态、通讯、整合或细胞身份修订。后续先评估本轮组织区域与候选证据；争议区域仍需病理/细胞边界层面的针对性核验。','']
# Lightweight review contains selected PNGs; complete PDF/legend/panel sets remain on server.
common_core={'spatial_clusters_historical','spatial_coarse_references','spatial_interpreted_three_lambda','spatial_recommended_regions_HE','markers_l02'}
core_by_section={
 'Cervilla':common_core|{'ROI_types_far_UMI_hotspot','ROI_clusters_small_nuclei','ROI_raw_markers_far_UMI_hotspot'},
 'D1_Johnson':common_core|{'D1_reference_conflict_location','D1_conflict_region_composition','D1_conflict_raw_marker_evidence','ROI_types_conflict_focus','ROI_clusters_conflict_focus','ROI_clusters_small_nuclei'}
}
def files():
 for p in r.rglob('*'):
  if not p.is_file():continue
  rel=p.relative_to(r)
  if rel.parts[0] in ['tmp','logs','Rlib','sources']:continue
  if p.name in ['review.zip','metadata.tsv.gz','type_metadata.tsv.gz','spatial_metadata.tsv.gz','neighbor_statistics.tsv.gz','coordinates.tsv.gz']:continue
  if 'tmp' in rel.parts or p.suffix=='.rds':continue
  if 'figures' in rel.parts and not (p.suffix=='.png' and p.stem in core_by_section.get(rel.parts[0],set())):continue
  if rel.parts[0] in ['Cervilla','D1_Johnson'] and 'roi_assets' in rel.parts:continue
  yield p
guide=[]
for section in ['Cervilla','D1_Johnson']:
 for p in sorted((r/section/'figures').glob('*.png')):
  stem=p.stem
  meaning=('Same-coordinate H&E, BLCA and HPCA coarse candidates; no consensus' if stem.startswith('ROI_types_') else 'Historical sketch and matched lambda branches; numeric IDs are branch-specific' if stem.startswith('ROI_clusters_') else 'Unsmoothed original RNA panel counts; square marks represent bin centers' if stem.startswith('ROI_raw_markers_') else 'Reused previous ROI RCTD first candidate; not independent truth or a new fit' if stem.startswith('ROI_RCTD_') else 'Original RNA cluster marker mean and detection; largest 25 clusters shown, all in tables' if stem.startswith('markers_') else 'Reviewed cluster interpretation, separate from candidate bin identity' if 'interpreted' in stem or 'recommended' in stem else 'Full-QC result; see filename and REPORT for branch/reference; no plotting subsample')
  guide.append(dict(section=section,file=str(p.relative_to(r)),description=meaning,format='white 300 dpi PNG; PDF + legend/panel versions where applicable'))
pd.DataFrame(guide).to_csv(r/'tables/figure_guide.tsv',sep='\t',index=False)
lines += ['review.zip 精简收录关键 PNG、必要表格和代码；完整 PDF（含主体/图例版本）及其余图件保留在服务器结果目录。', '']
report=r/'REPORT.md';report.write_text('\n'.join(lines)+'\n')
with zipfile.ZipFile(r/'review.zip','w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in files():z.write(p,p.relative_to(r))
size=sum(p.stat().st_size for p in r.rglob('*') if p.is_file());archive=(r/'review.zip').stat().st_size
report.write_text(report.read_text()+f'\n本轮目录实际新增文件总量（含 review.zip）：{size/1024**3:.3f} GiB；review.zip：{archive/1024**2:.1f} MiB。未复制 counts、完整 Seurat/h5ad、参考及全片图像。\n')
# Replace the report member without duplicate archive records.
with zipfile.ZipFile(r/'review.zip','w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in files():z.write(p,p.relative_to(r))
print('REPORT_AND_REVIEW_COMPLETE',size,archive)
