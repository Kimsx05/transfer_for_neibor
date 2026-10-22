#!/usr/bin/env python3
"""Create asset manifest, readiness decisions, lightweight report and integrity checks."""
from pathlib import Path
import os,json,re,hashlib,csv,subprocess
import pandas as pd,numpy as np,h5py
from PIL import Image
ROOT=Path(__file__).resolve().parents[1]
from audit_inputs import specifications,sha
GEO='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='
ZEN='https://zenodo.org/records/17999961/files/visiumhd_Bladder.zip'
FTP='https://ftp.ncbi.nlm.nih.gov/geo/samples/'
SAMPLE={'Cervilla':('visiumhd_Bladder','visiumhd_Bladder'),'GSE300692':('GSM9066380','D1_Johnson'),'GSE299573':('GSM9040999','UC_SCC')}

def write_manifest():
 rows=[];images=[];checksum_cache=ROOT/'metadata/asset_checksums.json';cache=json.loads(checksum_cache.read_text()) if checksum_cache.exists() else {}
 z={x['name']:x for x in json.loads((ROOT/'metadata/cervilla_zip_members.json').read_text())}
 tar={x['name']:x for x in json.loads((ROOT/'metadata/thor_archive_members.json').read_text())}
 transfers=[json.loads(l) for l in (ROOT/'logs/transfers.jsonl').read_text().splitlines()]
 downloaded={x['local']:x for x in transfers if 'local' in x}
 for p in sorted((ROOT/'data').rglob('*')):
  if not p.is_file() or p.is_symlink():continue
  rel=p.relative_to(ROOT/'data');dataset=rel.parts[0]
  sample,section=SAMPLE.get(dataset,('GSM9505917','MPBC1') if '9505917' in p.name else ('GSM9505918','MPBC2'))
  member='/'.join(rel.parts[1:]);expected='unknown';source='';integrity='local structural read + SHA256; no remote checksum supplied';status='reused_moved'
  if dataset=='Cervilla':source=ZEN;expected=z.get(member,{}).get('bytes','unknown');status='new_download';integrity='HTTP Content-Range + ZIP member CRC32 + local SHA256; whole ZIP MD5 NOT verified'
  elif dataset=='GSE300692':
   source=FTP+'GSM9066nnn/GSM9066380/suppl/GSM9066380_D1_Johnson.tar.gz';expected=tar.get(member,{}).get('bytes','unknown')
   if '/outs/spatial/' in member:status='recovered_from_local_archive'
   integrity='archive gzip CRC verified to EOF; member size + local SHA256'
  elif dataset=='GSE299573':
   if p.name.endswith('.h5'):source=FTP+'GSM9040nnn/GSM9040999/suppl/GSM9040999_square_008um_filtered_feature_bc_matrix.h5';status='new_download';expected=int(json.loads((ROOT/'metadata/geo_h5_headers.json').read_text()).get('Content-Length',p.stat().st_size))
   else:source=FTP+'GSM9040nnn/GSM9040999/suppl/GSM9040999_square_008um_spatial.tar.gz';member='spatial/'+p.name
  else:
   source=FTP+'GSM9505nnn/'+sample+'/suppl/'+p.name
   filelist=pd.read_csv(ROOT/'metadata/geo_318916_filelist.txt',sep='\t');matches=filelist[filelist.Name==p.name];expected=int(matches.Size.iloc[0]) if len(matches) else 'unknown'
  role='metadata';version='unknown';note=''
  if p.suffix=='.h5':
   role='filtered_count_matrix' if 'filtered_feature' in p.name else 'raw_probe_matrix' if 'probe' in p.name else 'raw_count_matrix'
   with h5py.File(p,'r') as h:
    version=str(h.attrs.get('software_version','unknown'));note='matrix shape='+str(h['matrix/shape'][:].tolist()) if 'matrix/shape' in h else 'header readable'
  elif p.suffix=='.parquet':role='spatial_positions'
  elif p.name=='scalefactors_json.json':role='spatial_scale';json.loads(p.read_text())
  elif p.suffix.lower() in ['.png','.jpg','.tif','.tiff','.jpeg']:
   role='cytassist_image' if 'cytassist' in p.name else 'tissue_hires_image' if 'hires' in p.name else 'tissue_lowres_image' if 'lowres' in p.name else 'alignment_preview'
   with Image.open(p) as im:
    rec={'dataset':dataset,'sample':sample,'section':section,'path':str(p),'role':role,'width':im.width,'height':im.height,'format':im.format,'mode':im.mode,'physical_scale_source':'section scalefactors for tissue_hires; unknown for CytAssist/alignment previews'};images.append(rec);note=f'{im.width}x{im.height}; {im.format}; original full-resolution H&E not established';im.verify()
  elif p.name=='web_summary.html':role='web_summary'
  elif p.name=='probe_set.csv':role='probe_set'
  elif 'metrics' in p.name:role='metrics'
  fingerprint=f'{p.stat().st_size}:{p.stat().st_mtime_ns}'
  if cache.get(str(p),{}).get('fingerprint')==fingerprint:digest=cache[str(p)]['sha256']
  else:digest=sha(p);cache[str(p)]={'fingerprint':fingerprint,'sha256':digest}
  rows.append(dict(dataset=dataset,sample=sample,section=section,role=role,source_version=version,remote_url=source,archive_member=member,local_path=str(p),expected_bytes=expected,actual_bytes=p.stat().st_size,sha256=digest,remote_checksum=('CRC32:'+z[member]['crc32']) if dataset=='Cervilla' and member in z else 'not supplied',status=status,integrity=integrity,notes=note))
 checksum_cache.write_text(json.dumps(cache,indent=2));(ROOT/'metadata/image_inventory.json').write_text(json.dumps(images,indent=2))
 archive=ROOT.parent/'highres_ST/highres_spatial_audit_v3_GSE300692/00_input_audit/extracted_source/GSM9066380_D1_Johnson.tar.gz'
 rows.append(dict(dataset='GSE300692',sample='GSM9066380',section='D1_Johnson',role='source_archive',remote_url=FTP+'GSM9066nnn/GSM9066380/suppl/'+archive.name,local_path=str(archive),actual_bytes=archive.stat().st_size,expected_bytes=27503110814,status='reused_in_place',integrity='gzip CRC verified through EOF in this run',notes='27.503 GB / 25.614 GiB; not copied or deleted; BAM and large unused members not extracted'))
 rows.append(dict(dataset='Cervilla',sample='visiumhd_Bladder',section='visiumhd_Bladder',role='remote_source_archive',remote_url=ZEN,expected_bytes=3362702746,actual_bytes=0,remote_checksum='MD5:e3f5e6dacae3b4c53f8007dabc64d9a6',status='remote_indexed_only',integrity='whole archive MD5 NOT verified; selected members validated separately'))
 for name,role in [('visiumhd_Bladder/probe_set.csv','optional_probe_set'),('visiumhd_Bladder/binned_outputs/square_008um/filtered_feature_bc_matrix.h5','optional_8um_count_matrix'),('visiumhd_Bladder/binned_outputs/square_016um/filtered_feature_bc_matrix.h5','optional_16um_count_matrix')]:
  if name in z and not (ROOT/'data/Cervilla'/name).exists():
   rows.append(dict(dataset='Cervilla',sample='visiumhd_Bladder',section='visiumhd_Bladder',role=role,remote_url=ZEN,archive_member=name,expected_bytes=z[name]['bytes'],actual_bytes=0,remote_checksum='CRC32:'+z[name]['crc32'],status='not_downloaded',notes='probe file optional request hit TLS EOF; version already recovered from web summary' if role=='optional_probe_set' else 'indexed only to minimize storage; not claimed as locally available resolution'))
 for dataset,sample,section in [('Cervilla','visiumhd_Bladder','visiumhd_Bladder'),('GSE300692','GSM9066380','D1_Johnson'),('GSE299573','GSM9040999','UC_SCC'),('GSE318916','GSM9505917','MPBC1'),('GSE318916','GSM9505918','MPBC2')]:
  for role in ['original_fullres_HE','cell_matrix','nucleus_or_cell_boundaries','bin_to_cell_map']+(['spatial_positions','scalefactors'] if dataset=='GSE318916' else []):
   rows.append(dict(dataset=dataset,sample=sample,section=section,role=role,remote_url=GEO+sample if sample.startswith('GSM') else 'https://zenodo.org/records/17999961',status='not_located_in_checked_sources',actual_bytes=0,notes='absence limited to local inventory, current official files and checked author publications/repositories; not a universal nonexistence claim'))
 pd.DataFrame(rows).to_csv(ROOT/'asset_manifest.tsv',sep='\t',index=False)
 return pd.DataFrame(rows),images

def repair_legacy_links():
 old=ROOT.parent/'highres_ST/VisiumHD_MP10_GSE299573_independent_validation_v1/02_resolution_construction_audit/input_binned_outputs/square_008um'
 for rel in ['filtered_feature_bc_matrix.h5','spatial/tissue_positions.parquet','spatial/scalefactors_json.json','spatial/tissue_hires_image.png']:
  p=old/rel;target=ROOT/'data/GSE299573/GSM9040999/square_008um'/rel
  if p.is_symlink() and not p.exists() and target.exists():
   original=os.readlink(p);p.unlink();p.symlink_to(target)
   with (ROOT/'logs/link_repairs.jsonl').open('a') as f:f.write(json.dumps({'path':str(p),'old_target':original,'new_target':str(target)})+'\n')

def finalize():
 repair_legacy_links();assets,images=write_manifest();a=pd.read_csv(ROOT/'matrix_spatial_audit.tsv',sep='\t');g=pd.read_csv(ROOT/'mp10_feature_detection.tsv',sep='\t')
 # Explicitly distinguish observed header versions from manuscript/GEO descriptions.
 a['reference']='unknown';a['probe_set_version']='unknown';a['reported_software_version']='unknown';a['patient_relationship']='unknown across datasets; no adjacent-section inference'
 for dataset,ref,probe,version in [('Cervilla','GRCh38-2020-A','Visium Human Transcriptome Probe Set v2.0','web summary 3.0.0'),('GSE300692','GRCh38-2020-A','Visium Human Transcriptome Probe Set v2.0','GEO 3.1.0; actual H5 and web summary 3.0.1'),('GSE299573','GRCh38-2020-A','unknown','GEO 3.1.1'),('GSE318916','hg38 (GEO); GRCh38 (paper)','unknown','paper 3.1.3; H5 MPBC1=3.0.0, MPBC2=4.0.1')]:
  mask=a.dataset==dataset;a.loc[mask,['reference','probe_set_version','reported_software_version']]=[ref,probe,version]
 # Paper explicitly identifies MPBC1/2 as two additional MIBC cases, despite HGT1 study title.
 a['pathology_context']='unknown';a.loc[a.dataset=='GSE318916','pathology_context']='MIBC mixed micropapillary/urothelial; primary paper PMC13232474';a.loc[a.dataset=='GSE300692','pathology_context']='MIBC; GSM9066380 extraction protocol';a.loc[a.dataset=='GSE299573','pathology_context']='urothelial carcinoma with squamous differentiation; exact stage unknown'
 a['in_tissue_bins_outside_hires']=pd.NA
 for idx,row in a.iterrows():
  if pd.notna(row.get('in_tissue_image_bounds_rate')):
   a.loc[idx,'in_tissue_bins_outside_hires']=int(round(row.matrix_in_tissue_bins*(1-row.in_tissue_image_bounds_rate)))
 a.to_csv(ROOT/'matrix_spatial_audit.tsv',sep='\t',index=False)
 ready=[]
 for dataset,sample,section in [('Cervilla','visiumhd_Bladder','visiumhd_Bladder'),('GSE300692','GSM9066380','D1_Johnson'),('GSE299573','GSM9040999','UC_SCC'),('GSE318916','GSM9505917','MPBC1'),('GSE318916','GSM9505918','MPBC2')]:
  sub=a[(a.dataset==dataset)&(a.section==section)];base=dict(dataset=dataset,sample=sample,section=section,available_resolutions_um=','.join(map(str,sub.resolution_um.astype(int))),reusable_cell_results='not located',official_Space_Ranger_count_ready='NO',cell_id_consistency='not assessable; no public cell outputs located',raw_counts_available='YES' if len(sub) and (sub.status=='audited').all() else 'NO',patient_id='unknown')
  if dataset=='Cervilla':base.update(bin_analysis='YES',reconstruction_pilot='CONDITIONAL',official_nuclei_then_aggregate='conditional on adequate image validation; no cell boundaries yet',recommended_route='Prioritize original registered H&E; conditional limited Bin2cell morphology pilot at native hires scale after review; retain unchanged-count aggregation separately from optional correction',missing_inputs='original full-resolution H&E; validated nuclear detail; full rerun FASTQ; independent registration JSON not located',evidence='2um counts+positions+scale; hires approximately 1.64 um/pixel; 1274 tissue bins beyond preview lower edge (0.02697%); bin overlay cannot validate nucleus-level accuracy')
  elif dataset=='GSE300692':base.update(bin_analysis='YES',reconstruction_pilot='NO_CURRENTLY',official_nuclei_then_aggregate='conditional after original H&E recovery and registration',recommended_route='Bin analysis now; recover JohnsonA1_D1_Scan1.qptiff and image-registration JSON before Bin2cell or official nuclei segmentation followed by actual-count aggregation',missing_inputs='original H&E; source image-registration JSON; FASTQ not publicly submitted',evidence='2/8/16um observed Space Ranger counts; hires approximately 5.59 um/pixel; no cell segmentation or Thor inferred output in archive')
  elif dataset=='GSE299573':base.update(bin_analysis='YES_8UM_ONLY',reconstruction_pilot='NO_CURRENTLY',official_nuclei_then_aggregate='not ready for faithful 2um-to-cell aggregation',recommended_route='Use supplied 8um bins; request 2um counts and original registered H&E; future count rerun only after complete inputs obtained',missing_inputs='2um counts; original H&E; full rerun registration/slide/probe metadata',evidence='only public 8um H5 and spatial; SRR33944517 publicly listed but FASTQ intentionally not downloaded; old derived16um is not original measurement')
  else:base.update(bin_analysis='EXPRESSION_AND_RELATIVE_GRID_ONLY',reconstruction_pilot='NO_CURRENTLY',official_nuclei_then_aggregate='not ready; missing image registration',recommended_route='Candidate: request author Bin2cell cell h5ad + boundaries + mapping or recover H&E/positions/scales; keep expression-only audit',missing_inputs='H&E; tissue positions; scalefactors; registered coordinate transform; author cell output',evidence='barcode-encoded grid retained; paper describes bin2cell mpp=0.5 prob_thresh=0.1 nms_thresh=0.5; no downloadable matching cell result located')
  ready.append(base)
 r=pd.DataFrame(ready);r.to_csv(ROOT/'segmentation_readiness.tsv',sep='\t',index=False)
 # Report budget and actual retained bytes without counting compatibility symlinks twice.
 transfers=[json.loads(l) for l in (ROOT/'logs/transfers.jsonl').read_text().splitlines()];traffic=sum(x.get('bytes',0) for x in transfers if 'url' in x and x.get('status')==206)
 perdata={d:sum(p.stat().st_size for p in (ROOT/'data'/d).rglob('*') if p.is_file() and not p.is_symlink()) for d in ['Cervilla','GSE300692','GSE299573','GSE318916']}
 for cmd,name in [(['df','-h',str(ROOT),'/tmp'],'disk_final.txt'),(['quota','-s'],'quota_final.txt'),(['free','-h'],'memory_final.txt')]:
  t=subprocess.run(cmd,capture_output=True,text=True);(ROOT/'logs'/name).write_text(t.stdout+t.stderr)
 manifest_sizes=sum(perdata.values());new_network=sum(p.stat().st_size for p in (ROOT/'data/Cervilla').rglob('*') if p.is_file())+(ROOT/'data/GSE299573/GSM9040999/square_008um/filtered_feature_bc_matrix.h5').stat().st_size
 resources={'initial_user_quota_used_GiB_rounded':458,'soft_quota_GiB':500,'hard_quota_GiB':505,'initial_filesystem_available':'12 TiB rounded','initial_memory_available':'1.1 TiB rounded','download_payload_bytes':traffic,'data_bytes_by_dataset':perdata,'final_data_bytes':manifest_sizes,'new_network_final_bytes':new_network,'download_scope':'selected ZIP members plus one GEO H5','temporary_peak_budget_bytes':1_000_000_000,'audit_seed':123,'threads':2,'resource_blocked':False}
 (ROOT/'logs/resources.json').write_text(json.dumps(resources,indent=2))
 checks={'audited_rows':len(a),'sections':len(r),'mp10_rows':len(g),'all_matrices_audited':bool((a.status=='audited').all()),'no_negative_nonfinite_noninteger':bool((a[['negative_values','nonfinite_values','noninteger_values']].fillna(-1)==0).all().all()),'mp10_counts_bounded':bool((g.total_counts_all_supplied_bins>=0).all()),'no_remaining_part_files':not list((ROOT/'data').rglob('*.part')),'frozen_MP10_exact_match':json.loads((ROOT/'metadata/mp10_lock.json').read_text())['exact_match']}
 assert checks['all_matrices_audited'] and checks['no_negative_nonfinite_noninteger'] and checks['frozen_MP10_exact_match'];assert len(g)==len(a)*20
 for _,row in a.iterrows():
  subset=g[(g.dataset==row.dataset)&(g.section==row.section)&(g.resolution_um==row.resolution_um)];assert subset.gene.nunique()==20
  assert (subset.total_counts_all_supplied_bins<=row.total_gene_expression_counts).all()
  if pd.notna(row.get('matrix_in_tissue_bins')):
   assert (subset.in_tissue_nonzero_bins<=subset.in_tissue_bin_denominator).all()
 for transfer in transfers:
  if 'local' in transfer and 'sha256' in transfer:
   match=assets[assets.local_path==transfer['local']]
   assert len(match)==1 and match.sha256.iloc[0]==transfer['sha256']
 for line in (ROOT/'logs/relocations.jsonl').read_text().splitlines():
  move=json.loads(line);old=Path(move['old']);new=Path(move['new']);assert old.is_symlink() and old.resolve()==new and new.stat().st_ino==move['inode']
 checks['downloaded_SHA256_crosscheck']='PASS';checks['relocation_inode_and_symlink_check']='PASS'
 checks['per_sample_gene_denominators_and_totals']='PASS';(ROOT/'logs/validation.json').write_text(json.dumps(checks,indent=2))
 lines=['# Visium HD 恢复、下载与输入审计','', '审计日期：2026-10-21（Asia/Shanghai）。工作目录：`'+str(ROOT)+'`。本轮只做输入恢复、计数/空间/影像审计和可测性描述；未运行分割、细胞注释、MP10 高低分组、差异分析或整合。','', '## 关键结论','', '| 数据 / section | 已审计分辨率 | 当前路线 |','|---|---|---|','| Cervilla / visiumhd_Bladder | 2 μm | bin 可用；原始 H&E 未定位，核级 pilot 仅有条件 |','| GSE300692 / D1_Johnson | 2、8、16 μm（同一 section） | bin 可用；需补原始高分辨 H&E 后再细胞重建 |','| GSE299573 / UC_SCC | 8 μm | 当前仅 bin；旧 derived 16 μm 不视为新的实测分辨率 |','| GSE318916 / MPBC1、MPBC2 | 各 2 μm | 表达与相对阵列网格可用；缺组织位置/影像/配准 |','', '所有 5 个 section 的 Top20 均需分别看 feature 与非零表达。逐基因精确计数和组织内分母见 `mp10_feature_detection.tsv`。MPBC1/2 的 in_tissue 未知，绝不以 filtered 文件名或相对网格替代组织掩膜。','', '## 恢复、迁移与预算','', '原始解压输入已实体移动到 `data/`；同文件系统 rename 保留 inode，不复制数据。旧入口保留兼容软链接；GSE299573 的四条失效输入链接已修复。迁移和修复日志为 `logs/relocations.jsonl`、`logs/link_repairs.jsonl`。历史分析结果未移动、删除或重算。','', '本地盘点先按文件名与旧 manifest 定位，并扩大到用户项目、生物数据及临时目录；完整搜索范围见脚本/日志。未找到项目根或本轮路径适用的 AGENTS.md。阅读了既有高分辨率审计的输入、结构、冻结与来源记录；其中历史“禁止联网”是旧任务范围，本轮按用户明确授权下载。','', '启动时 quota 约 458/500 GiB（硬上限 505 GiB），约 42 GiB 软配额余量；文件系统空闲约 12 TiB，/tmp 约 526 GiB；可用内存约 1.1 TiB。共享服务器 swap 已占满，因此使用最多 2 个数值线程、逐 bin 块稀疏读取，不生成密集全矩阵或完整 AnnData 副本。Cervilla 完整审计实测峰值 RSS 约 6.37 GiB（`logs/audit_complete_inputs.log` 的 Maximum resident set size），逐样本释放内存。本轮没有因磁盘资源阻断的数据集。','', '| 数据 | 本轮 data 内保留大小（MiB） | 获取方式 |','|---|---:|---|']
 for d,b in perdata.items():lines.append(f'| {d} | {b/2**20:.2f} | '+{'Cervilla':'ZIP Range 选择性下载','GSE300692':'移动既有输入 + 一次本地 TAR.GZ 遍历恢复小影像','GSE299573':'移动既有 spatial + 补下载 8 μm H5','GSE318916':'移动两份既有 H5；不重复下载'}[d]+' |')
 lines += ['',f'本轮记录到的大文件 HTTP 有效载荷总计 **{traffic/2**20:.2f} MiB**（含 Range 探测和重复预读；不含 TLS/HTTP 头）。小型元数据另不足数 MiB；精确成功访问字节见网络日志及 metadata 文件。新下载最终保留 **{new_network/2**20:.2f} MiB**。工作目录 data 合计 **{manifest_sizes/2**30:.3f} GiB**，大部分为已有文件原地迁移。','', 'D1_Johnson 的 27,503,110,814-byte TAR.GZ（25.614 GiB）保留旧位置引用；本轮仅从本地连续读取一次并消费到 gzip EOF，CRC 通过。未重新下载、未提取 BAM/大型可视化包。该归档不计入新工作目录的新增数据体积。临时下载峰值预算控制在约 1 GB 新增量；最终无压缩包副本或 `.part` 遗留。','', 'Cervilla 官方 ZIP 为 3,362,702,746 bytes，官方 MD5 `e3f5e6dacae3b4c53f8007dabc64d9a6` 与用户线索一致。每个范围请求要求 HTTP 206 且 Content-Range 精确匹配；只提取必要成员，逐成员 CRC 和本地 SHA256 通过。**未下载完整 ZIP，因此没有验证整包 MD5。** 不获取 CytAssist_Bladder 数据集。','', '## 表达、空间与影像','', '| section / 分辨率 | GEX features × bins | counts 中位数 | genes 中位数 | barcode 匹配 | hires μm/像素 |','|---|---:|---:|---:|---:|---:|']
 for _,row in a.iterrows():
  cq=json.loads(row.counts_quantiles);gq=json.loads(row.genes_quantiles);match='NA' if pd.isna(row.barcode_match_rate) else f'{100*row.barcode_match_rate:.1f}%';mpp='NA' if pd.isna(row.get('hires_microns_per_pixel')) else f'{row.hires_microns_per_pixel:.3f}'
  lines.append(f'| {row.section} / {int(row.resolution_um)} μm | {int(row.gene_expression_features):,} × {int(row.bins):,} | {cq[2]:g} | {gq[2]:g} | {match} | {mpp} |')
 lines += ['', '按 100,000 bins/块读取 CSC；检查全体存储值、CSC 索引和 pointer，Gene Expression 单独统计；未发现负值、非有限值或非整数值。整数性仅为数值检查，实验 counts 归属由 GEO 说明、原始 Space Ranger 输出路径、H5 header 和报告共同支持。未使用 Thor 推断/扩散/插补表达。重复 symbol/ID、重复 barcode/坐标、零计数 bins、完整分位数与仿射系数均保存在矩阵审计表。','', '空间表按 barcode 重排连接，不依赖行顺序。array_row/array_col 为阵列索引，pxl_row/col_in_fullres 为 source-image 像素，叠图使用 fullres × tissue_hires_scalef。独立按 seed=123 取样拟合 grid→pixel 仿射，数值残差与镜像/旋转方向记录在表中；仿射行列式本身不作为错误判据。可视化只能排查明显错配，不能证明亚核级配准准确。','', '所有已定位 H&E 风格图均为处理后的 hires/lowres 或 alignment preview，未确认原始全分辨 H&E。D1_Johnson 的 web summary 明确引用 `JohnsonA1_D1_Scan1.qptiff` 和 `H1-PF3RWFN-D1-fiducials-image-registration.json`，但它们没有包含在当前归档中；CytAssist TIFF 是另一类影像，不能替代该原始 H&E。图像实际格式/宽高见 asset manifest 与 `metadata/image_inventory.json`。','', 'Cervilla web summary 记录 Space Ranger 3.0.0、GRCh38-2020-A、Probe Set v2.0；原始输入路径指向 `Bladder_HD_he.tif` 与 `Bladder.json`，但这些文件未包含在 ZIP。Cervilla 的 6000×4711 hires，依据 scalefactors 为约 1.64 μm/像素，属于下采样图。仅在明确核形态仍可辨、配准通过人工复核的局部区域考虑有条件 pilot；优先找回原始注册 H&E。另有 1,274 个 in_tissue bins（0.02697%）位于提供的 hires 图下边界之外：图高 4711，组织内 y 最大 4879.343 pixels，约 168.343 pixels 超界。它们保留在 counts/MP10 统计中，未静默过滤；依赖影像的 pilot 应避开这些缺失视野并先核查原始影像。其余 99.9730% 的组织内 bins 在图像边界内。不得插值放大并声称恢复核细节。D1_Johnson 的约 5.59 μm/像素预览不作为合格分割输入。','', 'ROI 用 ROI1/2/3 命名，固定 seed=123，覆盖不同坐标位置与边缘；不能确认上皮/淋巴/间质身份，因此未赋细胞标签。所有图仅作描述；MPBC1/2 只画 barcode 相对网格和 counts/genes 分布，不能画出可信 H&E 叠图。','', '## MP10 冻结与检测','', '冻结来源：`../MP10_spatial_multisection_top20_vs_core19_v1/01_signatures/01_MP10_Top20_vs_Core19_gene_table.csv` 的 DirectTop20，按 MP10_rank 排序；SHA256、期望名单和匹配结果见 `metadata/mp10_lock.json`。另用旧正式 `highres_ST/highres_spatial_audit_v3_GSE300692/01_MP10_lock/frozen_EpiMP10.tsv` 的 top20 交叉核对，其来源为正式 `ST/mapping/07_EpiMP_gene_basis.qs`；两处与用户给出的 20 基因完全一致，未替换基因、未修改定义。','', '| section / 分辨率 | feature 存在 | 所有输入 bins 实际非零基因 | 组织内实际非零基因 |','|---|---:|---:|---:|']
 for keys,sub in g.groupby(['section','resolution_um'],sort=False):
  tissue='不可判定' if sub.in_tissue_nonzero_bins.isna().all() else str(int((sub.in_tissue_nonzero_bins>0).sum()))+'/20'
  lines.append(f'| {keys[0]} / {keys[1]} μm | {int(sub.feature_present.sum())}/20 | {int(sub.any_nonzero_all_supplied_bins.sum())}/20 | {tissue} |')
 lines += ['', '组织内比例的分母是 supplied matrix 中按 barcode 匹配且 in_tissue=1 的 bins；空间表中所有组织内 bins 另列，避免把 filtered-matrix 分母当作整张 capture array。MPBC1/2 仅报告 all supplied bins 的检测和总计数，组织内列留 NA；候补数据不能因此宣布完成组织内可测性验证。','', '## 来源冲突、细胞结果和下一阶段','', 'D1_Johnson：GEO 声称 Space Ranger 3.1.0；实际 H5 与 web summary 为 3.0.1。本报告保留两种来源，按实际输入版本审计。probe_set.csv 为 Visium Human Transcriptome Probe Set v2.0、GRCh38-2020-A。2/8/16 μm 的 filtered 集合略有不同，总计数差异不能直接解读为数据丢失。','', 'MPBC1/2：GSE318916 标题是 HGT1 研究，但原始论文明确将两例 MPBC 标为额外的 MIBC 混合 micropapillary/urothelial 手术病例，不能沿用旧记录的 NMIBC 归类。论文称 Space Ranger 3.1.3；实际 H5 分别为 MPBC1=3.0.0、MPBC2=4.0.1，版本不一致须向作者核实。样本到患者的独立匿名 ID 未公开；跨数据集患者重叠/相邻切片关系 unknown。','', '论文报道 MPBC 使用 Bin2cell（mpp=0.5、prob_thresh=0.1、nms_thresh=0.5），随后做归一化/log/scaling。当前 GEO 仍只提供两份 2 μm H5；本轮未在数据开放说明定位作者 cell h5ad、核/细胞边界或 bin-to-cell map。若以后获取 h5ad，必须区分 X 的变换值与未经改变的聚合 counts，并检查 cell ID 一致性；扩张距离/软件版本未能从公开结果确认，记 unknown。','', '当前五个 section 均没有可直接复用且来源/映射可核验的细胞矩阵。Thor 论文使用推断表达，不等同实验 counts；当前 D1_Johnson 归档中没有该推断矩阵/细胞边界。checked 作者仓库仅有分析代码，不以“软件能生成”当作“已有结果可复用”。','', '官方 `spaceranger segment` 可独立对合格 H&E 做核分割，不需要先下载 FASTQ；但其图像核结果仍需与 bins 配准并聚合 counts。`spaceranger count` 重跑是另一路线，需要测序及完整实验/影像输入。本轮五个 section 均不具备本地完整 count 重跑输入。D1_Johnson 的人类原始测序未公开；GSE299573 已核实 SRX29144058 / SRR33944517，ENA 列出两端 FASTQ 合计 22,692,009,528 bytes，仅登记不下载。','', '建议下一阶段优先找回 Cervilla 和 D1_Johnson 的原始注册 H&E；取得合格输入后再决定 Bin2cell 或官方核分割+实际 counts 汇总。GSE299573 保持 8 μm bin 层，不能反推原始 2 μm。MPBC1/2 保留候补状态，优先请求作者 cell counts、边界和映射或原始配准输入。','', '## 访问受阻与交付检查','', '服务器默认大小写代理配置冲突导致最初 HTTPS EOF；随后显式使用现有本地代理成功获取 Zenodo/GEO。CORA Dataverse 镜像返回 HTTP 200 的反机器人 HTML，而非 JSON，记录为不可验证；未试图规避其人机验证。此限制不影响 Zenodo 原始 record 的核心成员获取。核心成员完成后，额外的 6,134,282-byte probe_set.csv 请求连续遇到 TLS EOF（包括官方文件直链复核），因此仅登记远端成员；版本已由 web summary 核实。没有将该小文件故障作为阻断其他数据的理由。Github 个别初次 TLS/分支错误已重试，成功和失败证据保留。未以某一个失效来源终止其他审计。','', '交付物：四个 TSV（asset_manifest、matrix_spatial_audit、mp10_feature_detection、segmentation_readiness）、本报告、figures 下轻量 PNG、code 下可复跑脚本，以及 metadata/logs 中来源快照、校验值和验证结果。没有大型报告包。`logs/validation.json` 记录表间交叉检查；`logs/software_versions.json` 记录复用 Python 3.12 环境。','', '复跑：在本目录执行 `bash code/run_audit.sh`。下载脚本只取本轮最小输入；ZIP 支持按已压缩成员断点恢复、CRC 后清理自身临时片段，GEO H5 按范围块续传；已有文件和审计检查点均按状态复用。若输入发生变化，审计签名失效后重查。图示抽样 seed=123。','', '## 核验来源','', '- [Cervilla Zenodo record](https://zenodo.org/records/17999961)；[作者代码仓库](https://github.com/scervilla/SpatialBenchmarking)。','- [D1_Johnson GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9066380)；[Thor 原始论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC12325965/)。','- [UC_SCC GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9040999)；[ENA accession](https://www.ebi.ac.uk/ena/browser/view/SRR33944517)。','- [GSE318916](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318916)；[MPBC 原始论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC13232474/)。','- [10x 官方 segment 输入](https://www.10xgenomics.com/support/software/space-ranger/4.0/analysis/running-pipelines/space-ranger-segment)；[Bin2cell 官方示例](https://bin2cell.readthedocs.io/en/latest/notebooks/demo.html)。','']
 (ROOT/'AUDIT_REPORT.md').write_text('\n'.join(lines))
 print(json.dumps({'validation':checks,'resources':resources},indent=2))
if __name__=='__main__':finalize()
