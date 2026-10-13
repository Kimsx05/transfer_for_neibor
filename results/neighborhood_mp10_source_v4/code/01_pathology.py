from base import *
import urllib.request,shutil,html
P=OUT/'01_pathology_search';S=P/'sources';U=P/'user_review'
m=pd.read_parquet(V3/'01_features/metadata.parquet');man=pd.read_csv(ROOT/'cell2location_full_run/v2_30type/01_manifest/01_spatial_manifest_frozen.csv')
audit=ROOT/'spatial_ecology_top20_v2/patient_identity_audit_20261002/sources'
sources=[]
for dataset,g in man.groupby('dataset'):
 for r in g.itertuples():
  status='unknown';note='No traceable regional pathology labels found; H&E and registration evidence available; Epi/CNV are not pathology.'
  if 'GSE285715' in dataset:note+=' GEO source tissue describes adjacent normal for N1-N3 and tumor source for T1-T3; not a spot-level pathology determination.'
  if 'GSE238145' in dataset:note+=' Author tissue field only splits B12/B13 specimens, not tumor/normal histology.'
  if 'GSE319536' in dataset:note+=' Author code references Official_Annot, Level1_Annot and Level2_Annot; corresponding Seurat RDS/metadata not present in public code repository.'
  sources.append(dict(source_id=f'local_{r.section_id}',dataset_id=dataset,section_id=r.section_id,annotation_level='none',evidence_status=status,limitation=note,local_path=str(V3/'05_user_review/evidence'/r.section_id),exact_mapping='section matched existing registered H&E',coordinate_system='existing validated image_coordinates.csv',registration_evidence='registration.json',author_date='2026-10-13 audit'))
urls={'GSE319536_code':'https://raw.githubusercontent.com/VividKai/SpatialBC/main/Fig1.R','GSE319536_tree':'https://api.github.com/repos/VividKai/SpatialBC/git/trees/main?recursive=1','GSE171351_tree':'https://api.github.com/repos/KnottLab/bladder-snSeq/git/trees/master?recursive=1','GSE238145_paper':'https://pmc.ncbi.nlm.nih.gov/articles/PMC8361097/','CanDis_Zenodo':'https://zenodo.org/api/records/18237374','NPJ_paper':'https://www.nature.com/articles/s41698-026-01417-x','GSE285715_GEO':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285715','GSE319536_paper':'https://aacrjournals.org/cancerdiscovery/article/16/10/2037/788515/A-Spatial-Atlas-of-Muscle-Invasive-Bladder-Cancer'}
# Correct GSE238145 is documented locally; Gouin PMC8361097 belongs to GSE171351.
urls['GSE171351_paper']=urls.pop('GSE238145_paper')
logs=[]
for k,url in urls.items():
 try:
  r=urllib.request.urlopen(url,timeout=20);p=S/(k+'.txt');p.write_bytes(r.read());logs.append(dict(source_id=k,url=url,status=r.status,local_path=str(p),sha256=sha(p),date='2026-10-13'));log(f'pathology {k}: {r.status}')
 except Exception as e:logs.append(dict(source_id=k,url=url,status='request_failed',reason=str(e),date='2026-10-13'))
for p in audit.iterdir():
 if any(x in p.name for x in ['GSE238145_code_summary.R','GSE238145_article.txt','GSE171351_family.soft','GSE285715_family.soft','GSE319536_family.soft','NPJ_table1_retry','ZFP36']):
  if p.is_file() and p.stat().st_size<15e6:
   shutil.copy2(p,S/p.name);logs.append(dict(source_id=p.stem,url='local authoritative source archive',status='READ_LOCAL',local_path=str(p),sha256=sha(p),date='2026-10-13'))
logs.append(dict(source_id='Inhouse',status='NO_PUBLIC_ACCESSION',reason='Private specimens: searched project histology/annotation/metadata manifests; no author/public accession to resolve; do not infer mapping.'))
logs.append(dict(source_id='GSE238145_metadata',status='NOT_PATHOLOGY',local_path=str(audit/'GSE238145_metadata_retry.tsv.gz'),reason='81586 rows; tissue: NA74533, B13 4991, B12 2062; author code lines105-113 specimen split.'))
tab(P/'public_and_local_search_log.tsv',logs);tab(P/'pathology_source_inventory.tsv',sources)
tab(P/'pathology_resolved.tsv',pd.DataFrame(columns=['section_id','barcode','pathology_label','source_id','evidence_level']))
tab(P/'pathology_conflicts.tsv',[{'source':'GSE238145 tissue','reason':'Specimen split B12/B13 cannot be recoded as malignant/normal'},{'source':'GSE285715 N/T source tissue','reason':'Whole specimen origin cannot establish every central spot pathology'},{'source':'CanDis/NPJ','reason':'Local aliases lack verified public sample crosswalk'}])
ov=m[['spot_uid','section_id','barcode']].copy();ov['pathology_label']='unknown';ov['pathology_evidence_level']='unknown';ov['pathology_source_id']='';tab(P/'spot_pathology_overlay.tsv',ov)
roi=pd.read_csv(V3/'05_user_review/ROI_review_to_fill.tsv',sep='\t');assert roi.reviewed_label.isna().all(),'Existing review values require explicit import, do not discard'
E=U/'evidence';E.mkdir(exist_ok=True)
for i,r in roi.iterrows():
 for col in ['evidence_figure','raw_HE']:
  p=Path(r[col]);dest=E/(r.roi_id+'.png' if col=='evidence_figure' else r.section_id+'_HE.png')
  if not dest.exists():shutil.copy2(p,dest)
  roi.loc[i,col]='evidence/'+dest.name
 roi.loc[i,'question']='请结合 H&E 指认肿瘤上皮/非肿瘤上皮/其他；混合区请拆分并提供 polygon 或条码，记录坐标变换。CNV diploid 不能判正常。'
for c in ['polygon_coordinates','coordinate_transform','spot_members_file']:roi[c]=''
tab(U/'ROI_review_to_fill.tsv',roi)
body=['<!doctype html><meta charset="utf-8"><title>病理 ROI 待办</title><h1>133 个尚无可复用区域病理判定的 ROI</h1><p>顺序沿用预先确定的优先级；不根据 MP10 调整。<a href="ROI_review_to_fill.tsv">填写表</a></p>']
for r in roi.itertuples():body.append(f'<h2>{html.escape(r.roi_id)} — {r.n_common_anchors} anchors</h2><p>{html.escape(r.question)}</p><a href="{r.raw_HE}">整张 H&E</a><br><img loading="lazy" style="max-width:1100px" src="{r.evidence_figure}">')
(U/'review_index.html').write_text('\n'.join(body))
(U/'WHAT_USER_NEEDS_TO_DO.md').write_text('请打开 review_index.html，按表顺序复核尚无可信病理标签的 ROI。最先五项：CanDis Sample8、GSE319536 BC11、CanDis Sample1、CanDis Sample3、GSE285715 T1。\n\n在 ROI_review_to_fill.tsv 填 reviewed_label、review_comment、reviewer、review_date。允许值 tumor_epithelium / non_tumor_epithelium / other / mixed_needs_split / uncertain。混合区域需 parent_roi_id、polygon_coordinates 与 coordinate_transform 或 spot_members_file。\n\nGSE319536：请向作者索取 Fig1.R 引用的每切片 SeuratObj/*.seurat.umap.spatial.rds 或 V2503.figure1.metainfo.summary.txt，以及 Official_Annot/Level1_Annot/Level2_Annot 字典和坐标配准。CanDis 与 NPJ：需要本地 Sample 名到公开标本的原始交叉表。GSE238145 的 tissue 只是 B12/B13 标本拆分，不能当病理。GSE285715 的 N 是邻癌来源而非全切片病理正常。\n\n未知病理不影响本轮评分。所有图仅为证据，未替代用户判定。\n')
log('pathology first pass complete')
