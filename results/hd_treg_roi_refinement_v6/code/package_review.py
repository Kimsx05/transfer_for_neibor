from pathlib import Path
import pandas as pd,zipfile
from PIL import Image
root=Path(__file__).resolve().parents[1];r=root/'proseg_observed'
rows=['# 图像和表格入口','','正式主线在 proseg_observed。目录顶层 figures 中的空间比较和切片顶层 effects.tsv 是保守归属对照。','', '## 核修订与原图校准','']
for s in ['Cervilla','D1_Johnson']:
 for z,t in [('nuclear_changes','拆分、合并、新增与旧核保留'),('calibration','12 个校准窗口')]:
  p=root/'figures'/f'{s}_{z}.jpg';rows.append(f'- [{s} · {t}]({p.relative_to(root)})')
rows+=['','## 已实际打开的 MP10 空间窗口','','图内标注数量是矩形视窗内数量，正式统计取同组件累计圆，分母见 neighborhoods.parquet。']
nn=pd.read_csv(r/'tables/neighbor_visual.tsv',sep='\t')
for z in nn.itertuples():
 p=r/'figures'/f'{z.opened_view}.png';assert p.exists();rows+=['',f'- [{z.section} · {z.opened_view}]({p.relative_to(root)})　{z.observation}']
rows+=['','## 已实际打开的逐对象及未定位来源卡','','绿色方格代表最终唯一归属的实测捕获格，橙线为核，黄线为 Proseg 推定胞体。']
notes=pd.read_csv(r/'tables/visual_notes.tsv',sep='\t');notes['key']='cell_'+notes.id.astype(str);un=pd.read_csv(r/'tables/unlocalized_visual.tsv',sep='\t');un['key']=un.review_id.str.replace(':','_');allnotes=pd.concat([notes,un],ignore_index=True)
pics=[]
for z in allnotes.itertuples():
 p=r/z.section/'cards'/(z.key+'.jpg');assert p.is_file();pics.append(p);rows.append(f'- [{z.section} {z.key}]({p.relative_to(root)})　{z.observation}')
rows+=['','每个来源到细胞的多对多映射见各片 source_links.parquet。全队列审阅状态见 review_index.parquet，打开的卡与尚未打开的卡均可按 review_id 取出。','']
(root/'FIGURE_INDEX.md').write_text('\n'.join(rows))
files=[root/x for x in ['REPORT.md','METHODS.md','FIGURE_INDEX.md','frozen_config.yaml','input_manifest.tsv','run_status.tsv','rerun.sh']]+list((root/'code').glob('*.py'))+list((root/'code').glob('*.sh'))+list((r/'tables').glob('*.tsv'))+[root/'tables'/t for t in ['v5_funnel.tsv','v5_v6_changes.tsv','software.tsv','image_software.tsv','calibration_fields.tsv','nuclear_examples.tsv']]+list((root/'figures').glob('*nuclear_changes.jpg'))+list((root/'figures').glob('*calibration.jpg'))+list((r/'figures').glob('*.png'))+pics
for s in ['Cervilla','D1_Johnson']:
 for n in ['effects.tsv','roi_manifest.tsv','thresholds.tsv','priority_dominant_successor.tsv','display_windows.tsv','historical_score_check.tsv']:files.append(r/s/n)
with zipfile.ZipFile(root/'review.zip','w',zipfile.ZIP_DEFLATED,compresslevel=3) as archive:
 for f in files:archive.write(f,f.relative_to(root))
with zipfile.ZipFile(root/'review.zip') as archive:assert archive.testzip() is None
