from base import *
def index():
 rows=[];header=['# 图件与逐样本入口','','主矩阵为 50 µm 的全片亚群矩阵。HD 状态图为冻结 Epi64 H05，Xenium 状态图为 8 基因合成 UCell 连续分数。','', '| 样本 | 大类 | 亚群 | 状态 | 邻域 Z-score | high/other 构成 |','|---|---|---|---|---|---|']
 for s in HD+XE:
  d=O/'figures'/s
  for p in sorted(d.glob('*.png')):
   rows.append((s,p.stem,str(p.relative_to(O)),str(p.with_suffix('.pdf').relative_to(O)),p.with_suffix('.pdf').exists()))
  cols=[('大类','major'),('亚群','supported'),('状态','Epi64_H05' if s in HD else 'state_Chol_biosynthesis_8'),('50 µm','subtype_50_matrix'),('主半径','neighborhood_composition_'+('200' if s in HD else '50'))]
  header.append('| '+s+' | '+' | '.join(f'[{name}](../figures/{s}/{k}.png)' if (d/(k+'.png')).exists() else '等待计算' for name,k in cols)+' |')
 header+=['','[全部独立 PNG/PDF 图件清单](figure_manifest.tsv)']
 (O/'reports/FIGURES.md').write_text('\n'.join(header)+'\n');put(pd.DataFrame(rows,columns=['sample','figure','png','pdf','pdf_present']),O/'reports/figure_manifest.tsv')
if __name__=='__main__':index()

