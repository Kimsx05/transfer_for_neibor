from lib import *
import argparse
p=argparse.ArgumentParser();p.add_argument('--version',default='primary');arg=p.parse_args();out=O if arg.version=='primary' else O/'versions'/arg.version
f=read(out/'FIGURE_INDEX.tsv')
if (out/'figures/EIS_architecture.pdf').exists() and not f.figure.eq('EIS_architecture').any():
    f.loc[len(f)]=['EIS_architecture','E/I/S composition by architecture','plot_data/architecture_composition.tsv','Denominator is all classified tissue; each unit and then each available source gets equal weight.'];save(f,out/'FIGURE_INDEX.tsv')
for n in f.figure:
    for ext in ['pdf','png']:assert (out/'figures'/(n+'.'+ext)).stat().st_size>1000
(out/'FIGURE_INDEX.md').write_text('# 正式单图\n\nPDF 和 300 dpi PNG 各自独立保存。原始表路径相对于主分析目录，版本相关组别表相对于该版本。\n\n|图|文件|绘图表及图注|\n|---|---|---|\n'+''.join(f'|{r.title}|[PDF](figures/{r.figure}.pdf) · [PNG](figures/{r.figure}.png)|{r.plot_data} · {r.caption}|\n' for r in f.itertuples())+'\n所有指标的分来源散点、趋势和分布见 [离线指标图](metric_explorer.html)。')
review=read(out/'review/index.tsv')
(out/'review/INDEX.md').write_text('# 切片审图\n\n每张图分别显示 H&E 原图、E/I/S、T/I、局部 T Moran、获空间支持的窗口状态、双边界、MP10 和旧新覆盖。H&E 使用原像素方向，其他空间图使用 Visium 阵列坐标，并列图之间的方向可能不同。\n\n'+''.join(f'- [{r.section_id}]({r.section_id}.jpg) · {r.group} · E 域 {r.E_scored} 点\n' for r in review.itertuples()))
print(len(f),'formal figures;',len(review),'reviews')

