from common3 import *
from urllib.parse import quote
cl=read(OUT/'classification/section_classification.tsv').set_index('section_id')
sh=read(O/'tables/section_shapes_near.tsv');fig=read(O/'audit/section_figure_index.tsv');status=read(O/'audit/section_boundary_status.tsv')
txt=['# 逐切片双边界图与免疫地图','', '每行一个实际切片。四个链接依次为原始近窗、原始全域、中心化近窗、中心化全域。横轴从正距离向负距离，色带为箱内 IQR。近窗为 ±1000 μm。形状列是固定描述规则的结果。','', '|切片|队列、簇、归属频率|全部 E 图|ENN+ENS 图|负侧近窗形状（依前两列顺序）|免疫地图|','|---|---|---|---|---|---|'];rows=[]
for s,r in cl.iterrows():
    parts=[];shps=[]
    for b in B:
        v=fig[(fig.section_id==s)&(fig.boundary==b)];ss=status[(status.section_id==s)&(status.boundary==b)].iloc[0]
        if not ss.display_n:parts.append(ss.plot_status)
        else:
            pp=[]
            for scale,win,tag in [('raw','near','原近'),('raw','full','原全'),('centered','near','中近'),('centered','full','中全')]:
                path=v[(v.scale==scale)&(v.view==win)].file.iloc[0];dest=Path(os.path.relpath(O/path,OUT));pp.append(f'[{tag}]({dest})')
            parts.append(' · '.join(pp)+f' ({ss.display_n:,})')
        shp=sh[(sh.section_id==s)&(sh.boundary==b)].near_shape.iloc[0];shps.append(shp)
        rows.append(dict(section_id=s,boundary=b,source=r.dataset_id,unit=r.unit,tumor=r.tumor,cluster=r.exploratory_cluster,stability=r.stability,old_label=r.provisional_architecture_group,shape=shp,status=ss.plot_status,n=ss.display_n))
    txt.append('|'+s+'|'+('肿瘤来源' if r.tumor else '邻近正常投影')+' · '+r.exploratory_cluster+f' · {r.stability:.2f}|'+parts[0]+'|'+parts[1]+'|'+' / '.join(shps)+f'|[审图](review/{s}.jpg)|')
save(rows,OUT/'audit/sample_index.tsv');(OUT/'SAMPLES.md').write_text('\n'.join(txt)+'\n')
idx=read(OUT/'audit/group_figure_index.tsv')
for kind in ['region','R12']:
    key='composition_'+kind
    if (OUT/'figures/groups'/(key+'.pdf')).exists():
        idx.loc[len(idx)]=[key,'figures/groups/'+key+'.pdf','figures/groups/'+key+'.png','tables/composition_plot_data.tsv','All three groups share the contributing sources; equal sections, units and sources.']
save(idx,OUT/'audit/all_group_figure_index.tsv')
(OUT/'FIGURE_INDEX.md').write_text('# 图像索引\n\n[按切片查询 488 张距离图和 63 张免疫地图](SAMPLES.md)。缺界面组合也有状态。独立组图的区间、分母及直接数据见下表。\n\n|图|文件|数据|图注|\n|---|---|---|---|\n'+''.join(f'|{r.figure}|[PDF]({r.pdf}) · [PNG]({r.png})|{r.data}|{r.caption}|\n' for r in idx.itertuples()))

