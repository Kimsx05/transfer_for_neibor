from pathlib import Path
import csv
B=Path('mp10_spatial_communication_v1');O=B/'spatialcellchat'
def read(p):
 with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
mp={x['analysis_fine_label']:x['major'] for x in read(B/'input_review/spatial_label_crosswalk.tsv')}
secs=['GSE171351_A1','Inhousedata_T4469479']; ds={s:{(x['source'],x['target']):x for x in read(O/s/'main/fine_edges.tsv')} for s in secs}
rows=[];treg=[]
for sec,d in ds.items():
 for dr in ['Epi→TME','TME→Epi']:
  key=lambda node,state:(f'Epi_MP10_{state}',node) if dr=='Epi→TME' else (node,f'Epi_MP10_{state}')
  allnodes=sorted({s for s,t in d if not s.startswith('Epi_MP10_')})
  tn=[x for x in allnodes if mp.get(x)=='T/NK']
  w=lambda node,state:float(d[key(node,state)]['total_mean_strength'])
  f=lambda node,state:float(d[key(node,state)]['total_flow'])
  vals={st:w('Treg',st) for st in ['low','high']}
  rank={st:1+sum(w(x,st)>vals[st] for x in tn) for st in vals}
  sl={st:vals[st]/sum(w(x,st) for x in tn) for st in vals};sa={st:vals[st]/sum(w(x,st) for x in allnodes) for st in vals}
  treg.append(dict(section=sec,direction=dr,mean_low=vals['low'],mean_high=vals['high'],delta=vals['high']-vals['low'],flow_low=f('Treg','low'),flow_high=f('Treg','high'),rank_low=rank['low'],rank_high=rank['high'],n_TNK=len(tn),share_TNK_low=sl['low'],share_TNK_high=sl['high'],share_all_TME_low=sa['low'],share_all_TME_high=sa['high']))
  for node in allnodes:
   if mp.get(node) not in ['Endothelial','Fibroblast','Mural']:continue
   lineage=mp[node];nn=[x for x in allnodes if mp.get(x)==lineage];lo=w(node,'low');hi=w(node,'high')
   rows.append(dict(section=sec,lineage=lineage,node=node,direction=dr,low=lo,high=hi,delta=hi-lo,share_lineage_low=lo/sum(w(x,'low') for x in nn),share_lineage_high=hi/sum(w(x,'high') for x in nn)))
for rr,name in [(rows,'focus_descriptive_summary.tsv'),(treg,'treg_mean_flow_summary.tsv')]:
 with open(O/name,'w') as f:
  out=csv.DictWriter(f,fieldnames=rr[0],delimiter='\t');out.writeheader();out.writerows(rr)
lines=['# SpatialCellChat 两片试点重点解读','',
'本说明仅比较各片内部的 H/L 方向和谱系排序；两片不是原 scRNA 五队列的复刻，不比较不同切片的绝对数值。GSE171351_A1 与 Inhousedata_T4469479 各一片，重复数不足以作患者差异检验。以下 mean 均为总 mean strength：所有纳入且经位置 QC 的 LR 逐个 F_l/K_l 后相加。','',
'## Treg：主 mean 两片均下降，谱系份额与 flow 不完全同向','',
'|切片|方向|Low mean|High mean|Low flow|High flow|T/NK排名 L→H|T/NK内份额 L→H|全TME份额 L→H|',
'|---|---|---:|---:|---:|---:|---|---|---|']
for x in treg:
 lines.append(f"|{x['section']}|{x['direction']}|{x['mean_low']:.6f}|{x['mean_high']:.6f}|{x['flow_low']:.3f}|{x['flow_high']:.3f}|{x['rank_low']}/{x['n_TNK']}→{x['rank_high']}/{x['n_TNK']}|{x['share_TNK_low']:.2%}→{x['share_TNK_high']:.2%}|{x['share_all_TME_low']:.2%}→{x['share_all_TME_high']:.2%}|")
lines += ['',
'两片两向 Treg 的主 mean 都是 High<Low。A1 中 Treg 在可辨 T/NK 层级排名第二，但不属于 high 状态下增量突出者；其 T/NK 内份额下降。Inhouse 的 high 状态 Treg mean 排名5/5，尽管绝对 mean 下降，T/NK 内份额反而略升——该片其他 T/NK 成分下降幅度不同。因此不应把“份额上升”写成“Treg绝对通讯增强”。T_NK_unresolved 汇合了三个不可独立辨识标签，保留于分母和排名只是组成事实，不能命名为新亚型。','',
'Inhouse Treg 的 flow 与 mean 反向：Epi→Treg flow 106.585→120.582，而 mean 0.204514→0.147416；Treg→Epi flow 115.411→133.793，而 mean 0.215975→0.156735。这是 F 与每条LR的支持链接分母 K 口径不同所致，不是两套独立H/L网络。总 mean 是 sum_l(F_l/K_l)，一般不等于 sum_l F_l / sum_l K_l；表里的 summed_LR_supported_pair_denominators 仅作覆盖汇总，不能当主 mean 的单一分母。A1 的 mean 和 flow 两向都下降，故仅按 mean 的两片同向不能升级为总量的重复下降。','',
'## EC、Fibroblast 与 Mural：哪些方向重复，哪些相反','',
'在两个方向上都重复的主 mean 方向：Arterial EC、Venous EC、Lymphatic EC 两片均 High<Low；iCAF 两片均 High>Low；Pericyte 两片均 High<Low。','',
'在两个方向上均出现切片反例：Capillary EC、vCAF、Matrix_CAF、rCAF、SMC 在 A1 中 High>Low，在 Inhouse 中 High<Low。没有把这些不一致条目剔除或称为稳定高状态候选。','',
'EC 的绝对排序也不固定：A1 中 Venous EC 最大，Inhouse 中 Lymphatic EC 最大；Capillary EC 的正向结果只在首片出现。内皮内份额中，Venous EC 两片两向均下降；Arterial/Lymphatic 在 A1 下降而 Inhouse 上升，Capillary相反。因此绝对 mean 与谱系内份额必须并列看。参考中未独立建模 Tip/Cycling EC，本分支不制造其列。','',
'Fibroblast 中 vCAF 在两片两向均为绝对 mean 最大，但其 H−L 方向不重复；iCAF 是这两片唯一两向均 mean 增加的 CAF 标签，且两向谱系内份额均上升，但只可作为两片试点的探索性重复。Matrix_CAF 的谱系份额在 A1 下降、Inhouse 上升，即使其绝对方向相反。rCAF 的份额变化很小或因方向不同，不据此升级生物学结论。','',
'Mural 按原参考保留独立 major。A1 中 SMC 高于 Pericyte，Inhouse 中 Pericyte 高于 SMC；SMC 的 mean 方向跨片相反，Pericyte 两片均下降，但 Pericyte 在 Mural 内份额 A1 下降、Inhouse 上升。它们没有被强行改名为 CAF。','',
'|谱系/亚群|方向|A1 Low→High mean|A1 Δ|Inhouse Low→High mean|Inhouse Δ|','|---|---|---|---:|---|---:|']
for node in [n for n in mp if mp[n] in ['Endothelial','Fibroblast','Mural']]:
 for dr in ['Epi→TME','TME→Epi']:
  rr=[x for x in rows if x['node']==node and x['direction']==dr]
  if len(rr)!=2:continue
  a,b=sorted(rr,key=lambda x:secs.index(x['section']))
  lines.append(f"|{a['lineage']} / {node}|{dr}|{a['low']:.5f}→{a['high']:.5f}|{a['delta']:+.5f}|{b['low']:.5f}→{b['high']:.5f}|{b['delta']:+.5f}|")
lines += ['','## 首片有限范围敏感性','',
'固定1037个LR、相同Q与参数，仅将 native diffusible 范围200→300μm；该版本ECM与分泌共用此范围，接触保持105μm。70个major/fine方向比较68项符号保留。两项接近零的翻转为 Epi→Mural（−0.004674→+0.000984）和 Epi→Classical Mono（+0.000334→−0.000174）。Treg两向mean仍下降；其余可辨EC、Fibroblast、Mural fine方向保留。每条Contact位置概率与主运行完全一致，已作独立数值核对。','',
'所有结论均为 c2l 比例加权的跨spot通讯归属，依赖去卷积与混合表达；不能确定某个受体转录本一定来自同spot中的Treg，不能将这些数值称为免疫抑制强度、功能验证或因果机制。保留参考identifiability caution。另见 METHOD_NOTES.md、LR_category_summary.tsv、run_summary.tsv、numerical_checks.tsv。']
(O/'SUMMARY_NOTES.md').write_text('\n'.join(lines)+'\n')
print('wrote SUMMARY_NOTES.md and 2 focus tables')
