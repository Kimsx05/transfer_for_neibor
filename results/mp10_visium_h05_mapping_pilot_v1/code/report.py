from common import *
from html import escape
import markdown
def table(header,rows):
 def cell(v):return str(v).replace('|','&#124;')
 return '\n'+'| '+' | '.join(header)+' |\n| '+' | '.join(['---']*len(header))+' |\n'+''.join('| '+' | '.join(map(cell,r))+' |\n' for r in rows)+'\n'
def f(x):return f'{x:.3f}' if pd.notna(x) else '—'
metrics=read(O/'section_metrics.tsv')
metrics=metrics[metrics.branch=='main']
ab=read(O/'feature_ablation.tsv');pv=read(O/'pseudo_spot_validation.tsv')
m=['# 普通 Visium 四片 H05 映射与冻结 MP10 比较','**结论。** B17 的原始支持度与 MP10 最一致，但很大部分伴随 Epi 含量和测序深度变化。T1 的弱正相关在调整后保留，在移除 Top20 后接近零。B1 效应弱，BC01 的条件支持度与 MP10 呈弱负相关。H05 映射可以帮助定位参考谱系支持与原评分不一致的位置；目前证据不足以用它判定 H05 阳性 spot 或定量其细胞含量。','实际结果目录为 `'+str(O)+'`。四片全片共 12,979 spots，冻结 query 共 8,187。表达、旧评分、坐标和新预测逐键一一匹配，原始 UMI 列和与旧 nCount_RNA 的最大差均为 0。','## 各片主结果','相关在各自冻结 Epi query 上计算。r 为 Pearson，ρ 为 Spearman。调整项为 c2l Epi fraction 与 log1p(nCount)，残差相关与对应未调整相关使用同一组点。低分母为 Epi_support < 0.10，包含零。']
m.append(table(['切片','支持度','有效 n / query','缺失','低分母 n（有效中）','r','ρ','残差 r','残差 ρ'],[(r.section_id,'原始' if r.score=='H05_support' else '条件',f'{r.n} / {r.n_query}',r.missing,r.low_den_effective,f(r.pearson),f(r.spearman),f(r.resid_pearson),f(r.resid_spearman)) for r in metrics.itertuples()]))
m.extend(['- **B1。** 原始支持度的全片最大值仅 0.112。支持度随 Epi fraction 变化（r = 0.722），但 MP10 与该比例的片内相关约为 0.007，调整后原始支持度的相关略升。把 Epi 分母纳入解释后，仅见很弱的一致性。','- **B17。** 原始支持度与 Epi fraction 的 r = 0.827，与 log 深度的 r = 0.451。原始支持度与 Epi_support 几乎同步（r = 0.994）。即使限定在同一批 1,203 个分母非零的点，MP10 对原始和条件支持度的 r 仍分别为 0.223、0.043。69 个条件支持度 ≥ 0.9 的 query 中，68 个处于低分母区。','- **BC01。** 原始支持度同时伴随组成（r = 0.757）与 log 深度（r = 0.667）增加，MP10 与两项协变量的关系方向相反。query 中原 MP10 有 43.3% 为零，条件支持度的负向关系在调整后接近零。','- **T1。** 92.8% 的 query 具有 Epi_support > 0.99，全部 query 均超过低分母界限。这里原始与条件读数较接近，组成和深度调整对结果影响小。原评分有 39.5% 为零，仍保留在散点与统计中。','## 去 Top20 的敏感性','四片主映射实际纳入全部 20 个定义成员（含 ACAT2），消融后均为 0，重建了 PCA 和 anchors。下表在每行两个分支共同有效的点上比较，条件支持度行会与上表单分支分母不同。'])
m.append(table(['切片','支持度','成对 n','主 / 去除特征数','主 r','去除 r','两分支支持度 ρ'],[(r.section_id,'原始' if r.score=='H05_support' else '条件',r.n_paired,f'{r.main_features} / {r.masked_features}',f(r.main_r),f(r.masked_r),f(r.support_rank_agreement)) for r in ab.itertuples()]))
m.extend(['T1 原始支持度的 r 从 0.153 降至 0.027，两个版本的支持度排序相关仅 0.555。B17 排序保留较多（ρ = 0.990），但其较强的组成依赖也保留。BC01 在两个版本中的原始相关均接近零，条件相关为负。定义基因和谱系背景都影响结果，影响程度随片变化。','## 不一致区域与原始表达','按片内 query 的相反四分位数定位；BC01 和 T1 的 MP10 下四分位数均为 0，包含所有边界并列。表中 UMI 为实测 Top20 成员的原始计数之和，位置为已有组织内外侧标注。按原阵列六邻接划分的每个区域均可回查 barcode。'])
dc=read(O/'discordant_expression.tsv');dc=dc[dc.discordance!='other']
reg=read(O/'discordant_regions.tsv').sort_values(['n','component'],ascending=[False,True]).drop_duplicates(['section_id','discordance'])
dc=dc.merge(reg[['section_id','discordance','n','x_px_center','y_px_center']],on=['section_id','discordance'],suffixes=('','_largest'))
m.append(table(['切片','MP10 / 支持度','n','inside n','中位 nCount','中位 Top20 UMI','最大连通区 n @ 像素质心'],[(r.section_id,'高 / 低' if r.discordance.startswith('MP10 high') else '低 / 高',r.n,r.n_inside,round(r.median_depth),round(r.median_top20_UMI),f'{r.n_largest} @ ({r.x_px_center:.0f}, {r.y_px_center:.0f})') for r in dc.itertuples()]))
m.extend(['B1 的高 MP10 / 低支持度点多数在外侧（49/56），相反方向多数在内侧（33/37）。B17 两个方向的 Top20 UMI 中位数均为 14，但 Epi_support 中位数分别为 0 和 0.938。BC01 的低 MP10 / 高支持度点深度约为相反方向的 3.8 倍，其原始 Top20 UMI 中位数也更高；AUCell 描述相对排序，测序深度影响原始计数的解释。T1 的高 MP10 / 低支持度区域反而保留更高 Top20 UMI，中位数为 109 对 70。逐基因 UMI、EPCAM、KRT8 等可从各片主表核对。','## 全片非上皮位置','预定诊断区使用 c2l Epi fraction < 0.10 且其他大类占优。B1 无符合点；补查原 0.25 门槛以下且非 Epi 占优的 38 点，H05 支持度全部为零，中位 COL1A1 为 25 UMI、EPCAM 为 1。B17 的 452 个预定诊断点全部为零，以 Mural 为主（311 点）。BC01 的 2,439 点中位数为 0、Q95 为 0.012、最大值为 0.260，高尾见于 T/NK 占优位置。T1 仅 36 点，中位数为 0、最大值达 0.316；COL1A1 与 KRT8 的中位数分别为 126 与 9 UMI，可见混合组织信号。低 Epi 丰度可用于谱系定位检查，无法充当胆固醇生物学的绝对阴性真值。','## 小批来源留出伪 spot','旧 TEST 中选取完整来源留出的三个 specimen，共 240 个不同 mixture（两组来源，a/b 同患者），每个 mixture 对应两个深度。按原始 H05 成员逐细胞重标，38 个 mixture 的真实 H05 为零。为单独比较深度，最终两深度固定为共同检出的 859 个特征和完全相同的参考 PCA 输入，再分别寻找 anchors。下表的 ρ 为真实 H05 占全部组成细胞的份额与原始支持度的 Spearman。'])
m.append(table(['深度','留出 specimen','n','真零 n','真零支持度中位 / Q95','有 H05 的中位支持度','含量 ρ'],[(r.depth,r.source.split('::')[-1],r.n,r.n_true_zero,f'{f(r.zero_support_median)} / {f(r.zero_support_Q95)}',f(r.positive_support_median),f(r.fraction_vs_raw_spearman)) for r in pv.itertuples() if r.source!='all_reused_mixtures']))
m.extend(['三个来源样本内均见含量增加的方向性；真零和有 H05 的分数明显重叠，T4463209a 的真零中位数甚至更高。由 medium 到 low，三个 specimen 的配对支持度中位变化分别为 −0.178、−0.051、−0.037。两深度总体真零分布的 Q95 分别为 0.742 和 0.669。旧 TEST 的组成设计及细胞重复使用限制了外推，不能据此宣称新独立验证或细胞数量校准。','## 如何使用这一层','保留原 mixed-spot MP10 作为主读数。标签映射增加了两类信息：一个 spot 的表达整体更支持哪类参考谱系，以及该谱系支持中有多少相对偏向冻结 high 状态。这两层在 B17 的分离揭示了组成混淆，在 T1 的接近说明该片主要需要关注定义基因共享。两种方法共享 RNA；逐片相关和消融衡量方法一致性。','**扩展判断。** 现阶段适用于这四片的描述性定位、筛查不一致区域及回看组织。片间异质性、低分母、深度效应和真零伪 spot 的较高支持度共同限制了阳性判定。暂不建议把该流程作为 H05 分类或定量工具铺开到 63 片。','## 交付与复核','44 张独立 PNG/PDF 已分片保存。每片包括 H&E、MP10、两种 H05 支持度、Epi_support、c2l Epi fraction、nCount、低分母图、两张带趋势的散点和一张不一致区域图。图片预览见 HTML；原始表、消融和方法见 [文件目录](README.md)、[统计明细](section_metrics.tsv)、[图件索引](figure_index.tsv) 和 [方法](METHODS.md)。','全流程使用现有 Seurat 5.5.1。GNU time 记录主映射进程峰值 RSS 为 5.19 GiB，四片及两深度初次映射耗时 5 分 37 秒；固定配对特征后的伪 spot 重跑耗时 1 分 29 秒。本轮未运行 c2l 或扩增切片。'])
body='\n\n'.join(m)+'\n'
(O/'REPORT.md').write_text(body)
gallery=''
for s in S:
 gallery+='<h3>'+s+'</h3><div class="grid">'
 labels={'tissue':'组织背景','EpiMP10_Top20_AUCell':'冻结 query 的 MP10','H05_support':'Epi high 原始支持度','H05_conditional_support':'条件支持度','discordance':'不一致区域','Epi_support':'Epi 支持度','Epithelial_fraction':'c2l Epi fraction','nCount_RNA':'测序深度','low_denominator':'低分母','scatter_H05_support':'原始支持度散点','scatter_H05_conditional_support':'条件支持度散点'}
 for k,label in labels.items():
  p=f'figures/{s}/{k}';gallery+=f'<figure><a href="{p}.pdf"><img loading="lazy" src="{p}.png" alt="{escape(label)}"></a><figcaption>{label}</figcaption></figure>'
 gallery+='</div>'
html='<!DOCTYPE html><html lang="zh-CN"><meta charset="utf-8"><title>四片 H05 映射与 MP10</title><style>body{font:16px/1.65 Arial,sans-serif;max-width:1200px;margin:3em auto;padding:0 2em;color:#202830}table{border-collapse:collapse;font-size:14px}td,th{padding:.38em .7em;border-bottom:1px solid #ccc}th{background:#eff3f5}a{color:#135579}code{overflow-wrap:anywhere}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:12px}figure{margin:0}img{width:100%;height:210px;object-fit:contain}figcaption{text-align:center}h2{margin-top:2em}</style><main>'+markdown.markdown(body,extensions=['tables'])+'<h2>逐片图件（点击打开 PDF）</h2>'+gallery+'</main></html>'
(O/'REPORT.html').write_text(html)

