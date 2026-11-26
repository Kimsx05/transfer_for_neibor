from shared import *
from html import escape
import markdown,re,platform,scipy,matplotlib,pyarrow,h5py
fig=pd.concat([read(O/x) for x in ['figure_index_primary.tsv','figure_index_violin.tsv']],ignore_index=True)
save(fig,'figure_index.tsv')
def table(head,rows):
 def fmt(x):
  if isinstance(x,(float,np.floating)):return '—' if not np.isfinite(x) else f'{x:.3f}'
  return str(x).replace('|','&#124;')
 return '\n| '+' | '.join(fmt(h) for h in head)+' |\n| '+' | '.join(['---']*len(head))+' |\n'+''.join('| '+' | '.join(fmt(x) for x in row)+' |\n' for row in rows)
def paired(x,y):return f'{x:.3f} / {y:.3f}'
def stem(s):return s.split('_')[-1]
mods={'Synthesis':'胆固醇合成','Homeostasis':'胆固醇稳态','Fatty_acid':'脂肪酸代谢','Efflux':'胆固醇外排（补充）'}
q=read(T/'join_audit.tsv');fc=read(T/'functional_correlations.tsv');ct=read(T/'functional_contrasts.tsv')
f=read(T/'four_way.tsv');f=f[f.level=='subtype'];ft=read(T/'functional_TME.tsv.gz').query('level=="subtype"')
nc=read(T/'niche_scores.tsv');sw=read(T/'endpoint_swap_summary.tsv');old=read(T/'old_vs_new.tsv');low=read(T/'low_TME.tsv')
def record(df):return df.iloc[0]
def get(df,**kw):
 for k,v in kw.items():df=df[df[k]==v]
 return record(df)
body=['# 普通 Visium 四片的脂质转录、TME 与外围 niche','[打开全部独立图件](figures.html) · [四片 H05 功能对照图](figures/overview/four_sections_function_high.pdf) · [四片 MP10 功能对照图](figures/overview/four_sections_function_MP10.pdf)',
      (O/'code/interpretation.md').read_text(),
      '## 范围、参考与读数',
      '绝对路径为 `'+str(O)+'`。方案 A 分别把固定的二类上皮参考与独立的 29 类 TME 参考映射到相同的全片 RNA。连续高状态支持度与每类 TME 支持度在冻结合格点上比较。方案 B 读取同一位置原实测混合 spot MP10 和 log1p(c2l posterior mean absolute abundance)。两个交叉组合每次只替换其中一端。所有相关以切片为单位，原始和残差 Pearson 为主，另存 Spearman 及残差 Spearman。',
      table(['切片','全片','冻结 query','niche 已分配 / 未分配','上皮参考实际特征（Top20）','TME 参考实际特征（Top20）'],[(s.section_id,s.full,s.query,f'{s.assigned_query} / {s.unassigned_query}', *[f'{len(z)} ({sum(z.Top20)})' for z in [read(O/f'features/{s.section_id}.{model}.tsv') for model in ['Epi','TME']]]) for s in q.itertuples()]),
      '2,449 个上皮细胞逐 ID 继承旧 pilot，19,269 个 TME 细胞由相同未留出来源的细注释重新分层抽样。T/NK、Myeloid 等大类只在单独的大类表中汇总，各亚群原名见 annotation_mapping.tsv。旧节点表在 B17、BC01、T1 另有 49、3,209、8 个全片分析对象外的节点，四片对象内均无坐标或旧 query 连接丢失。',
      '## 连续功能关系和两端区域',
      '每个格子按“完整基因集原始 r / 去 Top20 原始 r / 去 Top20 调整 r”排列。调整为 Epi fraction 与 log1p(nCount) 的双端 OLS 残差。最后一列为去重叠合成分数在片内 high 与 low 区域的均值差。high/low 划分保留 Q75/Q25 边界并列，真实组大小列在同一行。']
rows=[]
for s in S:
 for x in [H,MP]:
  values=[]
  for k in mods:
   a=get(fc,section_id=s,exposure=x,module_score=k+'__full');b=get(fc,section_id=s,exposure=x,module_score=k+'__minus20')
   values.append(f'{a.pearson:.3f} / {b.pearson:.3f} / {b.resid_pearson:.3f}')
  t=get(ct,section_id=s,exposure=x,module_score='Synthesis__minus20')
  rows.append((stem(s),'H05 映射' if x==H else '实测 MP10',*values,f'{t.n_high} / {t.n_low}',t.high_minus_low_mean))
body+=[table(['片','读数',*mods.values(),'high / low n','合成 Δ均值'],rows),
       '合成、稳态、脂肪酸和补充外排的完整成员数分别为 27、74、158、82，去重叠后分别为 18、63、152、81。各片去重叠后的实际检出数分别为 17、59–60、135–141、51–62。四个过程均超过预设最低覆盖。空间 RNA 秩积分反映混合组织的转录特征，无法直接推断胆固醇含量、通量、细胞内在活动或细胞间代谢交换。',
       '## 同 spot 的 TME 关系与方法替换','每片列出两方案调整后最正和最负的共同细亚群。顺序在全部 29 类中确定。每个关系以该片共同有效的 query 计算；完整原始/调整矩阵和四组合细表均保留。']
rows=[]
for s in S:
 z=f[f.section_id==s]
 for combo in ['A','B']:
  g=z[z.combination==combo].dropna(subset=['resid_pearson']).sort_values('resid_pearson')
  b,a=g.iloc[0],g.iloc[-1]
  rows.append((stem(s),combo,a['n'],a.cell_type,paired(a.pearson,a.resid_pearson),b.cell_type,paired(b.pearson,b.resid_pearson)))
body.append(table(['片','方案','共同 n','最正亚群','原始 / 调整 r','最负亚群','原始 / 调整 r'],rows))
abund=read(T/'abundance_sensitivity.tsv.gz')
rows=[]
for s in S:
 r=f[(f.section_id==s)&(f.combination=='B')].nlargest(1,'resid_pearson').iloc[0]
 z=abund[(abund.section_id==s)&(abund.level=='subtype')&(abund.exposure==MP)&(abund.cell_type==r.cell_type)].set_index('TME_units')
 rows.append((stem(s),r.cell_type,r.resid_pearson,*[z.loc[k,'resid_pearson'] for k in ['mean_absolute','mean_all30_fraction','mean_TME29_fraction']]))
body+=['对上述 B 最正亚群，换用不同丰度单位后的调整 Pearson 如下。不同分母会改变问题，各单位结果均完整保存。',table(['片','B 亚群','log1p 均值','原始均值','全 spot 比例','TME 内比例'],rows),
       '下表在调整 Pearson 上比较两端替换。绝对变化先在固定另一端的两种口径间平均，再在 29 类取中位数。类型谱 ρ 在亚群间计算，比较本片 A 和 B 对各类关系的排序。',
       table(['片','可比类型','A/B 同号数','A/B 类型谱 ρ','替换功能端 |Δr|','替换 TME 端 |Δr|'],[(stem(r.section_id),r.n_types,r.same_sign,r.A_B_type_profile_spearman,r.median_function_swap,r.median_TME_swap) for r in sw[(sw.level=='subtype')&(sw.metric=='resid_pearson')].itertuples()]),
       '## 哪些细胞与实测功能共现',
       '直接将去 Top20 后的实测功能分数与两种 TME 读数配对。下表列出每片、每个主要模块的调整后最正亚群，便于区分 H05 映射关系与实测功能关系。完整正负结果保存在 functional_TME.tsv.gz。']
rows=[]
for s in S:
 for k in list(mods)[:3]:
  z=ft[(ft.section_id==s)&(ft.module_score==k+'__minus20')]
  a=z[z.TME_method=='transfer'].nlargest(1,'resid_pearson').iloc[0];b=z[z.TME_method=='log1p_mean'].nlargest(1,'resid_pearson').iloc[0]
  rows.append((stem(s),mods[k],a.cell_type,a.resid_pearson,b.cell_type,b.resid_pearson))
body+=[table(['片','实测功能','TME 映射最正','调整 r','c2l 丰度最正','调整 r'],rows),'## 不一致区域与旧联合模型',
       '表中位置和表达均来自本片相反分位尾部。Top20 UMI 和非 Top20 的 DHCR24 UMI 是原始计数中位数，功能值是去重叠合成分数中位数。所有 spot 和更多基因均可回查。']
ds=read(T/'discordant_expression.tsv')
body+=[table(['片','不一致方向','n / inside n','中位 nCount','Epi fraction','Top20 UMI','DHCR24 UMI','合成分数'],[(stem(r.section_id),r.group,f'{r.n} / {r.n_inside}',r.median_nCount_RNA,r.median_Epithelial_fraction,r.median_top20_UMI,r.median_UMI_DHCR24,r.median_Synthesis__minus20) for r in ds[ds.group!='other'].itertuples()])]
rows=[]
for s in S:
 a=get(old,section_id=s,reference_readout='H05_support');b=get(old,section_id=s,reference_readout='H05_conditional_support');c=get(old,section_id=s,reference_readout=MP)
 z=low[(low.section_id==s)&(low.population=='query')&(low.stratum=='bottom_quartile')].sort_values('n_top').iloc[-1]
 rows.append((stem(s),paired(a.pearson,b.pearson),paired(c.pearson,c.resid_pearson),z.Q25_absolute_TME,z.cell_type,z.median_winner_score))
body+=[table(['片','新 vs 旧原始 / 旧条件 r','新 vs MP10 原始 / 调整 r','低 TME Q25','低 TME 最常居首类','最高类支持度中位数'],rows),
       '旧条件支持度分母为零的点沿用 NA，仅影响对应对照。旧联合模型使用大类竞争、LogNormalize 与强制纳入可测 Top20。新模型同时改变了上皮与 TME 参考、归一化、特征和 PC 数，差异不能全部归于分开参考。低 TME 位置也会有一个预测最高的亚群，这是 TME-only 竞争性评分的重要限制。',
       '## 原外围 niche',
       '以下比较仅使用原已分配 query。列出至少 20 个中心的原 niche 中，H05、MP10 和去重叠合成分数平均值最高的 niche；完整数据保留所有规模的 niche。该门槛仅用于挑选可读的报告示例。']
rows=[]
for s in S:
 a=nc[(nc.section_id==s)&(nc.n>=20)]
 vals=[]
 for x in [H,MP,'Synthesis__minus20','Homeostasis__minus20','Fatty_acid__minus20']:
  z=a[a.readout==x].nlargest(1,'mean').iloc[0]
  vals.append(f'N{int(z.niche)} (n={int(z.n)}, {z["mean"]:.3f})')
 rows.append((stem(s),*vals))
body.append(table(['片','H05','MP10','合成','稳态','脂肪酸'],rows))
profile=read(R/'mp10_paper_distance_full_v1/tables/source_R12_profiles.tsv').copy()
profile['niche']=profile.cluster+1
save(profile,'tables/original_niche_profiles.tsv')
w=profile.sort_values('log2_enrichment',ascending=False).drop_duplicates(['niche','ring']).pivot(index='niche',columns='ring',values='cell_type')
body+=['原冻结 R12 跨队列画像中，按相对全局平均的 log2 富集排序，每个 niche 在 R1 / R2 居首的类型如下。该表只解释原标签的组成定义。',table(['原 niche','R1 富集首位','R2 富集首位'],[(int(i),r.R1,r.R2) for i,r in w.iterrows()])]
ring=read(T/'peripheral_vs_center.tsv.gz')
rows=[]
for s in S:
 a=ring[(ring.section_id==s)&(ring.exposure=='Synthesis__minus20')&(ring.cell_type=='SPP1 TAM')].set_index('scale')
 rows.append((stem(s),int(a.n.iloc[0]),*[a.loc[k,'resid_pearson'] for k in ['center_q05_all30','R1_q05_all30','R2_q05_all30']]))
body+=['以 SPP1 TAM 作固定尺度示例，比较完全相同的分配中心上，去重叠胆固醇合成分数与中心、R1 和 R2 原 q05/all30 比例的残差相关。R1/R2 均排除中心；用于定义 niche 的梯度及细胞画像不构成独立验证。',table(['片','相同中心 n','同 spot','外围 R1','外围 R2'],rows),
       '## 复核入口','四片冻结 query 共 8,187 个，其中 7,213 个连接到原 niche，974 个未分配。共交付 288 张独立图（每片 71 张，加 4 张跨片效应图）。测得单阶段最大 RSS 为 7.52 GiB。部分 SCT 基因离散度拟合出现迭代上限提示。八组预测均成功写出，映射进程尾部另有脚本解析报错；已重新按完整原轴复读八组结果并核验数值，见 [运行说明](METHODS.md) 和 [映射复核日志](logs/map_verify.log)。',
       '[独立图件浏览](figures.html) · [机器可读图件索引](figure_index.tsv) · [四组合共同样本表](tables/four_way.tsv) · [实测功能表](tables/functional_correlations.tsv) · [参考及输入说明](README.md) · [详细方法](METHODS.md)']
md='\n\n'.join(body)+'\n';(O/'REPORT.md').write_text(md)
css='body{font:16px/1.6 Arial,sans-serif;color:#182732;max-width:1320px;margin:2em auto;padding:0 1.4em}td,th{padding:.35em .6em;border-bottom:1px solid #ccd5db}table{border-collapse:collapse;font-size:.9em}th{background:#eef3f6}a{color:#175883}h2{margin-top:2em}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:1em}figure{margin:0}img{width:100%;height:210px;object-fit:contain}figcaption{overflow-wrap:anywhere}code{overflow-wrap:anywhere}details{margin:1em 0}'
def page(title,content):return '<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+title+'</title><style>'+css+'</style><main>'+content+'</main></html>'
(O/'REPORT.html').write_text(page('四片 Visium 功能定位',markdown.markdown(md,extensions=['tables'])))
gallery=['<h1>四片 Visium 独立图件</h1><p>展开切片下的图型，点击预览打开矢量 PDF。</p>']
for s in ['four_sections',*S]:
 gallery.append('<h2>'+s+'</h2>')
 for kind,g in fig[fig.section_id==s].groupby('kind',sort=False):
  gallery.append(f'<details><summary>{escape(kind)} ({len(g)})</summary><div class="grid">')
  for row in g.itertuples():gallery.append(f'<figure><a href="{escape(row.pdf)}"><img loading="lazy" src="{escape(row.png)}" alt="{escape(row.figure)}"></a><figcaption>{escape(row.figure)}</figcaption></figure>')
  gallery.append('</div></details>')
(O/'figures.html').write_text(page('独立图件',''.join(gallery)))
(O/'versions_python.txt').write_text('\n'.join([platform.python_version()]+[x.__name__+' '+x.__version__ for x in [pd,np,scipy,matplotlib,pyarrow,h5py]])+'\n')

