from base import *
from datetime import datetime
from zoneinfo import ZoneInfo
def read(p):return pd.read_csv(p,sep='\t')
def md(t):
 cols=t.columns;return '| '+' | '.join(cols)+' |\n|'+'|'.join(['---']*len(cols))+'|\n'+'\n'.join('| '+' | '.join('NA' if pd.isna(v) else f'{v:.3f}' if isinstance(v,(float,np.floating)) else str(v) for v in row)+' |' for row in t.itertuples(index=False,name=None))
stamp=datetime.now(ZoneInfo('Asia/Shanghai')).strftime('%F %T %Z')
status=[];audit=[];counts=[];test=[];effects=[];rels=[];sums=[];cv=[]
for plat,ss in [('HD',HD),('Xenium',XE)]:
 for s in ss:
  d=O/plat.lower()/s;is_hd=plat=='HD';labels=d/'labels.parquet'
  facts={'sample':s,'platform':plat,'unit':'16 µm bin' if is_hd else 'segmented cell','identity':'pending','graph':'ready' if (d/'graph_50.npz').exists() else 'pending','matrix':'pending','state comparison':'pending'}
  if is_hd:
   y=json.loads((d/'conservation.json').read_text());e=pd.read_parquet(d/'epi64_centers.parquet')
   audit.append((s,y['bins8'],y['bins16'],y['umi'],len(e),int(e.H05.sum()),int(e.shared_removed_200.sum()),int(e.independent_200.sum()),int(((e.coverage_200>=.5)&(e.independent_200>0)).sum())))
   if not (d/'fit_complete.txt').exists():facts['identity']=f'RCTD running ({len(list((d/"weights").glob("[0-9][0-9][0-9][0-9][0-9].tsv.gz")))} chunks saved)'
  if labels.exists():
   z=pd.read_parquet(labels);facts['identity']='computed';t=read(d/'label_counts.tsv');t.insert(0,'sample',s);counts.append(t)
   if not is_hd:
    u=read(d/'threshold.tsv').iloc[0];corr=read(d/'correlations.tsv');vs={tuple([r.a,r.b]):r.Spearman for r in corr.itertuples()}
    sums.append((s,len(z),int(z.epi_eligible.sum()),u.q95,u.strict_high,bool(u.group_evaluable),vs['Chol_biosynthesis_8','MP10_covered'],vs['Chol_biosynthesis_8','H05_transfer'],vs['MP10_covered','H05_transfer'],vs['Chol_biosynthesis_8','nGeneTranscripts']))
  want=[d/f'{k}_{r}_matrix.tsv.gz' for k in (['subtype','major'] if is_hd else ['subtype','major','chol']) for r in ([50] if is_hd else [50,200]) if k!='chol' or not (d/'threshold.tsv').exists() or read(d/'threshold.tsv').iloc[0].group_evaluable]
  if all(p.exists() for p in want):facts['matrix']='computed'
  if all((d/f'effects_{r}.tsv').exists() for r in [50,200]):facts['state comparison']='computed'
  for r in [50,200]:
   p=d/f'effects_{r}.tsv'
   if p.exists():
    v=read(p);v['platform']=plat;effects.append(v)
   p=d/f'subtype_{r}_matrix.tsv.gz'
   if p.exists():
    v=read(p);v=v[(v['from']<v.to)&~v['from'].isin(['Mixed','Unknown'])&~v.to.isin(['Mixed','Unknown'])].copy();v['sample']=s;v['platform']=plat;v['radius']=r;rels.append(v)
  status.append(facts)
status=pd.DataFrame(status);put(status,O/'reports/branch_status.tsv')
h=pd.DataFrame(audit,columns=['HD','8 µm 输入','16 µm 单元','UMI','合格 Epi64','H05','200 µm 共享边剔除','独立邻接边','覆盖合格区域']);put(h,O/'reports/hd_grid_and_regions.tsv')
ss=pd.DataFrame(sums,columns=['样本','细胞','合格上皮','Q95','严格 high','可分组','合成–MP10 ρ','合成–H05 ρ','MP10–H05 ρ','合成–深度 ρ']);put(ss,O/'reports/xenium_states.tsv')
if counts:
 ct=pd.concat(counts);put(ct,O/'reports/identity_coverage.tsv')
 flags=ct.pivot_table(index='sample',columns='flag',values='units',aggfunc='sum',fill_value=0)
 flags['all_units']=flags.sum(axis=1)
 for key in ['major_only','subtype_supported','mixed','uncertain']:
  if key not in flags:flags[key]=0
  flags[key+'_fraction']=flags[key]/flags.all_units
 put(flags.reset_index(),O/'reports/per_sample_coverage.tsv')
if effects:
 eff=pd.concat(effects);put(eff,O/'reports/per_sample_state_effects.tsv')
 rows=[]
 for (r,k),sub in eff.groupby(['radius','neighbor']):
  v={}
  for plat,sl in [('HD',HD),('Xenium',XE)]:
   df=sub[sub.platform==plat].set_index('sample').reindex(sl)
   valid=(df.n_high>=5)&(df.n_other>=5)&(df.section_nodes>=10)&df.delta.notna()
   sign=np.sign(df.delta[valid]);v[plat]=list(sign)
   for sample,row in df.iterrows():
    cv.append((sample,plat,r,k,row.delta,row.delta_density_mm2,row.n_high,row.n_other,row.section_nodes,bool(valid.loc[sample])))
  def direction(a):
   if len(a)<2:return '信息不足'
   if min(a)>0:return '正向'
   if max(a)<0:return '负向'
   return '异质或含零'
  hd,xe=map(direction,[v['HD'],v['Xenium']]);cross='信息不足' if '信息不足' in [hd,xe] else '方向一致' if hd==xe and hd in ['正向','负向'] else '异质'
  rows.append((r,k,len(v['HD']),hd,len(v['Xenium']),xe,cross))
 cross=pd.DataFrame(rows,columns=['radius','neighbor','HD可评片','HD方向','Xenium可评片','Xenium方向','跨平台']);put(cross,O/'reports/cross_platform_state_direction.tsv')
 hprim=cross[cross.radius==200][['neighbor','HD可评片','HD方向']];xprim=cross[cross.radius==50][['neighbor','Xenium可评片','Xenium方向']]
 primary=hprim.merge(xprim,on='neighbor',how='outer')
 primary['主尺度对照']=['信息不足' if pd.isna(hh) or pd.isna(xx) or '信息不足' in [hh,xx] else '方向一致' if hh==xx and hh in ['正向','负向'] else '异质' for hh,xx in zip(primary.HD方向,primary.Xenium方向)]
 put(primary,O/'reports/primary_scales_direction.tsv')
 put(pd.DataFrame(cv,columns=['sample','platform','radius','neighbor','delta_fraction','delta_density','n_high','n_other','section_nodes','descriptive_support']),O/'reports/cross_platform_direction_inputs.tsv')
if rels:
 allrel=pd.concat(rels,ignore_index=True);put(allrel,O/'reports/sample_subtype_pairs.tsv.gz')
 out=[]
 fulltax=read(O/'shared/taxonomy.tsv');focus=sorted(fulltax.loc[fulltax.major.isin(['T/NK','Myeloid','Endothelial','Fibroblast']),'label'].unique())
 from itertools import combinations
 for a,b in combinations(focus,2):
  v=allrel[(allrel['from']==a)&(allrel.to==b)&(allrel.radius==50)&(allrel.n_from>=10)&(allrel.n_to>=10)&allrel.z_score.notna()]
  dirs=[]
  for pl in ['HD','Xenium']:
   w=v[v.platform==pl].z_score;dirs.extend([len(w),'信息不足' if len(w)<2 else '富集' if w.min()>0 else '减少' if w.max()<0 else '异质或含零'])
  out.append((a,b,*dirs,'信息不足' if '信息不足' in dirs else '方向一致' if dirs[1]==dirs[3] and dirs[1] in ['富集','减少'] else '异质'))
 paired=pd.DataFrame(out,columns=['a','b','hd_n','hd','xenium_n','xenium','cross']);put(paired,O/'reports/cross_platform_subtype_pairs.tsv')
val=O/'reports/heldout_validation.tsv'
note=['# HD 与 Xenium 首轮空间分析','',f'更新时间 {stamp}。下表为实际落盘进度。', '',md(status),'','## HD 的表达来源与区域分母','',md(h),'','200 µm 边数按合格 Epi64 计，每个区域仅计一次。独立邻接边为剔除共享 8 µm 来源后的候选边，区域还需达到 0.5 实测覆盖。','', '## Xenium 三种状态','',md(ss) if len(ss) else '等待细胞级评分。','','五片均实测 5,001 个生物特征。8/8 合成、9/20 MP10、10/27 Reactome 和 35/74 Hallmark 的覆盖基因在各样本均有检出。每个细胞的实际检出另见对应样本的 epi_detection 表。分组使用每片合格上皮严格大于 Q95。未达分组支持的切片保留连续结果。']
if val.exists():
 v=read(val);note+=['','## 留整来源验证','',md(v.groupby(['mode','score'])[['Spearman','AUROC','average_precision','prevalence']].median().reset_index()),'','表内是逐留出来源指标的中位数。逐来源样本量、原始阳性率、Pearson 以及 PR 曲线分别保存于 heldout_validation.tsv 和 shared/validation/PR_curves.tsv.gz。']
if counts:
 t=ct[ct.flag=='subtype_supported'].pivot_table(index='supported',columns='sample',values='units',aggfunc='sum').reindex(columns=HD+XE)
 for s in ct['sample'].unique():t[s]=t[s].fillna(0)
 note+=['','## 有转移或拟合及原始 marker 支持的亚群','',md(t.rename_axis('亚群').reset_index()),'','计数单位为各平台内的 bin 或分割细胞。0 表示该片通过本轮支持规则的单元数为零；缺测写作 NA。细胞身份置信和总体丰度的推断还受 panel 与深度制约。']
if effects:
 note+=['','## 状态周围的组成与密度','',md(cross[(cross.radius.isin([50,200]))&cross.neighbor.isin(['Treg','Endothelial (all)','Fibroblast (all)','Myeloid (all)'])]),'','方向取 high 减同片其他合格上皮的邻域份额。流程同时保留单位覆盖面积密度和邻居数。逐片至少 5 个可评 high 和 other、目标类全片至少 10 个节点才进入方向摘要。至少 2 张可评片同号称为描述性同向；各平台具体分母见表。']
if rels:
 cnt=paired.cross.value_counts();note+=['','## 亚群间 50 µm 邻域','',f'限定重点亚群且剔除对角线。方向一致 {cnt.get("方向一致",0)} 对，异质 {cnt.get("异质",0)} 对，信息不足 {cnt.get("信息不足",0)} 对。完整名称和逐片方向见 cross_platform_subtype_pairs.tsv 与 sample_subtype_pairs.tsv.gz。']
note+=['','## 解读与优先级','','先检查亚群的可检测 marker 和 major-only 分母，再判断状态测量的可比性。H05_transfer 的量纲是参考支持度，实测程序反映转录状态。重复邻近只能用于选择待验证的位置和关系。全局置换无法充分控制空间自相关。未把切片或成千上万个 bin/cell 作为独立患者做显著性推断。','','图件按样本放在 figures，每张单独导出 PNG 和 PDF。状态相关散点有平滑线。']
(O/'REPORT.md').write_text('\n'.join(note)+'\n')


if val.exists() and len(ss)==5 and rels:
 from interpret import interpret
 interpret()
 with (O/'REPORT.md').open('a') as f:f.write('\n[亚群支持、三状态对照、逐片空间差异与验证优先级](reports/INTERPRETATION.md)\n')

from figure_index import index
index()
with (O/"REPORT.md").open("a") as fh:fh.write("\n[逐样本图件入口与全部 PNG/PDF](reports/FIGURES.md)\n")

from completion import write as completion_table
completion_table()
