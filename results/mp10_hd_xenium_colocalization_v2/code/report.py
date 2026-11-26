from common import *
from scipy.stats import false_discovery_control
import datetime,html
def read(p):
 return pd.read_csv(p,sep='\t') if p.exists() and p.stat().st_size>1 else pd.DataFrame()
def merged(name,sub=False):
 frames=[]
 for s in HD+XE:
  p=(OUT/'niche'/s) if sub else folders(s)[1];a=read(p/name)
  if len(a):
   if 'sample' not in a:a['sample']=s
   frames.append(a)
 return pd.concat(frames,ignore_index=True) if frames else pd.DataFrame()
def bh(a,group):
 a=a.copy();a['FDR']=np.nan
 for _,t in a.groupby(group,dropna=False):
  ix=t.index[t.P.notna()];a.loc[ix,'FDR']=false_discovery_control(a.loc[ix,'P'].values) if len(ix) else np.nan
 return a
def md(z):
 def val(x):return '—' if pd.isna(x) else (f'{x:.4g}' if isinstance(x,float) else str(x)).replace('|','/')
 cols=list(z.columns);return '| '+' | '.join(cols)+' |\n|'+'|'.join(['---']*len(cols))+'|\n'+'\n'.join('| '+' | '.join(val(k) for k in row)+' |' for row in z.itertuples(index=False,name=None))
def report():
 d=OUT/'reports';d.mkdir(exist_ok=True)
 a=merged('effects.tsv.gz');base=a.copy();rev=merged('reverse_effects.tsv.gz');nich=merged('within_niche_effects.tsv',True)
 for tab in [a,rev,nich]:
  if len(tab):tab['unit_resolution']=tab.unit_resolution.astype(str)
 a.loc[a.endpoint.str.startswith('CDF_'),'P']=np.nan
 a.loc[a.endpoint.str.startswith('CDF_'),'reason']='point on globally tested distance curve'
 alltab=bh(pd.concat([a,rev,nich],ignore_index=True),['platform','family'])
 alltab['computation_status']=np.where(alltab.reason.eq('estimable'),'completed','not_estimable')
 alltab.loc[alltab.endpoint.str.startswith('CDF_'),'computation_status']='completed'
 alltab['annotation_status']=np.where(alltab.human_review_status.eq('pending'),'Codex-reviewed','partially human-reviewed');alltab['graph_version']='physical-visibility-v2';alltab['score_version']=np.where(alltab.platform.eq('HD'),'AUCell-historical-18085-K905','v1-UCell/frozen-H05-transfer')
 put(alltab,d/'endpoints.tsv.gz');pri=alltab[alltab.family.eq('primary')].copy();put(pri,d/'primary.tsv')
 fam=alltab.groupby(['platform','family']).agg(defined=('P','count'),rows=('P','size')).reset_index();put(fam,d/'test_families.tsv')
 sc=merged('score_QC.tsv').rename(columns={'high_single_gene_fraction':'high_mean_max_gene_share'});put(sc,d/'scores.tsv');ov=merged('high_overlap.tsv');put(ov,d/'high_overlap.tsv');put(merged('score_correlations.tsv'),d/'score_correlations.tsv')
 tools=[];checks=[];stable=[];results=[]
 for s in HD+XE:
  v,dd=folders(s);n=OUT/'niche'/s
  for name,evidence in [('Squidpy',dd/'squidpy_native_distance.npz'),('CRAWDAD',dd/'CRAWDAD_formal.tsv'),('spatstat',dd/('spatstat_envelope_'+('MP10' if s in HD else 'H05_transfer')+'.tsv')),('MENDER',n/'MENDER.parquet'),('BuildNicheAssay',n/'BuildNicheAssay.tsv.gz'),('CellCharter',n/'CellCharter_stability.tsv')]:
   tools.append((s,name,'completed' if evidence.exists() else 'pending',str(evidence.relative_to(OUT)),'pending'))
  for name,t in [('MENDER',read(n/'MENDER_stability.tsv')),('BuildNicheAssay',read(n/'BuildNicheAssay_stability.tsv'))]:
   if len(t):stable.append((s,name,int(t.loc[t.median_seed_ARI.idxmax(),'K']) if name=='MENDER' else 4,t.median_seed_ARI.max() if name=='MENDER' else t.ARI.median()))
  for k in ['all','masked']:
   p=n/f'CellCharter_{k}_ARI.txt'
   if p.exists():stable.append((s,'CellCharter '+k,4,float(p.read_text())))
  for branch,files in [('identity',['identity.parquet']),('state scores',['score_QC.tsv']),('center inference',['effects.tsv.gz']),('reverse inference',['reverse_effects.tsv.gz']),('review cards',[])]:
   ok=all((dd/f).exists() for f in files) if files else (paths(s)/'index.html').exists()
   checks.append((s,branch,'completed' if ok else 'pending','human-pending'))
  ff=OUT/'figures'/s;figs=sorted(ff.glob('*.png'));(ff/'index.html').write_text('<meta charset="utf-8"><h2>'+s+' · provisional</h2><ul>'+''.join('<li><a href="'+f.name+'">'+f.stem+'</a> · <a href="'+f.with_suffix('.pdf').name+'">PDF</a></li>' for f in figs)+'</ul>') if ff.exists() and not ff.is_symlink() else None
  results.extend((s,str(f.relative_to(OUT)),f.with_suffix('.pdf').is_file()) for f in figs)
 tt=pd.DataFrame(tools,columns=['sample','tool','compute_status','evidence','human_review_status']);put(tt,d/'tools.tsv');put(pd.DataFrame(checks,columns=['sample','branch','computation_status','annotation_validation']),d/'completion.tsv');st=pd.DataFrame(stable,columns=['sample','method','K','seed_ARI']);put(st,d/'niche_stability.tsv');put(pd.DataFrame(results,columns=['sample','PNG','PDF_present']),d/'figures.tsv')
 cs=merged('composition.tsv',True);no=merged('niche_overlay.tsv',True)
 if len(cs):put(cs,d/'niche_composition.tsv');put(no,d/'niche_overlay.tsv')
 native=merged('CRAWDAD_formal.tsv')
 if len(native):native['platform']=np.where(native['sample'].isin(HD),'HD','Xenium');put(bh(native,['platform']),d/'CRAWDAD_formal.tsv')
 globalcurve=merged('distance_envelope.tsv')
 if len(globalcurve):
  meta=['sample','res','definition','state','block'];gt=globalcurve.drop_duplicates(meta).rename(columns={'global_P':'P'});gt['platform']=np.where(gt['sample'].isin(HD),'HD','Xenium');put(globalcurve.merge(bh(gt, ['platform'])[meta+['FDR']],on=meta),d/'distance_global.tsv')
 resource=pd.read_csv(OUT/'audit/resources.tsv',sep='\t',header=None,names=['time','stage','exit','seconds','cgroup_peak_bytes','pids_peak']);put(resource,d/'resources.tsv')
 links=['<li><a href="'+s+'/index.html">'+s+'</a> · <a href="'+s+'/human_review.tsv">TSV</a> · <a href="'+s+'/atlas.html">全部窗口</a></li>' for s in HD+XE]
 (OUT/'review/index.html').write_text('<meta charset="utf-8"><h1>MP10–Treg 人工复核</h1><p>Provisional · 请独立填写，所有窗口隐去代谢状态。</p><ul>'+''.join(links)+'</ul>')
 summaries=read(d/'review_summary.tsv');loss=read(d/'gate_loss.tsv');text=['# Visium HD 与 Xenium 的 MP10–Treg 空间关系','Provisional / Codex-reviewed，人工待签署。以下按切片报告。','## 执行与输入',f'原 8 µm 表达、成员、全 RCTD 权重、Xenium 原分割、参考转移及兼容 UCell 分数均已复用。运行核查包含 {len(sc)} 个状态/尺度组合。记录到的计算 scope 峰值为 {resource.cgroup_peak_bytes.max()/1e9:.2f} GB，实际最高 {resource.pids_peak.max()} 个进程或线程。任务使用单 CPU 亲和性和单线程数值库，硬上限为每次 7 tasks / 165 GB，重计算阶段不重叠。峰值计入子进程与 cgroup 文件缓存。各阶段秒数及历史报错见 [资源表](reports/resources.tsv)。','## Treg 初审']
 if len(summaries):
  text += [md(summaries[['sample','candidates','supported','possible','unresolved','unsupported','reliable_T','old_supported','new_supported','human_requested']]),'三个 HD 数量为成分 bin，五个 Xenium 数量为分割细胞。低支持切片的个别目标可改变符号。[逐片留一表入口](reports/leave_one_summary.tsv)']
  if len(loss):
   t=loss[loss['sample'].isin(HD)&loss.old_gate.isin(['old_parent_not_TNK','old_Treg_not_top1','old_missing_anchor','old_old_combination'])];text+=[md(t.pivot(index='sample',columns='old_gate',values='new_supported').reset_index()),'旧门控原因可以重叠。[全部入选入口](reports/candidate_routes.tsv) · [全部流失原因](reports/gate_loss.tsv)']
 text+=['## 新 16/32 µm MP10 与 Xenium 三状态',md(sc[['sample','res','state','n','high','q95','high_multigene','QC']]),'Cervilla 和 D1_Johnson 的 16 µm high 在 80% 稀释后的集合重叠低于冻结界限。连续排名和 high 稳定性分别解释。D1_Johnson 的 32 µm 分组达到稳定性门槛。三套 Xenium 状态在同片共享合格上皮。逐片重叠及相关见 [Jaccard](reports/high_overlap.tsv) 与 [得分相关](reports/score_correlations.tsv)。','## 50 µm 主检验','每行的密度与 Treg/T 差值都按 high 减 low，后者为比例单位。q 是平台内 BH 校正结果。']
 pp=pri.pivot(index=['sample','state_definition','targets','score_QC'],columns='endpoint',values=['effect','FDR']);pp.columns=['_'.join(vv) for vv in pp.columns];text.append(md(pp.reset_index()))
 sig=pri[pri.FDR.lt(.05)]
 text+=['通过平台内主检验 FDR 0.05 的项目：'+('；'.join(str(r['sample'])+' '+r.state_definition+' '+r.endpoint+f' Δ={r.effect:.3g}，q={r.FDR:.3g}' for _,r in sig.iterrows()) if len(sig) else '无')+'。D1_Johnson 的二分类稳定性限制其解释。Xenium 的 A1R5 三定义原始密度差为负，A2R2 三定义为正，A1R6 随定义变化。逐片信号不能建立跨队列普遍正向富集。','## 暴露、距离与分割边界近接触']
 # Combine direction-corrected descriptors side by side without selecting the most favorable state.
 cols=['sample','state_definition'];defmain=a[(a.identity_definition=='supported')&(a.block_um==400)&(a.unit_resolution.astype(str)!='32')]
 companion=pri[pri.endpoint=='density'][cols+['effect']].rename(columns={'effect':'density Δ'})
 def add(df,t,label,col='effect'):return df.merge(t[cols+[col]].rename(columns={col:label}),on=cols,how='left')
 companion=add(companion,pri[pri.endpoint=='Treg/T'],'Treg/T Δ')
 for met,r,title in [('restricted distance um',200,'Epi→Treg Δ µm'),('any contact',2,'≤2 µm 接触概率 Δ')]:companion=add(companion,defmain[(defmain.endpoint==met)&(defmain.radius_or_contact_um==r)],title)
 if len(rev):
  rr=rev[(rev.identity_definition=='supported')&(rev.unit_resolution.astype(str)!='32')]
  for met,r,title in [('reverse high exposure',50,'Treg中心 high 暴露−零期望'),('reverse restricted distance to high um',200,'Treg→high 距离−零期望')]:companion=add(companion,rr[(rr.endpoint==met)&(rr.radius_or_contact_um==r)],title)
 text+=[md(companion),'负距离差表示更近。反向暴露与最近 high 距离的零模型保留每层 high 数量。近接触仅限 Xenium 的有效原边界。不同几何读数有不同有效分母；零分母及无可达目标见逐中心文件。[完整终点与区间](reports/endpoints.tsv.gz) · [有效覆盖与 T 分母](reports/denominators.tsv)']
 text+=['## 三个空间工具','Squidpy 全局标签打乱衡量固定图上的互斥身份聚集；HD 的独立成分属性通过专项统计保留。CRAWDAD 固定焦点、改变局部标签，spatstat 在上皮内局部分层移动状态。各工具对应不同零假设。方向与校验见 [逐片工具对照](reports/tool_comparison.tsv)。','## 三个 niche 方法']
 if len(st):
  text.append(md(st.pivot(index='sample',columns='method',values='seed_ARI').reset_index()))
 if len(cs):
  pct=cs[cs.method=='MENDER'].sort_values('fraction',ascending=False).drop_duplicates(['sample','niche']);rows=[]
  for s,ss in pct.groupby('sample'):
   dominant=', '.join(str(r.niche)+' '+r.identity+' '+f'{r.fraction:.0%}' for _,r in ss.sort_values('niche').iterrows());rows.append((s,dominant))
  text+= [md(pd.DataFrame(rows,columns=['sample','MENDER 各区的首位身份及其占比'])),'Mixed/Unknown 主导的区首先反映身份覆盖；深度与组织边缘也随 niche 变化。niche 在揭示代谢状态之前冻结，随后逐片叠加 high 份额、密度和 Treg/T。可在 [niche 覆盖与叠加](reports/niche_overlay.tsv)、[组成相似性](reports/niche_similarity.tsv) 与统一长表的 within MENDER niche 家族审阅。']
 text+=['## 人工复核与计算状态','首批 [浏览入口](review/index.html) 含三片 HD 各 30 个和五片 Xenium 各 12 个，共 150 个。优先复核 HD D1_Johnson 的非 T/NK 主要身份但有 Treg 成分者，以及 A1R4、A1R6 稀疏目标的边界和竞争线索。Xenium 的 Venous EC、SPP1 TAM、SMC 缺少原 anchor，各片按 panel-limited 保留。H05 为转移支持。现有数据支持优先校准 Treg 身份及高分尾部稳定性，之后检验候选组织区域中的具体关系。','[计算状态](reports/completion.tsv) · [六工具逐片状态](reports/tools.tsv) · [图件索引](reports/figures.tsv) · [方法与复现](METHODS.md)']
 (OUT/'REPORT.md').write_text('\n\n'.join(text)+'\n')
if __name__=='__main__':report()
