from common import *
from report import md,read
def fmt(x):return f'{x:.3g}' if np.isfinite(x) else 'NA'
def interpret():
 d=OUT/'reports';pri=read(d/'primary.tsv');ct=read(d/'tool_comparison.tsv');ni=read(d/'niche_overlay.tsv');nu=read(d/'niche_quality.tsv');cv=read(d/'CellCharter_convergence.tsv');a=read(d/'endpoints.tsv.gz');st=read(d/'niche_stability.tsv')
 pr=pri[(pri['sample']=='D1_Johnson')&(pri.endpoint=='density')].iloc[0];test=[f'D1_Johnson 的 50 µm high−low 密度差为 {pr.effect:.2f} 个支持 bin/mm²，区组重采样 95% 区间 [{pr.CI_low:.2f}, {pr.CI_high:.2f}]。受限置换的中心为 {pr.null_mean:.2f}，双侧 P={pr.P:.4g}，HD 族 q={pr.FDR:.4g}。置换检验衡量相对分层零模型的偏离，区间描述切片内区组重采样的原始差值。上皮 high 的降深稳定性为低支持。']
 xp=pri[pri.platform=='Xenium'];test.append(f'Xenium 平台 {len(xp)} 个主检验中最小 BH q={xp.FDR.min():.4g}。A1R5 的三种 high 在 50 µm 均呈较低 Treg 密度，200 µm 内到固定 Treg 集合的限制平均距离均较长。A2R2 三种 high 的密度与 Treg/T 都较高、距离较短、≤2 µm 分割边界近接触比例较高，但全片主定义目标仅 13 个。A1R4 有 4 个，A1R6 有 5 个；A1R6 的 H05 高转移支持有局部正向信号，合成读数的密度方向可随单个 Treg 的剔除翻转。')
 xe=read(d/'high_overlap.tsv');corr=read(d/'score_correlations.tsv');rows=[]
 for _,t in xe.groupby(['state1','state2']):
  a1=t.state1.iloc[0];b=t.state2.iloc[0];cr=corr[(corr.state==a1)&(corr.variable==b)];rows.append((a1,b,t.Jaccard.min(),t.Jaccard.max(),cr.Spearman.min(),cr.Spearman.max()))
 text=['# 结果判读与方法校验','## 主检验与跨片方向',*test,'## 三个 Xenium 状态','下表为五片共同合格上皮上的范围。连续秩相关与极端尾部交集分别报告。',md(pd.DataFrame(rows,columns=['状态 A','状态 B','Jaccard 最小','Jaccard 最大','Spearman 最小','Spearman 最大'])),'## 原生工具、零模型与校验']
 if len(ct):
  c=ct[ct['sample']=='A1R5'];col=[k for k in ['sample','state','density_delta','density_q','Z','OE','CRAWDAD_Z_100','CRAWDAD_Z_800','curve_q','spatstat_outside_global_band'] if k in c];text.append(md(c[col]))
  text.append('Squidpy 的目标 Z 对比切片全局身份打乱。CRAWDAD 在去重的焦点缓冲区内计算上皮状态相对局部打乱的 Z，允许不同类型坐标组成共同变化。冻结图的主检验只交换合格上皮内状态。数值和显著性不能互换。')
 native=[];cl=[]
 for s in HD+XE:
  x=read(folders(s)[1]/'spatstat_primary.tsv')
  if len(x):native.append(x)
  y=read(folders(s)[1]/'CRAWDAD_native_check.tsv')
  if len(y):y['sample']=s;cl.append(y)
 if native:
  x=pd.concat(native,ignore_index=True);x['platform']=np.where(x['sample'].isin(HD),'HD','Xenium');from report import bh
  x=bh(x,['platform']);x=x.merge(pri[['sample','state_definition','endpoint','effect','P','FDR']].rename(columns={'P':'project_P','FDR':'project_FDR','state_definition':'state'}),on=['sample','state','endpoint']);assert np.allclose(x.observed,x.effect,atol=1e-8)
  put(x,d/'spatstat_primary_check.tsv');text.append(f'spatstat::rlabel 的独立 4,999 次联合随机标记对 {len(x)} 个主读数进行重算。观测差值最大绝对误差 {abs(x.observed-x.effect).max():.3g}。新随机化样本的 P 和 BH 另表保留。[原生受限置换核对](spatstat_primary_check.tsv)')
 if cl:
  cl=pd.concat(cl,ignore_index=True);text.append(f'CRAWDAD 的原生窗口与稀疏全片实现比较了 {len(cl)} 个类别/尺度，{cl.absolute_error.notna().sum()} 个有有限 Z。有限结果的最大绝对误差为 {fmt(cl.absolute_error.max())}。原生 Z 遇到任一置换的合并比例为 0 或 1 时缺失，另列可定义置换次数和按有限值条件计算的 Z。[完整校验](CRAWDAD_kernel_checks.tsv)');put(cl,d/'CRAWDAD_kernel_checks.tsv')
 text+=['## Niche 的生物学可解释范围','Cervilla 与 D1_Johnson 均有 Matrix_CAF 占首位的 MENDER 区，也有上皮、Mixed 和 Unknown 占首位的区。UC_SCC 及五片 Xenium 的首位类别常为上皮、Unknown 或 Mixed。身份未解析率与原始深度在不同区内变化。可复现的是这些组成差异，当前不足以给所有同号区赋予统一生物学命名。']
 if len(ni):
  n=ni[(ni.method=='MENDER')&(ni.endpoint=='density')];desc=n.groupby(['sample','state']).agg(niches=('niche','size'),Epi_high_min=('high_fraction','min'),Epi_high_max=('high_fraction','max'),supported_Treg_min=('targets','min'),supported_Treg_max=('targets','max')).reset_index()
  text.append(md(desc));f=a[a.family=='within MENDER niche'];text.append(f'niche 内共有 {f.P.notna().sum()} 个可评主定义次要检验，BH q<0.05 的项目为 {f.FDR.lt(.05).sum()} 个。high 分布和 Treg 所在 niche 的重叠仍属描述。')
 if len(cv):text.append(f'CellCharter 按固定预算迭代，{cv.converged.sum()}/{cv.fits.sum()} 次拟合在 40 轮内满足原生收敛判据。未收敛的结果作为有限配置的补充，解释还需结合种子 ARI 和屏蔽目标基因的敏感性。[收敛表](CellCharter_convergence.tsv)')
 text+=['## 最少必要人工任务','首批 150 个对象均保留原始表、独立坐标和局部图。16 个对象有 Codex 的逐图记录，其余候选为按冻结规则逐行完成的 table_only 判读。人工签署数为零。优先判断 D1_Johnson 混合 bin 中 FOXP3 与 T 谱系的归属，以及 A1R4/A1R6 稀少目标的边界与活化 T 竞争。Xenium 的六个 pan-T 候选仅 CD3E、CD3G 实测。完成这些身份核对后，以独立测量复核胆固醇合成程序；A1R5 和 A2R2 可作为方向相反的组织区域进行后续验证。[便携复核包](../review/first_batch.zip)']
 (d/'INTERPRETATION.md').write_text('\n\n'.join(text)+'\n')
 rp=OUT/'REPORT.md';main=rp.read_text()
 if 'reports/INTERPRETATION.md' not in main:rp.write_text(main+'\n[效应区间、三个状态、工具核验与 niche 的逐片判读](reports/INTERPRETATION.md) · [首批便携复核包](review/first_batch.zip)\n')
if __name__=='__main__':interpret()
