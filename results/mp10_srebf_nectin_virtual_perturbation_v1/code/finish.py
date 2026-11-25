from common import *
import re,subprocess
from scipy.stats import spearmanr
S=rd('tables/virtual_perturbation_summary.tsv');R=rd('tables/virtual_perturbation_robustness.tsv');B=rd('tables/refit_ranges.tsv')
K=['SREBF1','SREBF2'];x=S[(S.analysis_id=='main')&(S.endpoint_gene=='NECTIN2')].set_index('perturbed_gene')
res=rd('tables/resource_usage.tsv');inp=rd('logs/stages.tsv').iloc[-1];scope=rd('tables/eligible_metadata.tsv.gz');m=scope[scope.main];ids=scope.patient_id.nunique()
c=rd('plot_data/F04_controls.tsv');cov=rd('plot_data/F01_coverage.tsv').set_index('gene')
state=S[(S.analysis_id.str.startswith('state_'))&(S.endpoint_gene=='NECTIN2')]
data=rd('qc/jobs/main_effects.tsv.gz')
h=data.pivot(index='gene',columns='perturbed_gene',values='distance')
comp=[]
for k in K:
 group=c[(c.matched_to==k)&(c.role=='Control')]
 for ctl in group.perturbed_gene:
  sub=h.drop(index=[k,ctl]);comp.append(dict(ko=k,control=ctl,profile_distance_spearman=spearmanr(sub[k],sub[ctl]).statistic,target_distance=x.loc[k,'effect_value'],control_distance=x.loc[ctl,'effect_value'],target_rank=x.loc[k,'rank'],control_rank=x.loc[ctl,'rank']))
save(pd.DataFrame(comp),'tables/control_specificity.tsv')
pe=rd('tables/patient_conditioned.tsv')
pat=pe[pe.stratum=='All'].groupby('perturbed_gene').effect_value.agg(['count','median','min','max'])
prog=rd('tables/program_response.tsv'); prog=prog[prog.analysis_id=='main']
def fmt(x):return f'{x:.5g}'
lines=['# 上皮 SREBF1/2 单基因虚拟扰动','',
'**整体判断**　主网对两个单 KO 都给出 NECTIN2 向下方向，但该基因的调控距离在两次检验中均未通过模型内多重检验。对照扰动也能产生相近排名，支持有限。原簇 C7 的两个方向均与合并主网相反。','',
f'原始缓存有 300,033 个细胞、90,011 个特征。冻结合格 Epi 共 {len(scope):,} 个，来自 {scope.sample_uid.nunique()} 份标本、{ids} 个记录患者及 5 个来源。与元数据核对的 UMI、检出基因数完全一致。已有簇覆盖 {(scope.state!="Unassigned").sum():,} 个细胞。正式主网为 {len(m):,} 个细胞、{m.sample_uid.nunique()} 份标本、{m.patient_id.nunique()} 个记录患者、2,000 个基因。T4450398 缺原簇标签，T4463209a/b 共用同一患者键。','',
'| 扰动 | NECTIN2 距离 | 距离排名 / 1999 | 方向评分 | 模型内 P | 模型内 BH q |','|---|---:|---:|---:|---:|---:|']
for k in K:
 a=x.loc[k];lines.append(f'| {k} | {fmt(a.effect_value)} | {int(a["rank"])} | {fmt(a.direction_score)} | {a.p_value:.3f} | {a.q_value:.3f} |')
lines+=['','## 稳定性、患者与对照','']
for k in K:
 sub=B[(B.perturbed_gene==k)&(B.analysis_type=='patient_subsample')].iloc[0];l=B[(B.perturbed_gene==k)&(B.analysis_type=='LODO')].iloc[0];alt=R[(R.perturbed_gene==k)&(R.analysis_id=='alt')].iloc[0]
 ctl=c[(c.matched_to==k)&(c.role=='Control')]
 low=ctl.effect_value.min();high=ctl.effect_value.max();pp=pat.loc[k];ss=state[state.perturbed_gene==k]
 lines.append(f'- **{k}**　患者子集 {int(sub.n_down)}/{int(sub.n_evaluable)} 次向下，距离 {fmt(sub.distance_min)}–{fmt(sub.distance_max)}，排名 {int(sub.rank_min)}–{int(sub.rank_max)}。逐一剔除来源后 {int(l.n_down)}/{int(l.n_evaluable)} 次向下，距离 {fmt(l.distance_min)}–{fmt(l.distance_max)}；固定种子 456 的方向为 {alt.predicted_direction}。患者条件距离中位数 {fmt(pp["median"])}，观察范围 {fmt(pp["min"])}–{fmt(pp["max"])}，共 {int(pp["count"])} 位。五个匹配对照的 NECTIN2 距离为 {fmt(low)}–{fmt(high)}，其中 {(ctl.effect_value<x.loc[k,"effect_value"]).sum()}/5 小于主 KO。')
lines+=['',f'两个 KO 的深度调整以及再除去 PC1 后，主网方向评分均保持负号。剔除 InhouseData 时，SREBF1 的原始评分仍为负，深度回归和再去 PC1 后转正；SREBF2 在五种剔除设置下均保持负号。可评估的九个既有簇中，八个向下、C7 向上；C9 只抽入 29 个细胞、2 位患者，依据预设门槛标 NE。固定 WT 条件下的患者间差异以患者为单位描述。仅 {pe[(pe.stratum=="H05")&(pe.perturbed_gene=="SREBF1")].effect_value.notna().sum()} 位患者在 H05 层达到患者条件距离门槛，详细分母见表。','',
f'两组五对匹配实际涉及 {c[c.role=="Control"].perturbed_gene.nunique()} 个不同对照。全基因距离排序与配对对照的 Spearman 相关中位数为 {np.median([v["profile_distance_spearman"] for v in comp]):.3f}。SREBF2 虽有较大的绝对距离，其 NECTIN2 在全基因中的位置与 SREBF1 几乎相同。SREBF2 与若干对照在表达量、检出率上匹配较松；结果只能用于描述特异性。','',
'## 检出与校准','']
for g in ['SREBF1','SREBF2','NECTIN2']:
 a=cov.loc[g];lines.append(f'- {g} 有精确原始 symbol，{int(a.n_positive):,}/{int(a.n_cells):,} 个合格上皮细胞检出（{a.detect_frac:.1%}），涉及 {int(a.n_positive_patients)} 个记录患者。该基因进入主模型。')
lines+=['','两种 KO 的 HMGCR、SQLE、LDLR 方向均向下，各单基因的模型内 BH q 均大于 0.05。校准程序的方向一致性应与网络扰动幅度分开解读。','',
'| 主 KO | 程序 | 口径 | 可用 / 原成员 | 覆盖成员的距离中位数 |','|---|---|---|---:|---:|']
for a in prog[prog.program.isin(['MP10','CHOL_SYNTH','SCENIC_SREBF1','SCENIC_SREBF2'])].itertuples():lines.append(f'| {a.perturbed_gene} | {a.program} | {a.mask} | {a.n_covered}/{a.n_frozen} | {fmt(a.effect_value) if pd.notna(a.effect_value) else "NE"} |')
lines+=['','NECTIN1/3/4、PVR、CollecTRI 和 DoRothEA 的完整结果、覆盖和 NE 原因见主表及程序表。MP10 的 M20 口径无剩余成员。所有程序统计量均为覆盖成员的距离摘要，不能据此推断胆固醇含量变化。','',
'## 计算与交付','']
usage=res.set_index('analysis_id');rss=rd('logs/simultaneous_memory.tsv').total_RSS_bytes.max()/2**30
fitted=res[res.reason=='fitted'];rr=usage.loc['main']
lines+= [f'输入抽取 {inp.wall_s:.1f} 秒、CPU {inp.cpu_s:.1f} 秒、峰值 {inp.process_highwater_RSS_KiB/2**20:.2f} GiB。小试 2,500 个细胞从 2,000 个候选保留 {int(usage.loc["pilot","n_genes"])} 个，{usage.loc["pilot","wall_s"]:.1f} 秒。正式主作业（含两个单 KO、匹配对照及患者导出）{rr.wall_s:.1f} 秒、CPU {rr.cpu_s:.1f} 秒，峰值 {rr.peak_RSS_KiB/2**20:.2f} GiB。并行作业定时观测到的合计 RSS 最高 {rss:.2f} GiB。成功建网 {len(fitted)} 个分析配置，含小试 1、主网 1、替代种子 1、留一来源 5、状态 9、患者子集 20，另 1 个状态按门槛跳过。总 CPU 和各作业用时见 resource_usage.tsv；并行作业用时不能直接相加当作墙钟时间。','',
'首次试运行曾遇到包内部函数导出名问题；首次主网与替代种子作业在写出表后遇到脚本读取错误。两项均按相同配置重跑，发布验证要求最终作业成功退出。分阶段耗时、未成功作业的记录、参数与输入指纹分开保留。','',
'CellOracle 的状态位移与转移概率没有生成。本次没有测量真实 KO 后的表达、SREBF 结合、NECTIN2 蛋白或 Treg 的跨细胞功能。就既定模型而言，整体上可重复的负方向与有限的终点特异性可以同时成立。','',
'只读入口为 config/virtual_perturbation_manifest.json；三个 TSV 的用途分别为主估计、逐次重拟合与覆盖审计。FIGURE_INDEX.tsv 对应每张图、图注及直接数据。']
Path('RESULTS.md').write_text('\n'.join(lines)+'\n')

