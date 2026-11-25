from common import *
T=Path('tables'); R=lambda n:rd(T/n);co=R('celloracle_all_genes.tsv.gz');cp=R('celloracle_program_response.tsv');cross=R('cross_method_summary.tsv');main=cross[cross.analysis=='main'];fit=R('celloracle_node_audit.tsv');q=R('celloracle_simulation_qc.tsv');jobs=R('celloracle_job_status.tsv');sc=R('control_specificity.tsv');controls=R('control_matching.tsv')
pr=R('program_response.tsv');save(pd.concat([pr.assign(method='scTenifoldKnk'),cp.assign(method='CellOracle')],ignore_index=True),T/'program_response_all_methods.tsv.gz');mt=R('mtor_module_summary.tsv');program=R('program_members.tsv');ps=R('program_stability.tsv');e=R('endpoint_stability.tsv')
edges=R('celloracle_edges_and_paths.tsv');paths=edges[(edges.record_type=='path')&(edges.graph=='native')&(edges.steps==3)];ed=edges[edges.record_type=='fitted']
# Predefined five controls: count only legal, reachable, nonself simulations at that target.
parts=[];rr=[]
for (job,graph,steps),group in co.groupby(['analysis','graph','steps']):
 if graph!='native' or steps!=3:continue
 for k in K[:2]:
  a=group[(group.ko==k)&~group.is_KO].set_index('gene');c=group[group.ko.isin(controls[controls.ko==k].control)&~group.is_KO].pivot(index='gene',columns='ko',values='absolute_delta').reindex(a.index)
  a=a.assign(n_fixed_controls=5,n_evaluable_controls=c.count(1),median_control_absolute_delta=c.median(1),n_controls_lower=c.lt(a.absolute_delta,axis=0).sum(1).where(a.absolute_delta.notna()))
  parts.append(a.reset_index())
  p=cp[(cp.analysis==job)&(cp.graph==graph)&(cp.steps==steps)]
  y=p[p.ko.isin(controls[controls.ko==k].control)&p.status.eq('EVALUABLE')].groupby(['program','mask']).agg(n_evaluable_controls=('ko','nunique'),control_neg=('negative_fraction','median'),control_rank=('median_rank_fraction','median')).reset_index()
  rr.append(p[p.ko==k].merge(y,on=['program','mask'],how='left'))
oc=pd.concat(parts);save(oc,T/'celloracle_specificity.tsv.gz');pc=pd.concat(rr);save(pc,T/'celloracle_program_controls.tsv')
# Human-reader long tables and quantitative stability across independent refits.
st=cross[cross.ko.isin(K[:2])&cross.gene.isin(E+list(program.gene))].copy();st['type']=np.select([st.analysis.str.startswith('drop_'),st.analysis.eq('seed456')],['LODO','alternate seed'],default='main')
br=st[st.analysis=='main'].set_index(['ko','gene']).signed_delta
st['same_sign_as_main']=[(np.sign(v)==np.sign(br.get((k,g)))) if pd.notna(v) and pd.notna(br.get((k,g))) else np.nan for k,g,v in zip(st.ko,st.gene,st.signed_delta)]
save(st,T/'cross_method_endpoint_stability.tsv.gz')
# Full and M20 recovery. Include empty MP10 M20. Primary path support separate from native distance.
pco=cp[(cp.graph=='native')&(cp.steps==3)&cp.ko.isin(K[:2])].copy();pco['refit_type']=np.where(pco.analysis.str.startswith('drop_'),'LODO',pco.analysis)
save(pco.groupby(['ko','program','mask','refit_type']).agg(n_models=('analysis','size'),n_evaluable=('negative_fraction','count'),negative_min=('negative_fraction','min'),negative_max=('negative_fraction','max'),rank_fraction_min=('median_rank_fraction','min'),rank_fraction_max=('median_rank_fraction','max')).reset_index(),T/'celloracle_program_stability.tsv')
match=pd.concat([R('knk_job_status.tsv').assign(method='scTenifoldKnk',state=''),jobs.rename(columns={'cells':'n_cells','patients':'n_patients'}).assign(method='CellOracle')]);save(match,T/'task_status.tsv')
native=co[(co.analysis=='main')&(co.graph=='native')&(co.steps==3)&~co.is_KO]
def fmt(v,d=3):
 if pd.isna(v):return 'NE'
 return f'{float(v):.{d}f}'
def pct(v):return fmt(v*100,1)+'%' if pd.notna(v) else 'NE'
def table(headers,rows):
 return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,row))+' |' for row in rows])
def pv(t,k,n,mask='Full'):
 return t[(t.ko==k)&(t.program==n)&(t['mask']==mask)]
def pcrow(k,n):return pv(pc[pc.analysis=='main'],k,n).iloc[0]
def xrow(k,g):return main[(main.ko==k)&(main.gene==g)].iloc[0]
def end(k,g):return sc[(sc.ko==k)&(sc.gene==g)].iloc[0]
def range2(v,p=False):
 v=v.dropna()
 return 'NE' if not len(v) else (pct(v.min())+'–'+pct(v.max()) if p else fmt(v.min())+'–'+fmt(v.max()))
s=['# MP10–SREBF–mTOR–NECTIN 虚拟扰动 v2','8,000 个原上皮细胞、38 位记录患者、39 个标本、5 个来源。两种模型使用同一观察性队列。scTenifoldKnk 完成 38 个正式 WT，C9 因仅 29 个细胞、2 位患者而记 NE。CellOracle 已完成小试、主种子、替代种子和五组按来源留一的簇级模型。正式簇任务共 %d 个成功、%d 个 NE。'%(sum(jobs.status=='complete'),sum(jobs.status=='NE'))]
s+=['## 去除核糖体和线粒体模块后的前排基因']
g=R('ranking_comparisons.tsv');other=g[(g.top_n==50)&(g.right_model=='other_main_KO')];shared=R('shared_response.tsv')
s+=['在 P_full，五个主扰动的前 50 名中有 21–22 个细胞质核糖体结构蛋白及 0–8 个线粒体基因组基因。P_clean 精确去掉 39 个面板成员并重新建网。五个主扰动仍有 %d 个基因全部进入前 50，任意两者的全排名 Spearman 为 %s，前 50 两两重合 %d–%d 个。'%(sum(shared.n_main_top50==5),range2(other.spearman),other.top_intersection.min(),other.top_intersection.max())]
rows=[]
for k in K:
 z=sc[sc.ko==k].nsmallest(5,'rank');y=native[native.ko==k].nsmallest(5,'rank')
 rows.append([k,'、'.join(z.gene),'、'.join(y.gene) if len(y) else 'NE（非 TF）'])
s+=[table(['扰动源','scTenifold 前五（距离）','CellOracle 前五（患者等权绝对 Δ）'],rows),'PRDX1、NME2 均进入五个主源和全部 25 个匹配对照的前 50。S100A14 和 CALR 同时进入 SREBF2 与 RICTOR 的前 50；在此阈值未见仅某一主源独有的前排候选。scTenifold 适合提示共有响应，单靠其排名难以选出主源特异靶。']
s+=['## A　SREBF 与 MP10/脂质转录程序','两种模型均给出 MP10 成员总体向下的预测，支持程序方向；特异性和成员层面证据有差别。scTenifold 主模型中，两源均为 20/20 向下，20 次患者子抽样均维持这一比例，位次会漂移。两组匹配对照的程序负向比例中位数也都是 100%。']
rows=[]
for k in K[:2]:
 for n in ['MP10','CHOL_SYNTH','SCENIC_SREBF1','SCENIC_SREBF2']:
  y=pcrow(k,n);rep=pv(pc[pc.analysis.str.startswith('drop_')],k,n)
  rows.append([k,n,f'{y.n_evaluable}/{y.n_frozen}',pct(y.negative_fraction),pct(y.median_rank_fraction),pct(y.control_neg),range2(rep.negative_fraction,True)])
s+=[table(['源','程序 Full','CellOracle 可评估/冻结','负向成员','排名中位分位','可评估对照负向中位','五次留一负向范围'],rows),'覆盖门槛为至少 3 个成员且达到冻结集合的一半。NE 程序仍公开每个可达成员。M20 结果及成员见程序表；MP10 的 M20 无成员。']
for k,t in [('SREBF1','SCD'),('SREBF2','HMGCS1')]:
 x=xrow(k,t);de=ed[(ed.analysis=='main')&(ed.source==k)&(ed.target==t)]
 s.append(f'{k} 对 {t} 的 CellOracle 患者等权有符号 Δ 为 {fmt(x.signed_delta,4)}，绝对变化排名 {int(x["rank"])}/{int(x.rank_denominator)}。{len(de)} 个主模型簇保留该直接推断边，其中 {sum(de.coefficient>0)} 个系数为正。')
s+=['## B　RPTOR、MTOR、RICTOR 与脂质程序','三者与 SREBF 相关脂质程序的负向转录轮廓相容。对 C1 特异性的支持不足。下表在五源共同可评估的非源 MP10 成员上比较，低排名分位代表更靠前。']
p=mt[(mt.kind=='program')&(mt.readout=='MP10')&(mt['mask']=='Full')]
s+=[table(['源','MP10 负向','成员排名中位分位','五对照排名中位分位','20次重拟合排名中位分位范围'],[[r.ko,pct(r.negative_fraction),pct(r.median_rank_fraction),pct(r.control_rank_fraction),pct(r.refit_rank_fraction_min)+'–'+pct(r.refit_rank_fraction_max)] for r in p.itertuples()]),'RPTOR 和 MTOR 未显示区别于 RICTOR 的一致程序选择性。mTORC1 Hallmark 含 200 个成员，各主源与各自对照的排名也接近。MTOR 同时参与两种复合体。LPIN1、RPS6KB1/2、EIF4EBP1 的基因行均在模块表；它们的 RNA 读出不能确定底物磷酸化。CellOracle 标准 TF 模型不评估这三个非 TF 源。']
s+=['## C　SREBF 与 NECTIN/PVR','NECTIN2 在两种主模型中的有符号预测均向下。scTenifold 距离未通过 BH，CellOracle 的路径、幅度、对照和留一结果需要分开看。']
rows=[]
for k in K[:2]:
 for gene in ['NECTIN2','NECTIN1','NECTIN3','NECTIN4','PVR']:
  x=xrow(k,gene);z=oc[(oc.analysis=='main')&(oc.ko==k)&(oc.gene==gene)].iloc[0];pa=paths[(paths.analysis=='main')&(paths.source==k)&(paths.target==gene)&(paths.reachable==True)];l=st[(st.ko==k)&(st.gene==gene)&(st.type=='LODO')]
  rows.append([k,gene,int(x.sc_rank),fmt(x.sc_q),fmt(x.signed_delta,4),'NE' if pd.isna(x['rank']) else str(int(x['rank'])),f'{len(pa)}/9',str(int(z.n_controls_lower))+'/'+str(int(z.n_evaluable_controls)) if pd.notna(z.n_controls_lower) else 'NE',f'{sum(l.signed_delta<0)}/{l.signed_delta.count()}'])
s+=[table(['源','终点','sc 排名/2193','sc q','CellOracle Δ','CO 排名','可达簇','低于主源的有效 TF 对照','留一 Δ<0/可评估'],rows)]
for k in K[:2]:
 r=end(k,'NECTIN2');z=xrow(k,'NECTIN2');subset=ed[(ed.analysis=='main')&(ed.source==k)&(ed.target=='NECTIN2')];pa=paths[(paths.analysis=='main')&(paths.source==k)&(paths.target=='NECTIN2')&(paths.reachable==True)]
 s.append(f'{k} 的 NECTIN2 在 scTenifold 原生、回归测序深度、再去 PC1 后的方向分数为 {fmt(r.directionScore,5)}、{fmt(r.depth_score,5)}、{fmt(r.pc1_score,5)}；距离超过 {int(r.n_controls_lower)}/5 个匹配对照。CellOracle 在 {int(z.n_patients)} 位患者有可评估细胞，{len(subset)} 个主簇保留直接边，可达簇最短路径为 '+('、'.join(str(int(i)) for i in sorted(pa.path_length.dropna().unique())) or 'NE')+' 步。')
s+=['映射至面板的官方启动子先验未列出 SREBF1/2→NECTIN2 直接边。两源的 NECTIN2 在 20 次 scTenifold 患者子抽样和 5 次来源留一中均保留原生负方向，均无 q<0.05。SREBF1 去除 PC1 后方向翻转。合并考虑，NECTIN2 获得弱的计算方向线索，直接调控和源特异性未获确证。SREBF2–NECTIN1 在两种主模型间方向相反；表中其余成员按各自覆盖与数值解释。']
s+=['## 诊断、分层和可复核边界']
for k in K[:2]:
 h=R('endpoint_model_results.tsv.gz');y=h[(h.analysis=='clean')&(h.ko==k)&(h.gene=='NECTIN2')&(h.branch=='manifold')].iloc[0];c7=h[(h.analysis=='cluster_C7')&(h.ko==k)&(h.gene=='NECTIN2')&(h.branch=='heat')].iloc[0]
 z=co[(co.analysis=='main')&(co.ko==k)&(co.gene=='NECTIN2')];vals=', '.join([str(int(v.steps))+' 步 '+fmt(v.signed_delta,4) for v in z[z.graph=='native'].sort_values('steps').itertuples()]);ran=z[z.graph=='randomized'].iloc[0]
 s.append(f'{k} 的经典 manifold 对 NECTIN2 排名 {int(y["rank"])}/{int(y.rank_denominator)}，q={fmt(y["p.adj"])}；C7 独立 heat 模型（698 个细胞、28 位患者）排名 {int(c7["rank"])}/{int(c7.rank_denominator)}，方向 {c7.direction}。CellOracle 主簇合并后，{vals}；官方随机图 3 步为 {fmt(ran.signed_delta,4)}。各步和随机图的患者分母可不同。')
s+=['冻结分层的 SREBF1/2–NECTIN2 配对 High 5%、10%、50% 与 Low 50% 分别有 23、34、35 位共同合格患者。两源各分层 High−Low 距离差的患者中位数均为正，个别患者为负。该距离反映固定 WT 加上层内表达分布的条件读出。',f'CellOracle 空扰动最大 |Δ| 为 {q.no_perturbation_max_abs.max():.2g}。主运行保留未额外裁剪的输出，实际越出 WT 范围和取零计数见 celloracle_simulation_qc。人工添加边为零，P/q 空白表示该算法未提供该统计量。', '## 后续验证优先级','1. 优先检查 SREBF1–SCD、SREBF2–HMGCS1。在独立材料中区分转录响应与蛋白/脂质输出。这两对有较靠前的 CellOracle 位次及跨簇同号拟合边。','2. 若比较 RPTOR 与 RICTOR，应同步检查核内成熟 SREBP、复合体底物磷酸化及脂质读出。当前 RNA 模型无法区分 C1 与 C2 特异性。','3. NECTIN2 优先级低于上述脂质读出。若继续，分别测 RNA 与细胞表面蛋白，并检验直接调控证据。本数据未检验 NECTIN2–TIGIT 结合、Treg 抑制、胆固醇通量或药敏。','## 文件导航','主结果在 tables/。knk_all_gene_index.tsv 和 celloracle_cluster_index.tsv 定位每次真实拟合的全基因记录；celloracle_all_genes.tsv.gz 为跨簇患者等权排名。control_specificity、celloracle_specificity、cross_method_summary、mtor_module_summary 以及两种方法的 program 与 stability 表分别回答特异性、跨法、程序和重拟合问题。图见 FIGURE_INDEX.tsv，输入、范围和复现入口见 METHODS.md。']
Path('RESULTS.md').write_text('\n\n'.join(s)+'\n')

