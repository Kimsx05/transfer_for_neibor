from pathlib import Path
import csv,gzip,math,json,hashlib
from refresh_virtual import D,run as refresh
T=D/'tables'; ext=['GSE32894','GSE31684','GSE13507']
def read(p):
 with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return list(csv.DictReader(f,delimiter='\t'))
def put(p,rows):
 with p.open('w') as f:
  w=csv.DictWriter(f,list(dict.fromkeys(k for row in rows for k in row)),delimiter='\t');w.writeheader();w.writerows(rows)
def ff(v,d=2):
 try:return f'{float(v):.{d}f}' if math.isfinite(float(v)) else 'NE'
 except:return 'NE'
def qq(v):
 try:return f'{float(v):.2g}'
 except:return 'NE'
def tab(header,rows):return '\n'.join(['| '+' | '.join(header)+' |','|'+'|'.join(['---']*len(header))+'|']+['| '+' | '.join(str(t) for t in row)+' |' for row in rows])
def interval(u,key='r',low='lo',high='hi'):return ff(u[key])+' ['+ff(u[low])+', '+ff(u[high])+']'
z=read(T/'all_cohorts.tsv');meta=read(T/'meta_external.tsv')
desc={'TCGA':'PanCancer Atlas BLCA，每个记录患者保留一份原发标本。原始 Pearson 表的 Full 层为 407 人，纯度子集 397 人。原图主展示 Spearman，本轮从原统计文件恢复 Pearson。','GSE32894':'308 份肿瘤均纳入，包括 213 NMIBC、93 MIBC 和 2 份 Tx。随访排除标记不参与横断面纳入，首次诊断状态未逐例核实。按来源重建未中心化的强度背景，额外正常与技术重复仅用于归一化。辅助表中 UC_0169_1.P5 与强度列 UC_0169_1.P4 不匹配，未将二者配对。13,849 个可与发布结果比对的基因中，跨样本表达相关的中位数为 0.99946。这项核对检验样本间表达轮廓，无法单凭中心化结果验证样本内排序。','GSE31684':'93 份膀胱切除标本中 78 份为 MIBC、15 份为 NMIBC。按官方临床补表核实 3 人接受术前化疗，排除后剩 90 人。分型敏感性使用原研究 cluster。','GSE13507':'165 份原发样本均与临床工作表匹配，按研究原侵袭性标签为 103 NMIBC、62 MIBC。BT142 的 TaN3M1 与研究的 MIBC 标记有分歧，另列纯 T 分类的 104/61 敏感性。排除 23 份复发、58 份癌旁、9 份可识别正常和 1 份标签冲突样本 GSM340542。系列设计称复发来自 14 人，GSM 标题含 15 个不同 BT 前缀，无法据此重建可靠配对。全部复发均排除。BT002（GSM340606）保留在主分析中；仅在影响诊断中留出它时，mTORC1–NECTIN2 的 r 从 0.226 变为 0.020，MP10–NECTIN2 从 0.173 变为 0.009。'}
intro={'TCGA':'四个上游指标与 NECTIN2 的相关均接近零。','GSE32894':'四个上游指标与 NECTIN2 同为弱正相关。','GSE31684':'mTORC1 与 NECTIN2 的点估计为负，SREBF1 与 NECTIN2 弱正相关且区间跨零。','GSE13507':'NECTIN2 与 SREBF1 的正向相关较清楚，与 SREBF2 的区间跨零。'}
secs=[]
for a in ['TCGA']+ext:
 td=T/'bulk'/a;u=read(td/'core.tsv');cov=read(td/'score_coverage.tsv');aux=read(td/'Treg.tsv');fm=read(td/'family.tsv');ss=read(td/'sensitivities.tsv')
 sec='## '+a+'\n\n'+desc[a]+' '+intro[a]+'\n\n'+tab(['配对','n','Pearson r [95% CI]','BH q'],[(v['x']+' – '+v['y'],v['n'],interval(v),qq(v['q'])) for v in u])
 used=[v for v in cov if v['mask']=='Full' and v['variable'] in ['MP10','mTORC1','SREBF1','SREBF2','Treg','Tlineage']]
 sec+='\n\n背景 '+used[0]['background']+' 个基因。实际覆盖为'+'、'.join(v['variable']+' '+v['used']+'/'+v['expected'] for v in used)+'。'
 sec+=' 家族比较中，SREBF1 与各成员的 r 为'+'、'.join(v['y']+' '+ff(v['r']) for v in fm if v['x']=='SREBF1')+'。Treg 辅助层的相关为'+'、'.join(v['x']+'–'+v['y']+' '+interval(v) for v in aux)+'。'
 par=[v for v in ss if v['family_id']=='core' and v['x']=='SREBF1' and v['y']=='NECTIN2' and v['model'] in ['MIBC','NMIBC','subtype_or_study_cluster','T_stage','DNA_purity','without_preRC_chemo','T_lineage','CD3D_composition_proxy','M20','remove_shared_members']]
 sec+='\n\n'+tab(['SREBF1–NECTIN2 敏感性','有效 n','相关 [95% CI]'],[(v['model'],v['n'],interval(v,'partial_r' if v['model'] in ['subtype_or_study_cluster','T_stage','DNA_purity','T_lineage','CD3D_composition_proxy'] else 'r')) for v in par])
 rob={v['model']:v for v in ss if v['family_id']=='core' and v['x']=='MP10' and v['y']=='SREBF2'}
 lead=max(u,key=lambda v:float(v['max_loo_delta']))
 sec+='\n\nMP10–SREBF2 的原始 r 为 '+ff(rob['raw']['r'])+'，M20 后为 '+ff(rob['M20']['r'])+'，两侧去共享成员后为 '+ff(rob['remove_shared_members']['r'])+'。九项中删一变化最大的是 '+lead['x']+'–'+lead['y']+'，'+lead['most_influential']+' 对应 |Δr|='+ff(lead['max_loo_delta'])+'；该样本保留。'
 if a in ['GSE13507','GSE32894']:sec+='\n\nT-lineage 仅 4/5 个可用，低于原最少 5 个成员的门槛。另列单基因 CD3D 成分代理。'
 if a!='TCGA':(D/'reports'/(a+'.md')).write_text(sec+'\n\n核心 q 在三队列 27 项中校正；家族、辅助、敏感性各有独立校正范围。完整分层、影响诊断及每个标本的测量见本队列表。\n')
 secs.append(sec)
m=[v for v in meta if v['family_id']=='core']
cross='## 外部汇总\n\n'+tab(['配对','研究数','REML r [95% CI]','I² (%)','τ²'],[(v['x']+' – '+v['y'],v['k'],interval(v),ff(v['I2_pct'],1),ff(v['tau2'],3)) for v in m])+'\n\n正向重现最清楚的是上游程序间共变。SREBF1–NECTIN2 是下游较一致但较弱的关联。mTORC1–NECTIN2 在 GSE31684 方向相反，合并均值和宽区间会掩盖差异。三个队列对异质性的估计不精确，TCGA 单列参考。家族成员间强弱因队列而异，转录关联不能指定配体竞争或受体占据。'
pert=[];evidence=[]
def ev(link,source,status,detail):evidence.append(dict(link=link,source=source,status=status,result=detail))
ev('上游代谢程序共变','TCGA + 3 BLCA bulk','支持组织水平共变','见各队列 n、Pearson r 与重叠去除敏感性')
ev('SREBF1 程序与 NECTIN2','外部3队列','支持弱正相关','REML '+interval(next(v for v in m if v['x']=='SREBF1' and v['y']=='NECTIN2')))
ev('mTORC1 程序与 NECTIN2','外部3队列','方向不一致','切除队列为负，其他两队列为正')
for acc,lines in [('GSE236742',['HCV29']),('GSE143803',['KYSE150','TE5']),('GSE267018',['HeLa'])]:
 sec='## '+acc+' 真实扰动\n\n'
 for line in lines:
  td=T/'perturbation'/acc/line;p=read(td/'programs.tsv');fg=[v for v in read(td/'family.tsv') if v['feature']=='NECTIN2'];glob=read(td/'global.tsv')
  sec+='**'+line+'**　'+{'HCV29':'人膀胱来源 TSC1-null，24 h。两个药物各与匹配 DMSO 比较，每组 3 份 RNA。','KYSE150':'食管细胞系，siSREBF1 对 siNC，各 2 份 RNA。','TE5':'食管细胞系，siSREBF1 对 siNC，各 2 份 RNA。两份重复的设计信息不足以把它们称作两种独立 siRNA。','HeLa':'人宫颈细胞系，每个基因型和培养条件各 3 份 RNA。公开 HeLa 补表为 FPKM。交互表示耗竭相对 FBS 的 KO 效应变化。'}[line]+'\n\n'
  sec+=tab(['干预或交互','SREBF1 程序 Δ / q','SREBF2 程序 Δ / q','NECTIN2 效应 [95% CI]','家族 q','全基因 q<.05'],[(v['contrast'],*[next(ff(w['effect'],0)+' / '+qq(w['q']) for w in p if w['feature']==n and w['contrast']==v['contrast']) for n in ['SREBF1','SREBF2']],interval(v,'effect','low','high'),qq(v['q']),next(g['n_BH05']+'/'+g['n_tested'] for g in glob if g['contrast']==v['contrast'])) for v in fg])+'\n\n'
 sec+= {'GSE236742':'两药均降低 mTORC1、两个 SREBF 程序和 G2M，凋亡转录程序升高；NECTIN2 在 rapamycin 下升高、RMC-6272 下降。前半段有干预支持，后一段受药物或细胞状态影响。','GSE143803':'两条细胞系的 SREBF1 程序均下降。NECTIN2 在 TE5 下降，在 KYSE150 点估计上升且区间跨零。GEO 未在独立小型补表列出 SREBF1 peaks，保留直接占据为未评估。','GSE267018':'SREBP2 KO 在两种培养条件下均降低甾醇稳态程序与 NECTIN2，两个 KO2 条件效应的交互区间跨零。SREBP1 KO 未显示同向下降。表达效应以 Δlog2(FPKM+1) 表示。该干预支持跨癌种的 SREBF2 候选环节，BLCA 中的直接转录调控仍待检验。'}[acc]
 pert.append(sec)
ev('mTOR 干预到 SREBF 程序','HCV29','支持','两种药物的 mTORC1、SREBF1、SREBF2 程序下降')
ev('mTOR 干预到 NECTIN2','HCV29','方向不一致','rapamycin 上升，RMC-6272 下降')
ev('SREBF1 干预到 NECTIN2','TE5 / KYSE150 / HeLa','背景依赖','TE5 下降；KYSE150、HeLa 未见一致下降')
ev('SREBF2 干预到 NECTIN2','HeLa','跨癌种支持','两种培养条件下下降，交互不精确')
cc=read(T/'competition/models.tsv')
bc='## 配对 Treg 与有限竞争模型\n\n原 Treg OXPHOS 固定为结局。每次在 dataset 模型中加入一个去 OXPHOS 成员的受者读数。系数为每个分析分母内 Epi 暴露增加 1 SD 时结局的变化。\n\n'+tab(['暴露 / 范围','控制','记录患者','β [HC3 95% CI]'],[(v['exposure']+' / '+v['population'],v['control'],v['patients'],interval(v,'beta')) for v in cc if not(v['exposure']=='Q50_z' and v['population']=='replicated_datasets')])+'\n\nQ50 的负向点估计在加入活化读数后减弱，控制自主合成时近似原值，各区间均宽。18 人 High 匹配集中有 singleton dataset，杠杆为 1；表里只有点估计。预设有重复标本的数据集范围为 17 个记录患者，区间同样跨零。既有 Treg 的普通 TCR 投影 β=0.199 [−0.019, 0.417]；SCAP、FASN、CD36 与 SREBF2 依赖未显示稳定增强。既有 G2M 上升、OXPHOS 下降及 noncycling/cycling 拆分均来自同批标本。普通活化是应保留的竞争解释，自主代谢不能被排除。\n\n既有通讯中，Treg→Epi 和 Epi→Treg 的 High/Low 权重比经各方向背景调整后分别为 1.482 和 0.991。NECTIN2_TIGIT 和 TIGIT_NECTIN2 是镜像记录。这些结果无法给出绝对作用强度、膜蛋白结合或单位 Treg 的抑制能力。'
ev('NECTIN2 到 Treg 功能','bulk / 配对 / CellChat','功能不可评估','组织成分、受者 OXPHOS 转录及通讯候选均未直接测量单位 Treg 功能')
ev('共同活化与自主合成','32 记录患者与18人 High 交集','效应不精确','TCR 调整减弱 Q50 点估计，自主合成调整影响小')
put(T/'evidence_by_link.tsv',evidence)
final='## 机制判断\n\n优先保留膀胱来源模型中的 mTOR 干预到 SREBF 程序，以及跨癌种模型中的 SREBF2 到 NECTIN2 两段证据。BLCA bulk 的 SREBF1–NECTIN2 弱正相关可用于选择组织验证指标。两种 mTOR 抑制剂对 NECTIN2 方向相反，且配对 Treg 效应不精确，现有数据不足以连成 mTORC1→SREBF→NECTIN2→Treg 功能的因果链。要区分共同活化、受体表达与功能改变，仍缺同组织中 Epi 膜 NECTIN2、Treg 表面受体和独立功能读出的对应测量。'
(D/'RESULTS.md').write_text('# BLCA 外部验证与真实上皮扰动\n\n'+ '\n\n'.join(secs+[cross]+pert+[bc,'## 虚拟接口\n\n<!--virtual-->\n等待读取接口状态。\n<!--/virtual-->',final])+'\n')
refresh()

