from common import *
g=read('global_incoming_outgoing.tsv');c=read('epi_allTME_total_contrasts.tsv');obj=read('object_index.tsv');tr=read('treg_total_summary.tsv');pc=read('treg_partner_categories.tsv');lr=read('treg_LR_contribution_decomposition.tsv.gz');cross=read('cross_dataset_direction_summary.tsv');supp=read('supplementary_NECTIN2_TIGIT_contribution_context.tsv');expr=read('supplementary_NECTIN2_TIGIT_dataset_expression.tsv');u=read('epi_state_totals.tsv')
def md(df):
 def fmt(v):
  if pd.isna(v):return 'NA'
  if isinstance(v,(float,np.floating)):return f'{v:.5g}'
  return str(v).replace('|',' / ').replace('\n',' ')
 return '| '+' | '.join(map(str,df.columns))+' |\n|'+'|'.join(['---']*len(df.columns))+'|\n'+'\n'.join('| '+' | '.join(fmt(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))+'\n'
major=g[g.scope.eq('ALL')&g.network_level.eq('major')&~g.node.isin(EPI)][['node','incoming_all_W','rank_incoming_all_W','outgoing_all_W','rank_outgoing_all_W','self_W']]
fine=g[g.scope.eq('ALL')&g.network_level.eq('fine')&~g.node.isin(EPI)];fin=[fine.nlargest(8,'incoming_all_W')[['node','incoming_all_W','rank_incoming_all_W','n_matched_recorded_IDs']],fine.nlargest(10,'outgoing_all_W')[['node','outgoing_all_W','rank_outgoing_all_W','n_matched_recorded_IDs']]]
ct=c[c.network_level.eq('fine')&c.node.eq('Treg')].copy();ct['low_share_pct']=ct.low_share*100;ct['high_share_pct']=ct.high_share*100
rank=tr[['scope','incoming_all_W','outgoing_all_W','rank_incoming_all_W','rank_outgoing_all_W','n_TME_ranked','self_W','rank_incoming_other_W','rank_outgoing_other_W','n_matched_specimens','n_matched_recorded_IDs']]
sets=c[c.network_level.eq('fine')].groupby(['scope','direction']).agg(n_computed=('computed','sum'),n_high_greater=('delta_W',lambda x:x.gt(0).sum()),n_low_greater=('delta_W',lambda x:x.lt(0).sum()),n_zero=('delta_W',lambda x:x.eq(0).sum())).reset_index();save(sets,'whole_TME_state_direction_counts.tsv')
selected=['Treg','CD8 Exhausted','SPP1 TAM','C1QC TAM','FOLR2 TAM','Capillary EC','Tip EC','Matrix_CAF','iCAF','Pericyte','dataset_bias_fibro'];cr=cross[cross.network_level.eq('fine')&cross.node.isin(selected)][['node','direction','n_computable_datasets','n_high_greater','n_low_greater','n_zero_difference','n_datasets_at_least_two_matched_IDs','n_high_share','n_low_share']]
cont=lr[lr.scope.eq('ALL')];cont_top=pd.concat([cont[cont.direction.eq(dr)].nlargest(7,'contribution') for dr in ['epi_to_TME','TME_to_epi']])[['direction','LR_id','supported_low','supported_high','contribution','raw_difference']];h=lr[lr.scope.eq('HRA003620')&lr.direction.eq('epi_to_TME')].nsmallest(5,'contribution')[['LR_id','contribution']]
selfpart=pc[pc.scope.eq('ALL')].copy();selfpart['percent_of_direction_total']=selfpart.W/selfpart.groupby('direction').W.transform('sum')*100
ss=read('supplementary_NECTIN2_TIGIT.tsv');ss=ss[ss.LR_id.eq('NECTIN2_TIGIT')][['scope','epi_state','raw_probability','internal_P','display_strength','status']]
# Important standard-record dependence is quantified without altering main W/C definitions.
rel=supp[supp.LR_id.eq('TIGIT_NECTIN2')&supp.direction.eq('TME_to_epi')].copy();rel['fraction_of_net_W_difference']=rel.contribution/rel.total_delta_W.replace(0,np.nan);save(rel,'reverse_DB_record_contribution_context.tsv')
# Summarize actual complete result for all groups, not a new weighted candidate score.
pooled=g[g.scope.eq('ALL')][['network_level','node','incoming_all_W','outgoing_all_W','rank_incoming_all_W','rank_outgoing_all_W']].rename(columns={k:'pooled_'+k for k in ['incoming_all_W','outgoing_all_W','rank_incoming_all_W','rank_outgoing_all_W']});allgroup=cross.merge(pooled,on=['network_level','node'],how='left');save(allgroup,'all_group_reading_summary.tsv')
text=f'''# MP10 round3：全TME总通讯、上皮状态比较与Treg定位

## 1. 全TME的大群与亚群总通讯格局

本轮完成模块0–4，主结果覆盖全部TME，不以B/Plasma、T/NK或旧AX列表限定。复用v2全部12个可用共享状态对象：五队列与ALL、major/fine各一套。实际纳入40个肿瘤标本、39个derived/provisional记录患者ID；Inhouse为12标本、11个记录ID。两类epi是每个标本内部的相对MP10状态，共享同一份TME。未重新评分、推断网络或跑丰度模型，v1/v2保持只读。

**全局总输出最突出的谱系是Fibroblast；Treg并不是全局总收发最强的节点。** 全局incoming、outgoing、与epi的强度、MP10状态差和通讯份额是不同问题，下面分开回答。W为全部内部P<0.05的已计算LR概率之和，C为支持LR记录数量，R为不按P筛选的原始和。这些总量没有患者差异显著性P。

### 1.1 读取对象与实际覆盖

{md(obj[['scope','network_level','n_nodes','n_computed_LR','n_specimens','n_recorded_patient_IDs']])}

major固定9个显示节点（两个epi＋七TME），fine固定47个显示节点（两个epi＋45个TME标签）。队列中不可计算节点显示NE，排名不把它们当0纳入。对象只覆盖其已计算LR维度，不代表所有可能分子联系。完整表保存TME↔TME、epi↔TME、epi↔epi及自环。

### 1.2 大群主格局

ALL的七个TME大群结果如下；排名仅在七个可计算TME内进行：

{md(major)}

Fibroblast的ALL outgoing W为10.8828，Endothelial为6.1562，Myeloid为1.5691。**这三个大群在五个dataset中均分别排outgoing第1、2、3**，比单一pooled图更有重复覆盖。incoming方面，Endothelial在四队列第1，HRA由T/NK居第1。B/Plasma/Mast等仍保留完整结果，没有因为主线改变而删除。

自环不能忽视：ALL Fibroblast的incoming中约51.7%来自大群自环；这是同一大群内不同细胞/亚型也被合并后的联系，不能称为单细胞自我作用。所有IO另给去自环版本。major和fine各自独立推断，不能把细群和当作大群，也不能把两层绝对W直接比较。

### 1.3 细群格局

ALL incoming前列的代表项（完整45群见表/图）：

{md(fin[0])}

ALL outgoing前列的代表项：

{md(fin[1])}

CD8 Exhausted的ALL incoming W为48.9758，排1/45；SPP1 TAM和多种内皮亚群也处于前列。outgoing则主要由成纤维/基质类标签占据。dataset_bias_fibro虽数值最高，但其实际细群网络只在Inhouse可计算、配对覆盖3个ID，且有既有偏倚标签，不能把这个排名升级为跨队列生物学候选。Cycling stromal、vCAF等也存在稀少覆盖；Matrix_CAF、iCAF、Pericyte更适合结合各自队列范围继续阅读。

全量circle保留所有正边与实际自环；fine图较密，精读时使用配套完整source×target热图。主图入口是聚合网络，不是单LR图。

## 2. 两类上皮对全部TME的总输出、总输入与群体差异

对同一对象内的同一组可计算TME，U_out只汇总epi→TME、U_in只汇总TME→epi，均排除epi↔epi及epi自环。以下fine总量每一行只在对应对象内比较：

{md(u[u.network_level.eq('fine')][['scope','direction','n_common_TME','U_low','U_high','U_difference']])}

ALL中，epi输出总W由11.3465增至14.3185（+26.2%），从TME接收的总W由17.7825增至26.0795（+46.7%）。五队列里CNP、GSE、Inhouse、PRJNA两个总方向均偏high；HRA两个方向的总体W均略偏low。major也呈同样的队列级方向格局，但绝对值不能跨层比较。

全部TME节点的对象内方向计数：

{md(sets)}

**ALL细群的两个方向均为45/45节点high增强。** 因而“某群high联系更强”在pooled结果中本身缺乏相对特异性，需要再看份额、状态差排名和队列覆盖。不同对象的原始W没有被直接平均成患者效应，ALL也不计作第六个独立队列。

ALL epi→TME的最大W差出现在Tip EC（+0.20844）、Capillary EC（+0.20380）、Cycling EC（+0.19325）、Arterial EC（+0.17044）等；Treg差为+0.03993，排37/45。反方向最大差包括Matrix_CAF（+0.85653）、dataset_bias_fibro（+0.80091）、Cycling stromal（+0.78106）、Pericyte（+0.69100）；Treg为+0.01887，排27/45。Tip EC只有Inhouse一个细群网络可计算，Cycling EC主要单队列覆盖，不能仅凭ALL差最大就称为最稳健。

## 3. Treg的全网络地位、与两类上皮的聚合联系和份额

### 3.1 Treg总体收发是否强？

同一shared对象只有一个Treg节点，下面全局值只有一套；没有high组Treg/low组Treg两套总量。排名是当前对象可计算TME中的名次：

{md(rank)}

ALL Treg incoming W=18.1792（36/45）、outgoing W=8.0240（34/45）；去自环后incoming=18.1128（34/45）、outgoing=7.9576（34/45）。自环W=0.06634，只占incoming的0.365%、outgoing的0.827%，所以较低排名主要不是自环造成的。其global C分别为823和430，是跨伙伴细胞群对累计的支持LR记录数，不是823/430个独立机制。

Inhouse Treg incoming=25.0988（27/39）、outgoing=11.8206（31/39）；去自环incoming排23/39、outgoing仍31/39。各队列都不能据此称Treg为“全局通讯最强”或“免疫抑制最强”。

### 3.2 Epi-high与Treg的总联系是否更强？

{md(ct[['scope','direction','low_W','high_W','delta_W','low_C','high_C','delta_R','W_R_direction_agree','rank_delta_W']])}

epi→Treg在CNP、GSE、Inhouse、PRJNA偏high，HRA偏low；反方向Treg→epi在五队列均偏high。Treg这12个对象/方向比较的W和辅助R符号一致，但内部支持阈值会改变增量大小，不能把W的增加自动当作统计显著。

ALL epi→Treg：0.18372→0.22366，支持LR数14→18；Treg→epi：0.01956→0.03843，支持LR数5→7。Inhouse分别为0.28556→0.33465（10→14条）和0.01487→0.04175（3→5条）。HRA epi→Treg的W下降，但C从20增至21，直接说明“数量更多”和“总强度更高”不能混用。

### 3.3 增强是否相对偏向Treg？

{md(ct[['scope','direction','low_share_pct','high_share_pct','rank_low_W','rank_high_W','rank_delta_share']])}

ALL epi分配给Treg的输出份额由1.6192%降至1.5620%，强度排名由31/45变为33/45；所以虽然绝对联系增强，**并没有相对更偏向Treg**。该输出份额只有CNP、PRJNA微升，GSE、HRA、Inhouse下降。Inhouse对应2.5760%→2.3997%。

反方向，Treg占epi全部TME输入来源的份额由0.1100%增至0.1474%，份额差排名8/45；五队列有四个上升、GSE下降。它在这个组成指标上比全局排名更突出，但绝对份额仍很小，而且与数据库反向分子记录相关，不能直接称为Treg特异反馈。

### 3.4 Treg总量主要来自谁、流向谁？

ALL四类来源/去向分解：

{md(selfpart[['direction','category','W','C','percent_of_direction_total']])}

Treg incoming的97.39%、outgoing的98.45%来自/流向其他TME；两类epi合计只占其incoming的2.24%、outgoing的0.72%。因此全局Treg总量与它的epi状态边不是同一个量。所有单独伙伴均在`treg_partner_breakdown.tsv`和`Treg_all_partners.pdf`中保留，没有仅展示强伙伴。

配对覆盖为34个标本、33个记录ID（CNP7、GSE4、HRA6、Inhouse11、PRJNA5个ID）；网络整体仍是40标本/39 ID。覆盖差异没有被当作删除全谱系网络的新门槛，也没有用这个数冒充患者级总W的重复。

## 4. 其他值得讨论的群体与反向/不一致队列

下表从全量结果中摘出便于对照的群体，不是另建加权候选分数；全部群仍在完整表和固定顺序方向矩阵中：

{md(cr)}

可以分别提出以下具体问题：

- **Fibroblast/CAF是全局主要输出谱系。** 大群outgoing五队列均第1；Matrix_CAF反向TME→epi的ALL差最大，但细群只在3队列计算、其中2队列有至少2个匹配ID。iCAF更广泛，HRA存在反向；Pericyte在4个可计算队列中两个方向均high增强，配对覆盖较完整的有3队列。它们值得研究的理由与Treg不同，不能统称免疫抑制。
- **CD8 Exhausted是强接收节点。** ALL global incoming第1、epi→该群在五队列均high增强，五队列均有至少两个匹配ID；但输出份额只在两队列增加、三队列下降，ALL份额也下降，仍需区分绝对增量和相对偏好。
- **C1QC TAM/FOLR2 TAM是覆盖较完整的双向对照。** 五队列都可计算且有至少两个匹配ID；其→epi总量五队列均high增强，而epi→它们在HRA反向。这些是可以与Treg共同比较的细胞群，而非沿用v1的单LR shortlist。
- **Capillary EC有较大epi状态差。** 四个可计算队列中均偏high，但只有CNP、PRJNA达到至少两个同标本匹配ID。Tip EC、Cycling EC及部分CAF/dataset_bias标签更偏单队列或低覆盖，保留数值而不宣称普遍性。
- **B、Plasma、Mast、其他T/NK及未解析标签完整保留。** 例如B/Plasma并未成为本轮入口限制；其队列内结果与所有群使用相同W/C/R、排名和覆盖口径。

HRA是结构性不一致队列：epi总输出/输入均略低，但部分节点（例如CD8 Exhausted的epi接收方向、C1QC/FOLR2的→epi方向）仍偏high。不能把HRA整队列简单删除，也不应把缺失内皮/CAF节点记成“反向”。

## 5. 单条LR如何解释Treg的总量变化

每条贡献是支持阈值后的高低差，不是未过滤概率差；所有已计算LR都进入分解。ALL的正贡献代表项：

{md(cont_top)}

ALL epi→Treg有18条正贡献、无负贡献，总差+0.03993；NECTIN2_TIGIT贡献+0.01099，约27.51%，随后有APP_CD74 +0.00702、LAMA5_CD44 +0.00368、MIF_CD74_CD44 +0.00300等。NECTIN2在这个pooled对象是最大单条贡献，但在CNP、GSE、Inhouse的正贡献排名分别是5、6、3，不是每队列都第一。

ALL Treg→epi的TIGIT_NECTIN2贡献也是+0.01099，约占该方向净差58.21%。数据库的正反记录共享分子对，当前原始概率也相同；主总量按标准DB维度保留，但这不能当作两条独立重复或已证实双向功能反馈。另一个ADGRE5_CD55贡献+0.00512，其中low内部P=0.9未支持、high支持，故W贡献大于其原始概率差；这正是必须区分W差与R差的例子。

HRA的负贡献并未丢弃：

{md(h)}

在HRA，NECTIN2_TIGIT仍增加约0.00149，但MDK_NCL、COL4A5_CD44、APP_CD74等负贡献更大，合计epi→Treg总W下降0.02183。**一条轴增强不能替代总通讯结论。** 完整LR分解和“其余正/负贡献”图已交付，贡献和与总差的误差处于浮点精度范围。

## 6. NECTIN2–TIGIT小型补充

此补充放在`supplementary/NECTIN2_TIGIT`，不替代全网络主图。真实单LR值如下：

{md(ss)}

所有五队列和ALL的low/high两端均有内部支持，圈图两条边都画，没有制造“high有、low无”。内部P=0是100次置换下的经验输出，不是患者差异P，也不说明真实概率精确为0。

同标本表达复核：

{md(expr)}

共34个合格标本/33个记录ID；先把同ID重复标本的表达差平均后，31个ID的epi NECTIN2高状态均值较高，2个较低（均在HRA）。五队列的留一ID平均差均保持正向。检测支持分别记录：34标本中low端32、high端33个满足该状态NECTIN2及共享Treg TIGIT各≥10%，未要求两端同时通过才保留。保留真实零和低表达，不以细胞作为患者重复。

CD8 Exhausted及其他实际TIGIT阳性T/NK接收群的同口径表达/网络简洁对照在`supplementary_TIGIT_positive_receiver_comparison.tsv`。这说明分子轴可供下一轮具体核查；**没有验证Treg全网络总强度，也没有证明Treg特异接触或抑制。**

## 7. 当前结论与下一步具体问题

本轮支持的描述是：全网络主要输出来自成纤维/基质谱系；部分髓系、内皮和CD8 Exhausted具有较强接收或上皮联系；Treg整体收发排名偏后，但epi状态边在多数队列增强。这个增强通常伴随广泛TME增强，并不等于epi更偏向Treg。Treg输出份额和输入来源份额可出现不同方向，HRA保留反向。内部支持总和没有自动获得患者级差异显著性。

下一步应把问题拆开：首先在明确标本/患者单位下，检验CAF/基质输出、内皮/髓系与epi的广泛联系；若聚焦Treg，检验其相对于CD8 Exhausted、C1QC/FOLR2等完整对照的邻近性与边界位置，以及NECTIN2/TIGIT是否在对应细胞和空间关系中匹配。HRA的总量下降但单轴增加是优先解释的不一致，不宜通过改阈值或删队列消除。

本轮没有从pooled对象反推出患者总W，没有把LR表达乘积命名为患者CellChat总通讯，没有新空间/注释/CNV/预后分析。实际主结论已有五队列Treg覆盖说明，不需要为绘图或单一预期重推网络；稀少其他节点按覆盖限制描述。

阅读入口为[逐队列图索引](FIGURE_INDEX.md)。主circle是全LR聚合图A/B/C，W与C分别输出；完整热图、IO含/不含自环、角色散点、全群高低比较与份额、Treg面板均按major/fine及五队列/ALL输出。透明PDF和300dpi PNG保存在完整目录，`review.zip`仅含核心PDF、主要表、README、方法和本报告；不含旧对象、表达矩阵或PNG。
'''
(O/'round3_report.md').write_text(text)
print('REPORT WRITTEN',len(text),'characters')
