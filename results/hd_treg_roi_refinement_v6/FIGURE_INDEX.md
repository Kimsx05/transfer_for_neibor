# 图像和表格入口

正式主线在 proseg_observed。目录顶层 figures 中的空间比较和切片顶层 effects.tsv 是保守归属对照。

## 核修订与原图校准

- [Cervilla · 拆分、合并、新增与旧核保留](figures/Cervilla_nuclear_changes.jpg)
- [Cervilla · 12 个校准窗口](figures/Cervilla_calibration.jpg)
- [D1_Johnson · 拆分、合并、新增与旧核保留](figures/D1_Johnson_nuclear_changes.jpg)
- [D1_Johnson · 12 个校准窗口](figures/D1_Johnson_calibration.jpg)

## 已实际打开的 MP10 空间窗口

图内标注数量是矩形视窗内数量，正式统计取同组件累计圆，分母见 neighborhoods.parquet。

- [Cervilla · Cervilla_ROI1_MP10](proseg_observed/figures/Cervilla_ROI1_MP10.png)　旧106665近处有两例T对照，可评分上皮主要在200微米圈右缘。

- [Cervilla · Cervilla_ROI2_MP10](proseg_observed/figures/Cervilla_ROI2_MP10.png)　小孤立覆盖岛未见可评分上皮，无比较机会。

- [Cervilla · Cervilla_ROI3_MP10](proseg_observed/figures/Cervilla_ROI3_MP10.png)　中心有数个灰色对象，无上皮分数或T中心。

- [Cervilla · Cervilla_ROI4_MP10](proseg_observed/figures/Cervilla_ROI4_MP10.png)　两个圈旁散在极少对象，无可评分上皮。

- [Cervilla · Cervilla_ROI5_MP10](proseg_observed/figures/Cervilla_ROI5_MP10.png)　上皮从西北和南侧进入视野，多数分值低；未见主Treg。

- [Cervilla · Cervilla_ROI6_MP10](proseg_observed/figures/Cervilla_ROI6_MP10.png)　组织灰色对象沿左下方向分布，展示窗口内无可评分上皮。

- [Cervilla · Cervilla_ROI7_MP10](proseg_observed/figures/Cervilla_ROI7_MP10.png)　对象沿边界成笔直的两行；原图黑边提示捕获边缘。

- [D1_Johnson · D1_Johnson_ROI1_MP10](proseg_observed/figures/D1_Johnson_ROI1_MP10.png)　旧49654周围有多例T对照，可评分上皮仅出现在200微米左下边缘。

- [D1_Johnson · D1_Johnson_ROI2_MP10](proseg_observed/figures/D1_Johnson_ROI2_MP10.png)　捕获边缘密度向左增加，圈内无可评分上皮。

- [D1_Johnson · D1_Johnson_ROI1_cell3391_MP10](proseg_observed/figures/D1_Johnson_ROI1_cell3391_MP10.png)　A路径中心附近有对照，200微米内无可评分上皮。

- [D1_Johnson · D1_Johnson_ROI1_cell62224_MP10](proseg_observed/figures/D1_Johnson_ROI1_cell62224_MP10.png)　严格个体附近同时有A与B及多个T对照，圈内缺可评分上皮。

- [D1_Johnson · D1_Johnson_ROI1_cell32362_MP10](proseg_observed/figures/D1_Johnson_ROI1_cell32362_MP10.png)　A中心200微米内左下仅见一例低分上皮，50微米空缺。

- [D1_Johnson · D1_Johnson_ROI1_cell138645_MP10](proseg_observed/figures/D1_Johnson_ROI1_cell138645_MP10.png)　B中心的上皮主要在西南外圈，H05位于较远一侧，50微米内无可评分上皮。

## 已实际打开的逐对象及未定位来源卡

绿色方格代表最终唯一归属的实测捕获格，橙线为核，黄线为 Proseg 推定胞体。
- [Cervilla cell_101487](proseg_observed/Cervilla/cards/cell_101487.jpg)　旧106665的核有两种直接T标记及CTLA4，最终未见归属FOXP3；细分模型对Treg的支持不足。
- [Cervilla cell_102236](proseg_observed/Cervilla/cards/cell_102236.jpg)　直接T位于核周，普通CD4竞争明确，可作为对照。
- [Cervilla cell_102345](proseg_observed/Cervilla/cards/cell_102345.jpg)　直接T在右侧核缘，可见其他位置的FOXP3，未归给此核。
- [Cervilla cell_33002](proseg_observed/Cervilla/cards/cell_33002.jpg)　旧25104的现核周围有间质表达，FOXP3在最终推定胞体之外。
- [Cervilla cell_39197](proseg_observed/Cervilla/cards/cell_39197.jpg)　旧59121可见FOXP3捕获格，但直接T为零，邻近间质表达丰富。
- [Cervilla cell_56997](proseg_observed/Cervilla/cards/cell_56997.jpg)　相邻两核间有T及IKZF2表达，当前归属缺FOXP3，模型更支持非T。
- [Cervilla cell_72058](proseg_observed/Cervilla/cards/cell_72058.jpg)　胶原丰富组织中可见单格FOXP3及直接T，整体主群仍偏间质。
- [Cervilla cell_78199](proseg_observed/Cervilla/cards/cell_78199.jpg)　旧92281的旁侧有FOXP3，核支持对象自身无直接T，髓系竞争强。
- [Cervilla cell_11311](proseg_observed/Cervilla/cards/cell_11311.jpg)　视野周围为空白，深染小物体只有6 UMI，无直接T或核心捕获格。
- [Cervilla cell_11322](proseg_observed/Cervilla/cards/cell_11322.jpg)　密集多边形核间存在大量上皮表达，IKZF2单独检出且未见直接T。
- [Cervilla cell_18722](proseg_observed/Cervilla/cards/cell_18722.jpg)　分散的深染点状结构中只有3 UMI，缺支持身份的转录本。
- [Cervilla cell_26757](proseg_observed/Cervilla/cards/cell_26757.jpg)　背景近乎空白，轮廓内5 UMI。不能从单个深染结构判定T身份。
- [Cervilla cell_109894](proseg_observed/Cervilla/cards/cell_109894.jpg)　上皮组织中的小核，核边界与上下邻核可分，周边上皮信号较多，自身未见T标记。
- [Cervilla cell_131167](proseg_observed/Cervilla/cards/cell_131167.jpg)　靠组织空隙的椭圆核只分到11 UMI，零星周边上皮点不能解释为T身份。
- [Cervilla cell_135749](proseg_observed/Cervilla/cards/cell_135749.jpg)　原图黑色边界缺乏可辨核，1 UMI。标记为图像或捕获边缘，不能用作身份判定。
- [D1_Johnson cell_3391](proseg_observed/D1_Johnson/cards/cell_3391.jpg)　原图中为较宽的独立核，两个 FOXP3 信号与直接 T 在核内或紧贴轮廓。黄色 Proseg 范围没有跨入左邻核。
- [D1_Johnson cell_9377](proseg_observed/D1_Johnson/cards/cell_9377.jpg)　小横椭圆核贴近空隙，核心信号在核右半，T 在左缘；右侧邻核轮廓仍能分开。
- [D1_Johnson cell_9913](proseg_observed/D1_Johnson/cards/cell_9913.jpg)　FOXP3 位于较小的梨形核下部，同核及核周有多个 T 捕获格；竞争点少量。
- [D1_Johnson cell_12265](proseg_observed/D1_Johnson/cards/cell_12265.jpg)　直接 T 靠下侧另一狭长核，核心在左下缘。Proseg 核周格不能解决此接缝来源。
- [D1_Johnson cell_12516](proseg_observed/D1_Johnson/cards/cell_12516.jpg)　孤立窄核，直接 T 与 FOXP3 在周边，邻核距离可辨；pT 靠近门槛。
- [D1_Johnson cell_26817](proseg_observed/D1_Johnson/cards/cell_26817.jpg)　原图小圆核内分见 T 与 FOXP3，70 UMI；细亚型仍以常规 CD4 竞争。
- [D1_Johnson cell_32063](proseg_observed/D1_Johnson/cards/cell_32063.jpg)　FOXP3 格与下方贴靠的核轮廓交界，Proseg 边缘量化后仍不能明确其属于哪一核。
- [D1_Johnson cell_32362](proseg_observed/D1_Johnson/cards/cell_32362.jpg)　较孤立的细长核，T 与 FOXP3 在上核缘两侧，背景胶原可辨；旧 CAF 解释保留。
- [D1_Johnson cell_34418](proseg_observed/D1_Johnson/cards/cell_34418.jpg)　较圆核左缘 FOXP3、上缘直接 T 与明确 T 主群；细分类明显支持 CD8，按 A 路径保留冲突。
- [D1_Johnson cell_37080](proseg_observed/D1_Johnson/cards/cell_37080.jpg)　弯形小核上缘有直接 T 和 IKZF2，黄色区域局限，不见第二个可分核中心。
- [D1_Johnson cell_39383](proseg_observed/D1_Johnson/cards/cell_39383.jpg)　小核下缘呈 FOXP3/T 支持，周围相邻 T 捕获格单独可见；未扩大其归属。
- [D1_Johnson cell_39660](proseg_observed/D1_Johnson/cards/cell_39660.jpg)　核内见 IL2RA，但对应 T 信号靠左侧相贴的另一核，黄色延伸仅少量格宽。
- [D1_Johnson cell_41343](proseg_observed/D1_Johnson/cards/cell_41343.jpg)　同核内 FOXP3、上核缘直接 T；黄色轮廓没有吞并右上相邻核，胶原背景保留。
- [D1_Johnson cell_42872](proseg_observed/D1_Johnson/cards/cell_42872.jpg)　较长形核右侧 FOXP3 位于橙线内，T 位于两侧；CD4 及背景竞争仍在。
- [D1_Johnson cell_43869](proseg_observed/D1_Johnson/cards/cell_43869.jpg)　中心位于密集相触核间，T 点在左上接缝，CTLA4 在右侧贴邻核边，归属不稳定。
- [D1_Johnson cell_56118](proseg_observed/D1_Johnson/cards/cell_56118.jpg)　中心小核右上方三个核紧靠，FOXP3 捕获格落在右上交界，无法仅凭推定胞体消除竞争。
- [D1_Johnson cell_57539](proseg_observed/D1_Johnson/cards/cell_57539.jpg)　大核右上与另一大核相贴，下缘散在 T、IL2RA；旧 DC 背景与边缘支持均保留。
- [D1_Johnson cell_58642](proseg_observed/D1_Johnson/cards/cell_58642.jpg)　大轮廓内部多个深染焦点，核心与 T 集中在不同边缘；不凭 H&E 强拆，严格候选降级。
- [D1_Johnson cell_59265](proseg_observed/D1_Johnson/cards/cell_59265.jpg)　宽大核有明显核内深染结构，主要 T 点近左下接触面，IKZF2 在下界，胞体来源有疑问。
- [D1_Johnson cell_61533](proseg_observed/D1_Johnson/cards/cell_61533.jpg)　T 位于核上缘，IL2RA 在下缘，独立圆核上下左右轮廓可分。
- [D1_Johnson cell_61994](proseg_observed/D1_Johnson/cards/cell_61994.jpg)　直接 T 在核内，IL2RA 位于左上缘，旧 DC 与新 T 分类的分歧保留。
- [D1_Johnson cell_62224](proseg_observed/D1_Johnson/cards/cell_62224.jpg)　核内 FOXP3 与 IL2RA，周边直接 T；附近胶原/髓系点密集但有清晰核支持，保留旁源标志。
- [D1_Johnson cell_78835](proseg_observed/D1_Johnson/cards/cell_78835.jpg)　FOXP3 位于中心核下半，T 位于底缘，邻核另有独立 T 捕获格。
- [D1_Johnson cell_79078](proseg_observed/D1_Johnson/cards/cell_79078.jpg)　核下缘有两个直接 T，IL2RA 支持在左缘；核间分界清楚，亚型支持接近门槛。
- [D1_Johnson cell_79577](proseg_observed/D1_Johnson/cards/cell_79577.jpg)　细长大核内部同时有胶原/髓系/上皮信号，T 与 FOXP3 靠下缘，归属到同一 T 细胞不足。
- [D1_Johnson cell_138645](proseg_observed/D1_Johnson/cards/cell_138645.jpg)　顶部接触处新核成凹口，T 在左下缘，IKZF2 在右下核旁，保留 B 及余留边界不确定性。
- [D1_Johnson cell_161473](proseg_observed/D1_Johnson/cards/cell_161473.jpg)　旧49654的核确有图像支持。较紧的黄色范围未包括左下核周 FOXP3 捕获格，新计数为0，主群偏髓系。
- [D1_Johnson cell_2](proseg_observed/D1_Johnson/cards/cell_2.jpg)　组织边缘核，黄色边界部分覆盖明亮背景，只有22 UMI，缺直接 T。
- [D1_Johnson cell_151583](proseg_observed/D1_Johnson/cards/cell_151583.jpg)　边界修复后仅左半核有捕获格，三条 UMI 均无直接T；切片外侧为空白。
- [D1_Johnson cell_138409](proseg_observed/D1_Johnson/cards/cell_138409.jpg)　夹在两核下侧的小横核，T在顶缘、IKZF2在右边缘，旧Tfh仍是竞争解释；仅扩展层。
- [D1_Johnson cell_151693](proseg_observed/D1_Johnson/cards/cell_151693.jpg)　上皮丰富背景，右下有T和IKZF2捕获格，大量上皮信号分布于邻核。不升级扩展标签。
- [D1_Johnson cell_138643](proseg_observed/D1_Johnson/cards/cell_138643.jpg)　T与FOXP3在同一橙色小核中，主群支持未过0.5；维持扩展并记录上皮竞争。
- [D1_Johnson cell_138619](proseg_observed/D1_Johnson/cards/cell_138619.jpg)　接触纵向核底部的小核，直接T在远离接触线的一侧，30 UMI及CD8支持。保留对照的低深度标记。
- [D1_Johnson cell_138676](proseg_observed/D1_Johnson/cards/cell_138676.jpg)　不规则凹形核沿左上邻核延伸，直接T位于橙线内右侧，相邻另有多个T格；保留常规CD4及旁源标记。
- [D1_Johnson cell_57143](proseg_observed/D1_Johnson/cards/cell_57143.jpg)　密集区小核底部直接T可辨，左下邻核有独立T格，细分类为Tfh。
- [D1_Johnson cell_32462](proseg_observed/D1_Johnson/cards/cell_32462.jpg)　横椭圆核下缘直接T与右下方邻核可分，旧新均支持常规CD4。
- [D1_Johnson cell_57553](proseg_observed/D1_Johnson/cards/cell_57553.jpg)　斜长核内直接T，左下另一核的紫色核心点未归给目标。维持普通CD4对照。
- [D1_Johnson cell_56863](proseg_observed/D1_Johnson/cards/cell_56863.jpg)　细小核右上有T格、下方相邻核也有独立T；细亚型支持CD4，保留旁源说明。
- [D1_Johnson cell_159794](proseg_observed/D1_Johnson/cards/cell_159794.jpg)　纵向核内及顶端直接T，FOXP3和其他核心未见，旧CD8与新CD4竞争不影响非Treg T对照。
- [D1_Johnson cell_4905](proseg_observed/D1_Johnson/cards/cell_4905.jpg)　旧2233有核内FOXP3及IL2RA，但髓系竞争明显，周边有数个T格，pT很低。
- [D1_Johnson cell_12602](proseg_observed/D1_Johnson/cards/cell_12602.jpg)　旧14187靠血管样空隙，狭长核顶端T、IL2RA在下侧边缘，未达到主B路径。
- [D1_Johnson cell_29555](proseg_observed/D1_Johnson/cards/cell_29555.jpg)　旧40701小核周围多个T和FOXP3格，模型强烈偏髓系。保留定位线索。
- [D1_Johnson cell_12780](proseg_observed/D1_Johnson/cards/cell_12780.jpg)　旧48458的大椭圆核分见多个T格，IKZF2在核下缘，FOXP3点在推定胞体外。T内模型只支持扩展。
- [D1_Johnson cell_56470](proseg_observed/D1_Johnson/cards/cell_56470.jpg)　旧86643的核下方有T与FOXP3，黄色胞体未包含这些格；间质主群与周边胶原一致。
- [D1_Johnson cell_78997](proseg_observed/D1_Johnson/cards/cell_78997.jpg)　旧90099大核内有多个T和IL2RA，FOXP3在左下外侧；pReg|T仅0.31。
- [D1_Johnson cell_4335](proseg_observed/D1_Johnson/cards/cell_4335.jpg)　T见于核下缘，CTLA4在旁侧；原内皮与新上皮主群竞争，pT只过扩展门槛。
- [D1_Johnson cell_5172](proseg_observed/D1_Johnson/cards/cell_5172.jpg)　相贴两小核的右核中有IKZF2、下缘有T，T内模型0.49，维持扩展。
- [D1_Johnson cell_10826](proseg_observed/D1_Johnson/cards/cell_10826.jpg)　凹入的较大推定细胞与上方大核紧贴，唯一FOXP3正落向右下相邻核的接缝，无法排除旁源。
- [D1_Johnson cell_11443](proseg_observed/D1_Johnson/cards/cell_11443.jpg)　小核右上缘直接T与IKZF2，27 UMI；旧Treg分数高但新主群支持有限。
- [D1_Johnson cell_12088](proseg_observed/D1_Johnson/cards/cell_12088.jpg)　核内有T与IKZF2，周围多有T，未见FOXP3，两个模型支持均不足主集合。
- [D1_Johnson cell_12941](proseg_observed/D1_Johnson/cards/cell_12941.jpg)　橙色凹形核上下均有邻核，直接T及IL2RA落在同一左下接缝。核归属难辨。
- [D1_Johnson cell_13055](proseg_observed/D1_Johnson/cards/cell_13055.jpg)　新检出核沿右侧大核呈极窄凹月形；虽有T及IL2RA，难从本图分辨独立胞体。
- [D1_Johnson cell_39290](proseg_observed/D1_Johnson/cards/cell_39290.jpg)　背景较稀疏，核右缘见IL2RA，T在上下边；主群支持接近门槛，保留扩展。
- [D1_Johnson cell_39968](proseg_observed/D1_Johnson/cards/cell_39968.jpg)　核内T清楚，单个IL2RA格在下方紧贴的另一核交界；主群另以上皮为首选。
- [D1_Johnson cell_56524](proseg_observed/D1_Johnson/cards/cell_56524.jpg)　大核下缘有数格T，中心IL2RA，浆细胞竞争且深度高；保留为扩展冲突对象。
- [D1_Johnson cell_57656](proseg_observed/D1_Johnson/cards/cell_57656.jpg)　分界清楚的小核中直接T和FOXP3同在，旧常规CD4、新模型更偏上皮与Tfh，仍仅扩展。
- [D1_Johnson cell_60550](proseg_observed/D1_Johnson/cards/cell_60550.jpg)　下半核右缘两个T及IL2RA，左上与另一核接壤，标记位于另一侧；pT未达主门槛。
- [D1_Johnson cell_60816](proseg_observed/D1_Johnson/cards/cell_60816.jpg)　直接T集中于右核缘，IL2RA与CTLA4支持在两侧，FOXP3属左邻核。T内模型略低于0.5。
- [D1_Johnson cell_62600](proseg_observed/D1_Johnson/cards/cell_62600.jpg)　直接T有多格，IL2RA在顶缘近邻核；上皮主群竞争，不升级为主集合。
- [D1_Johnson cell_79326](proseg_observed/D1_Johnson/cards/cell_79326.jpg)　椭圆核旁分见直接T与CTLA4，旧Treg及新T内分类均支持，但主群pT为0.46。
- [D1_Johnson cell_79555](proseg_observed/D1_Johnson/cards/cell_79555.jpg)　哑铃形双叶轮廓，直接T集中左上、IL2RA在凹口；不能确认同一细胞的转录来源。
- [D1_Johnson cell_11496](proseg_observed/D1_Johnson/cards/cell_11496.jpg)　小核底缘直接T，28 UMI；FOXP3在左下远于黄色范围。
- [D1_Johnson cell_9577](proseg_observed/D1_Johnson/cards/cell_9577.jpg)　细长核为浆细胞首选，FOXP3及直接T散在右侧，未归到本核。
- [D1_Johnson cell_12430](proseg_observed/D1_Johnson/cards/cell_12430.jpg)　核内多格T与顶缘FOXP3可见，旧树突样、新上皮分类相冲，主群pT仅0.04。
- [D1_Johnson cell_13019](proseg_observed/D1_Johnson/cards/cell_13019.jpg)　核左下缘单个T，FOXP3在左上方，主要证据仍为普通T线索。
- [D1_Johnson cell_10722](proseg_observed/D1_Johnson/cards/cell_10722.jpg)　空隙边缘核顶凹入，T与FOXP3均在周围，中心缺直接T。
- [D1_Johnson cell_11629](proseg_observed/D1_Johnson/cards/cell_11629.jpg)　密集区T及FOXP3与核接近，旧CD8、新髓系；T条件概率高但主群T概率低。
- [D1_Johnson cell_11749](proseg_observed/D1_Johnson/cards/cell_11749.jpg)　狭长核右缘有T，FOXP3在核右侧远处，新髓系竞争明显。
- [D1_Johnson cell_5207](proseg_observed/D1_Johnson/cards/cell_5207.jpg)　直接T与FOXP3靠圆核底边，主群T支持接近零；保留分子线索。
- [D1_Johnson cell_54972](proseg_observed/D1_Johnson/cards/cell_54972.jpg)　血管样空隙旁核上缘见FOXP3，原树突与新髓系一致，若干T格位于外围。
- [D1_Johnson cell_34393](proseg_observed/D1_Johnson/cards/cell_34393.jpg)　橙核左下的T与FOXP3均在圈外，43 UMI，新模型偏B。
- [D1_Johnson cell_32061](proseg_observed/D1_Johnson/cards/cell_32061.jpg)　孤立核上侧的FOXP3格在推定范围外，胶原背景丰富，旧周细胞与新间质相符。
- [D1_Johnson cell_32298](proseg_observed/D1_Johnson/cards/cell_32298.jpg)　宽而长的核边有T和FOXP3，肥大细胞竞争强，pT仅0.13。
- [D1_Johnson cell_60017](proseg_observed/D1_Johnson/cards/cell_60017.jpg)　小核右下缘单个FOXP3，T点在邻核周围，中心没有直接T。
- [D1_Johnson cell_43407](proseg_observed/D1_Johnson/cards/cell_43407.jpg)　核内T及FOXP3空间可辨，59 UMI且主群偏髓系，保留冲突。
- [D1_Johnson cell_41316](proseg_observed/D1_Johnson/cards/cell_41316.jpg)　近圆核左下有T，FOXP3位于右侧胞体外；新浆细胞与旧树突竞争。
- [D1_Johnson cell_55068](proseg_observed/D1_Johnson/cards/cell_55068.jpg)　靠空隙长核的下缘见FOXP3，T在周围多个格；最终中心直接T为0。
- [D1_Johnson cell_62896](proseg_observed/D1_Johnson/cards/cell_62896.jpg)　密集免疫区的T和FOXP3常在本核外，新pT大于0.5，但中心直接T为0。
- [D1_Johnson cell_62943](proseg_observed/D1_Johnson/cards/cell_62943.jpg)　右核缘一格T，FOXP3在下方远离核。常规CD4支持强于Treg。
- [D1_Johnson cell_79794](proseg_observed/D1_Johnson/cards/cell_79794.jpg)　长核左缘FOXP3和T共存，模型主群偏B，旧Tfh解释仍保留。
- [D1_Johnson cell_63262](proseg_observed/D1_Johnson/cards/cell_63262.jpg)　浆细胞首选，中心一格T，FOXP3在下缘邻接处，未落入本次胞体。
- [D1_Johnson cell_79884](proseg_observed/D1_Johnson/cards/cell_79884.jpg)　不规则小核内一格T，FOXP3留在左下方，pT仅0.25。
- [D1_Johnson cell_80172](proseg_observed/D1_Johnson/cards/cell_80172.jpg)　大圆核上缘贴近数个T格，旧树突与新髓系一致，FOXP3未归属。
- [D1_Johnson cell_60890](proseg_observed/D1_Johnson/cards/cell_60890.jpg)　单FOXP3贴核右侧，背景胶原/髓系密集，新旧均优先浆细胞。
- [D1_Johnson cell_60913](proseg_observed/D1_Johnson/cards/cell_60913.jpg)　楔形小核内有多格直接T，FOXP3在外侧，常规CD4优先。
- [D1_Johnson cell_103894](proseg_observed/D1_Johnson/cards/cell_103894.jpg)　上皮丰富区中窄核，强新旧上皮证据；T仅一格，FOXP3未归属。
- [D1_Johnson cell_130418](proseg_observed/D1_Johnson/cards/cell_130418.jpg)　胶原点在视野内密集，旧周细胞与新间质一致。T及FOXP3空间存在但主群不支持T。
- [D1_Johnson cell_136311](proseg_observed/D1_Johnson/cards/cell_136311.jpg)　上皮长核左尖端有T，FOXP3在左上圈外，新旧上皮一致。
- [Cervilla n_9571](proseg_observed/Cervilla/cards/n_9571.jpg)　上皮/间质密集，T点和单个核心点分开，中心无归属格。
- [Cervilla n_121947](proseg_observed/Cervilla/cards/n_121947.jpg)　核列旁的空白夹隙中没有中心归属，周围见上皮点，无法定位T。
- [Cervilla n_152794](proseg_observed/Cervilla/cards/n_152794.jpg)　组织边缘疏松间质及狭长核可见，中心无指向某个核的T或核心捕获格。
- [D1_Johnson n_42](proseg_observed/D1_Johnson/cards/n_42.jpg)　组织突起与明亮腔隙相接，重检后核数有变化，中心没有可归属的T或核心格。
- [D1_Johnson n_193524](proseg_observed/D1_Johnson/cards/n_193524.jpg)　组织小梁间空隙，旧中心大核未保留，未见T或FOXP3捕获格。维持未定位。
- [D1_Johnson n_15182](proseg_observed/D1_Johnson/cards/n_15182.jpg)　旧小核附近有T、FOXP3及另一核心，但重检中央未保留单一实例。不能排除漏核。
- [D1_Johnson n_37297](proseg_observed/D1_Johnson/cards/n_37297.jpg)　旧中央轮廓消失，T在几个不同核周，单个FOXP3靠右下腔隙。
- [D1_Johnson n_43448](proseg_observed/D1_Johnson/cards/n_43448.jpg)　强红染间质内旧小核未保留，T与FOXP3格互相分散。
- [D1_Johnson n_64753](proseg_observed/D1_Johnson/cards/n_64753.jpg)　几个新核围绕中央未归属区，FOXP3在偏左方，无法定位到同一核。
- [D1_Johnson n_65237](proseg_observed/D1_Johnson/cards/n_65237.jpg)　密集免疫背景中中央旧轮廓消失，单个FOXP3邻近多条离散T信号。
- [D1_Johnson n_65710](proseg_observed/D1_Johnson/cards/n_65710.jpg)　旧核未能与新实例对上，FOXP3靠上方两个核间隙，T多在其下。
- [D1_Johnson n_94968](proseg_observed/D1_Johnson/cards/n_94968.jpg)　旧辅助位点旁有多格T和FOXP3；中央未见可可靠分割的对应新核，保留局部核验价值。
- [D1_Johnson n_95364](proseg_observed/D1_Johnson/cards/n_95364.jpg)　一个FOXP3在上方邻核边，另一个靠近中央空隙，不能按最近中心任意归并。
- [D1_Johnson n_95509](proseg_observed/D1_Johnson/cards/n_95509.jpg)　旧中央核未保留；FOXP3在其左下、其他核心在上方邻核，证据分离。
- [D1_Johnson n_170561](proseg_observed/D1_Johnson/cards/n_170561.jpg)　旧辅助位点处深红组织无对应新核，T和FOXP3捕获格靠近；有线索但归属未定。
- [D1_Johnson n_175312](proseg_observed/D1_Johnson/cards/n_175312.jpg)　强红染区旧实例消失，多格T落在无核归属区；单格FOXP3靠下方。
- [Cervilla n_59711](proseg_observed/Cervilla/cards/n_59711.jpg)　旧非mixed FOXP3位点的中心旧核未保留；原图低分辨率，中央粉色格、两侧直接T与另一核心信号分散，无法归属于同一核。

每个来源到细胞的多对多映射见各片 source_links.parquet。全队列审阅状态见 review_index.parquet，打开的卡与尚未打开的卡均可按 review_id 取出。
