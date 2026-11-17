# 实际图像复核

已打开检查Cervilla candidate_1、D1 candidate_4、UC_SCC candidate_1、MPBC2 candidate_1完整图册、全部24个ROI的原始marker热图及主森林图。其他ROI同样交付，可按索引逐图查阅；不声称对每幅图作病理专家判读。

- Cervilla：注册H&E可见密集核；T-lineage信号成片，但核心Treg基因低UMI、空间稀疏。候选还可见CXCL13/PDCD1，不能只靠FOXP3或形态升级为Treg细胞。
- D1：展示区域存在组织内T信号、IL2RA/IKZF2/CCR8等组合支持，附近竞争状态不完全消失。原始面板支持区域Treg-like线索，不证明每个组成bin为Treg。
- UC_SCC：下采样H&E只能支持粗组织位置；原始核心Treg信号零散、Tfh可见。不能据图推断细胞边界。
- MPBC2：纯相对阵列图、无H&E，明确保留组织未知；弱Treg候选与强T-lineage对照并不等价。
- 多个候选ROI中Epi评分面板为空：当前窗口没有满足资格的Epi区域，不是MP10=0，不扩大50µm来填补。
- 森林图原来图例靠近首行，已移出主体；仅描述性终点改为空心点，避免把没有CI的点解释为高度确定。

所有ROI黑线为64µm、红圆为50µm，实际坐标范围为256µm。未采用平滑、插值或空间多数票填标签。

补充支持范围：图册内已加入橙色中心T区域成员bin与青色Epi评分成员bin空心标记，用实际成员位置展示评分支持，未用凸包跨越空缺。
