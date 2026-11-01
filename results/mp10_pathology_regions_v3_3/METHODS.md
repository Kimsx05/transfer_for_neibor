# v3.3 方法：候选位置、支持程度、距离分开

## 输入与主规则

沿用63张/7个dataset的既有manifest，逐点输入为完整v3.1，v3.2提供历史状态、局部/共享分数、小ROI审阅及作者材料。没有重训上游模型，也没有新MP10/TME统计。原始barcode、坐标、全部marker/组成/MP10数值与缺失状态保持不变。

**主位置从v3.1初始化。inside与outside都表示计算候选估计，结论标准相同。** v3.2仅因狭窄outside条件撤掉的位置恢复为候选outside，支持程度另记；原unlocalized/QC不自动恢复。没有按颜色比例或下游结果调参。v3.2的新增inside也不自动保留，尤其不因取消同域收缩，就把局部分数过线的孤点升为新ROI。

主要字段：`lesion_location`、`position_support`、`reference_status`、`assignment_basis`、`change_class`；历史位置分别在`v3_1_lesion_location`和`v3_2_lesion_location`。`boundary_interpretation`与`reference_boundary_support`解释边界含义；它们不决定点是否有几何坐标，也不作为候选距离的准入条件。

## 小ROI与B16

逐一处理v3.2识别的57个弱借域ROI、87点。复用局部及共享分数，不重算上游恶性模型。只有对象内各点的六邻位全部观测、外邻域在不变的v3.1候选图上全为outside、自身局部分数最大值仍<0.42，才撤销独立ROI并继承outside。缺测、混合邻域或其他依据不足时保留v3.1位置，标为不稳健小候选。保留位置不等于确认恶性；作者粗参考可在后续步骤校准，逐对象表分别记录局部决定和最终决定。

B16的`AGATAACACGGTATAA-1_4`（89,175）局部分数0.463318，共享分数0.415650，六邻均为v3.1 outside。它回到基线outside，较高局部分数保留为属性。上一轮将它解释为“inside邻域支持”不符合实际，已明确纠正。

本轮不对其他小岛按大小删除，不做多数滤波、凸包或全局填孔。所有邻域判断都基于实际六邻接，不跨缺测区/组织块。

## 作者粗区域参考：同时检查两侧

复用[SpatialBC论文](https://doi.org/10.1158/2159-8290.CD-26-0099)、[作者代码](https://github.com/VividKai/SpatialBC)及[已下载S1补充图](https://ndownloader.figshare.com/files/69527807)。未取得作者逐barcode标注；S9相对距离区域不替代本项目距离带。

对16张同名BC切片，以**原始H&E到作者图配对H&E**做相似变换配准，完全不使用位置/表达标签拟合：SIFT、比值检验0.75、RANSAC误差2作者像素；至少8个内点、内点比例≥0.35、空间覆盖≥0.12、95%内点误差≤2.5像素、H&E灰度相关≥0.5才允许粗参考转移。保存矩阵、旋转、尺度、误差及H&E叠图。相似变换不能拟合镜像；不能获得可靠正常方向对应的切片拒绝转移，而非强配准。BC01/BC13/BC18经过图像检查；其余也保留同格式诊断。

从原图图例采样颜色，仅提取明确色块。Tumor Cell+Tumor Stroma构成粗inside参考；Normal Stroma+Urothelium构成粗outside参考。肌肉、坏死、免疫、血管、Others等不自动二分。原图颜色欧氏差<28，原同宏观色类连通区至少60像素。为抑制印刷/抗锯齿小孔影响，在作者像素邻域中要求同类占比≥0.70、相反宏观类≤0.02，形成粗区域共识。

**误差排除半径 r=max(4, 2×H&E配准误差P95+2) 作者像素。** 点到粗共识区域外的距离必须≥r，且中心本身仍有明确原作者对应色类，才接入主位置校准；不是只用邻域平均来代替边缘排除。其他位置保留计算候选并记录参考边缘/潜在冲突。对Others不会通过平滑继承二元标签。这里的像素处理只用于低分辨率参考的不确定带，不改变组织网格或填充病灶ROI。

可靠粗inside/outside对主位置的校准是对称的：既可修正漏掉的Tumor Stroma，也可撤回与Normal Stroma核心冲突的红区。**称为“作者图派生的粗区域参考”，不是作者barcode真值，更不是agent新增病理真值。** 改变点数是本流程输出数量，不是分割准确率；不报告Dice/IoU。参与修订的参考是校准信息。

BC18不再沿用“D8、两个域、60%接触”的特例；按同样的双侧粗参考规则复核全部可对应区域。T3没有同片作者宏观轮廓，恢复可检查的候选分区，明确边界可能为内部组成/状态分界；广泛EPI和肿瘤来源都不被用于将整片强行归inside。NPJ等未核验的图像/样本桥接不转移标签。

## 候选几何、距离与分带

复制并修改已有v3.2几何函数，全部切片从v3.3位置重建ROI/边段/参考ID。主边段是**计算候选ROI的可观测外缘**：观测inside/outside直接邻接、outside观测连通部分到达组织外侧；仅经QC缺测通道连通空白不算外部观测路径。内部区室接触、unlocalized接触、QC缺位、缺损、组织/捕获截断分别记录。原小巢接触对象单独保留，更新位置成员关系，不混成主外缘。

`candidate_observed_outer`直接进入主距离，**不再先筛选`supported_macro_outer`**。作者支持仅在几何计算完成后追加为边段注释（单侧/双侧350μm邻域粗背景；不能据此声称轮廓精确匹配）。v3.2旧candidate字段未被直接复用。

普通Visium标称中心间距100μm，`x=50×array_col`，`y=50√3×array_row`；图像方向及尺度复用已验证仿射。零点是共享hex边，内部负、外部正。inside只参考自身新ROI，outside只参考同一组织块，直线路径不得跨空白；未改成测地距离。

重建已知几何K：所有当前inside/outside观测点，不因支持较弱删掉。合格内部技术缺位需原始组织网格记录、每点图像前景≥0.65、同一周边块/位置（inside还需同ROI）及合格坐标映射。只用于几何，无新增分子观测。真实/未知缺口、未定位及组织截断仍限制距离。

沿用v3.1保守距离范围：g为到可观测候选外缘距离，q为到K外全部未知区域边界距离；L=min(g,q)，直线可见时U=g。q>g+1e−6时精确距离可估；相等仍删失。该未知模型包含outside一侧，不只看ROI的缺损边。路径每25μm采样，距同块支持几何中心≤58μm。没有候选界面时距离和范围NA，位置仍可保留。

同一几何下，100/200/300μm分别按U≤w判margin、L>w判超过阈值，否则分带未定。仅分带可判定不会填补精确距离。字段兼容保留`Interior (>w confirmed)`等名称；**confirmed仅指当前候选几何下距离超过阈值，不是病理core确认**。整块inside、无外缘时不生成core/margin。

## 图形、身份、复现

各片：三版本位置对照；当前位置/H&E/支持程度；分带/距离可靠性；EPI多marker/SpaCET/CNA；作者参考与冲突。GSE319536另有配准与原作者图对照。表达色标原样复用v3.2全队列色标，缺失marker仍unavailable；EPI是共享表达证据检查，不是独立验证。身份表原样复用v3.1已校正版本，B5/B18同P1，T1/T2/T3保留未核验独立切片单位，不新增配对。

先完成7张代表片，再两片并发全队列；没有上游重训、MP10/TME新分析或大对象缓存。原版数据和脚本不改。OpenCV已在项目环境可用，其余几何/绘图依赖沿用；小型版本表见`software_versions.json`。

项目根目录复现：

```bash
python3.12 results/mp10_pathology_regions_v3_3/scripts/register_author.py
python3.12 results/mp10_pathology_regions_v3_3/scripts/run_v3_3.py --pilot
python3.12 results/mp10_pathology_regions_v3_3/scripts/run_all.py
python3.12 results/mp10_pathology_regions_v3_3/scripts/summarize.py
python3.12 results/mp10_pathology_regions_v3_3/scripts/write_results.py
python3.12 results/mp10_pathology_regions_v3_3/scripts/check_results.py
```

已完成片自动跳过；单片修改后使用`run_v3_3.py SECTION --force`。若改变作者配准/粗参考，需要重跑相应GSE319536片。完整点表、已执行代码、配置和参考变换均在本目录，原S1图从v3.2本地来源复用。
