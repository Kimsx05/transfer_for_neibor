# v3.2：宏观病灶位置与外缘修订

## 对象与输入

输入为 v2 manifest 中的 63 张 ordinary Visium 切片、7 个 dataset。逐点数值来自完整 `results/mp10_pathology_regions_v3_1/by_dataset/`。保留 v3/v3.1，不重跑 BANKSY、SpaCET、CopyKAT、cell2location 或 MP10。不开展新的 MP10/TME 比较。原始 barcode、坐标、分子数值和缺失状态逐字段校验；完整组成保留在各片 parquet。

宏观位置、恶性/细胞组成证据、界面可信程度、精确距离和距离分带分别保存。`inside` 是计算支持的病灶内位置，**不是该点细胞全为恶性，也不是 deep core**。`extent_unresolved` 表示宏观范围不能确定，不代表新的组织类型、正常组织或组成混合。

## 从旧分数到新位置

复算旧点分数：0.65 × SpaCET malignant fraction + 0.20 × CNV 支持编码 + 0.15 × clip(EPI/2.5, 0, 1) − 0.45 × 平滑肌背景指示量。CNV 编码沿用 strong=1、moderate=0.55、baseline caution=0.35、native-diploid conflict=0.30、weak/unavailable=0；该编码仅用于复算旧规则，不将 unavailable 宣称为阴性。平滑肌指示量沿用 SMC≥0.5、检出平滑肌 marker≥3、EPI<0.8。

在真实六邻接、同组织块、同 BANKSY 域中计算局部均值；不足 5 点的域连通部分另复算旧 `0.25×local + 0.75×同块同域均值`。保存 `domain_local_evidence_score`、`domain_shared_evidence_score` 及仅共享后跨过 0.42 的标记。新旧共享分数必须逐点一致。分数不是经病理校准的概率，旧 Interface-related domain 不是 200 μm margin。

固定规则如下：

1. 局部、未借用远处信息的分数≥0.42，作为病灶包含组织的计算候选种子；seq_depth 回退不能单独升级。沿用 N1/N2/N3 已有正常来源与尿路上皮/肌组织背景裁决及局部冲突，不把模型高值自动视为肿瘤。
2. 保留与上述局部种子处在**同一原 ROI 观测连通组**内的病灶主体及原有内部区室。通过实际观测邻接继承位置，不转移恶性身份，不跨缺位连接。独立的小 ROI 没有局部种子时，不能仅由同域共享分数保留。
3. 对现存完整观测邻域内、周围只有 inside 的未定组织区室，可以赋予内部位置；要求相邻位置完整（每点六个邻居）、已有组织图像前景中位数≥0.35、nFeature 中位数≥500。最多两次同步更新，不填表达、不填未观测洞、不做多数滤波或凸包。
4. 正向 outside 支持限于现有正常来源背景裁决，或连通域中的明确平滑肌程序：平均 SMC≥0.5、平均检出平滑肌基因≥3、平均 EPI<0.8、平均恶性比例<0.2、strong CNA 比例<0.1。它是计算的非病灶组织候选，不能排除显微浸润。普通低恶性基质/免疫域不因此成为 outside；不满足 inside/outside 正向依据的范围保留未定。
5. BC18 另采用作者 S1 的同片 Tumor Cell/Tumor Stroma 共存背景作**有限校准**：一个未定域必须直接接触至少两个不同的局部肿瘤支持域、跨域观测接触边至少60%接向这些域，且本域平均 malignant≥0.2、EPI≥0.5、SMC<0.5，才能扩展为内部候选。不是图片转 barcode，不能将 Others 全部吞入。该规则仅在有此同片背景的 BC18 使用，明确记录于 config；不据下游结果选择。
6. 仅共享后跨阈值的小域，若完整的直接邻域一致为已支持的位置，可继承位置；否则缺乏同 ROI 局部支撑的独立候选改为范围未定。大小本身不是删除标准，强局部支持的小岛保留为候选，不称卫星灶。

不把没有 outside 种子当作全块 inside，不把没有肿瘤种子当作正常。全块 inside 要求每个保留位置都有上述正向或同一观测病灶结构内部支持；全块 outside 同理。组织块状态和点级位置分开汇总，微小断片与较大单一区域分别报告。

## 外缘、几何、距离

全部切片从新标签重新生成 ROI ID（`:v3_2:L*`）、hex 多边形、边段 ID、参考 ROI、连续距离与三种宽度分带。实际六邻接为 `(0,±2),(±1,±1)`。同块内 inside 与有正向依据的 outside 直接相邻，且 outside 组达到组织外侧（必须接向原始全组织网格以外；仅因 QC 缺测通道连通空白不算），才记为 `supported_macro_outer`。未知位置接触、内部非病灶区室接触、QC/缺位边、真实或未知缺损边、采样/组织切缘分别保存。**这里的 supported 是有限计算支持，不是病理医师确认，也不自动等于侵袭前沿。**

无支持外缘时，inside 保持红色位置标签，距离 NA，显示“Inside; outer edge unobserved”，不生成 margin/core。有内部红区但宏观范围未定时，不能用红区候选轮廓冒充生物外缘。

距离模型延续 v3.1。普通 Visium 标称中心间距100 μm，坐标为 `x=50×array_col`、`y=50√3×array_row`；使用已有验证的图像仿射，不冒称独立组织尺寸测量。距离零点为 hex 共享边而非 spot 中心，内部负、外部正。

重新检查内部缺孔：只有原始组织网格确实存在、图像前景每点≥0.65、周围同一组织块/位置且 inside 时同一新 ROI、坐标残差合格，才作为几何路径的技术缺位支持。无新增观测行、无分子插补。空白图像缺孔为 cavity/defect 候选，不自动叫坏死；其余未知缺口保留。路径仍为直线路径，25 μm 采样、距同块已知几何中心≤58 μm，不改测地距离。

令 K 为观测且位置已知的 hex 区域加合格技术缺位；所有 K 外（包括范围未定点与采样外）允许有未观察前沿。g 为到同 ROI（inside）/同组织块（outside）支持外缘的距离，q 为到 K 的未知区域边界距离。该明确未知模型下，绝对距离下界 L=min(g,q)；直线可见时上界 U=g。**q 包含 outside 一侧全部未知区域，不只是 ROI 缺损边**。q>g+1e−6 且路径可见时才保存精确距离；相等仍删失。没有支持外缘则不构造数值界。

同一几何下 w=100/200/300 μm：精确距离可估时按距离；U≤w 判相应侧 margin，L>w 判超过阈值，其他分带未定。仅分带恢复不会将候选 g 写成精确距离。

保留原 v2 肿瘤小巢接触边作为次级几何对照，仅重连其 v3.2 内部位置成员关系；不与宏观外缘混用，不新增接触分析。

## 分子与作者分区对照

EPI 是原先可用 EPCAM/KRT8/KRT18/KRT19/KRT7 的 log1p(CP10k) 均值。保持原分数；每片可用基因见 `marker_availability.tsv`。缺失基因显示 unavailable，绝不填0。各基因共用全队列有限数值99百分位色标；SpaCET用0–1，距离用±2000 μm。截断说明写入图中；这些表达数据来源相同，不能声称独立验证。

作者专业注释只采用已公开材料：

- [SpatialBC / GSE319536](https://doi.org/10.1158/2159-8290.CD-26-0099)：[作者代码](https://github.com/VividKai/SpatialBC)读取未随仓库提供的 `*.seurat.umap.spatial.rds` 与 Level1_Annot/Level2_Annot。Level1 是 Tumor/Normal/Other Region，Level2 含 Tumor Cell、Tumor Stromas、Normal Stromas、Others 等。实际下载的[补充文件](https://ndownloader.figshare.com/files/69527807)是 DOCX，S1/S9 原图已提取留存。S1 仅作同片定性比较；S9 是相对距离分区，其大范围蓝区在部分切片与 S1 Tumor Stroma 重叠，二者不可直接视作本项目同一个宏观对象，更不可替代固定 μm 分带。未取得 barcode 级作者标签，不计算 Dice/IoU。
- [GSE285715](https://doi.org/10.1073/pnas.2505125123)：GEO 支持 T/N 原始样本来源；论文 marker 细胞组成图不证明 fibroblast/immune region 位于病灶外。没有已核验 Tumor-419a↔T3 映射，不作伪同片验证。
- [NPJ](https://doi.org/10.1038/s41698-026-01417-x)：未恢复本地 sample1/sample2 与作者编号的可靠桥接，不能转移作者肿瘤/侵袭前沿标签；同用 SpaCET/表达也不是独立实验验证。
- 其他数据集只继承现有来源；没有可用同片专业分区时明确写未获得参考，不编造 Inhouse 出处。

完整比较简表逐片记录来源、可用形式、层级及限制。BC18 的作者信息用于校准，因此仅称校准后的定性对照，不称独立验证。未建立用户或 agent 的新病理真值。

## 执行与复现

先运行 BC18、T3、B16、BC01、N1、B15，检查主体保留、局部借分伪影、图像方向与无外缘行为，随后固定规则全队列两片并发。单进程 BLAS/OpenMP 线程为1，无全队列稠密距离矩阵。Python 3.12；沿用项目 NumPy/Pandas/SciPy/Shapely/Matplotlib/Pillow/PyArrow 环境；补充 DOCX 的图像使用标准库 zipfile 提取；临时阅读软件已清理，不依赖新增分析环境。

在项目根目录执行：

```bash
python3.12 results/mp10_pathology_regions_v3_2/scripts/prepare.py
python3.12 results/mp10_pathology_regions_v3_2/scripts/run_v3_2.py --pilot
python3.12 results/mp10_pathology_regions_v3_2/scripts/run_all.py
python3.12 results/mp10_pathology_regions_v3_2/scripts/references.py
python3.12 results/mp10_pathology_regions_v3_2/scripts/summarize.py
python3.12 results/mp10_pathology_regions_v3_2/scripts/finalize_metadata.py
python3.12 results/mp10_pathology_regions_v3_2/scripts/check_results.py
```

已完成切片自动跳过；单片重做用 `run_v3_2.py SECTION --force`。代码依赖同项目已有 v3 只读几何帮助函数和 manifest。未修改上游。身份/分析单位沿用 v3.1 校正：B5/B18 不因状态文字不同拆单位；T1/T2/T3 仍为未核验独立切片单位，不新增患者或T/N配对。当前无新生物学聚合结果。
