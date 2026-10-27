# 空间通讯输入与预先固定的设计

本目录只读复用既有 63 张普通 Visium 的 hard-QC 后组织内完整切片、冻结 Top20 MP10、accepted cell2location v2_30type 后验均值、真实阵列索引，以及已修订患者身份表。没有新训练、CNV、niche 或区域分析。

在查看任何通讯结果前，以完整表达/坐标/后验均值、20/20 MP10 基因可测和计算规模为标准，选择最小公共队列切片 GSE171351_A1（768 spots）及最小 Inhouse 切片 Inhousedata_T4469479（3917 spots）。全片原 hard-QC 后位置均保留。`candidate_sections.tsv` 列出63片候选，不使用Treg通讯大小选择。

## 文件接口

每片 `sections/{section_id}/`：

- `counts.mtx.gz`：gene × spot，非负实测整数counts；`counts_spot_gene.npz` 是同一矩阵的Python稀疏CSR格式。
- `log1p.mtx.gz`：gene × spot，使用全部原始测量基因进行每spot总量1万归一化，再自然log1p；`log1p_spot_gene.npz` 同一数据的CSR格式。不得再次对这一层log1p或用c2l数值代替表达。
- `genes.tsv` 顺序对应矩阵行；`meta.tsv` / `spots.tsv` 顺序对应矩阵列。`spot_id = section_id::barcode` 为唯一位置键，原barcode单列保留。
- `coords.tsv` 与meta同序，含标称平台阵列μm坐标及原lowres图像像素坐标、array_row/col。
- `Q_original_30.tsv` 是完整30个accepted后验均值除以每spot全部30类均值之和，不使用q05。`abundance_mean_original_30.tsv` 保存同一原始后验均值。
- `Q_fine.tsv`、`Q_major.tsv` 是用于软聚合的完整组成背景；第一列spot_id，其他列均为Q。每行和为1，major/fine独立从同一原始30类构造，不能将fine均值通讯相加代替major。
- `input_summary.json` 提供来源、规模、分组和尺度信息。

## MP10与共享组成

MP10原值复用 `EpiMP10_Top20_AUCell`，两片均完整20基因，未重新计算。历史冻结合格Epi anchor为posterior mean Epithelial / 全部30类型mean总和 ≥0.25。此定义不同于历史RCTD anchor，不替换。

本轮没有找到现成可靠的该空间通讯H/L分组，故新增 `epi_state`：在每片既有合格anchors内计算连续冻结MP10中位数，严格大于为high、严格小于为low，等于或不合格为other。保存阈值而不覆盖原字段。Q_EpiH=Q_Epi×high指示，Q_EpiL=Q_Epi×low指示，其余质量为Epi_MP10_other。三个上皮列相加等于原Epithelial。TME仅有同一份Q，从未按H/L复制。

节点实际含义是“MP10-high/low关联spot中的上皮成分”，并非同spot内的上皮状态已被分辨。Visium结果是“细胞类型比例加权的跨spot通讯归属”，依赖去卷积与实测混合表达，不能确定每条LR由同spot哪类细胞表达。

## 可辨识层级及支持

使用 `cell2location_final_validation_gate_v1/09_final_status/DOWNSTREAM_USAGE_TABLE.csv`。三个原列 CD4 Helper/Conventional、CD4 Naive/Memory-like、CD8 GZMK/Effector-memory-like 被原审查标记LOW_IDENTIFIABILITY / DO_NOT_INTERPRET_SEPARATELY，本轮将其Q质量相加为 `T_NK_unresolved` 背景而不重新分配给其他命名细胞。不把该合并称为新生物亚型；原30列仍完整保存。Tip EC、Cycling EC未在空间参考独立建模，明确NE，不制造。

Treg有独立Q列，保留且标记INTERPRET_WITH_CAUTION，即使从不占dominant也不删除。其余多数可辨EC/CAF等同样需要谨慎。各类型支持按原验证已有规则记录：至少10个spot后验均值abundance ≥0.01、全片总abundance >1。此规则用于类型可评估性，而非对位置Q截断。所有positive Q均保留，不使用0.1比例截断；未新增噪声筛选。本次所有正式Q列均达到该覆盖规则。

空间原major有独立 `Mural`（Pericyte、SMC），保持原映射。与scRNA的Fibroblast/基质可联合文字解读，但不能改名CAF或默默合并major。`T_NK`仅显示规范化为 `T/NK`。

## 坐标与距离

原对象是普通Visium，有经过旧审计恢复的真实交错array_row/array_col。依据[10x官方平台技术说明](https://cdn.10xgenomics.com/image/upload/v1722285090/support/documents/CG000686_VisiumHDSpatialGeneExpressionPerformance_TechnicalNote_RevA.pdf)所给标准Visium 100 μm中心间距/55 μm捕获点径，建立标称阵列坐标：x_um=50×array_col，y_um=50×sqrt(3)×array_row。以完整切片实测最近邻中位数核对为100 μm；阵列与原图像位置的仿射拟合RMS均<0.01像素，确认位置几何相符。

这是以明确平台规格标定的阵列物理坐标，不声称恢复了显微镜扫描μm/pixel。图像中的 `scale_factors/spot` 不可盲用：Inhouse值0.06860122与hires倍率相同，不是可信的捕获点像素径。原图像像素仍保存；两片实测图像最近邻间距与阵列映射均在表中。没有照搬教程的spot size/contact.range。

在看H/L通讯前固定分泌范围2×d_NN=200 μm，接触范围1.05×d_NN=105 μm（最近邻+5%容差），ECM为200 μm（依SpatialCellChat原生API与分泌共用interaction.range；在首片结果出现前统一固定）；算法原生类别/API限制另由各方法记录。第一片仅附加一次较宽范围敏感性：SpatialCellChat因原生API同时将ECM与分泌改为3×d_NN=300 μm，COMMOT只将分泌改为300 μm、ECM保持200 μm。此区别单列，不宣称完全匹配。邻spot接触类信号不证明真实膜接触；位置距离是采样分辨率下约束。

## 身份与外推

GSE171351_A1来自既有来源核验患者Gouin2021::1204-1A（VERIFIED_WITHIN_STUDY）。Inhousedata_T4469479患者身份没有核验，保留USER_ASSUMED_SAMPLE_UNIT_PATIENT_ID_UNVERIFIED及匿名记录单位；不能把两片称为两个已核验患者。每片单独模型，不拼坐标、不以spot/LR数充当患者。

目标目录High/Low图只是裁剪同一网络显示，不重新推断。原空间网络的同位置对角贡献须另存；跨位置聚合不意味着重新求解了禁止同spot通讯模型。

## 高分辨率现状

既有高分辨率审计记录Visium HD 8μm bin及Xenium原始包，但未发现可直接复用且有可靠最终细胞注释的完整表达+分割+坐标交付。故本轮不额外开展可选cell-level ROI；HD bin不被称为真实单细胞，不为此重新分割或注释。
