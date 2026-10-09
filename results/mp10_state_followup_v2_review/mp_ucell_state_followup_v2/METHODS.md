# 方法与解释范围

## 冻结输入
所有原NMF、UCell、聚类和降维均只读。FULL严格沿用冻结01_input_state.qs$cells和旧cell_metadata；123699细胞、57个specimen、5个dataset，保留22436个Normal。按cell_id连接标签/分数/坐标。C8来自r0.4，C12/C16为r0.6完整群，两套resolution分别处理。sample=dataset|orig.ident。37个sample的患者映射为derived、20个为provisional；41个recorded patient key不是41位已核实独立患者，其中16个key各对应两个sample。

## 原usage衔接
只使用usage_raw，未变成行和比例。每群14MP的均值、中位数、Q25/Q75、零比例同时保存usage及UCell。六项比较均正向=左侧高；各组至少30细胞的sample可配对，使用组内全部细胞，不重采样。完整57样本资格及文库大小另存。pooled细胞加权结果与sample等权结果分开。样本均值双侧配对Wilcoxon，exact=FALSE，至少3对；全零/不足样本返回NA及原因。BH在每个comparison×scale内14MP单独校正；dataset/source分层表同样在各自层内14MP校正，不声称整个多表筛选的家族错误受控。冻结MP10-high复用10_high_usage_global_cutoffs.csv的原全局90% quantile阈值和usage>=阈值定义；此轮没有重新选择阈值。

## RNA与pseudobulk
原Seurat RNA单一counts layer，90011 features×123699 cells，388983646非零项，全部非负整数。未把完整表达矩阵转dense。已有RNA data在固定25个细胞的所有features上与log1p(counts/library×10000)一致，最大差8.88e-16，因此复用普通LogNormalize data做展示。先用稀疏指示矩阵聚合sample×原cluster的count sum、LogNormalize sum、detected cell sum，再按固定比较合成sample×group。未复制完整Seurat；原对象仅一次加载。

edgeR 4.8.2/limma 3.66.0，先保存并检查本地官方help。filterByExpr(y,group=group)默认min.count=10,min.total.count=15,large.n=10,min.prop=0.7；过滤按两大group，不用含sample阻断的高杠杆design放宽检测要求。过滤后重新计算文库大小、TMM归一化；estimateDisp(robust=TRUE)、glmQLFit(robust=TRUE,legacy=FALSE)、glmQLFTest。设计~sample+group，Right reference，检验groupLeft系数；满秩、文库非零、残余自由度>0均核对。dataset/source在每个sample内唯一，已被sample固定效应吸收，不再加共线主效应。各比较在全部通过过滤的基因内BH，完整结果/过滤记录保留。没有基于细胞P值的DE。即使全部模型运行成功，也不意味着统计独立性、混杂或分群后选择偏倚已经消除。

## marker审阅
完整基因ID不变；gene_symbol仅对已知SYMBOL或唯一ENSEMBL映射填写，歧义保留NA，不合并行。表达与检测率先按sample×group求均值，再对该比较合格样本等权。检测率即RNA counts>0的细胞比例，不等同蛋白阳性。短名单门槛为FDR<0.05、|log2FC|≥0.5、目标平均检测率≥0.10、样本内表达方向支持≥0.70；排序依次为dataset均值差方向支持、sample方向支持、朝目标方向的检测率差、|log2FC|、gene_id。每comparison至多20，C12 vs C16每方向至多10，其余只取左侧；这是一套透明的实用筛选规则，不是生物学验证阈值或综合验证分。

展示marker每state至多3个，在固定短名单与冻结基因诊断结果中人工复核局部对照及功能多样性，不更改DE或群。C8的ELOVL5来自冻结诊断/符合条件但未进入该比较前20的候选，其加入为展示脂质支持，而非改变短名单排序。输出其余原群中最接近/更高表达群，不新增全群两两DE。同分辨率所有非空sample×cluster（含小群）用于描述性特异性汇总，n_samples/n_cells均保留；因此此汇总不能替代严格配对局部对照。MP07-only/MP10-only/shared=38/38/12；逐基因表达与检测率，不重算任何UCell/删共享后分数。功能仅用org.Hs.eg.db 3.22.0及GO.db本地注释的有限关键词记录（含GO ID/term/evidence），没有通路富集网格。

## 图形与复现
冻结FULL_Z UMAP及原Seurat umap_allgene按cell_id取子集，不重拟合。点序seed42；精确零黑色并先画，正值用plasma；NA灰色。颜色范围跨dataset固定并单独记录；usage连续图为明确标注的log1p显示，统计仍是raw。RNA DotPlot为sample等权普通log-normalized RNA，不是usage；大小为平均检测率。QC counts/features展示log1p(mean)，percent.mt保留百分数；无缺失填0。PDF和300dpi透明PNG，每图保留实际数据和说明。GSE326225保留所有细胞/图及零资格记录，但六项比较无合格样本，不伪造配对或DotPlot。最多2线程，原Seurat读取/模型顺序执行。没有重新运行NMF/UCell/UMAP/聚类、空间/CNV/预后或其他任务。

## 参考
- [edgeR官方入口](https://bioconductor.org/packages/edgeR/)；执行API以00_audit中edgeR4.8.2官方help为准。
- [OSCA配对多样本比较](https://bioconductor.org/books/3.17/OSCA.multisample/multi-sample-comparisons.html)
- [生物学重复与单细胞DE](https://www.nature.com/articles/s41467-021-25960-2)
- [EGR1原始功能研究](https://pubmed.ncbi.nlm.nih.gov/8336701/)
