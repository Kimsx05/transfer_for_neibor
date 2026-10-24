# CellCharter 两片联合探索 v1：待用户验收

本轮已完成 PRIMARY、EXPANDED_LOWINFO 和 NO_CURRENT_CONFLICT 的实际计算。只到空间簇发现、全谱系解释和两个敏感性版本；未计算 MP10/胆固醇、spottedpy、通讯、CNV 或正式 TLS。原标签不是人工确认，空间簇也不是机制性 niche 的证明。

## 输入与身份解释

全部保留308,549，输入对象、counts、分割、scVI和工作注释不修改。按cell_id对齐，实测latent20所有维度及空间坐标有限才有建模资格。三集合与给定预期完全一致：

| section_id | PRIMARY_included | EXPANDED_LOWINFO_included | NO_CURRENT_CONFLICT_included | review_strong_competition | morphology_identity_uncertain |
| --- | --- | --- | --- | --- | --- |
| Cervilla | 85154 | 95979 | 83207 | 1947 | 6 |
| D1_Johnson | 179803 | 206292 | 157793 | 22010 | 2 |

PRIMARY=264,957，EXPANDED=302,271，NO_CURRENT_CONFLICT=241,000；6,278缺有效latent者未补零。PRIMARY不是旧UMAP/core集合，noncore不一刀切排除。37,314个有latent的low_info仅进入expanded。三集合都保留Unknown；模型输入未使用细胞身份标签。

当前强竞争23,957（1,947/22,010）全部标为Disputed，不是只处理其中2,402个仍有已知大类的子集。此子集中的1,038个继承细类解释也回退Unknown。8个形态身份疑问单独标Disputed（7个PRIMARY），不再新增拟合版本。历史artifact/mixed标志原样保存，不断言污染或doublet。

## 构图、尺度、孤立对象

每片分别以真实spatial_um进行SciPy Delaunay（三角剖分等价于Squidpy generic Delaunay底层），显式按section_id分组；聚合sample_key=section_id。逐版本断言跨切片边为0，重复坐标检查通过。预先声明PRIMARY每片99百分位与40µm取小值，后续两版沿用同一物理上限。

| version | section_id | n | threshold_um | edges | removed_edges | edge_median_um | edge_p95_um | neighbors_median | isolated | components | largest_component | largest_components |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PRIMARY | Cervilla | 85154 | 40.0 | 250324 | 5101 | 12.925 | 25.762 | 6.0 | 71 | 128 | 67276 | [67276, 14194, 2893, 168, 51, 47, 42, 38, 38, 38] |
| PRIMARY | D1_Johnson | 179803 | 38.322 | 533972 | 5394 | 10.227 | 21.639 | 6.0 | 101 | 184 | 177570 | [177570, 959, 244, 129, 114, 113, 52, 52, 37, 37] |
| EXPANDED_LOWINFO | Cervilla | 95979 | 40.0 | 282685 | 5210 | 12.968 | 26.202 | 6.0 | 44 | 87 | 76119 | [76119, 14368, 4454, 339, 229, 117, 109, 50, 21, 18] |
| EXPANDED_LOWINFO | D1_Johnson | 206292 | 38.322 | 615383 | 3452 | 10.229 | 21.663 | 6.0 | 33 | 65 | 204584 | [204584, 962, 305, 113, 60, 39, 34, 29, 19, 15] |
| NO_CURRENT_CONFLICT | Cervilla | 83207 | 40.0 | 244257 | 5327 | 12.98 | 26.327 | 6.0 | 71 | 128 | 65462 | [65462, 14065, 2893, 168, 47, 47, 42, 38, 38, 38] |
| NO_CURRENT_CONFLICT | D1_Johnson | 157793 | 38.322 | 465766 | 7570 | 10.181 | 24.71 | 6.0 | 101 | 184 | 155614 | [155614, 905, 244, 129, 114, 113, 52, 52, 37, 37] |

没有为了连接碎片而补边。主分析的最长100条保留边图已审阅：大片空白与分离组织块未被长连接跨越；距离裁剪无法保证识别所有小裂隙，仍是局限。Cervilla有多个大组织分量，D1主体更连续。孤立点保留拟合并标 `isolated_self_expression_only`：其邻居层无邻居，不能把分类称为完整邻域证据。

CellCharter官方aggregate_neighbors：n_layers=3、mean，中心20维+第1/2/3跳各20维，实际80维。以下为固定seed123每片300个中心的抽样，3-hop不是固定微米半径，也不是普通Visium ring：

| version | section_id | hop | sampled_cells | median_neighbors | p95_neighbors | median_max_um | p95_max_um |
| --- | --- | --- | --- | --- | --- | --- | --- |
| EXPANDED_LOWINFO | Cervilla | 1 | 300 | 6.0 | 8.0 | 18.51 | 34.72 |
| EXPANDED_LOWINFO | Cervilla | 2 | 300 | 13.0 | 16.0 | 36.61 | 58.22 |
| EXPANDED_LOWINFO | Cervilla | 3 | 300 | 20.0 | 26.0 | 53.89 | 82.29 |
| EXPANDED_LOWINFO | D1_Johnson | 1 | 300 | 6.0 | 8.0 | 14.85 | 27.03 |
| EXPANDED_LOWINFO | D1_Johnson | 2 | 300 | 13.0 | 17.0 | 29.43 | 47.56 |
| EXPANDED_LOWINFO | D1_Johnson | 3 | 300 | 21.0 | 26.0 | 45.0 | 67.64 |
| NO_CURRENT_CONFLICT | Cervilla | 1 | 300 | 6.0 | 8.0 | 18.32 | 33.71 |
| NO_CURRENT_CONFLICT | Cervilla | 2 | 300 | 12.0 | 16.0 | 36.12 | 57.86 |
| NO_CURRENT_CONFLICT | Cervilla | 3 | 300 | 20.0 | 26.0 | 55.57 | 78.9 |
| NO_CURRENT_CONFLICT | D1_Johnson | 1 | 300 | 6.0 | 8.0 | 14.82 | 32.95 |
| NO_CURRENT_CONFLICT | D1_Johnson | 2 | 300 | 13.0 | 17.0 | 30.45 | 56.23 |
| NO_CURRENT_CONFLICT | D1_Johnson | 3 | 300 | 20.0 | 27.0 | 48.2 | 78.45 |
| PRIMARY | Cervilla | 1 | 300 | 6.0 | 7.05 | 18.05 | 33.49 |
| PRIMARY | Cervilla | 2 | 300 | 13.0 | 17.0 | 35.35 | 55.91 |
| PRIMARY | Cervilla | 3 | 300 | 21.0 | 25.0 | 54.51 | 78.64 |
| PRIMARY | D1_Johnson | 1 | 300 | 6.0 | 8.0 | 14.61 | 26.0 |
| PRIMARY | D1_Johnson | 2 | 300 | 13.0 | 17.0 | 28.57 | 48.44 |
| PRIMARY | D1_Johnson | 3 | 300 | 21.0 | 26.0 | 44.48 | 68.55 |

## K选择与计算方法

CellCharter0.3.8、Squidpy1.8.3；实际接口已核对。GPU复用现有CUDA PyTorch，8 CPU线程、GPU进程上限40%，batch16384，full covariance GMM，协方差正则1e-5，最多100epochs。AutoK基础seed123；选定PRIMARY模型的实际seed1514383052，两项敏感性均固定为同一实际seed。参数完整见 `metadata/parameters.json`。AutoK初始K2–12、max_runs5、官方稳定性收敛阈值0.01。该版本为评估相邻K稳定性还实际拟合边界邻居K1/13，非额外生物学候选。

| stage | K | stability | sd | runs | min_cluster_n | min_cluster_fraction | cluster_sizes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| initial | 2 | 0.7026 | 0.0111 | 5 | 126562 | 0.4777 | [138395, 126562] |
| initial | 3 | 0.7321 | 0.0373 | 5 | 75995 | 0.2868 | [75995, 101719, 87243] |
| initial | 4 | 0.7638 | 0.0601 | 5 | 52967 | 0.1999 | [52967, 72882, 60212, 78896] |
| initial | 5 | 0.7361 | 0.0794 | 5 | 25092 | 0.0947 | [25092, 70441, 63273, 46437, 59714] |
| initial | 6 | 0.7272 | 0.0908 | 5 | 24058 | 0.0908 | [59160, 38142, 31065, 67525, 24058, 45007] |
| initial | 7 | 0.7337 | 0.0946 | 5 | 23511 | 0.0887 | [43498, 31479, 23511, 28930, 30577, 57201, 49761] |
| initial | 8 | 0.7369 | 0.0908 | 5 | 20562 | 0.0776 | [25106, 21089, 40504, 20562, 25141, 48393, 45217, 38945] |
| initial | 9 | 0.7173 | 0.0851 | 5 | 18601 | 0.0702 | [20608, 18601, 39725, 20356, 21638, 47269, 21833, 35242, 39685] |
| initial | 10 | 0.7227 | 0.059 | 5 | 9409 | 0.0355 | [20992, 35393, 27553, 20989, 17968, 39992, 23704, 29579, 39378, 9409] |
| initial | 11 | 0.7773 | 0.0501 | 5 | 10303 | 0.0389 | [12972, 16256, 34383, 21137, 18581, 39915, 17576, 34669, 38821, 10303, 20344] |
| initial | 12 | 0.772 | 0.0619 | 5 | 8083 | 0.0305 | [8083, 19324, 16923, 8214, 38133, 18950, 20757, 26507, 25166, 31012, 38545, 13343] |

主K=**11**。初始内部峰=[4, 8, 11]；是否仅扩展一次到16：False。优先内部稳定性峰，并检查最小簇≥0.5%以排查退化小簇，再比较官方稳定性；没有用身份标签、预期阳性或MP10选K。来源/深度/组织分布在下表逐簇解释；不是以完全混合为验收标准。如多个K接近，不把本K称为唯一真实数量。实际重复次数见runs列，独立运行种子见K_selection.json。保存的是选定K的最佳似然模型 `models/primary_gmm.pkl`，不是所有候选模型。

[AutoK曲线](figures/autok_stability.png)；[主分析两片并列](figures/PRIMARY_two_sections.png)。两片均出现的簇：['niche0', 'niche1', 'niche2', 'niche3', 'niche4', 'niche5', 'niche6', 'niche7', 'niche8', 'niche9']；只在一片有成员：['niche10']。出现不等于同等比例或同一机制。尤其Cervilla的niche7仅3个，不能据此宣称两片共享的可靠群体；niche10仅见D1。K11稳定性约0.777，K4约0.764，差距有限，当前K是探索性分辨率而非唯一最优真值。

## 每片每簇解释（分母包括Unknown/Disputed）

| section | niche | n | section_pct | interpretable_pct | composition_all | representative_measured_genes | median_genes | Unknown_pct | Disputed_pct | review |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Cervilla | niche0 | 790 | 0.93 | 67.72 | Fibroblast 46.3%; Unknown 19.4%; Disputed 12.9%; Epithelial 11.8%; Myeloid 5.3% | COL3A1(91%,Δ3.82); COL1A1(83%,Δ3.60); TAGLN(69%,Δ3.52); TIMP1(67%,Δ3.48); COL1A2(79%,Δ3.44); SPARC(71%,Δ3.20) | 326.0 | 19.37 | 12.91 | 组成/实测marker/位置共同审阅，暂不赋机制性名称 |
| Cervilla | niche1 | 15871 | 18.64 | 92.0 | Epithelial 89.0%; Unknown 6.3%; Fibroblast 1.9%; Disputed 1.7%; Myeloid 1.0% | S100P(83%,Δ1.51); ID1(65%,Δ1.48); PSCA(70%,Δ1.44); HPGD(49%,Δ1.39); WFDC2(44%,Δ1.22); FOS(41%,Δ1.20) | 517.0 | 6.26 | 1.73 | 主要由本来源贡献，技术/生物差异混杂 |
| Cervilla | niche2 | 14797 | 17.38 | 99.83 | Epithelial 99.8%; Unknown 0.2%; Fibroblast 0.1%; Disputed 0.0%; Myeloid 0.0% | RAB3IP(94%,Δ3.09); MDM2(94%,Δ3.06); CPSF6(85%,Δ2.73); ADIRF(95%,Δ2.67); WSB2(86%,Δ2.63); PEBP1(92%,Δ2.62) | 693.0 | 0.16 | 0.01 | 组成/实测marker/位置共同审阅，暂不赋机制性名称 |
| Cervilla | niche3 | 5562 | 6.53 | 45.09 | Unknown 53.9%; Fibroblast 15.0%; Myeloid 14.7%; Plasma 10.5%; Epithelial 2.5% | IGKC(93%,Δ8.06); IGHG1(71%,Δ5.25); IGHA1(52%,Δ3.73); JCHAIN(39%,Δ2.48); IGHG3(28%,Δ1.92); CD74(36%,Δ1.48) | 124.0 | 53.9 | 1.01 | 未知身份占比较高 |
| Cervilla | niche4 | 5134 | 6.03 | 57.36 | Unknown 38.3%; Fibroblast 29.2%; Myeloid 13.5%; Endothelial 6.9%; Disputed 4.4% | COL1A1(66%,Δ2.50); COL3A1(71%,Δ2.46); COL1A2(59%,Δ2.29); VIM(47%,Δ2.19); IGKC(43%,Δ2.03); IGFBP7(41%,Δ1.87) | 185.0 | 38.27 | 4.36 | 组成/实测marker/位置共同审阅，暂不赋机制性名称 |
| Cervilla | niche5 | 596 | 0.7 | 33.56 | Unknown 66.4%; T/NK 10.4%; Myeloid 9.6%; Fibroblast 9.2%; Epithelial 1.5% | IGKC(46%,Δ2.54); TMSB4X(70%,Δ1.20); CD74(29%,Δ1.09); TMSB10(27%,Δ0.87); PTGDS(16%,Δ0.86); CORO1A(13%,Δ0.76) | 87.0 | 66.44 | 0.0 | 未知身份占比较高；深度偏低 |
| Cervilla | niche6 | 16632 | 19.53 | 79.12 | Epithelial 66.0%; Unknown 13.6%; Fibroblast 9.8%; Disputed 7.2%; Myeloid 1.4% | COL1A1(75%,Δ2.59); COL3A1(80%,Δ2.58); COL1A2(64%,Δ2.12); MDM2(78%,Δ1.86); SPARC(51%,Δ1.75); RAB3IP(76%,Δ1.71) | 559.0 | 13.65 | 7.23 | 组成/实测marker/位置共同审阅，暂不赋机制性名称 |
| Cervilla | niche7 | 3 | 0.0 | 33.33 | Disputed 33.3%; Epithelial 33.3%; Unknown 33.3% | MPEG1(67%,Δ3.57); ANXA5(67%,Δ3.41); SMAD2(67%,Δ3.37); MDM2(100%,Δ3.34); DRG2(67%,Δ3.30); LGALS1(67%,Δ3.29) | 366.0 | 33.33 | 33.33 | 争议身份占比较高；成员极少，不作独立生物群解释 |
| Cervilla | niche8 | 24863 | 29.2 | 82.77 | Epithelial 81.3%; Unknown 17.2%; Fibroblast 1.2%; Mural 0.2%; Myeloid 0.0% | KRT13(63%,Δ0.30); KRT19(34%,Δ0.28); MYCL(16%,Δ0.06); S100A2(19%,Δ0.03); FGFR3(31%,Δ0.02); KRT17(31%,Δ0.01) | 78.0 | 17.21 | 0.02 | 深度偏低 |
| Cervilla | niche9 | 906 | 1.06 | 60.6 | Unknown 30.1%; Epithelial 27.6%; Fibroblast 18.0%; Disputed 9.3%; Myeloid 8.7% | IGKC(88%,Δ6.21); IGHG1(72%,Δ4.41); IGHG3(49%,Δ3.12); TIMP1(65%,Δ2.85); COL1A1(74%,Δ2.44); COL3A1(79%,Δ2.40) | 299.0 | 30.13 | 9.27 | 组成/实测marker/位置共同审阅，暂不赋机制性名称 |
| D1_Johnson | niche0 | 12182 | 6.78 | 63.8 | Fibroblast 46.7%; Unknown 21.3%; Disputed 14.9%; Epithelial 10.8%; Myeloid 5.1% | COL1A1(91%,Δ3.08); COL1A2(84%,Δ3.02); SPARC(74%,Δ2.44); COL3A1(87%,Δ2.36); LUM(60%,Δ1.95); TAGLN(47%,Δ1.71) | 373.0 | 21.33 | 14.87 | 组成/实测marker/位置共同审阅，暂不赋机制性名称 |
| D1_Johnson | niche1 | 385 | 0.21 | 89.35 | Epithelial 88.3%; Unknown 9.4%; Disputed 1.3%; Fibroblast 0.5%; Myeloid 0.5% | SPINK1(74%,Δ4.15); UPK1A(62%,Δ3.09); ADIRF(65%,Δ3.05); HPGD(60%,Δ2.88); KRT7(59%,Δ2.52); BHMT(48%,Δ2.47) | 438.0 | 9.35 | 1.3 | 组成/实测marker/位置共同审阅，暂不赋机制性名称 |
| D1_Johnson | niche10 | 20344 | 11.31 | 68.1 | Epithelial 57.4%; Disputed 18.9%; Unknown 13.0%; Myeloid 8.8%; Fibroblast 1.5% | S100A9(95%,Δ4.17); S100A8(90%,Δ3.97); KRT16(76%,Δ2.80); KRT5(70%,Δ2.45); IFI6(73%,Δ2.38); KRT17(76%,Δ2.38) | 782.0 | 13.01 | 18.89 | 争议身份占比较高；主要由本来源贡献，技术/生物差异混杂 |
| D1_Johnson | niche2 | 19586 | 10.89 | 92.72 | Epithelial 89.7%; Unknown 6.1%; Myeloid 2.4%; Disputed 1.2%; Fibroblast 0.3% | S100A9(95%,Δ4.20); S100A8(86%,Δ3.29); S100P(66%,Δ2.68); SLC44A2(73%,Δ2.65); SMARCA4(72%,Δ2.57); HSPB1(87%,Δ2.53) | 688.5 | 6.07 | 1.21 | 组成/实测marker/位置共同审阅，暂不赋机制性名称 |
| D1_Johnson | niche3 | 15575 | 8.66 | 43.78 | Unknown 54.1%; Fibroblast 26.2%; Myeloid 8.4%; Plasma 3.7%; Disputed 2.2% | IGKC(78%,Δ4.99); IGHG1(64%,Δ4.15); IGHA1(36%,Δ2.35); IGHG3(29%,Δ1.59); IGLC1(22%,Δ1.36); JCHAIN(20%,Δ1.10) | 120.0 | 54.05 | 2.16 | 未知身份占比较高 |
| D1_Johnson | niche4 | 13447 | 7.48 | 49.78 | Unknown 42.5%; Fibroblast 29.1%; Disputed 7.8%; Endothelial 7.1%; Mural 5.3% | SPARC(57%,Δ1.60); COL4A1(31%,Δ1.39); VIM(59%,Δ1.38); COL3A1(68%,Δ1.30); A2M(31%,Δ1.20); COL1A1(68%,Δ1.16) | 173.0 | 42.46 | 7.76 | 未知身份占比较高 |
| D1_Johnson | niche5 | 39319 | 21.87 | 40.77 | Unknown 57.8%; Fibroblast 19.0%; T/NK 9.3%; Myeloid 7.4%; B 3.6% | IGHM(16%,Δ0.76); TRAC(18%,Δ0.68); CXCR4(15%,Δ0.66); PTGDS(17%,Δ0.66); IGKC(38%,Δ0.59); TRBC2(14%,Δ0.58) | 108.0 | 57.78 | 1.45 | 未知身份占比较高；主要由本来源贡献，技术/生物差异混杂 |
| D1_Johnson | niche6 | 944 | 0.53 | 50.0 | Unknown 46.9%; Epithelial 28.2%; Fibroblast 11.7%; Myeloid 8.7%; Disputed 3.1% | HBA2(13%,Δ0.83); KRT5(36%,Δ0.35); KRT14(35%,Δ0.34); COL1A1(60%,Δ0.25); S100A2(29%,Δ0.24); COL3A1(58%,Δ0.23) | 127.0 | 46.93 | 3.07 | 未知身份占比较高 |
| D1_Johnson | niche7 | 34666 | 19.28 | 55.6 | Epithelial 40.7%; Disputed 32.1%; Fibroblast 12.5%; Unknown 12.3%; Myeloid 2.3% | KRT14(84%,Δ3.64); KRT5(84%,Δ3.57); KRT6A(82%,Δ3.45); KRT17(84%,Δ3.29); S100A2(79%,Δ3.17); S100A10(82%,Δ3.08) | 701.0 | 12.27 | 32.12 | 争议身份占比较高；主要由本来源贡献，技术/生物差异混杂 |
| D1_Johnson | niche8 | 13958 | 7.76 | 39.25 | Unknown 60.6%; Epithelial 31.2%; Fibroblast 5.7%; Myeloid 1.3%; T/NK 0.4% |  | 52.0 | 60.63 | 0.11 | 未知身份占比较高；深度偏低；没有满足当前检出率/正表达差规则的代表基因，不等于无表达 |
| D1_Johnson | niche9 | 9397 | 5.23 | 51.67 | Fibroblast 33.5%; Disputed 31.8%; Unknown 16.6%; Epithelial 10.9%; Myeloid 6.3% | IGKC(83%,Δ4.25); IGHG1(74%,Δ3.82); IGHG3(59%,Δ3.45); COL1A1(96%,Δ3.23); COL3A1(95%,Δ3.04); COL1A2(89%,Δ2.97) | 462.0 | 16.58 | 31.75 | 争议身份占比较高 |

便于审阅的生物学描述如下，所有身份都属于工作解释：

| niche | Cervilla | D1_Johnson |
| --- | --- | --- |
| 0 | 成纤维占46.3%，COL3A1/COL1A1/TAGLN；小群790 | 成纤维占46.7%，COL1A1/COL1A2/LUM/TAGLN；仍含21.3% Unknown |
| 1 | 上皮占89.0%，S100P/PSCA/HPGD；主要来自本片 | 仅385个，上皮占88.3%，SPINK1/UPK1A/KRT7；敏感性较弱 |
| 2 | 上皮占99.8%，RAB3IP/MDM2/CPSF6 | 上皮占89.7%，S100A8/S100A9/S100P；与C片同编号不代表完全相同表达状态 |
| 3 | Ig/JCHAIN高，Unknown53.9%，并有成纤维/髓系/浆细胞 | Ig/JCHAIN高，Unknown54.1%，成纤维26.2%；两片均不宜命名纯浆细胞区或TLS |
| 4 | COL1A1/COL3A1/COL1A2，基质相关混合区；Unknown38.3% | SPARC/COL4A1/VIM，基质与血管相关成分；Unknown42.5% |
| 5 | 仅596个，Unknown66.4%，中位87genes，免疫marker伴低深度 | 39,319个，Unknown57.8%，TRAC/TRBC2/IGHM支持免疫混合成分；来源偏倚明显 |
| 6 | 上皮66.0%、成纤维9.8%，MDM2/RAB3IP与胶原共同支持混合邻域 | 仅944个，Unknown46.9%，KRT5/KRT14与胶原；不稳定小群 |
| 7 | 仅3个，不作生物群解释 | 上皮40.7%，KRT14/KRT5/KRT6A/KRT17；Disputed32.1%，不能视为纯上皮区 |
| 8 | 上皮81.3%，KRT13/KRT19；中位78genes，空间外围广泛分布 | Unknown60.6%，中位52genes；没有满足筛选规则的正向代表基因，深度解释优先 |
| 9 | Ig及胶原共同表达，含上皮/成纤维/未知身份 | 成纤维33.5%、Disputed31.8%，Ig及胶原；争议身份影响较大 |
| 10 | 无成员 | 上皮57.4%，S100A8/S100A9/KRT16/KRT5；Disputed18.9%，D1特有且敏感 |

空间图显示Cervilla上皮相关niche1/2/6占据不同连续区域，niche8多分布外围；D1的niche2主要位于图中下部，niche5偏上部/右部，niche7/10在主体中部互相邻接。上述位置是聚类图观察，尚无病理分区证实，不能据此确定侵袭边界或机制。niche8具有空间连续性，也不能排除深度驱动。全体逐片逐簇的数值证据见上表和marker文件。

本轮统一称niche0…niche10，表示中心及周围表达共同定义的空间簇。依据组成、实测marker和空间位置作描述，未强行命名TLS/GC/机制。全谱系均纳入：major_composition.tsv含全体建模分母与可解释大类分母两列；同时有interpretable_coverage，避免隐藏Unknown/Disputed。subtype_parent_composition.tsv明确沿用细类或父类，父类不冒充细类。

marker来自未改动v2整数counts，合并同名基因、log1p10k，逐片相对于其他PRIMARY空间簇计算平均log表达差和检出率。未使用插补、P值筛选或独立验证措辞；仅两张切片，细胞不等于患者重复。`representative_markers.tsv`有限代表基因，`measured_marker_evidence.parquet`同时保留全谱系已知marker。原feature缺失（如TPSAB1/TPSB2、KRT18等）不同于检测为零，缺失明细沿用上轮 `../joint_umap_major_review/metadata/review_panels.json`。

每片质量图分别显示Unknown/Disputed、n_genes、n_counts与强竞争，同原物理坐标；深度采用对数颜色。以上review列指出深度、争议身份、未知覆盖或来源偏倚，不自动删掉这些簇。

## 固定K的两个敏感性版本

两个版本均真正重新Delaunay构图、按主分析物理阈值裁剪、3层聚合、固定K重新拟合，实际seed1514383052与选定PRIMARY模型一致。早期seed123敏感性试算已被此统一种子的结果替换；以下为最终数值。相同种子约束初始化随机性，但不同样本集合仍会改变初始化与优化路径，不能把差异全归于生物学。不重新AutoK，也不是只投射新增细胞。

| version | section_id | common_cells | ARI | mapped_agreement |
| --- | --- | --- | --- | --- |
| EXPANDED_LOWINFO | joint | 264957 | 0.7455 | 0.8091 |
| EXPANDED_LOWINFO | Cervilla | 85154 | 0.6289 | 0.6964 |
| EXPANDED_LOWINFO | D1_Johnson | 179803 | 0.7823 | 0.8625 |
| NO_CURRENT_CONFLICT | joint | 241000 | 0.6949 | 0.7685 |
| NO_CURRENT_CONFLICT | Cervilla | 83207 | 0.5661 | 0.6214 |
| NO_CURRENT_CONFLICT | D1_Johnson | 157793 | 0.7615 | 0.8461 |

ARI及交叉表只在每版与PRIMARY共有ID上计算，排除者不当不一致。显示编号在两片联合交叉表上用Hungarian匹配一次，映射文件保留，每片不独立改名；raw_cluster仍保留。NO_CURRENT_CONFLICT不是无污染真值。

两版比较中的具体结构变化（百分比均以共有ID为分母）：

- 较稳定：两片niche3的Jaccard分别为0.909/0.837（expanded）与0.865/0.860（去冲突）；D1 niche5为0.893/0.911。Cervilla niche8为0.892/0.948，但稳定性不消除其深度偏倚。D1 niche2两版分别有99.3%/99.4%成员保留在对应簇。
- 纳入low_info后，Cervilla niche1保留63.3%，另有21.8%进入显示niche10、14.2%进入niche2；niche6的65.8%进入niche1、23.4%进入niche0，仅10.1%保持匹配niche6，显示混合邻域边界重组。D1 niche6有92.6%保留，但该匹配簇共同成员从944增至12,880，Jaccard仅0.067，属于显著扩张而非稳定边界。
- 去冲突后，Cervilla niche2有97.0%共同成员进入显示为niche6的新簇；旧niche1的50.5%及旧niche6的61.0%也进入该簇。旧niche2内部大部分保持共同归属，但与其他上皮/混合邻域的区分消失，属于合并，不能仅凭对角Jaccard接近0就说全部结构消失。D1 niche6有96.3%进入niche0，同样是小群并入更大区域。
- D1 niche10两版中保留到同一匹配簇的比例为65.1%/59.4%，Jaccard0.648/0.462；D1 niche7为89.1%/78.8%，Jaccard0.801/0.683。niche7比niche10更稳定，但高Disputed比例仍限制其身份解释。
- Cervilla niche7原本只有3个，去冲突后共有2个；不对这些极小分母的稳定性作生物解释。Cervilla niche0/4也会合并或吸收其他成员，需结合组成变化审阅。

[逐簇转移图](figures/sensitivity_transitions.png)同时显示非对角的大块归属；最大去向比例只说明共同归属，不能代替Jaccard证明簇边界稳定。变化最明显的簇（完整逐簇表见per_niche_sensitivity.tsv）：

| version | section_id | niche | primary_common_n | sensitivity_common_n | retained_same_mapped_n | retention | Jaccard | largest_destination | largest_destination_fraction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NO_CURRENT_CONFLICT | D1_Johnson | niche6 | 915 | 4 | 0 | 0.0 | 0.0 | niche0 | 0.9628 |
| NO_CURRENT_CONFLICT | Cervilla | niche2 | 14795 | 6 | 6 | 0.0004 | 0.0004 | niche6 | 0.9701 |
| NO_CURRENT_CONFLICT | Cervilla | niche0 | 688 | 2178 | 152 | 0.2209 | 0.056 | niche4 | 0.7485 |
| EXPANDED_LOWINFO | D1_Johnson | niche6 | 944 | 12880 | 874 | 0.9258 | 0.0675 | niche6 | 0.9258 |
| EXPANDED_LOWINFO | Cervilla | niche6 | 16632 | 3590 | 1680 | 0.101 | 0.0906 | niche1 | 0.6579 |
| NO_CURRENT_CONFLICT | Cervilla | niche7 | 2 | 8 | 1 | 0.5 | 0.1111 | niche4 | 0.5 |
| EXPANDED_LOWINFO | Cervilla | niche0 | 790 | 4855 | 696 | 0.881 | 0.1406 | niche0 | 0.881 |
| EXPANDED_LOWINFO | Cervilla | niche7 | 3 | 3 | 1 | 0.3333 | 0.2 | niche0 | 0.3333 |
| NO_CURRENT_CONFLICT | D1_Johnson | niche10 | 16502 | 14546 | 9806 | 0.5942 | 0.4616 | niche10 | 0.5942 |
| NO_CURRENT_CONFLICT | D1_Johnson | niche1 | 380 | 335 | 254 | 0.6684 | 0.551 | niche1 | 0.6684 |
| EXPANDED_LOWINFO | D1_Johnson | niche10 | 20344 | 13343 | 13243 | 0.651 | 0.6478 | niche10 | 0.651 |
| EXPANDED_LOWINFO | D1_Johnson | niche8 | 13958 | 10008 | 9739 | 0.6977 | 0.6845 | niche8 | 0.6977 |

[纳入low_info空间图](figures/EXPANDED_LOWINFO_two_sections.png)；[去当前强竞争空间图](figures/NO_CURRENT_CONFLICT_two_sections.png)。请结合逐簇retention/Jaccard、两版组成/规模和hop范围比较，整体ARI不能掩盖小群变化。低信息增点会改变邻接和聚合；去强竞争也会改变局部密度，因此改变不等于某版一定错误或正确。结果不稳定处保留待讨论，没有新增筛选版本。

## 是否进入下游

目前值得把这些结果作为下一步研究设计和核对材料，但**尚不宜直接作为确认性MP10/共定位结论的固定niche**。先验收空间结构及敏感性较稳定的簇，明确是否把Disputed/孤立点仅用于建模而限制解释；对高Unknown、强来源或逐簇Jaccard偏低的结构尤其保留不确定性。两片无独立患者重复，marker是同数据描述性支持。B/Plasma聚集不自动等于TLS，GC B未标出不证明组织无GC。任何后续评分/共定位均等待用户明确授权。

## 文件、复现和检查

主入口 `cellcharter_assignments.parquet` 覆盖全部308,549个ID，各版未纳入者niche为空且有排除原因；工作标签、解释标签、历史/current质量字段均保留。`code/prepare.py → code/fit.py → code/summarize.py → code/refine_marker_plots.py → code/review_sensitivity.py → code/report.py`；每步用 `code/run.sh code/脚本名.py`（以上路径相对于本结果目录）执行；运行科学流程用 `code/run.sh` 复用旧GPU环境和本轮隔离软件依赖。主聚合表示只保存一份 `metadata/primary_aggregate.npy`（约85MB），其余邻接和聚合只在内存。图件为透明PNG、真实µm尺度；没有新h5ad或全体两两距离。

必要交付检查11项通过（delivery_checks.tsv）。计数、ID、纳入掩码、原工作标签/空间坐标、完整强竞争排除、解释层级及每版标签空值核对通过；运行中验证所有版本跨片边为0和80维聚合有限。

官方依据：[教程](https://cellcharter.readthedocs.io/en/stable/notebooks/cosmx_human_nsclc.html)、[聚合API](https://cellcharter.readthedocs.io/en/stable/generated/cellcharter.gr.aggregate_neighbors.html)、[AutoK API](https://cellcharter.readthedocs.io/en/stable/generated/cellcharter.tl.ClusterAutoK.html)。安装环境版本可能较网页实现更新，实际参数以本地代码和日志为准。

本轮已停在CellCharter空间簇及敏感性结果验收，等待用户审阅。
