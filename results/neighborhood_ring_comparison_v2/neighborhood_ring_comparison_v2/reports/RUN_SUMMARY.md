# 本轮实际运行总结

已完成完整组织拓扑与共同特征：7 datasets、63切片、290,678组织节点、265,002有效组成spots、129,632共同anchors；完整6/12环108,345；参考20,000。全队列来源已知patient_uid 29，共同anchors覆盖27；未核实身份不算新患者。

区域全部为自动候选，未发现真实审定区域来源，未虚构复核人。已生成全63切片CNV/Epi/候选证据图、ROI完整表和优先队列。

## 六方案

| config | status | dimensions | n_pc | k | n_reference | stability_median | stability_min | stability_max | native_common_agreement | low_vote_fraction |
|---|---|---|---|---|---|---|---|---|---|---|
| R1_KM | BLIND_FIT_COMPLETE | 30 | 0 | 2 | 20000 | 0.9903 | 0.9873 | 0.9968 | 0.971 | 0.0195 |
| R1_LEIDEN | BLIND_FIT_COMPLETE | 30 | 18 | 45 | 20000 | 0.913 | 0.8902 | 0.9215 | 1.0 | 0.0506 |
| R2_KM | BLIND_FIT_COMPLETE | 30 | 0 | 2 | 20000 | 0.9898 | 0.9856 | 0.9913 | 0.9735 | 0.0177 |
| R2_LEIDEN | BLIND_FIT_COMPLETE | 30 | 17 | 43 | 20000 | 0.8931 | 0.8647 | 0.9005 | 1.0 | 0.0384 |
| R12_KM | BLIND_FIT_COMPLETE | 60 | 0 | 2 | 20000 | 0.9865 | 0.9838 | 0.99 | 0.972 | 0.0168 |
| R12_LEIDEN | BLIND_FIT_COMPLETE | 60 | 22 | 38 | 20000 | 0.8782 | 0.8666 | 0.9009 | 1.0 | 0.0305 |



## 直接回答科学问题

1. R1、R2、R12三种KM均推荐K=2，表现为上皮占主导环境中的相对较高/较低上皮组成两类（R1工作群的R1上皮组成约0.75 vs 0.52）。Leiden进一步区分CAF/免疫等次要组分，但所有工作群的平均首位成分仍是上皮。原始R1/R2 profile和逐群命名描述见environment_dictionary.tsv。不能把编号跨配置等同，也不宣称唯一赢家。

2. 同输入参考ARI如下；不同K时属于完整流程差异，不能单独归因PCA。K_L>20按规则不补同K，不为桥接扩展网格。

| config1 | config2 | ARI | k1 | k2 |
|---|---|---|---|---|
| R1_KM | R1_LEIDEN | 0.0362 | 2 | 45 |
| R2_KM | R2_LEIDEN | 0.0406 | 2 | 43 |
| R12_KM | R12_LEIDEN | 0.0433 | 2 | 38 |

3. KM各组均覆盖多个dataset及真实患者；Leiden的R1/R2/R12分别有31/45、28/43、27/38个群，超过90%的anchors来自同一dataset。虽部分群在其他来源也有少量spots，不能据此称普遍复现。每群真实患者/dataset/切片支持见cluster_patient_dataset_support.tsv；来源集中群保留并标记。拟合重采样稳定性不等于患者生物学复现。

4. 六方案中心MP10效应、真实患者方向一致性尚未检验：区域身份尚未审定，空间MP10保持盲态。没有MP10效应量/P值可报告。

5. 肿瘤支持身份、界面与Epi混杂关系尚待真实区域复核；当前深度无效NA，不能称病理界面。

6. 已保留完整环子集。六方案MP10完整环敏感性待解盲；本轮不提出近/远MP10更强结论，因此未触发R2降采样6邻居的条件敏感性。

7. 来源验证状态：PARTIAL_SOURCE_VALIDATION。已实际审计9600细胞/32类型/58来源样本，原Top20。一般混合模拟260 bags。MP10不一定是上皮专属；source_by_dataset_sample_type.tsv和gene_contribution_by_type.tsv展示其他谱系的得分/贡献。全参考RNA universe不直接等于空间量尺。恶性/正常独立标签、环境条件composition-only与真实高Epi敏感性仍未完成。

8. 文献复用及本项目设计/未核对Methods分开列在METHODS.md。

## 实际状态与未决事项

M00–M03完成；M04–M05以本表状态和frozen_manifest为准。M01候选材料完成、终稿待人工。M06门控未解盲。M07部分完成且不能删除。

- 最小人工输入：在ROI_review.tsv的对应ROI填写reviewer、真实日期、reviewed_label、证据comment；优先ROI_priority_review.tsv中的冲突/小巢。允许复核后保留uncertain，不自动将diploid判正常。

- 无核实patient_uid样本保留环境分析，但不参与真实患者独立重复统计；section_identity.tsv列明缺少的crosswalk。

- 原MP10基因可用性见MP10_gene_availability.tsv，不能跨缺基因集合直接合并分数。

- 完成M07关联归因还需环境条件composition-only、空间feature universe/深度匹配，以及有来源证据的恶性/正常单细胞身份；目前一般混合模拟不能替代。

## 续跑

```bash
bash '/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_comparison_v2/code/resume.sh'
```

收到真实复核表后：

```bash
bash '/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_comparison_v2/code/resume.sh' --review '/absolute/path/ROI_review_completed.tsv'
```

使用checkpoint，不从头重训；仅区域身份变更不重拟合环境标签。

## 实现修复记录

历史QC-only坐标遗漏25676个tissue==1节点，已实际重提取核验并重建；无效早期输出保留archive_qc_only_topology_INVALID，禁止用于结论。未改accepted结果。首次一般pseudo-spot池每来源上皮仅5–6个，未满足20细胞无放回需求，空模拟未视为完成；03b单独扩大有限来源池后才记录实际bags。

MP10 上皮来源验证未完成；现有证据只支持部分评分来源/混合敏感性审计，不支持“上皮特异”或因果诱导。
