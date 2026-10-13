# v5 实际完成结果

全部为结果知情后的探索/校准。v4、niche、Top20及真实模型保持只读。

完整上皮：123,699细胞匹配原存储raw usage，固定90011基因背景；总体rho=0.5756，来源等权中位rho=0.5752。49个来源可计算，8个不足/常量；这些来源并非已核实患者。

M03A完成1998个缓存bags回放。M03B完成540个控制bags、5个reference重叠来源、3种背景、3个深度及24个模型匹配基线。M03C完成270bags，状态CONVERGED_OBJECTIVE_RULE；两次新拟合尝试仅一个可用模型。

## 三个同切片直接比较

|对比|层|单位|均值差|探索95%区间|
|---|---|---:|---:|---|
|niche1_minus_niche2|spot_MP10_model_universe|30|+0.00809|[-0.00134,+0.01595]|
|niche1_minus_niche6|spot_MP10_model_universe|16|+0.01310|[+0.00591,+0.02025]|
|niche7_minus_niche6|spot_MP10_model_universe|15|+0.00182|[-0.00681,+0.01103]|
|niche1_minus_niche2|epi_MP10_model_universe|30|+0.00679|[-0.00078,+0.01268]|
|niche1_minus_niche6|epi_MP10_model_universe|16|+0.00782|[+0.00328,+0.01285]|
|niche7_minus_niche6|epi_MP10_model_universe|15|+0.00008|[-0.00704,+0.00763]|

## 每dataset效应

值为observed/model Epi，括号为有效单位；不足5单位只描述。

|dataset|1−2|1−6|7−6|
|---|---|---|---|
|CanDis__UNVERIFIED|-0.0025/-0.0036 (7)|+0.0118/+0.0037 (3)|+0.0096/+0.0033 (3)|
|GSE171351|-0.0019/+0.0012 (2)|+0.0097/+0.0130 (2)|NA：无合格共同覆盖|
|GSE238145|+0.0029/-0.0016 (4)|+0.0052/+0.0047 (3)|-0.0025/-0.0004 (4)|
|GSE285715|+0.0091/+0.0117 (3)|+0.0280/+0.0249 (1)|NA：无合格共同覆盖|
|GSE319536|+0.0191/+0.0130 (6)|+0.0228/+0.0073 (4)|+0.0007/-0.0028 (5)|
|Inhousedata__UNVERIFIED|-0.0045/+0.0095 (3)|+0.0043/+0.0023 (2)|-0.0414/-0.0324 (1)|
|NPJ__UNVERIFIED|+0.0247/+0.0182 (5)|+0.0117/+0.0153 (1)|+0.0232/+0.0197 (2)|

## 恢复校准

|分支|背景|状态对|平均真实Δ|平均估计Δ|rank筛查反转率|
|---|---|---:|---:|---:|---:|
|FROZEN_WEIGHT_REPLAY|actual_v4_center|99|0.0779|0.0507|19.2%|
|KNOWN_MIXTURE_FIXED_SIGNATURE|GC B|30|0.0489|0.0403|6.7%|
|KNOWN_MIXTURE_FIXED_SIGNATURE|Neutrophil|30|0.0549|0.0380|6.7%|
|KNOWN_MIXTURE_FIXED_SIGNATURE|SPP1 TAM|30|0.0499|0.0430|6.7%|
|SYNTHETIC_REFIT|GC B|15|0.0279|0.0302|13.3%|
|SYNTHETIC_REFIT|Neutrophil|15|0.0348|0.0289|13.3%|
|SYNTHETIC_REFIT|SPP1 TAM|15|0.0456|0.0397|6.7%|

固定真上皮状态的组成偏差见 composition_bias_summary.tsv；按背景/深度分层的恢复误差见 calibration_summary_by_background.tsv。两者不同于真实niche效应。

优先保留niche1−niche6作为下一步受限验证对象；niche1−niche2次之；niche7−niche6当前覆盖/方向/删一稳健性较弱。来源校准限制所有机制解释。

完整九群与三个pair见 ../05_evidence_synthesis/evidence_matrix_all9.tsv；已知/匿名单位数与七dataset完整结果见 paired_dataset_effects.tsv。

续跑入口 code/resume_mp10_calibration_v5.sh；--verify-only 仅核验，不重新训练。

合成端到端模型虽满足预设目标收敛规则，但真实上皮RNA占比平均49.165%，估计平均17.946%，存在明显低估。这是本合成设计中的来源恢复限制；不能外推为全部真实切片同幅度偏差，也不能把目标收敛或分数相关性称为表达归属验证通过。
