# 续跑与核验

环境和代码均位于本项目，无需修改旧成果。任意目录执行：

```bash
"/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_spottedpy_v8/code/run_spottedpy_v8.sh" --resume
```

无参数同样按指纹复用 checkpoint，补齐缺失阶段。输入、spec或锁定源码/依赖身份变化会拒绝复用，不自动覆盖旧输入。

仅核验（不会重新进行置换）：

```bash
"/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_spottedpy_v8/code/run_spottedpy_v8.sh" --verify-only
```

运行入口为code/run.py；环境env/bin/python；每次日志保留于logs。--resume需要本服务器完整数据，review包用于离线审阅，不包含表达矩阵或依赖环境。

点级表分布在02_hotspots/<section>/<variable>_k<k>.parquet，包含UID、原始值、G、Z、native P、双侧P、local q、原始/正式hot/cold和组件ID；逻辑清单为02_hotspots/point_hotspot_statistics_manifest.tsv。*_components.parquet保留小区域。计算状态JSON包含输入/spec/锁定源与依赖指纹。

不以修改阈值、减少9999置换、填最大距离、跨组织块或自动挑选显著组合来补结果。完成状态为V8_COMPLETE_AWAITING_REVIEW。
