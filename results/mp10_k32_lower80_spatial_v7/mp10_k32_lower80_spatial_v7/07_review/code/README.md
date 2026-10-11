# 执行与恢复

所有命令在项目根目录执行，Python入口为本目录的python.sh。00_prepare.py完成范围/标签与NB输入；01_reference.py训练唯一新NB；02_inputs.py和03_pseudospots.py可与NB并行，分别固定空间/展示基因与DEV bag。05_run.py等待reference冻结后管理默认三个映射，最多2个GPU进程，不管理其他任务。不要重复启动同一个模型或重复调度器。

输入脚本已有READY记录时核验hash并退出，避免重复抽样或覆盖已冻结输入。模型恢复核验counts/signature/基因/spot/标签及配置fingerprint；恢复模型参数和Pyro实际优化器状态。参考一旦供下游使用不得继续更换。发生输入冲突请使用新派生run，保留旧证据。

日常更新：`bash code/python.sh code/07_report.py incremental`（将code替换为本运行的绝对code目录）。最终核验：08_validate.py；最终报告打包：07_report.py final。最终生物学解释放07_review/INTERPRETATION.md，不能由技术成功自动推断。

MODEL_INDEX.tsv给出当前8小时窗口内的续训命令。预算结束后脚本默认不再训练；如以后明确决定继续某一映射，可在其命令末尾加 `--resume-seconds 2100`，每次最多3600秒，并生成显式手动恢复记录；该选项不会自动运行。epoch上限和数值停止规则不变。到时仅导出已有模型：`04_mapping.py MODEL_NAME final`。不要重新初始化模型冒称续训。

原始counts、模型.pt、Pyro优化器.pt和w_sf_draws.npy保留服务器，轻量包不包含它们。软件不在本轮升级。输入、种子、实际训练块、环境、检查点及可重现路径见00_manifest和各模型目录。
