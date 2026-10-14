# 可用续跑及验收

```bash
bash "/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring123_broad_v7/code/run_ring123_broad_v7.sh" --stage all --approval "/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring123_broad_v7/00_approval/approved_gate_v7.tsv"
bash "/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring123_broad_v7/code/run_ring123_broad_v7.sh" --verify-only
```

检查点绑定输入manifest和冻结spec哈希。完整性验证通过后跳过一致阶段；绘图只补缺失文件。命令实际执行全量三环与关联，不再调用v6 Gate-only接口。原始输入变化会报冲突，不能静默复用缓存。
