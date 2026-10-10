# scope_frozen_resume_v1

打开 REPORT.html。原v5范围阻塞已解除，本轮使用Epi_K32/Epi_other和其余29类。

当前不支持。两套完整31类NB reference与四套DEV/TEST解卷积均完成并达到预先冻结的技术收敛规则。TEST显示K32数量校准不足、低丰度检出很弱，并有来源特异的困难阴性误分配；去除P07T未改变这一结论。当前不支持将本方案直接用于真实Visium的K32定量试运行。本轮未拟合真实空间切片。

- 01_scope：完整scope manifest、18182排除细胞、实际分组与交叉表。
- 02_inputs：共享genes、稀疏counts、输入指纹。
- 03_reference：两套NB模型、signatures、训练与重建QC。
- 04_pseudospots：DEV/TEST counts、真值、成员、RNA贡献和场景可用性。
- 05_deconvolution：四套正式模型、w_sf draws及四种摘要和draw-wise比例。
- 06_evaluation：DEV冻结阈值、逐spot/来源/dataset/场景结果。
- 07_figures：透明PDF/300dpi PNG、plot_data、英文图注。
- 08_review：审阅回答、模型索引、轻量ZIP。

复现按 code/README.md。大体积counts、注册AnnData、模型和w_sf draws保留在服务器，ZIP不包含这些文件。无真实空间拟合。
