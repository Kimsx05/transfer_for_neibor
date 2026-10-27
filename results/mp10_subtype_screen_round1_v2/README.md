# MP10 round1 v2

全上皮MP10、标本内epi状态、B/Plasma重点分析；v1保持只读。本目录独立承载全部新增结果。

- [中文报告](round1_v2_report.md)：主要发现、反向/非显著结果、覆盖限制与空间承接。
- [方法](methods.md)：两个分组定义、评分、统计、CellChat及表达复核口径。
- [轻量阅读包](review.zip)：核心PDF、关键表、报告和方法；无表达矩阵或原始母对象。
- [逐队列图索引](review/FIGURE_INDEX.md)：含ALL与五个dataset，各含无支持和NE图。
- `tables/main_candidate_evidence_v2.tsv`：主要B/Plasma代表轴的分层证据；`all_groups_evidence_summary.tsv`保留其他谱系。
- `tables/epi_state_LR_contrasts_v2.tsv.gz`：主状态原始概率差；`secondary_condition_LR_contrasts_v2.tsv.gz`：次级独立对象条件差。
- `objects/networks/`：68个小型数值网络；`objects/networks_matched_B_Plasma/`：6个共同覆盖敏感性网络。
- `scripts/run_round1_v2.sh`：本地重现顺序；现有网络按冻结参数缓存复用。R库和Python依赖路径在脚本内明确记录。

标本级字段是 `sample_MP10_group_allEpi`；细胞状态字段是 `epi_MP10_state_within_specimen`。请勿把两者统称为同一种high/low。内部CellChat P仅用于网络展示支持，不是患者/条件差异P。患者ID全部为现有derived/provisional记录，不宣称已核实。

不自动开展空间模块，不发布或上传到外部服务。
