# 第一批人工复核

从每片的 index.html 打开抽查列表；cards/ 存独立 PNG/PDF，atlas.html 可选取全体候选 ID 和层内对照。每个窗口以蓝框标中心，方向 x 向右、y 向下，比例尺 20 µm。HD H&E 上覆盖 8 µm 单元和 16 µm 网格。marker 显示各自 8 µm 原始 UMI，候选表中的计数是所属 16 µm 单元的总和。Xenium 展示原分割边界，未获得同片配准染色图。紫灰为面板缺失，白色为实测零。邻格信号不作为中心单位内共表达。

先检查单个中心是否有可靠 T 线索，再比较 FOXP3、IL2RA、CTLA4、IKZF2、竞争身份、深度及边界。记录 supported / possible / unresolved / unsupported，指出有实质影响的证据。人工仅需完成每片 human_review.tsv 的待填区。保留 id；human_grade 填判读，human_reviewer 填署名，human_status 改为 reviewed；若需修订父类或上皮资格，可填 human_reliable_T 或 human_epi_eligible 为 true/false，comment 记理由。未查看条目保留 pending。复核表可复制到独立目录后按切片名命名。系统性问题应描述受影响层和证据，由整层修订生成新版本。

导入入口 code/import_review.py 验证署名及父子关系，以签署内容哈希建立新的 revisions 子目录并记录前后变化。初始批次仅用于发现系统性误判，不能估计总体真值灵敏度或特异度。
