# 仍需人工访问的输入



下列操作仅由你在自己的浏览器申请/联系；没有发送作者邮件，也不需要提供密码、cookie或令牌。先获取目录/小文件，再按目标成员下载。

| 优先级 / 对象 | 打开的完整入口；已尝试与实际结果 | 具体需要你取得的内容 | 服务器目标目录 / 能解锁什么 |
|---|---|---|---|
| 1 / MPBC1、MPBC2 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318916 ； https://www.thno.org/v16p6892.htm 。GEO、论文全文/补充入口、关联accession已查，公开GEO缺空间文件；不是网页403导致漏判 | 向作者索取两例 **Bin2cell原始细胞counts、cell坐标、边界/掩膜、2µm bin-to-cell映射**及样本对应；若不能提供，索取同片H&E、tissue_positions.parquet、scalefactors_json.json及变换。作者对象真实文件名未公开，请先取得目录，不猜h5ad名称 | `/home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/data/GSE318916/MPBC1/` 和 `MPBC2/`；可能同时补齐两张，最快达到≥3张 |
| 2 / UC_SCC | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9040999 ； https://link.springer.com/article/10.1007/s12672-026-04927-z 。GEO与既有完整spatial目录仅8µm；补充材料入口已查 | 先索取作者原始 `square_002um/` 目录的 counts＋positions＋scales，以及同片原生HE或含图像Loupe/feature_slice；或者完整作者细胞counts＋空间映射。没有现成结果时才讨论FASTQ重跑所需slide/probe/配准 | `/home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/data/GSE299573/UC_SCC/`；解锁该HD section的忠实2µm到细胞路线 |
| 3 / STT0000103 | https://db.cngb.org/stomics/project/STT0000103/spatialGeneExpression ； https://db.cngb.org/data_access/ 。直连超时后HTTPS代理已成功，官方文件表明确 `file_is_public=false` 且空file_link，属于受限文件而非不存在 | 页面 Spatial Gene Expression → Download 按提示进入CDA申请。指定人类 STSA0000790–STSA0000811 / STTS0001304–STTS0001325，先请求目录、患者/切片表、bin或cell counts粒度、图像、坐标/变换、cellbin边界/掩膜。**看不到真实文件名，需提供目录内容或可下载链接**，不要先拉全部82.96GB | `/home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/data/STT0000103/`；可争取18患者22张人类切片的全转录组高分辨输入，是否具备细胞级形态支持仍须验证 |
| 4 / Oh/Jung TMA | https://zenodo.org/records/18244722 ； https://zenodo.org/records/18328035 。两条官方API均restricted、未公开文件列表 | 浏览器 Request access；先取得两部分文件目录、12患者/2slide/core映射和panel，再挑counts＋cells/边界＋focused morphology。隐藏文件名称未知；若分卷压缩需先确认可单独获取小型输入，避免盲目下载全部分卷 | `/home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/data/Xenium_OhJung_12patient/`；解锁新患者的Xenium输入，但完整Top20是否测量需panel核对 |
| 低 / STOmics演示线索 | https://en.stomics.tech/CSPTestedData/list.html?id=2844&sTypes=%2C1470%2C 。官方索引出现Human Bladder Cancer/FFPE/ssDNA，直接页面访问未得到可验证目录 | 浏览器选择对应demo；提供真实目录、样本来源和公开链接。不要猜文件名；需确认是否与研究数据重复及是否有counts、形态图、坐标 | `/home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/data/STOmics_Bladder_demo/`；只有线索，未计入已确认新研究或切片 |


本轮没有联系作者、发送邮件或索取登录凭据；上述操作不阻塞当前两张HD。
