# HD Epi-associated bin niches: dual methods, scales and references

2026-11-12 10:41:33.855931+08:00

## Completion and denominator

Eight main combinations executed on all evaluable shared anchors. Main model outputs are not confined to training/ROI. Training is 80000 anchors per ref/radius (20000 each C/D/UC/MPBC1); MPBC2 is projection-only. Niche counts describe anchors, not cells. Original inputs and labels are unchanged. MPBC original tissue mask remains unknown.

| reference | radius | method | selected_parameter | n_clusters | seed_ARI | block_ARI | silhouette | shared_anchors | evaluable | native_train | external | unresolved | complete |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| blca | 50 | A | 2.0 | 2 | 0.9411693211620089 | 0.9327173106364337 | 0.4081923478869996 | 525552 | 525551 | 80000 | 445549 | 3 | True |
| blca | 50 | B | 0.1 | 182 | 0.639616955281987 | 0.48599288774484517 | 0.054762373346389566 | 525552 | 525551 | 80000 | 407910 | 37642 | True |
| blca | 200 | A | 2.0 | 2 | 0.9102074075676339 | 0.9785137043071513 | 0.3930208519845492 | 525552 | 525551 | 80000 | 445551 | 1 | True |
| blca | 200 | B | 0.8 | 67 | 0.6608970557976348 | 0.596490186331766 | 0.0915398435204939 | 525552 | 525551 | 80000 | 439216 | 6336 | True |
| hpca | 50 | A | 2.0 | 2 | 0.9509034476741607 | 0.9745139706456609 | 0.4693674550561619 | 525552 | 525550 | 80000 | 445544 | 8 | True |
| hpca | 50 | B | 0.6 | 716 | 0.612077054026865 | 0.589193426914143 | 0.01145098719290212 | 525552 | 525550 | 80000 | 267561 | 177991 | True |
| hpca | 200 | A | 2.0 | 2 | 0.8843890049595402 | 0.7644391329790738 | 0.3918256039781895 | 525552 | 525551 | 80000 | 445551 | 1 | True |
| hpca | 200 | B | 1.0 | 75 | 0.6736588746501868 | 0.5696821288911605 | 0.18457697093076603 | 525552 | 525551 | 80000 | 436043 | 9509 | True |

## Shared anchors and support

| section | radius | reference | anchors | known_min | known_q01 | known_median | known_fraction_q01 | coverage_q01 | under10 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Cervilla | 50 | blca | 72171 | 11.0 | 79.0 | 120.0 | 0.9833333333333333 | 0.6583333333333333 | 0 |
| Cervilla | 50 | hpca | 72171 | 11.0 | 79.0 | 120.0 | 0.9658119658119658 | 0.6583333333333333 | 0 |
| Cervilla | 200 | blca | 72171 | 183.0 | 1005.0 | 1960.0 | 0.979732468426736 | 0.5127551020408163 | 0 |
| Cervilla | 200 | hpca | 72171 | 176.0 | 1003.7 | 1960.0 | 0.9548983101120713 | 0.5127551020408163 | 0 |
| D1_Johnson | 50 | blca | 97858 | 12.0 | 104.0 | 120.0 | 0.9833333333333333 | 1.0 | 0 |
| D1_Johnson | 50 | hpca | 97858 | 9.0 | 103.0 | 120.0 | 0.9745762711864406 | 1.0 | 1 |
| D1_Johnson | 200 | blca | 97858 | 123.0 | 1282.0 | 1940.0 | 0.9716237519705728 | 0.8086734693877551 | 0 |
| D1_Johnson | 200 | hpca | 97858 | 125.0 | 1259.0 | 1940.0 | 0.9503636651846707 | 0.8086734693877551 | 0 |
| UC_SCC | 50 | blca | 263581 | 33.0 | 79.0 | 120.0 | 0.9915966386554622 | 0.6583333333333333 | 0 |
| UC_SCC | 50 | hpca | 263581 | 33.0 | 79.0 | 120.0 | 1.0 | 0.6583333333333333 | 0 |
| UC_SCC | 200 | blca | 263581 | 500.0 | 1054.0 | 1960.0 | 0.9764395036936097 | 0.5377551020408163 | 0 |
| UC_SCC | 200 | hpca | 263581 | 500.0 | 1054.0 | 1960.0 | 0.9903902621254462 | 0.5377551020408163 | 0 |
| MPBC1 | 50 | blca | 77767 | 38.0 | 90.0 | 120.0 | 0.975 | 0.75 | 0 |
| MPBC1 | 50 | hpca | 77767 | 34.0 | 89.0 | 120.0 | 0.9583333333333334 | 0.75 | 0 |
| MPBC1 | 200 | blca | 77767 | 603.0 | 1049.0 | 1958.0 | 0.9475904653716893 | 0.5627551020408164 | 0 |
| MPBC1 | 200 | hpca | 77767 | 553.0 | 1041.0 | 1955.0 | 0.9108895673666064 | 0.5627551020408164 | 0 |
| MPBC2 | 50 | blca | 14175 | 4.0 | 70.0 | 120.0 | 0.9745762711864406 | 0.6228333333333335 | 1 |
| MPBC2 | 50 | hpca | 14175 | 4.0 | 68.0 | 120.0 | 0.8607309941520468 | 0.6228333333333335 | 1 |
| MPBC2 | 200 | blca | 14175 | 4.0 | 659.22 | 1960.0 | 0.9358079264094302 | 0.43544897959183676 | 1 |
| MPBC2 | 200 | hpca | 14175 | 4.0 | 628.48 | 1960.0 | 0.828779768458804 | 0.43544897959183676 | 1 |

## Native versus external assignment coverage

| section | reference | radius | method | anchors | evaluable | training | external | unresolved |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Cervilla | blca | 50 | A | 72171 | 72171 | 20000 | 52171 | 0 |
| D1_Johnson | blca | 50 | A | 97858 | 97858 | 20000 | 77858 | 0 |
| MPBC1 | blca | 50 | A | 77767 | 77767 | 20000 | 57767 | 0 |
| MPBC2 | blca | 50 | A | 14175 | 14174 | 0 | 14174 | 1 |
| UC_SCC | blca | 50 | A | 263581 | 263581 | 20000 | 243579 | 2 |
| Cervilla | blca | 50 | B | 72171 | 72171 | 20000 | 42517 | 9654 |
| D1_Johnson | blca | 50 | B | 97858 | 97858 | 20000 | 76535 | 1323 |
| MPBC1 | blca | 50 | B | 77767 | 77767 | 20000 | 45456 | 12311 |
| MPBC2 | blca | 50 | B | 14175 | 14174 | 0 | 12885 | 1290 |
| UC_SCC | blca | 50 | B | 263581 | 263581 | 20000 | 230517 | 13064 |
| Cervilla | blca | 200 | A | 72171 | 72171 | 20000 | 52171 | 0 |
| D1_Johnson | blca | 200 | A | 97858 | 97858 | 20000 | 77858 | 0 |
| MPBC1 | blca | 200 | A | 77767 | 77767 | 20000 | 57767 | 0 |
| MPBC2 | blca | 200 | A | 14175 | 14174 | 0 | 14174 | 1 |
| UC_SCC | blca | 200 | A | 263581 | 263581 | 20000 | 243581 | 0 |
| Cervilla | blca | 200 | B | 72171 | 72171 | 20000 | 51871 | 300 |
| D1_Johnson | blca | 200 | B | 97858 | 97858 | 20000 | 77358 | 500 |
| MPBC1 | blca | 200 | B | 77767 | 77767 | 20000 | 57595 | 172 |
| MPBC2 | blca | 200 | B | 14175 | 14174 | 0 | 12976 | 1199 |
| UC_SCC | blca | 200 | B | 263581 | 263581 | 20000 | 239416 | 4165 |
| Cervilla | hpca | 50 | A | 72171 | 72171 | 20000 | 52171 | 0 |
| D1_Johnson | hpca | 50 | A | 97858 | 97857 | 20000 | 77857 | 1 |
| MPBC1 | hpca | 50 | A | 77767 | 77767 | 20000 | 57767 | 0 |
| MPBC2 | hpca | 50 | A | 14175 | 14174 | 0 | 14174 | 1 |
| UC_SCC | hpca | 50 | A | 263581 | 263581 | 20000 | 243575 | 6 |
| Cervilla | hpca | 50 | B | 72171 | 72171 | 20000 | 35435 | 16736 |
| D1_Johnson | hpca | 50 | B | 97858 | 97857 | 20000 | 71535 | 6323 |
| MPBC1 | hpca | 50 | B | 77767 | 77767 | 20000 | 39779 | 17988 |
| MPBC2 | hpca | 50 | B | 14175 | 14174 | 0 | 11087 | 3088 |
| UC_SCC | hpca | 50 | B | 263581 | 263581 | 20000 | 109725 | 133856 |
| Cervilla | hpca | 200 | A | 72171 | 72171 | 20000 | 52171 | 0 |
| D1_Johnson | hpca | 200 | A | 97858 | 97858 | 20000 | 77858 | 0 |
| MPBC1 | hpca | 200 | A | 77767 | 77767 | 20000 | 57767 | 0 |
| MPBC2 | hpca | 200 | A | 14175 | 14174 | 0 | 14174 | 1 |
| UC_SCC | hpca | 200 | A | 263581 | 263581 | 20000 | 243581 | 0 |
| Cervilla | hpca | 200 | B | 72171 | 72171 | 20000 | 51848 | 323 |
| D1_Johnson | hpca | 200 | B | 97858 | 97858 | 20000 | 77259 | 599 |
| MPBC1 | hpca | 200 | B | 77767 | 77767 | 20000 | 57416 | 351 |
| MPBC2 | hpca | 200 | B | 14175 | 14174 | 0 | 13315 | 860 |
| UC_SCC | hpca | 200 | B | 263581 | 263581 | 20000 | 236205 | 7376 |

## Native method comparison and full assignments

| reference | radius | scope | n | ARI | NMI | AMI |
| --- | --- | --- | --- | --- | --- | --- |
| blca | 50 | native_T | 80000 | 0.24584585433702208 | 0.3459709007285173 | 0.3455026757208551 |
| blca | 50 | full_S | 525551 | 0.22203434466366717 | 0.3605906335573754 | 0.3605159508295443 |
| blca | 50 | full_assigned_both | 487909 | 0.2534597579930537 | 0.36965125892772366 | 0.369573339811931 |
| blca | 200 | native_T | 80000 | 0.03713308574500039 | 0.2580505754251187 | 0.25791894059120424 |
| blca | 200 | full_S | 525551 | 0.050098461807440264 | 0.2696422111347309 | 0.26962140728464784 |
| blca | 200 | full_assigned_both | 519216 | 0.05132600958777423 | 0.27414278573078443 | 0.27412207824738066 |
| hpca | 50 | native_T | 80000 | 0.022325819200868915 | 0.20726907899672103 | 0.2060799765665532 |
| hpca | 50 | full_S | 525550 | 0.1657239643410281 | 0.23647751375706927 | 0.23623978042496818 |
| hpca | 50 | full_assigned_both | 347561 | 0.025707179230169328 | 0.22254333853905509 | 0.22226009859344642 |
| hpca | 200 | native_T | 80000 | 0.03279316173525839 | 0.25907139678016583 | 0.25892686328532766 |
| hpca | 200 | full_S | 525551 | 0.036546554778917675 | 0.2580829157348123 | 0.25805991415867513 |
| hpca | 200 | full_assigned_both | 516043 | 0.03794961275851152 | 0.26323226904659996 | 0.26320922635407784 |

The primary method comparison is native_T. Full_S contains A prediction and B external projection plus unresolved and is not a substitute. ARI/NMI/AMI measure partitions, not identity accuracy. A changes both representation (sqrt P) and clustering relative to B.

## Duplicate-composition diagnostic

| reference | radius | training_profiles | profiles_split_across_native_labels | training_bins_in_split_profiles | pure_epi_train | pure_epi_native_labels | external_exact_ambiguous_unresolved | native_labels_preserved |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| blca | 50 | 55691 | 6 | 5474 | 4608 | 89 | 37441 | True |
| blca | 200 | 79783 | 1 | 98 | 0 | 0 | 0 | True |
| hpca | 50 | 23865 | 56 | 20621 | 11424 | 245 | 176388 | True |
| hpca | 200 | 78666 | 2 | 356 | 258 | 5 | 5671 | True |

Identical compositions split into multiple native B labels cannot be interpreted as biologically distinct environments. Their native labels are retained and flagged; all exact-profile training ties participate in external transfer. No-majority ties stay unresolved. This is not evidence that pure-Epi environments are invalid; it is an implementation/partition sensitivity of the discrete B solution. B labels should be treated as exploratory wherever this occurs.

## Scale and reference comparisons

| comparison | reference_or_radius | method | scope | n | ARI | NMI |
| --- | --- | --- | --- | --- | --- | --- |
| scale | blca | A | common_evaluable_including_unresolved | 525551 | 0.46865014603671873 | 0.4079247108311696 |
| scale | blca | A | common_assigned_both | 525549 | 0.46865412177298826 | 0.40793711432917057 |
| scale | blca | B | common_evaluable_including_unresolved | 525551 | 0.03494710434212095 | 0.205094814937951 |
| scale | blca | B | common_assigned_both | 481886 | 0.03464916496732252 | 0.20310796429474728 |
| scale | hpca | A | common_evaluable_including_unresolved | 525550 | 0.36425489877285155 | 0.2623728731889248 |
| scale | hpca | A | common_assigned_both | 525544 | 0.3642605852000113 | 0.2624000298891322 |
| scale | hpca | B | common_evaluable_including_unresolved | 525550 | 0.031877585133575144 | 0.1626913806976469 |
| scale | hpca | B | common_assigned_both | 344760 | 0.02973737807411414 | 0.17160073508043222 |
| reference | 50 | A | common_evaluable_including_unresolved | 525550 | 0.44746212641527444 | 0.35525704799169 |
| reference | 50 | A | common_assigned_both | 525542 | 0.4474694815132002 | 0.35529156447133087 |
| reference | 50 | B | common_evaluable_including_unresolved | 525550 | -0.0005953686178165496 | 0.19908984089844822 |
| reference | 50 | B | common_assigned_both | 343086 | 0.03179050058869423 | 0.2054778849263446 |
| reference | 200 | A | common_evaluable_including_unresolved | 525551 | 0.3766605396770126 | 0.34956660639611065 |
| reference | 200 | A | common_assigned_both | 525551 | 0.3766605396770126 | 0.34956660639611065 |
| reference | 200 | B | common_evaluable_including_unresolved | 525551 | 0.14090151836472167 | 0.426898418760263 |
| reference | 200 | B | common_assigned_both | 509904 | 0.14358072635476335 | 0.43530597425400264 |

50/200 are nested complete physical circles, not disjoint rings or independent replicates. Continuous common-four-type differences are in reference_coarse_composition.tsv; missing reference classes are not encoded as biological absences. Comparisons use common evaluable centers and disclose separate coverage.

## Biological and technical review

## 本轮结论与使用范围

8个主组合均已执行全量可评价中心计算，非训练中心使用明示的A predict或B外层投影。共享中心525,552个；MPBC2的14,175个中心全部为投影候选，仅代表其局部上皮相关表达，不代表整片组织。五片加载原Seurat后逐一验证，counts、QC、ID保持原样，非中心的niche字段为空。

**主要可重复信号是宽的上皮候选占优—上皮/间质复合环境轴，而不是稳定的一组精细离散niche。** A在四个参考/尺度内部均选K=2，种子ARI0.884–0.951，空间区块诊断ARI0.764–0.979，silhouette0.392–0.469。这不证明只有两种真实生态位；它表明本轮完整候选类型组成更支持这一级宽分割。原生B产生182/67/716/75簇（BLCA50/200、HPCA50/200），silhouette约0.011–0.185、区块ARI约0.486–0.596；不少候选参数得分相近，主resolution并非唯一明确最优。

A/B共同T原生ARI分别0.246/0.037/0.022/0.033，不能声称两方法一致。部分B200分区沿组成梯度呈连续区域，但其细簇不能仅因空间连续就获得生物名称。B50的直接问题尤为明确：BLCA4,608个完全相同的纯Epi组成被分到89个原生标签；HPCA11,424个纯Epi组成被分到245个。输入P完全相同，因而这些区别不能由完整组成特征本身支持。保留原生结果和明确警示，不强行合并成新正式模型，也不把它们叫作数百种真实niche。

零距离外层转移改为汇总全部完全相同训练组成，避免任意取15个同距离样本造成假多数。B最终未定为37,642/6,336/177,991/9,509；其中完全相同组成的原生标签分裂解释37,441/0/176,388/5,671个非训练未定。故B50尤其HPCA50的未定主要是模型/投影歧义，不能误判成低表达、没有上皮环境或空间缺口。T原生标签没有被投影覆盖；其技术敏感性另有标记。A未定仅3/1/8/1，较少未定不是准确率优势。

## 原始邻居表达与局部复核

实际打开五片方法×半径全片图、D1 conflict_focus和UC typical固定ROI，并对照已有UC四ROI H&E、其14基因原始RNA图及前轮确认的D1角蛋白/胶原证据。新图只给中心赋值，淡背景为supplied非中心，深灰为未定中心，空白为没有所画输入；MPBC没有伪造H&E。各ROI的数量、范围与原始图路径在ROI_review_index.tsv。D1既有dense ROI没有共享上皮中心，如实记录，没有替换成更好看的窗口。

BLCA50的两个A环境中，周围原始胶原面板>=2基因检出比例分别约6.5%和46.4%，髓系面板约4.0%和15.8%；上皮面板分别67.4%和71.9%。因此“间质候选较高”环境仍有强上皮原始表达，尤其D1/UC，不应解释为纯成纤维或CAF环境。HPCA50对应胶原比例约9.7%与48.8%，同样支持较宽的胶原相关复合轴。原始marker是相关支持证据，不是独立于SingleR的身份真值。这里统计真实邻居，未仅检查中心自身；也没有使用BANKSY平滑表达。

Cervilla中A显示穿插的局部上皮相关环境与胶原相关环境，200um保留较大环境但改变边界；无supplied图像区域始终无niche。D1中大部分已选中心落在胶原相关复合环境，局部较少胶原的候选Epi占优小区在50um较多、200um更受周围环境影响。conflict_focus在两个参考的A中主要维持复合解释，不能用本轮niche倒推上游Fibroblast预测为正确细胞身份。

UC typical ROI的原始KRT5/KRT6A/KRT14/COL17A1支持右侧上皮样巢，左侧胶原信号及H&E纤维样间隔清楚。A50较能区分巢内与边界，A200将部分巢内中心归入较大复合环境，符合半径扩大覆盖间隔的定义，不是前后注释错误。B50在部分相同组成处表现为散点标签/未定，不能由H&E支持这些细分。

MPBC1上半部分较多复合环境、下半部分更多上皮候选占优，与前轮原始尿路上皮对照基因和胶原分布相符；A两尺度语义一致率约95.9%(BLCA)、91.5%(HPCA)。缺组织掩膜仍限制病理解释。MPBC2只在少量局部中心有支持，外围稀疏中心与缺口附近应保留覆盖警示，不能把投影成功推广至整片。

没有得到可跨两方法稳定命名的独立淋巴、髓系、血管离散niche；这些成分继续作为全组成中的连续量保留。未设置最少TME数量，充分支持的纯Epi环境有效。某些免疫面板与其他谱系共享基因，其检出不等于对应细胞比例。

## 尺度、参考和技术因素

在共同可赋值中心上，A两尺度按组成对齐后的宽环境一致率为Cervilla约79–81%、D1约76–88%、UC约79–80%、MPBC1约92–96%、MPBC2约78–85%。这是嵌套环境的相似性，不是独立重复验证。50/200全局分区ARI约0.469(BLCA)、0.364(HPCA)，不能把两尺度结果视为相同。

A跨参考语义一致率约74–97%，但绝对组成明显参考依赖：BLCA50复合组平均Epi约38%、Fibroblast44%、Mural11%；HPCA50复合组平均Epi约63%、Fibroblast33%，HPCA没有独立Mural/Mast/Plasma类别。两者不能当作同一组真实细胞比例。全共同中心四粗类组成差异与已知分母差别均有表格；没有用原始SingleR分数比较赢家。UC200跨参考宽标签一致率约73.8%，是需特别保留参考敏感性的例子。

所选中心周围已知候选比例中位数接近1，A两组低信息邻居平均约0.6–2.8%；宽组成差别不像主要由低检测量驱动。但“已知候选”不等于高置信身份。200um部分上皮占优区域的输入覆盖较低，边界不能忽略；MPBC组织范围未知，input coverage不是组织覆盖。A复合环境的参考冲突比例明显较高（如BLCA50约35%，另一组约9%），故其生物学名称要同时保留参考敏感说明。B50的相同P分裂则是明确的方法敏感性证据。

## 下一轮MP10共定位的适用性（本轮未运行）

可以优先进入预先声明的区域层面探索：A的宽环境及完整连续组成，保留50/200两尺度和BLCA/HPCA两个参考作为敏感性分析，不选显著性更好的一支。Cervilla、D1、UC以及MPBC1输入充分区域均有可用的原始表达支持；D1/UC的胶原相关复合区域不能作纯CAF或纯上皮细胞特异归因。MPBC2仅限明确支持的局部中心，保持组织范围未知。

B尤其50um的细标签、完全相同组成被拆分的native标签、未定/低支持/输入缺口位置，不适合直接作为独立生物类别开展正式MP10关联。B200可作为连续组成梯度的探索性分区对照，尚非稳定命名生态位。后续统计必须考虑重叠邻域和切片层级，不能把几十万个中心当成独立生物重复。没有运行MP10、通讯、恶性判定或精细亚群。

## 与普通Visium的只读对应
已读取冻结的129,632个普通Visium中心及正式R12标签，未重跑63张切片，也未载入MP10值。9个已接受niche在共同四类上同时包含Epi和不同程度免疫/间质成分；例如旧niche0的描述性等权环平均约81.4%Epi、13.0%Immune，部分其他旧类约36–39%Immune。本轮HD的A主要保留较宽的Epi—胶原相关轴，没有对应地恢复这些独立免疫细分，不能将旧niche编号套用HD。普通Visium数值来自每spot归一化c2l q05并对R1/R2等权平均，HD是已知类型bin等票份额，不能把绝对比例差直接解释为真实细胞生态差异。完整30→4映射和冻结旧组成已交付。


## Files and targeted verification

Each section contains anchor_manifest.parquet, r50/r200.parquet shared counts/proportions/support, geometric_coverage.parquet and original-neighbor marker support. Each blca_r50/blca_r200/hpca_r50/hpca_r200 directory contains labels.parquet (both methods, all shared centers), training_anchor_ids.tsv, A_model.joblib and B_training.npz. No eight copies of expression objects are made. tables/coverage_wide.tsv provides all 40 section/reference/radius/method coverage rows; tables/main_combinations.tsv provides eight completed main combinations.

| section | n_original | anchors | UMI | IDs_aligned | original_counts_unchanged | original_QC_unchanged | nonanchors_not_labeled | loader_executed |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Cervilla | 268924 | 72171 | 64443819 | True | True | True | True | True |
| D1_Johnson | 364805 | 97858 | 131440767 | True | True | True | True | True |
| UC_SCC | 692681 | 263581 | 126389649 | True | True | True | True | True |
| MPBC1 | 675380 | 77767 | 107299054 | True | True | True | True | True |
| MPBC2 | 422741 | 14175 | 39606795 | True | True | True | True | True |

Geometry, public API equivalence and native-label preservation records: geometry_validation.tsv, api_validation.json and validation.tsv. No computational main branch remains pending. Biological discreteness of B fine labels and MPBC tissue membership remain unresolved limitations.

## Reuse

Dependencies: sibling hd_bin8_banksy_extension_v1, hd_bin8_comparison_v1, hd_bin8_banksy_v1. No expression copies are needed.

```r
root <- "/home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/analysis/hd_epi_niche_dual_methods_v1"
source(file.path(root,"code/load_result.R"))
x <- attach_hd_niches(x, root, "UC_SCC") # existing original Seurat object
```

Python: code/load_result.py::attach_niches(adata, root, section), or load_niche_metadata(root,section). Non-anchor rows retain missing niche labels; epi_anchor=False. No circle members inherit their center label.

code/ plus rerun.sh reproduce definitions and completed steps; models and P/counts remain on server. Review zip includes only summaries, reports, code, and key PNGs. Ordinary Visium correspondence is detailed in METHODS; no ordinary Visium or MP10 rerun.

New allocated disk: 1.923 GiB (including isolated environment); review.zip 35.0 MiB.
