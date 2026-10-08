from pathlib import Path
import csv,html,statistics
R=Path('/home/data/t070721/codex_workspace/Bladder Metabolism/mp_ucell_state_v1')
def rd(p):
 with (R/p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def tab(rows,cols):
 return '<div class="scroll"><table><tr>'+''.join('<th>'+html.escape(k)+'</th>' for k in cols)+'</tr>'+''.join('<tr>'+''.join('<td>'+html.escape(str(r.get(k,'')))+'</td>' for k in cols)+'</tr>' for r in rows)+'</table></div>'
best=rd('06_comparisons/best_cluster_Jaccard.tsv');rows=[x for x in best if x['branch_A']=='FULL_Z_res0.4' and x['A']=='C8' and x['direction']=='A']
profiles=rd('07_mp10_review/all_cluster_14MP_profiles.tsv');p8=[x for x in profiles if x['instance']=='FULL_Z_res0.4' and x['cluster']=='C8']
for row in rows:
 for k in ['fraction_A','fraction_B','Jaccard']:row[k]=f'{float(row[k]):.4f}'
for row in p8:
 for k in ['raw_mean','raw_median','global_z_mean']:row[k]=f'{float(row[k]):.4f}'
body='''<p><strong>结论：14维MP活动组合可以形成预设图聚类，但目前不足以认定一个跨尺度、跨样本构成均稳定的唯一MP10细胞亚群。更符合“多个MP10相关活动组合/连续变化区域，边界随尺度、样本构成及resolution拆分或合并”的解释。</strong>没有独立MP10群不等于技术失败；有MP10高群也不等于已验证的新亚群。</p>
<p>Primary FULL_Z/res0.4中C8是最明确的MP10偏高候选：8,603细胞，占FULL的6.95%；MP10均值0.27365，中位数0.25023，FULL全局z均值1.413。它同时升高MP07(z=1.450)、MP03(1.341)、MP11(1.280)、MP04(1.236)和MP08(0.807)，而MP01、02、06、12、14不整体升高。因此应暂称“C8：MP07/10/03/11/04共同偏高组合”，不命名为脂肪酸代谢阳性或恶性亚群。MP10与MP07的top50共享12个基因，部分共同变化可能与signature重叠有关。</p>
<p>C8来自45个样本，最大样本CNP0000460_P01T贡献14.91%；45个有群内/群外细胞的样本中44个MP10均值差为正。但两大dataset CNP0000460与HRA003620合计88.67%，GSE326225仅3个细胞，不能把“五dataset有出现”解释成五dataset均有充分重复。另有InhouseData 442、PRJNA662018 530个C8细胞。C8包含Normal 2,862（33.27%）、NMIBC 1,860（21.62%）、MIBC 3,881（45.11%）；来源混合进一步说明不能直接命名为恶性亚群。患者层面未知。</p>
<p>其他高于FULL MP10均值的群为C6、C3、C13、C7、C12，完整轮廓与分布见主表。C6的MP10均值0.20750，主要伴随MP03/14/04，单一P7T7样本占69.05%；C7和C12仅轻度高于FULL MP10均值，且最大样本分别占80.48%、60.40%。C3（MP10均值0.19278）更突出MP08/06/13，C13（0.18837）更突出MP11/04；它们并非同一种“MP10高状态”的简单复本。</p>
<p>尺度比较：FULL_Z vs FULL_RAW全局ARI=0.671、NMI=0.700，但C8的最佳单群对应仅Jaccard=0.286（RAW C4）。C8的54.74%进入RAW C4，19.44%进入RAW C11；两个RAW群还合并了其他FULL_Z细胞。这是实质的拆分/合并与成员重排，不能仅称为小幅边界抖动。RAW中的MP10偏高组合仍存在（C4均值0.24379、C11均值0.22989），但群成员不等价。</p>
<p>样本上限抽样：三次BAL_Z都保留MP10偏高群（seed42 C5均值0.252、seed43 C3均值0.264、seed44 C3均值0.265），但C8最佳匹配Jaccard为0.379、0.289、0.462；共同细胞内FULL_Z与BAL_Z全局ARI为0.429、0.418、0.510。三次BAL_RAW同样保留MP10偏高组合，但更常分散在多个群，FULL_RAW与BAL_RAW ARI为0.379、0.448、0.421。seed42保持预先指定的主要展示，不按结果改选seed44。BAL细胞比例不用于估计原始丰度。</p>
<p>Resolution对照：FULL_Z在0.2/0.4/0.6分别为7/13/18群。0.2时C8主要并入大群；0.6时C8的63.28%和31.08%分别进入C12、C16，两群中来自原C8的比例分别为97.28%、96.22%。两者仍均为MP10高，但完整MP组合不同。这支持对C8内部异质性继续检查，不能据此挑选0.6作为新的主分析。</p>'''
body+='<p>C8在各固定对照中的最佳对应（size_A/B均为共同细胞内大小）：</p>'+tab(rows,['branch_B','B','overlap','size_A','size_B','fraction_A','fraction_B','Jaccard'])
body+='<p>C8完整14MP轮廓：</p>'+tab(p8,['MP','raw_mean','raw_median','global_z_mean'])
body+='''<p>原MP10 usage与UCell总体Spearman=0.4450（N=123,699）；57个样本中56个可评价，相关中位数0.4414，范围−0.1054至0.8485，55个为正。另1个样本CNP0000460_P02T仅有2个细胞，原MP10 usage均为0而为常数，保留NA；极小样本即使可计算相关也不作为重复验证。该关系说明有共享信息但并非可互换的尺度，也不是外部验证。</p>
<p>下一轮marker优先级：①以primary C8为首要候选，逐样本检查表达和检出，并检查其固定res0.6的C12/C16内部差异，及RAW C4/C11中的对应细胞；②以C6、C3作为MP10偏高但组合/来源不同的比较群，C13作为MP11偏高的补充比较。优先寻找能区分这些组合且跨多个生物样本保持表达/检出的marker，而不是只重复MP10 signature里的基因。应核实深度、增殖、分化状态、肿瘤/Normal来源及dataset影响；“脂肪酸代谢”“恶性”“新亚群”等功能或身份命名仍待独立证据。</p>
<p>本轮没有建立离散状态生成模型或外部验证，所以“连续梯度”也不是已被正式证明的替代机制。现有数据支持MP10在多个MP组合上变化，不能可靠区分真实离散边界与图分割连续结构。样本内方向一致、usage相关及UMAP分离均不能单独替代患者层面或外部验证。</p>'''
(R/'07_mp10_review/interpretation_final.html').write_text(body)
with (R/'07_mp10_review/C8_correspondence_summary.tsv').open('w') as f:
 cols=['branch_B','B','overlap','size_A','size_B','fraction_A','fraction_B','Jaccard'];w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
print('Interpretation saved')
