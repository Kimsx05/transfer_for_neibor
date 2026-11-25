from base import *
def write():
 rows=[]
 def add(s,branch,files,reason=''):
  ok=all(p.is_file() and p.stat().st_size>0 for p in files)
  rows.append((s,branch,'完成' if ok else '等待上游或执行中',reason if ok else '; '.join(str(p.relative_to(O)) for p in files if not p.exists()),' | '.join(str(p.relative_to(O)) for p in files)))
 for s in HD:
  d=O/'hd'/s
  add(s,'完整 8→16 µm 聚合及成员',[d/'conservation.json',d/'members.parquet'])
  add(s,'RCTD full 及未平滑表达支持',[d/'fit_complete.txt',d/'weights.parquet',d/'labels.parquet',d/'label_counts.tsv'])
  add(s,'50 µm 大类及亚群置换矩阵',[d/'major_50_matrix.tsv.gz',d/'subtype_50_matrix.tsv.gz',d/'matrix_nodes.tsv'])
  add(s,'原 Epi64 H05 成员分离和 200/50 µm 构成',[d/'epi64_centers.parquet',d/'effects_200.tsv',d/'effects_50.tsv',d/'center_composition_200.parquet'])
 for s in XE:
  d=O/'xenium'/s
  for b,ff in [('表达/细胞/原分割及 panel',['cells.parquet','panel.tsv.gz','symbol_audit.tsv.gz','boundary_geometry.parquet']),('大类、细类及独立表达簇',['labels.parquet','subtype_support.tsv','cluster_marker_review.tsv.gz']),('三状态和两个补充表达程序',['ucell.tsv.gz','state_transfer.tsv.gz','threshold.tsv','correlations.tsv']),('分歧区域原表达及来源诊断',['local_state_expression_80um.parquet','discordant_expression.parquet','reference_source_diagnostic.tsv.gz']),('50/200 µm 邻域矩阵',['subtype_50_matrix.tsv.gz','subtype_200_matrix.tsv.gz','chol_50_matrix.tsv.gz','chol_200_matrix.tsv.gz','matrix_nodes.tsv']),('high/other 与连续读数周围组成',['effects_50.tsv','effects_200.tsv','continuous_neighbors_50.tsv','continuous_neighbors_200.tsv'])]:
   add(s,b,[d/p for p in ff])
 add('shared','留整来源及按真实深度降采样验证',[O/'reports/heldout_validation.tsv',O/'shared/validation/PR_curves.tsv.gz'])
 add('Cervilla','同来源大类、同半径 8/16 µm 局部对照',[O/'hd/Cervilla/scale_direction_50.tsv'])
 add('shared','冻结 H05 全成员核验与抽样记录',[O/'shared/H05_membership_check.tsv',O/'shared/reference/membership.tsv.gz'])
 add('shared','区域面积和边界异常敏感性',[O/'reports/boundary_sensitivity.tsv'])
 for s in ['MPBC1','MPBC2']:rows.append((s,'HD 补充切片','仅登记','用户指定本轮三片','shared/samples.tsv'))
 rows.append(('GSE326226','作者最后细胞标签','未取得','GEO 逐文件清单含矩阵、坐标和分割；SpatialBC 代码引用未随仓库提供的 RDS。本轮使用参考转移。','shared/spatialbc_tree.json'))
 put(pd.DataFrame(rows,columns=['sample','branch','status','reason_or_missing','evidence']),O/'reports/completion.tsv')
if __name__=='__main__':write()

