"""Derive compact delivery summaries from completed authoritative v9 results."""
from pathlib import Path
import json,hashlib,datetime
import pandas as pd
OUT=Path(__file__).resolve().parent
p=pd.read_parquet(OUT/'primary_relations.parquet')
c=pd.read_csv(OUT/'combination_summary.tsv',sep='\t',dtype={'combo_id':str})
assert len(c)==99 and c.status.eq('COMPLETE').all() and p.global_status.eq('FINAL').all()
assert json.loads((OUT/'audit/final_validation.json').read_text())['status']=='PASSED'
masks={'within_BH':p.q_within_BH.lt(.05),'global_BH':p.q_global_BH.lt(.05),'within_only':p.q_within_BH.lt(.05)&~p.q_global_BH.lt(.05),'global_BY':p.q_global_BY.lt(.05),'inhouse_supported_global_BH':p.q_global_BH.lt(.05)&p.contains_inhouse&p.inhouse_actual_support}
cat=c[['combo_id','datasets','n_datasets','contains_inhouse','K']].copy()
for name,mask in masks.items():cat[name+'_relations']=cat.combo_id.map(p[mask].groupby('combo_id').size()).fillna(0).astype(int)
cat.to_csv(OUT/'final_combination_categories.tsv',sep='\t',index=False)
selected=pd.read_csv(OUT/'review/selected_relations.tsv',sep='\t',dtype={'combo_id':str});sens=[]
for cid in selected.combo_id.unique():
 s=pd.read_parquet(OUT/'review'/cid/'spotted_sensitivity.parquet')
 for hit in selected[selected.combo_id.eq(cid)].itertuples():
  a=s[s.target_id.eq(hit.target_id)&s.layer.eq(hit.layer)&s.metric.eq(hit.metric)&s.track.eq('SOURCE_P05')&s.eligible&s.effect.notna()]
  for (k,domain),g in a.groupby(['k','analysis_domain']):
   units=g.groupby('analysis_unit_uid')[['effect','stratified_effect','component_equal_effect']].mean()
   for statistic in units.columns:
    v=units[statistic].dropna();sens.append(dict(combo_id=cid,target_id=hit.target_id,layer=hit.layer,metric=hit.metric,k=k,analysis_domain=domain,statistic=statistic,n_units=len(v),mean_effect=v.mean(),same_direction=bool(v.mean()*hit.mean_effect>0),interpretation='descriptive sensitivity; no new primary P value'))
pd.DataFrame(sens).to_csv(OUT/'review/sensitivity_summary.tsv',sep='\t',index=False)
lines=['# v9 最终结果摘要','完成：99/99，含 inhouse 57/57，不含 42/42。2026-10-19 北京时间 07:34 主分析完成；08:47 图包与敏感性完成；08:49 最终核验完成（902 项通过）。','主分析：REFIT、SOURCE_P05、k10、ALL_POLICY_UNITS、两层、各组合全部 niche 的 O/D；planned family 2,876 条，其中 2,621 条可估计，255 条科学性不可估计均保留。','| 结果类别（q<0.05） | 关系数 | 涉及组合数 |','|---|---:|---:|']
labels={'within_BH':'组合内 BH 阳性（含全局阳性）','global_BH':'全局 BH 阳性','within_only':'仅组合内 BH 阳性','global_BY':'全局 BY 支持','inhouse_supported_global_BH':'全局 BH 阳性且 inhouse 有实际支持'}
for name,mask in masks.items():lines.append(f'| {labels[name]} | {int(mask.sum())} | {p[mask].combo_id.nunique()} |')
lines += ['', '完整的三个类别组合名单、真实 dataset 成员及各类别关系数见 [final_combination_categories.tsv](final_combination_categories.tsv)；各类别对应计数列大于 0 即为该类别名单。权威完整关系见 [primary_relations.parquet](primary_relations.parquet)，全局命中见 [global_primary_hits.tsv](global_primary_hits.tsv)，仅组合内探索性命中见 [within_combo_exploratory_hits.tsv](within_combo_exploratory_hits.tsv)。', '85 条全局 BH 阳性中，spot 层 80 条、epi 层 5 条；O 重叠减少 69 条、增加 1 条，D 距离增加 14 条、减少 1 条。含 inhouse 的全局阳性共 39 条/23 组合，其中 37 条/22 组合有可估计且同向的 inhouse 支持。实际支持不表示 inhouse 单独显著。', '## 重点组合', '| 组合 | 数据集范围 | 重点关系（spot 层） | 单位平均 O 效应 | 全局 BH q | 全局 BY q | inhouse 可估计单位数 |','|---|---|---|---:|---:|---:|---:|']
for cid,scope in [('1111111','全部七个数据集'),('1011111','除 GSE171351 外的六个数据集'),('1010111','除 GSE171351、GSE285715 外的五个数据集')]:
 r=p[p.combo_id.eq(cid)&p.target_id.eq('niche3')&p.metric.eq('O')&p.layer.eq('spot_MP10_model_universe')].iloc[0]
 lines.append(f'| {cid} | {scope} | niche3 的重叠减少 | {100*r.mean_effect:.2f} 个百分点 | {r.q_global_BH:.6g} | {r.q_global_BY:.6g} | {r.inhouse_n_units} |')
lines += ['','niche 编号属于各自重拟合，不能仅凭同编号认定跨组合身份完全相同；上述 niche3 的环组成均以上皮为主、Matrix_CAF 为主要非上皮成分之一。这是 MP10 热点空间重叠关系，不是代谢物浓度或主动排斥机制的证据。', '全七数据集的 spot 层 niche2、niche3、niche8 均通过全局 BH，q 分别为 0.02936、0.003924、0.01006；仅 niche3 通过全局 BY。epi 层 niche2 的组合内 q=0.02128，但全局 q=0.06029，属于仅组合内探索性结果。', '## 敏感性与限制','选中的 3 个组合共 8 条关系：逐一剔除 dataset 或 analysis unit 后均未出现平均效应方向逆转；组成匹配后的 REFIT/FROZEN_LABEL 方向一致 8/8。k6/k18、几何支持、既有技术分层和组件等权的单位平均效应也均同向，完整数值见 [review/sensitivity_summary.tsv](review/sensitivity_summary.tsv)。这些是方向和影响诊断，不等于剔除后依然显著或排除单个样本对 P 值的影响。','1011111 的 niche1 虽是全局 BH 阳性，但 inhouse 的平均效应与总体反向，未被计入实际支持。inhouse 七个样本患者身份尚未完整核实，沿用独立样本政策；不能称为七名已核实独立患者。99 组合共享数据，不能视为独立重复验证；BH/BY 不消除结果选择偏差。','全部组合统计、辅助 family、紧凑标签、切片/单位效应及复现记录均已保留。3 个重点组合共 164 张 PNG 及对应 PDF，另有全局 overview。无尚未完成的计算项；身份核实与独立验证仍是科学解释的限制。','资源：最终核验时实际分配约 1.13 GiB（预算 3 GiB），累计清理本轮可再生临时文件约 6.94 GiB。新增摘要为轻量衍生文件，最终字节数见 audit/delivery_summary_manifest.json。','完整执行报告与核验见 [RUN_SUMMARY.md](RUN_SUMMARY.md) 和 [audit/final_validation.json](audit/final_validation.json)。']
import re
(OUT/'FINAL_RESULTS.md').write_text(re.sub(r'(?<=\|)\n\n(?=\|)', '\n', '\n\n'.join(lines))+'\n')
used=sum(f.stat().st_blocks*512 for f in OUT.rglob('*') if f.is_file() and not f.is_symlink());budget=json.loads((OUT/'config.json').read_text())['budget_bytes'];assert used<budget
files=['final_combination_categories.tsv','review/sensitivity_summary.tsv','FINAL_RESULTS.md']
manifest=dict(status='PASSED',generated_at=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),authoritative_primary_sha256=hashlib.sha256((OUT/'primary_relations.parquet').read_bytes()).hexdigest(),files={f:hashlib.sha256((OUT/f).read_bytes()).hexdigest() for f in files},allocated_bytes=used,budget_bytes=budget,sensitivity_rows=len(sens),all_selected_sensitivity_directions_agree=all(r['same_direction'] for r in sens))
(OUT/'audit/delivery_summary_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps({k:v for k,v in manifest.items() if k!='files'},indent=2))
