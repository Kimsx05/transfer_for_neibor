from common import *
import numpy as np,pandas as pd
from programs import MAJOR,STATE
from expression_evidence import SUBTYPE
from scipy.spatial.distance import cdist

def tag(s):return s.replace('/','_')
o=pd.read_parquet(ROOT/'cell_annotations.parquet');p=pd.read_parquet(ROOT/'metadata/program_evidence.parquet');methods=json.loads((ROOT/'metadata/lineage_structure_methods.json').read_text());pc=np.load(ROOT/'metadata/lineage_pca.npz')['X_pca'];rows=[];nearest=[]
sets=[('global','global','global_cluster_v21',o)]+[(m,'lineage_'+tag(m),'lineage_cluster_v21',o[o.proposed_major.eq(m)]) for m in methods]
for lineage,prefix,field,ob in sets:
 markers=pd.read_parquet(ROOT/f'metadata/{prefix}_markers.parquet');summary=pd.read_csv(ROOT/f'metadata/{prefix}_summary.tsv',sep='\t');summary.cluster=summary.cluster.astype(str)
 neighbor={}
 if lineage!='global':
  means=[];cn=[]
  for cl,grp in ob.groupby(field,observed=True):
   a=pc[o.index.get_indexer(grp.index)];valid=np.isfinite(a[:,0]);
   if valid.any():means.append(np.nanmean(a[valid],axis=0));cn.append(cl)
  if len(cn)>1:
   d=cdist(np.nan_to_num(means),np.nan_to_num(means));np.fill_diagonal(d,np.inf);neighbor=dict(zip(cn,np.array(cn)[d.argmin(1)]))
 for _,s in summary.iterrows():
  cl=str(s.cluster);scope=s.scope;z=ob[ob[field].astype(str).eq(cl)];z=z if scope=='joint' else z[z.section_id.eq(scope)]
  if not len(z):continue
  m=markers[markers.cluster.astype(str).eq(cl)&markers.scope.eq(scope)].set_index('gene');major=z.proposed_major.value_counts();sub=z.proposed_subtype.value_counts();dominant=major.index[0];fine=sub.index[0];label=fine if fine!='Unknown' else dominant
  gs=SUBTYPE.get(fine,MAJOR.get(dominant,[]));positive=[g for g in gs if g in m.index and m.loc[g,'detection_fraction']>0];positive=sorted(positive,key=lambda g:-m.loc[g,'detection_fraction'])[:12]
  positives=';'.join(f'{g}:pct={m.loc[g,"detection_fraction"]:.3f},delta={m.loc[g,"mean_log2_difference"]:.3f}' for g in positive)
  allcomp=[]
  for other,gg in MAJOR.items():
   if other==dominant:continue
   for gene in gg:
    if gene in m.index:allcomp.append((gene,other,float(m.loc[gene,'detection_fraction'])))
  allcomp=sorted(allcomp,key=lambda q:-q[2])[:8];comp=';'.join(f'{g}({a}):{v:.3f}' for g,a,v in allcomp)
  states={c[6:]:round(float(z[c].mean()),3) for c in z if c.startswith('state_') and z[c].dtype==bool and z[c].any()};state=sorted(states,key=states.get,reverse=True)[:5]
  frac_source=z.section_id.value_counts(normalize=True).max() if scope=='joint' else np.nan;mtgenes=sum(g.startswith('MT-') for g in str(s.top_genes).split(';'));nonspecific=mtgenes>=5 or s.median_genes<60 or pd.isna(s.top_genes)
  hetero=1-float(major.iloc[0]/len(z));question=[]
  if mtgenes>=5 or s.median_genes<60:question.append('低深度/线粒体主导；保留可支持的大类，勿作为新身份')
  if pd.isna(s.top_genes):question.append('无满足描述性筛选的特异富集基因或样本太少；不由此创造新身份')
  if lineage=='Erythroid':question.append('仅2个参考外血红蛋白程序对象；优先核查非核/边界污染，不能认定稳定红系群体')
  if z.mixed_observed_programs.mean()>.05:question.append('同对象竞争谱系共表达；分割/污染与真实共表达尚不能区分')
  if z.secondary_lineage_signal.mean()>.05:question.append('主要身份外仍有成套竞争转录本；局部簇可能由背景/边界混入驱动，勿将其差异直接视为新身份')
  if hetero>.2:question.append('该簇包含不同大类成员；按逐细胞建议审阅，禁止整簇覆盖')
  if lineage=='Mural':question.append('核对 RGS5/PDGFRB 与 MYH11/CNN1 的 pericyte/SMC 分界及 fibroblast 反证')
  if fine=='Unknown' and lineage not in ['global','Epithelial','Plasma','Mast']:question.append('大类可支持，细亚群 marker 稀疏或相近程序竞争；保留父类')
  if frac_source>.95:question.append('来源偏倚；目前每来源一张，不能拆分技术与生物差异')
  if not question:question=['核对建议命名与所列 marker；当前无明显专门阻塞']
  neighborcl=neighbor.get(cl,'not_applicable')
  if neighborcl!='not_applicable':
   c=markers[markers.cluster.astype(str).eq(neighborcl)&markers.scope.eq(scope)].set_index('gene');shared=m.index.intersection(c.index)
   for gene in shared:
    nearest.append(dict(lineage=lineage,scope=scope,cluster=cl,nearest_expression_cluster=neighborcl,gene=gene,n=len(z),neighbor_n=int(c.n.iloc[0]),mean_log2_difference=float((m.loc[gene,'mean_log1p10k']-c.loc[gene,'mean_log1p10k'])/np.log(2)),detection_fraction=float(m.loc[gene,'detection_fraction']),neighbor_detection_fraction=float(c.loc[gene,'detection_fraction']),comparison='nearest local PCA centroid within corrected lineage; descriptive'))
  rows.append(dict(review_entry=('global::'+cl if lineage=='global' else cl),analysis_level='global' if lineage=='global' else 'lineage',lineage=lineage,scope=scope,cluster=cl,n=len(z),section_composition=s.section_composition,median_counts=s.median_counts,median_genes=s.median_genes,low_information_fraction=s.low_information_fraction,old_major_composition=s.old_major_composition,old_v2_labels=s.v2_labels,proposed_major_composition=json.dumps(major.to_dict()),proposed_subtype_composition=json.dumps(sub.to_dict()),suggested_dominant_label=label,supporting_markers=positives,competing_lineage_markers=comp,data_driven_top_markers=s.top_genes,SingleR_major_composition=json.dumps(z.SingleR_major_candidate.value_counts().head(5).to_dict()),SingleR_subtype_composition=json.dumps(z.SingleR_subtype_candidate.value_counts().head(5).to_dict()),CellTypist_SingleR_major_agreement=float((z.predicted_major.astype(str)==z.SingleR_major_candidate).mean()),SingleR_major_pruned_fraction=float(z.SingleR_major_pruned.mean()),measured_mixed_fraction=float(z.mixed_observed_programs.mean()),secondary_lineage_signal_fraction=float(z.secondary_lineage_signal.mean()),major_heterogeneity_fraction=hetero,depth_or_mito_driven_flag=nonspecific,max_section_fraction=frac_source,state_program_fractions=json.dumps(states),nearest_competitor_cluster=neighborcl,review_question='；'.join(question),review_status='computational_proposal_pending_user_review'))
r=pd.DataFrame(rows);r['review_priority']=np.where(r.measured_mixed_fraction>.1,'1_conflict',np.where(r.secondary_lineage_signal_fraction>.05,'2_secondary_signal',np.where(r.depth_or_mito_driven_flag,'3_low_information','4_naming')));r.sort_values(['review_priority','analysis_level','lineage','cluster','scope']).to_csv(ROOT/'ALL_CLUSTER_REVIEW.tsv',sep='\t',index=False);pd.DataFrame(nearest).to_parquet(ROOT/'metadata/nearest_cluster_marker_comparisons.parquet',compression='zstd')
for name,col in [('major','proposed_major'),('subtype','proposed_subtype'),('unknown_reasons','unknown_reason')]:pd.crosstab(o[col],o.section_id).to_csv(ROOT/f'metadata/counts_{name}.tsv',sep='\t')
o.groupby(['section_id','old_major_label','proposed_major'],observed=True).size().rename('n').reset_index().to_csv(ROOT/'metadata/major_transition_v2_to_v21.tsv',sep='\t',index=False)
o.groupby(['section_id','old_final_label','display_label'],observed=True).size().rename('n').reset_index().to_csv(ROOT/'metadata/label_transition_v2_to_v21.tsv',sep='\t',index=False)
print('REVIEW_TABLE_COMPLETE',len(r),'rows',r.loc[r.scope.eq('joint')].groupby('lineage').size().to_dict(),flush=True)
