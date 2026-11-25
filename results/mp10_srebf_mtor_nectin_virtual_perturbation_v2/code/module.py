from common import *
from itertools import combinations
T=Path('tables');a=rd('cache/jobs/clean/heat.tsv.gz');m=rd(T/'control_matching.tsv');members=rd(T/'program_members.tsv');stab=rd(T/'program_stability.tsv');es=rd(T/'endpoint_stability.tsv')
spec=rd(T/'control_specificity.tsv');rows=[];patterns=[]
classes=dict(zip(K,['SREBF1','SREBF2','mTORC1 component','pan-mTOR','mTORC2 component']))
for name in members.program.unique():
 for mask in ['Full','M20']:
  ids=set(members[(members.program==name)&(members['mask']==mask)].gene)
  common=set.intersection(*[set(a[(a.ko==k)&~a.is_KO].gene) for k in K]) & ids
  for k in K:
   subset=spec[(spec.ko==k)&spec.gene.isin(common)];ctl=a[a.ko.isin(m[m.ko==k].control)&a.gene.isin(common)&~a.is_KO];each=ctl.groupby('ko')
   ok=len(common)>=3 and 2*len(common)>=len(ids)
   r=dict(ko=k,role=classes[k],readout=name,kind='program',mask=mask,n_frozen=len(ids),n_common=len(common),status='EVALUABLE' if ok else 'NE')
   if ok:r.update(median_rank_fraction=subset.rank_fraction.median(),negative_fraction=(subset.directionScore<0).mean(),control_negative_fraction=each.directionScore.apply(lambda x:(x<0).mean()).median(),control_rank_fraction=each.rank_fraction.median().median(),median_direction=subset.directionScore.median())
   y=stab[(stab.ko==k)&(stab.program==name)&(stab['mask']==mask)&(stab.type=='patient_subsample')].iloc[0]
   r.update(refit_rank_fraction_min=y.rank_min,refit_rank_fraction_max=y.rank_max,refit_negative_fraction_min=y.fraction_negative_min,refit_negative_fraction_max=y.fraction_negative_max)
   rows.append(r)
  if len(common)>2:
   for left,right in combinations(K,2):
    x=a[(a.ko==left)&a.gene.isin(common)].set_index('gene');y=a[(a.ko==right)&a.gene.isin(common)].set_index('gene').loc[x.index]
    patterns.append(dict(program=name,mask=mask,left=left,right=right,n_common=len(x),spearman_rank=x.rank_fraction.corr(y.rank_fraction,method='spearman'),same_sign=(np.sign(x.directionScore)==np.sign(y.directionScore)).mean()))
for k in K:
 for g in E+K:
  if g==k:continue
  x=spec[(spec.ko==k)&(spec.gene==g)]
  row=dict(ko=k,role=classes[k],readout=g,kind='gene',mask='Full',status='NE' if x.empty else 'EVALUABLE',selection='v1 exploratory' if g in ['SCD','MARCHF6','RABGEF1'] else 'prespecified')
  if len(x):
   u=x.iloc[0];v=es[(es.ko==k)&(es.endpoint==g)&(es.type=='patient_subsample')].iloc[0]
   row.update(median_rank_fraction=u.rank_fraction,median_direction=u.directionScore,control_rank_fraction=u.control_median_rank,distance=u.distance,q=u['p.adj'],depth_score=u.depth_score,pc1_score=u.pc1_score,n_controls_lower=u.n_controls_lower,refit_rank_fraction_min=v.rank_fraction_min,refit_rank_fraction_max=v.rank_fraction_max,refit_same_direction=v.direction_same_as_clean)
  rows.append(row)
save(pd.DataFrame(rows),T/'mtor_module_summary.tsv');save(pd.DataFrame(patterns),T/'common_program_profiles.tsv')
# Supplement v1 control and candidate provenance and panel list for the exact clean contrast.
coverage=rd(T/'coverage.tsv');rm=rd(T/'removed_genes.tsv');rm['in_P_full']=rm.gene.isin(coverage[coverage.P_full==True].gene);rm['removal_reason']=rm.category;save(rm,T/'removed_genes.tsv')
hist=rd(V/'tables/matched_controls.tsv');hist['label']='旧匹配对照基因';save(hist,T/'legacy_controls.tsv')

