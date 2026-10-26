from pathlib import Path
import sys
P=Path('/home/data/t070721/codex_workspace/Bladder Metabolism');O=P/'mp10_subtype_screen_round1';T=O/'tables'
sys.path=[p for p in sys.path if 'dist-packages' not in p];sys.path[:0]=[str(P/'spatial_ecology_top20_v2/runtime'),str(P/'cell2location_setup_v1/cell2location_env/lib/python3.12/site-packages')]
import numpy as np,pandas as pd
axes=pd.read_csv(T/'candidate_union_axes.tsv',sep='\t');su=pd.read_csv(T/'candidate_union_sample_expression_support.tsv.gz',sep='\t');net=pd.read_csv(T/'CellChat_epithelial_comparisons.tsv.gz',sep='\t');props=pd.read_csv(T/'proportion_models.tsv',sep='\t');rows=[];details=[]
for _,a in axes.iterrows():
 z=su[su.LR_ID.eq(a.LR_ID)&su.source.eq(a.source)&su.target.eq(a.target)&su.MP10_group.isin(['low','high'])].copy();nn=net[net.LR_ID.eq(a.LR_ID)&net.source.eq(a.source)&net.target.eq(a.target)&net.comparison.eq('MP10_low_high')];z['expression_product']=z.ligand_mean*z.receptor_mean
 for ds,b in z.groupby('dataset'):
  n=nn[nn.scope.eq(ds)];r={'candidate_subtype':a.candidate_subtype,'LR_ID':a.LR_ID,'source':a.source,'target':a.target,'dataset':ds};r['network_difference']=n.difference.iloc[0] if len(n) else np.nan
  for g in ['low','high']:
   c=b[b.MP10_group.eq(g)];v=c[c.evaluable];support=v[v.expression_supported]
   r.update({f'{g}_n_specimens':len(c),f'{g}_n_evaluable':len(v),f'{g}_n_evaluable_patient_ids':v.patient_uid.nunique(),f'{g}_n_supported_specimens':len(support),f'{g}_n_supported_patient_ids':support.patient_uid.nunique(),f'{g}_ligand_mean':v.ligand_mean.mean(),f'{g}_receptor_mean':v.receptor_mean.mean(),f'{g}_expression_product':v.expression_product.mean()})
  r['replicated_expression_support']=r['low_n_supported_patient_ids']>=2 and r['high_n_supported_patient_ids']>=2
  r['expression_product_difference']=r['high_expression_product']-r['low_expression_product']
  signs=[]
  for pat in b.patient_uid.unique():
   c=b[b.patient_uid.ne(pat)&b.evaluable];avg=c.groupby('MP10_group').expression_product.mean()
   if {'low','high'}.issubset(avg.index):signs.append(np.sign(avg['high']-avg['low']))
  r['leave_one_recorded_patient_out_same_expression_direction_fraction']=np.mean(np.array(signs)==np.sign(r['expression_product_difference'])) if signs else np.nan
  details.append(r)
 d=pd.DataFrame([v for v in details if v['LR_ID']==a.LR_ID and v['source']==a.source and v['target']==a.target]);rep=d[d.replicated_expression_support&d.network_difference.notna()]
 o=a.to_dict();o.update(n_datasets_replicated_expression_support=len(rep),replicated_datasets=';'.join(rep.dataset),replicated_network_up=int(rep.network_difference.gt(0).sum()),replicated_network_down=int(rep.network_difference.lt(0).sum()),replicated_expression_up=int(rep.expression_product_difference.gt(0).sum()),replicated_expression_down=int(rep.expression_product_difference.lt(0).sum()),replicated_network_expression_same_sign=int((np.sign(rep.network_difference)==np.sign(rep.expression_product_difference)).sum()),min_leave_one_patient_expression_consistency=rep.leave_one_recorded_patient_out_same_expression_direction_fraction.min())
 for comp in ['Normal_Tumor','T1_T2','MP10_low_high','MP10_continuous','MP10_stage_complete_unadjusted','MP10_stage_adjusted']:
  pr=props[props.scope.eq('ALL')&props.subtype.eq(a.candidate_subtype)&props.level.eq('subtype')&props.denominator_type.eq('parent')&props.comparison.eq(comp)]
  if len(pr):
   for key in ['coefficient','CI_low','CI_high','P','q','raw_difference_pp','n_specimens','n_patient_ids_recorded','status']:o[comp+'_'+key]=pr.iloc[0][key]
 rows.append(o)
pd.DataFrame(details).to_csv(T/'candidate_axis_dataset_support.tsv',sep='\t',index=False,na_rep='NA')
r=pd.DataFrame(rows);r['replicated_direction_max']=r[['replicated_network_up','replicated_network_down']].max(axis=1)
r=r.sort_values(['n_datasets_replicated_expression_support','replicated_direction_max','replicated_network_expression_same_sign','median_abs_network_difference'],ascending=False)
r.to_csv(T/'candidate_axis_evidence.tsv',sep='\t',index=False,na_rep='NA')
# Every subtype retains its full abundance evidence, irrespective of shortlist selection.
mp=pd.read_csv(T/'subtype_parent.tsv',sep='\t');ab=props[props.level.eq('subtype')]
wide=ab.pivot(index=['subtype','major'],columns=['scope','denominator_type','comparison'],values=['coefficient','P','q','n_specimens','n_patient_ids_recorded','status']);wide.columns=['|'.join(v) for v in wide.columns];wide.reset_index().to_csv(T/'candidate_all_subtypes_abundance_evidence.tsv',sep='\t',index=False,na_rep='NA')
print(r[r.priority_interpretable][['candidate_subtype','source','target','LR_ID','n_datasets_replicated_expression_support','replicated_datasets','replicated_network_up','replicated_network_down','replicated_network_expression_same_sign','min_leave_one_patient_expression_consistency','MP10_continuous_q']].head(30).to_string(index=False))
