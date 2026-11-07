from common import *

def main():
 reg=read('endpoints.tsv.gz');c=read('combination_associations.tsv.gz');ds=read('dataset_endpoint_summary.tsv.gz');u=read('unit_endpoints.tsv.gz');unitref=combo_id([d for d in DS if d!='GSE285715']);allref=combo_id(DS)
 specs=[('B1_R0','HALLMARK_ANDROGEN_RESPONSE','none'),('B3_Epi_distance','HALLMARK_UNFOLDED_PROTEIN_RESPONSE','none'),('B3_Epi_distance','HALLMARK_MTORC1_SIGNALING','none'),('B3_tissue_distance','HALLMARK_MYC_TARGETS_V1','none'),('C1','HALLMARK_TNFA_SIGNALING_VIA_NFKB','Classical Mono'),('C1','HALLMARK_OXIDATIVE_PHOSPHORYLATION','Neutrophil'),('C2','HALLMARK_OXIDATIVE_PHOSPHORYLATION','Plasma'),('C3','distance','Treg')]
 selected=[]
 for mo,p,ca in specs:
  rr=reg[reg.module.eq(mo)&reg.canonical_program.eq(p)&reg.candidate_type.eq(ca)&reg.signature.isin(['full','not_applicable'])]
  if mo in ['C1','C2']:rr=rr[rr.score_method.eq('AddModuleScore_SCT_data')&rr.range_type.eq('common_support_0_1000')&rr.statistic.eq('r' if mo=='C1' else 'paired_delta_r')]
  elif mo=='C3':rr=rr[rr.statistic.eq('mean_difference_um')]
  else:rr=rr[rr.statistic.eq('r')]
  if len(rr)!=1:raise ValueError((mo,p,ca,len(rr)))
  eid=int(rr.endpoint_id.iloc[0]);v=c[c.endpoint_id.eq(eid)&c.n_min_rule.eq(10 if mo in ['C1','C2'] else 0)].copy();best=v[v.n_datasets.ge(2)&v.n_analysis_units.ge(15)&v.n_section_only.eq(0)].copy();best['selection_score']=best.estimate.abs()*best.unit_direction_fraction*best.dataset_direction_fraction*np.sqrt(best.n_analysis_units/best.n_analysis_units.max());best=best.sort_values('selection_score',ascending=False).head(1)
  for label,z in [('all_seven_descriptive',v[v.combination.eq(allref)]),('all_identified_units',v[v.combination.eq(unitref)]),('exploratory_selected_subset',best)]:
   if len(z):selected.append(z.iloc[0].to_dict()|dict(selection=label))
 selected=pd.DataFrame(selected).merge(reg,on='endpoint_id',validate='m:1');save(selected,'discussion_relations.tsv')
 # Full/noMP10, raw/partial, method, range and n>=30 counterparts stay distinct.
 sens=[]
 for mo,p,ca in specs:
  rr=reg[reg.module.eq(mo)&reg.canonical_program.eq(p)&reg.candidate_type.eq(ca)];q=c[c.endpoint_id.isin(rr.endpoint_id)&c.combination.isin([unitref,allref])].merge(rr,on='endpoint_id');sens.append(q)
 save(pd.concat(sens,ignore_index=True),'discussion_sensitivity.tsv.gz')
 # C4 low-query sensitivity explicitly uses total Epi n>=20, unrelated to MP10-P95 count.
 rr=reg[reg.module.eq('C4')];a=u[u.endpoint_id.isin(rr.endpoint_id)&u.n_min.ge(20)&u.value.notna()];rows=[]
 for eid,z in a.groupby('endpoint_id'):
  z=z.copy();z['z']=transform(z.value,'r');m=z.groupby('dataset_id').z.mean();vv=z.groupby('dataset_id').z.var();nn=z.groupby('dataset_id').size();known=z.groupby('dataset_id').identity_level.apply(lambda x:(x=='analysis_unit').all())
  for comb in combos(list(m.index)):
   k=len(comb);mu=m.loc[list(comb)].mean();ve=(vv.loc[list(comb)]/nn.loc[list(comb)]).to_numpy();df=ve.sum()**2/np.sum(ve**2/(nn.loc[list(comb)].to_numpy()-1));se=np.sqrt(ve.sum())/k;valid=known.loc[list(comb)].all() and np.isfinite(se) and se>0;lo=mu-stats.t.ppf(.975,df)*se if valid else np.nan;hi=mu+stats.t.ppf(.975,df)*se if valid else np.nan;pv=2*stats.t.sf(abs(mu/se),df) if valid else np.nan
   rows.append(dict(endpoint_id=eid,combination=combo_id(comb),total_Epi_min_n=20,n_records=int(nn.loc[list(comb)].sum()),estimate=np.tanh(mu),low=np.tanh(lo),high=np.tanh(hi),p_conditional=pv))
 out=pd.DataFrame(rows);out['q_across_C4_sensitivity_search']=bh(out.p_conditional);save(out,'C4_total_Epi_n20_sensitivity.tsv.gz')
 # Source q and spatial intervals remain available at their exact source files; no repurposing.
 idx=[dict(module='A/B',path=str(S/'tables/context85_effects.parquet'),columns='analysis,origin,r,partial_r,n,blocks; old beta P/CI not reused'),dict(module='C1',path=str(S/'tables/pearson_all.parquet'),columns='pearson_r,spatial_fisher_z_se,n_spots,partial_r_epi_depth,range_type,source q'),dict(module='C2',path=str(S/'tables/delta_all.parquet'),columns='r_P95,r_nonP95,delta_r,joint spatial CI,source spatial q'),dict(module='C3',path=str(S/'tables/closeness_all.tsv.gz'),columns='mean/median/zero/matched differences; matched source descriptive'),dict(module='C4',path=str(S/'tables/MP10_candidate_all.tsv.gz'),columns='all frozen Epi Pearson and partial r'),dict(module='niche',path=str(S/'tables/niche_counts_all.tsv.gz'),columns='scope,population,n_labelled,n_unlabelled,n_niche,fraction'),dict(module='state',path=str(OLD/'tables/<section>_rings.tsv.gz'),columns='R0 Epi P95/control means; new5_reference_rings appended'),dict(module='composition',path=str(R/'mp10_paper_distance_full_v1/tables/region_composition.tsv'),columns='disjoint lesion_location region mean fractions weighted by n_valid')];save(idx,'input_index.tsv')
 print('discussion selection complete',flush=True)
if __name__=='__main__':main()
