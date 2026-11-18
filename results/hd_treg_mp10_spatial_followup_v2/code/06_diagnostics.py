from fu_common import *
from statsmodels.stats.multitest import multipletests
# One unified auxiliary association family, including newly computed step B diagnostic.
r=pd.read_csv(R/'tables/paired_effects.tsv',sep='\t');b=pd.read_csv(R/'tables/stepwise_decomposition.tsv',sep='\t');reg=[]
for q in r.itertuples():
 if np.isfinite(q.P_delta):reg.append(dict(test=q.variant+f'_block{q.block_um}',family='primary'if q.test_family=='single_primary_delta'else'auxiliary_all',P=q.P_delta,replicates=999,null='paired standardized outcome difference restricted-null spatial wild bootstrap'))
q=b[b.step.eq('B_common_centers_full')].iloc[0];reg.append(dict(test='B_common_centers_original_bin_full',family='auxiliary_all',P=q.P,replicates=499,null='v1 restricted-null spatial wild bootstrap; diagnostic, not primary'))
reg=pd.DataFrame(reg);reg['BH_q_aux']=np.nan;ii=reg.index[reg.family.eq('auxiliary_all')&reg.P.notna()];reg.loc[ii,'BH_q_aux']=multipletests(reg.loc[ii,'P'],method='fdr_bh')[1];put(reg,'new_test_registry.tsv')
for i,q in r.iterrows():
 key=q['variant']+f'_block{q.block_um}';z=reg[reg.test.eq(key)]
 if len(z):r.loc[i,'BH_q_auxiliary']=z.BH_q_aux.iloc[0]
put(r,'paired_effects.tsv')
# Signature overlap and actual predictor collinearity, not effect-driven exclusions.
rows=[]
for a in ['Treg_main','T_lineage','Tfh','Activated_T']:
 for b in ['Treg_main','T_lineage','Tfh','Activated_T']:
  ga=set(cfg['signatures'][a]);gb=set(cfg['signatures'][b]);d=pd.read_parquet(R/'tables/Cervilla_common_centers_r200.parquet');rows.append(dict(signature_a=a,signature_b=b,overlap_n=len(ga&gb),overlap_genes=';'.join(sorted(ga&gb)),score_rho=spearmanr(d[a],d[b]).statistic))
put(pd.DataFrame(rows),'signature_overlap_collinearity.tsv')
# Profile-level activation diagnostics without changing a single reference cell, threshold or split.
parent=pd.read_parquet(R/'tables/mixture_parent_membership.parquet');rs=pd.read_parquet(V/'tables/reference_scores_full.parquet').set_index('cell_id');cp=parent[parent.role.isin(['CD4 Helper/Conventional','Tfh/CXCL13-like','CD8 Exhausted'])].copy();cp['activation']=rs.loc[cp.cell_id,'Activated_T'].to_numpy();act=cp.groupby('profile_id').activation.mean();mix=pd.read_parquet(R/'tables/mixture_rule_results.parquet');mix['competitor_activation']=mix.profile_id.map(act);cut=mix[(mix.partition=='train')&mix.competitor.eq('CD4 Helper/Conventional')].drop_duplicates('profile_id').competitor_activation.quantile(.75);mix['CD4_high_activation']=mix.competitor.eq('CD4 Helper/Conventional')&(mix.competitor_activation>=cut);rows=[]
for key,v in mix[mix.competitor.eq('CD4 Helper/Conventional')].groupby(['partition','CD4_high_activation','depth']):
 neg=v[v.Treg_RNA_fraction==0];pos=v[v.Treg_RNA_fraction>0];rows.append(dict(partition=key[0],CD4_high_activation=key[1],depth=key[2],activation_cut_training_q75=cut,negative_n=len(neg),positive_n=len(pos),old_FPR=neg.joint_old.mean(),candidate_FPR=neg.joint_candidate.mean(),old_sensitivity=pos.joint_old.mean(),candidate_sensitivity=pos.joint_candidate.mean(),interpretation='activation-score stratum of fixed conventional CD4 mixtures, not a new cell annotation or new threshold selection'))
put(pd.DataFrame(rows),'CD4_activation_competitor_diagnostic.tsv')
# Flow and zero-score diagnostics among qualified HD regions; all5 remain represented.
flow=[]
for s in SECTIONS:
 t=pd.read_parquet(V/s/'T32_regions.parquet');v=t[t.measurement_eligible];n=len(v);cond={'multigene':v.Treg_main_detected>=3,'core':v.Treg_core_detected>=1,'positive_rank_score':v.Treg_main>0,'old_high':v.Treg_main>cfg['Treg_like_threshold_training_competitor_q95']};keep=np.ones(n,bool)
 for name,mask in cond.items():
  flow.append(dict(section=s,qualified_n=n,stage=name,first_loss=int((keep&~mask.to_numpy()).sum()),failure_alone=int((~mask).sum()),survives=int((keep&mask.to_numpy()).sum())));keep&=mask.to_numpy()
put(pd.DataFrame(flow),'qualified_HD_flow.tsv');print('DIAGNOSTICS_DONE')
