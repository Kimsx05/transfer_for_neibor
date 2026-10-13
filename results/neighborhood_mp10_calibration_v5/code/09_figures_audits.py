from common import *
import matplotlib.pyplot as plt

def main(out):
 G=out/'02_gene_evidence';P=out/'04_fixed_niche_contrasts';F=out/'06_figures';z=pd.read_csv(G/'gene_by_section_niche.tsv.gz',sep='\t');gene=pd.read_csv(G/'gene_direct_contrasts.tsv',sep='\t');loo=pd.read_csv(G/'leave_one_gene_out_unit_effects.tsv',sep='\t');gl=gene.groupby(['pair','layer','gene','analysis_unit_uid']).auc_mean_difference.mean().reset_index();gs=gl.groupby(['pair','layer','gene']).agg(mean_unit_difference=('auc_mean_difference','mean'),positive_unit_fraction=('auc_mean_difference',lambda x:(x>0).mean()),n_units=('analysis_unit_uid','nunique')).reset_index();write_tsv(out,G/'gene_direction_unit_summary.tsv',gs);ls=loo.groupby(['pair','layer','removed_gene'])[['original_effect','renormalized19_effect','fixed_denominator_effect']].mean().reset_index();ls['sign_changed']=ls.original_effect*ls.renormalized19_effect<0;write_tsv(out,G/'leave_one_gene_out_summary.tsv',ls)
 for dataset,g in z.groupby('dataset_id'):
  for col in ['raw_detection_rate','raw_mean_per10k','epi_mean_counts','observed_auc_mean','epi_auc_mean']:
   h=g.groupby(['gene','niche_name'])[col].mean().unstack().reindex(index=TOP,columns=[f'niche{i}' for i in range(9)]);fig,ax=plt.subplots(figsize=(9,8));im=ax.imshow(h,aspect='auto',cmap='viridis');ax.set_yticks(range(20),TOP);ax.set_xticks(range(9),[f'n{i}' for i in range(9)]);ax.set_title(dataset+' '+col+'\nEqual section summaries; NA kept');fig.colorbar(im,ax=ax);fig.tight_layout()
   for ext in ['png','pdf']:fig.savefig(safe(out,F/dataset/(col+'.'+ext)),transparent=True,dpi=140)
   plt.close(fig)
 for layer,g in gs.groupby('layer'):
  h=g.pivot(index='gene',columns='pair',values='mean_unit_difference').reindex(TOP);fig,ax=plt.subplots(figsize=(8,8));mx=np.nanmax(np.abs(h));im=ax.imshow(h,aspect='auto',cmap='coolwarm',vmin=-mx,vmax=mx);ax.set_yticks(range(20),TOP);ax.set_xticks(range(3),h.columns,rotation=15);ax.set_title(layer+' per-gene AUC contribution direct contrasts');fig.colorbar(im,ax=ax);fig.tight_layout()
  for ext in ['png','pdf']:fig.savefig(safe(out,F/('gene_direction_'+layer+'.'+ext)),transparent=True,dpi=150)
  plt.close(fig)
 for layer,g in ls.groupby('layer'):
  fig,axes=plt.subplots(1,3,figsize=(13,7),sharey=True)
  for ax,(pair,gg) in zip(axes,g.groupby('pair')):
   gg=gg.set_index('removed_gene').reindex(TOP);ax.scatter(gg.renormalized19_effect,range(20));ax.axvline(gg.original_effect.iloc[0],c='tab:orange',label='formal20');ax.axvline(0,c='gray');ax.set_yticks(range(20),TOP);ax.set_title(pair,fontsize=9);ax.set_xlabel('mean direct effect after removing gene')
  axes[0].legend();fig.suptitle(layer+' fixed-background leave-one-out, all20');fig.tight_layout()
  for ext in ['png','pdf']:fig.savefig(safe(out,F/('leave_one_out_'+layer+'.'+ext)),transparent=True,dpi=150)
  plt.close(fig)
 # Tie strategy summary by predefined finite comparison types.
 f=out/'01_score_fidelity';t=pd.read_csv(f/'tie_and_implementation_audit.tsv',sep='\t');rs=[]
 for typ,g in t.groupby('input_type'):
  spread=g.groupby('cell_id').official_score.agg(lambda x:x.max()-x.min());rs.append(dict(input_type=typ,n_inputs=g.cell_id.nunique(),seeds='123,124,125',official_version=g.AUCell_version.iloc[0],mean_official_minus_historical=g.difference.mean(),max_abs_difference=g.difference.abs().max(),max_seed_spread=spread.max(),official_ties='random',historical_ties='average',decision='formal scores remain unchanged'))
 write_tsv(out,f/'implementation_sensitivity_summary.tsv',rs)
 m=pd.read_parquet(f/'cell_level_fidelity.parquet');u=pd.read_csv(f/'fidelity_by_patient_dataset.tsv',sep='\t');fs=json.loads((f/'fidelity_summary.json').read_text());fs['n_calculable_units']=int(u.status.eq('CALCULABLE').sum());fs['n_insufficient_or_constant']=int(u.status.ne('CALCULABLE').sum());fs['within_unit_spearman_min']=u.spearman.min();fs['within_unit_spearman_max']=u.spearman.max();fs['median_technical_partial_rank_correlation']=u.technical_partial_rank_correlation.median();fs['median_score_depth_correlation']=u.score_depth_spearman.median();write_json(out,f/'fidelity_summary.json',fs)
 # Retain technical-filter losses for all nine (same existing definitions).
 cv=[]
 for n,g in pd.read_parquet(G/'spatial_technical_metadata.parquet').query('accepted_label>=0').groupby('accepted_label'):
  for name,sel in [('MAIN',np.ones(len(g),bool)),('EPI_GE_0.5',g.epi_fraction_center_mean.ge(.5)),('EPI_GE_0.7',g.epi_fraction_center_mean.ge(.7)),('COMPLETE_RINGS',g.complete_rings)]:
   h=g[sel];cv.append(dict(niche_name=f'niche{n}',analysis_set=name,n=len(h),retained_fraction=len(h)/len(g),units=h.analysis_unit_uid.nunique(),sections=h.section_id.nunique()))
 write_tsv(out,P/'all9_sensitivity_retention.tsv',cv)
 # Model hashes checked against v4's before-v5 audit, not just a post-run snapshot.
 old=pd.read_csv(V4/'00_audit_and_freeze/model_checkpoint_hashes.tsv',sep='\t');old['current_sha256']=old.path.map(sha);assert old.sha256.eq(old.current_sha256).all();write_tsv(out,out/'00_audit_and_spec/real_model_integrity.tsv',old)
 checks=pd.read_csv(out/'00_audit_and_spec/input_checksums.tsv',sep='\t');checks['current_sha256']=checks.path.map(sha);assert checks.sha256.eq(checks.current_sha256).all();write_tsv(out,out/'00_audit_and_spec/input_integrity.tsv',checks);write_json(out,out/'00_audit_and_spec/specification_hash.json',dict(sha256=sha(out/'00_audit_and_spec/frozen_analysis_spec.yaml')));log('gene figures, all9 retention and immutable input audits complete')
if __name__=='__main__':main(OUT)
