from common import *
import anndata as ad
import matplotlib.pyplot as plt

def main(out):
 P=out/'03_attribution_calibration';axis=json.loads((P/'panel_axes.json').read_text());genes=axis['genes'];gi=pd.Index(genes).get_indexer(TOP);K=int(np.ceil(.05*len(genes)));step=1/sum(K-i for i in range(1,21));Y=sparse.load_npz(P/'panel_observed.npz').tocsr();E=sparse.load_npz(P/'panel_true_epi.npz').tocsr();W=np.load(P/'panel_known_weights.npy');md=pd.read_parquet(P/'panel_metadata.parquet');C=pd.read_csv(P/'frozen_reference.tsv',sep='\t',index_col=0);known=pd.read_parquet(P/'known_mixture_scores.parquet');known['rank_step']=step;generec=[];fit_path=P/'synthetic_refits/model1/posterior.h5ad';fit=None;abundance=[]
 for branch in ['KNOWN_MIXTURE_FIXED_SIGNATURE','SYNTHETIC_REFIT']:
  if branch=='SYNTHETIC_REFIT':
   if not fit_path.exists():continue
   a=ad.read_h5ad(fit_path);assert a.var_names.tolist()==genes;idx=pd.Index(md.bag_id).get_indexer(a.obs_names);assert (idx>=0).all();s=a.uns['mod']['post_sample_means'];cw=np.asarray(s['w_sf'],float);mg=s['m_g'];amb=np.asarray(s['s_g_gene_add'])[a.obs['_scvi_batch'].to_numpy(int)];status=json.loads((P/'synthetic_refits/model1/fit_status.json').read_text())['status']
  else:idx=np.arange(len(md));cw=W;mg=np.ones(len(genes));amb=np.zeros((len(md),len(genes)));status='CONDITIONAL_KNOWN_MIXTURE'
  sr=[]
  for start in range(0,len(idx),64):
   ii=idx[start:start+64];yy=Y[ii].toarray();ee=E[ii].toarray();ww=cw[start:start+64];eh=allocate(yy,C,ww,mg,amb[start:start+64]);ts,tc,_=auc_parts(ee,genes);hs,hc,_=auc_parts(eh,genes)
   for j,i in enumerate(ii):
    sr.append(dict(bag_id=md.bag_id.iloc[i],true_score=ts[j],estimated_score=hs[j],score_error=hs[j]-ts[j],rank_step=step,branch=branch,fit_status=status,estimated_epi_RNA_fraction=eh[j].sum()/yy[j].sum()))
    if branch=='SYNTHETIC_REFIT':
     truew=W[i]/W[i].sum();hatw=ww[j]/ww[j].sum();abundance.append(dict(bag_id=md.bag_id.iloc[i],signature_weight_fraction_MAE=float(np.abs(truew-hatw).mean()),true_epi_signature_weight_fraction=truew[8],estimated_epi_signature_weight_fraction=hatw[8],true_epi_RNA_fraction=ee[j].sum()/yy[j].sum(),estimated_epi_RNA_fraction=eh[j].sum()/yy[j].sum(),scope='signature weights are not unscaled cell fractions'))
    for gj,gene in enumerate(TOP):generec.append(dict(branch=branch,bag_id=md.bag_id.iloc[i],donor=md.donor.iloc[i],background=md.background.iloc[i],condition=md.condition.iloc[i],depth_band=md.depth_band.iloc[i],gene=gene,true_counts=ee[j,gi[gj]],estimated_counts=eh[j,gi[gj]],observed_counts=yy[j,gi[gj]],true_auc_contribution=tc[j,gj],estimated_auc_contribution=hc[j,gj]))
  scored=md.merge(pd.DataFrame(sr),on='bag_id',how='inner',validate='one_to_one');assert len(scored)==len(idx)
  if branch=='SYNTHETIC_REFIT':fit=scored;fit.to_parquet(safe(out,P/'end_to_end_scores.parquet'),index=False)
  else:known=scored;known.to_parquet(safe(out,P/'known_mixture_scores.parquet'),index=False)
 if abundance:write_tsv(out,P/'synthetic_abundance_recovery.tsv',abundance)
 replay=pd.read_parquet(P/'frozen_replay_scores.parquet');allscore=pd.concat([replay,known]+([fit] if fit is not None else []),ignore_index=True);allscore['depth_band']=allscore.depth_band.fillna(-1);summ=[]
 for keys,g in allscore.groupby(['branch','background','condition','depth_band'],dropna=False):
  z=g.dropna(subset=['true_score','estimated_score']);x=z.true_score.to_numpy();y=z.estimated_score.to_numpy();ok=len(z)>2 and np.std(x)>1e-10;summ.append(dict(zip(['branch','background','condition','depth_band'],keys),n=len(z),n_donors=z.donor.nunique(),spearman=spearmanr(x,y).statistic if ok else np.nan,MAE=np.mean(np.abs(y-x)),bias=np.mean(y-x),slope=np.polyfit(x,y,1)[0] if ok else np.nan,status='CALCULABLE' if ok else 'CONSTANT_OR_INSUFFICIENT'))
 write_tsv(out,P/'calibration_summary_by_background.tsv',summ)
 donor=[]
 for keys,g in allscore.groupby(['branch','donor','background'],dropna=False):
  z=g.dropna(subset=['true_score','estimated_score']);ok=len(z)>2 and z.true_score.std()>1e-10;donor.append(dict(zip(['branch','donor','background'],keys),n=len(z),spearman=spearmanr(z.true_score,z.estimated_score).statistic if ok else np.nan,MAE=(z.estimated_score-z.true_score).abs().mean(),bias=(z.estimated_score-z.true_score).mean()))
 write_tsv(out,P/'calibration_by_donor.tsv',donor);paired=[]
 for (branch,pair),g in allscore[allscore.condition.eq('FIXED_COMPOSITION_STATE')].groupby(['branch','pair_id']):
  if len(g)!=2 or set(g.state)!=set(['low','high']):raise ValueError('paired state axis mismatch '+str(pair))
  l=g[g.state.eq('low')].iloc[0];h=g[g.state.eq('high')].iloc[0];dt=h.true_score-l.true_score;dh=h.estimated_score-l.estimated_score;tol=max(h.rank_step,l.rank_step);paired.append(dict(branch=branch,pair_id=pair,donor=h.donor,background=h.background,depth_band=h.depth_band,rep=h.get('rep',np.nan),Delta_true=dt,Delta_hat=dh,delta_error=dh-dt,strict_sign_reversal=dt*dh<0,near_zero_true=abs(dt)<=tol,near_zero_hat=abs(dh)<=tol,rank_tolerance=tol,robust_rank_sign_reversal=(dt*dh<0 and abs(dt)>tol and abs(dh)>tol),abs_signal_change=abs(dh)-abs(dt),expected_high_but_realized_nonpositive=dt<=0))
 pa=pd.DataFrame(paired);write_tsv(out,P/'paired_state_recovery.tsv',pa)
 # State-pair replicate spread is an additional noise diagnostic, not a post-hoc pass threshold.
 noise=[]
 for keys,g in pa[pa.branch.ne('FROZEN_WEIGHT_REPLAY')].groupby(['branch','donor','background','depth_band']):
  noise.append(dict(zip(['branch','donor','background','depth_band'],keys),n_repeats=len(g),delta_true_repeat_sd=g.Delta_true.std(),delta_hat_repeat_sd=g.Delta_hat.std(),scope='repeat noise displayed; not independent patients'))
 write_tsv(out,P/'paired_repeat_noise.tsv',noise)
 noise_df=pd.DataFrame(noise);truth_noise=noise_df[noise_df.branch.eq('KNOWN_MIXTURE_FIXED_SIGNATURE')].set_index(['donor','background','depth_band']).delta_true_repeat_sd
 pa['true_delta_repeat_sd']=[truth_noise.get((r.donor,r.background,r.depth_band),np.nan) for r in pa.itertuples()]
 pa['near_zero_with_truth_repeat_noise']=np.where(pa.true_delta_repeat_sd.notna(),pa.Delta_true.abs()<=np.maximum(pa.rank_tolerance,pa.true_delta_repeat_sd),pd.NA)
 pa['truth_noise_screened_reversal']=np.where(pa.true_delta_repeat_sd.notna(),pa.robust_rank_sign_reversal & (pa.Delta_true.abs()>np.maximum(pa.rank_tolerance,pa.true_delta_repeat_sd)),pd.NA)
 write_tsv(out,P/'paired_state_recovery.tsv',pa);bias=[]
 for (branch,group,target),g in allscore[allscore.condition.eq('FIXED_STATE_COMPOSITION')].groupby(['branch','fixed_epi_group','background']):
  g=g.sort_values('target_fraction');assert len(g)==2;assert np.isclose(g.true_score.iloc[0],g.true_score.iloc[1],atol=1e-12);dx=g.target_fraction.iloc[1]-g.target_fraction.iloc[0];delta=g.estimated_score.iloc[1]-g.estimated_score.iloc[0];bias.append(dict(branch=branch,group=group,donor=g.donor.iloc[0],background=target,depth_band=g.depth_band.iloc[0],true_score=g.true_score.iloc[0],low_target_fraction=g.target_fraction.iloc[0],high_target_fraction=g.target_fraction.iloc[1],estimated_high_minus_low=delta,error_slope=delta/dx if dx>0 else np.nan,rank_step=step,near_zero=abs(delta)<=step,scope='two-point gradient, not proof of monotonicity over all compositions'))
 write_tsv(out,P/'composition_dependent_bias.tsv',bias);gr=pd.DataFrame(generec);write_tsv(out,P/'gene_source_recovery.tsv',gr.groupby(['branch','background','depth_band','gene'])[['true_counts','estimated_counts','observed_counts','true_auc_contribution','estimated_auc_contribution']].mean().reset_index());gr.to_parquet(safe(out,P/'gene_source_recovery_bags.parquet'),index=False)
 for branch,g in allscore.groupby('branch'):
  fig,axes=plt.subplots(1,3,figsize=(13,4));targets=['GC B','SPP1 TAM','Neutrophil'] if branch!='FROZEN_WEIGHT_REPLAY' else ['actual_v4_center']
  for ax,target in zip(axes,targets):
   z=g[g.background.eq(target)];ax.scatter(z.true_score,z.estimated_score,s=7,alpha=.3);mx=max(z.true_score.max(),z.estimated_score.max());ax.plot([0,mx],[0,mx],c='gray');ax.set(title=target,xlabel='true simulated Epi score',ylabel='estimated Epi score')
  fig.suptitle(branch);fig.tight_layout()
  for ext in ['png','pdf']:fig.savefig(safe(out,out/'06_figures'/(branch+'_recovery.'+ext)),transparent=True,dpi=150)
  plt.close(fig)
 for branch,g in pa.groupby('branch'):
  fig,ax=plt.subplots(figsize=(6,5));ax.scatter(g.Delta_true,g.Delta_hat,s=12,c=np.where(g.near_zero_true,'gray','tab:blue'),alpha=.5);ax.axhline(0,c='gray');ax.axvline(0,c='gray');ax.set(title=branch,xlabel='Delta true Epi high-low',ylabel='Delta estimated Epi high-low');fig.tight_layout()
  for ext in ['png','pdf']:fig.savefig(safe(out,out/'06_figures'/(branch+'_state_delta.'+ext)),transparent=True,dpi=150)
  plt.close(fig)
 bd=pd.DataFrame(bias)
 for branch,g in bd.groupby('branch'):
  fig,axes=plt.subplots(1,3,figsize=(12,4),sharey=True)
  for ax,target in zip(axes,['GC B','SPP1 TAM','Neutrophil']):
   z=g[g.background.eq(target)];ax.scatter(z.depth_band,z.estimated_high_minus_low,s=15);ax.axhline(0,c='gray');ax.set(title=target,xlabel='depth quantile band',ylabel='false score shift, fixed true Epi')
  fig.suptitle(branch);fig.tight_layout()
  for ext in ['png','pdf']:fig.savefig(safe(out,out/'06_figures'/(branch+'_composition_bias.'+ext)),transparent=True,dpi=150)
  plt.close(fig)
 write_json(out,P/'calibration_completion.json',dict(replay_bags=len(replay),known_mixture_bags=len(known),new_refit_bags=len(fit) if fit is not None else 0,new_models=1 if fit is not None else 0,end_to_end_status=json.loads((P/'synthetic_refits/model1/fit_status.json').read_text())['status'] if fit is not None else 'NOT_YET_COMPLETED',donors=known.donor.nunique(),backgrounds=['GC B','SPP1 TAM','Neutrophil'],held_out=False));safe(out,P/'calibration_report.md').write_text('Three distinct recovery branches: frozen real weights replayed on synthetic Y (platform/ambient mismatch stress test); known-mixture fixed reference with library-calibrated signature weights; synthetic c2l fit learns weights from synthetic Y alone. All truth used only after estimation. Five reference-overlapping donor sources; no independent donor validation. Numerical positive control is not donor generalization. Score and gene recovery, two-point composition shifts, state-pair sign reversals and noise are provided separately. No universal pass correlation threshold.\n');log('calibration summaries complete')
if __name__=='__main__':main(OUT)
