from base import *
import matplotlib.pyplot as plt
from patsy import dmatrix
P=OUT/'04_niche_comparison';F=OUT/'06_figures';rng=np.random.default_rng(20261013)
SCORES=['spot_MP10_frozen','spot_MP10_model_universe','epi_MP10_model_universe']
def effects(d,score,label='accepted_label'):
 rows=[]
 for section,g in d.groupby('section_id',observed=True):
  g=g[g[score].notna()];
  for n in range(9):
   a=g.loc[g[label].eq(n),score];b=g.loc[g[label].ne(n),score];r=g.iloc[0] if len(g) else None
   rows.append(dict(section_id=section,dataset_id=r.dataset_id if r is not None else '',analysis_unit_uid=r.analysis_unit_uid if r is not None else '',patient_uid=r.patient_uid if r is not None else '',niche_name=f'niche{n}',score=score,n=len(a),n_rest=len(b),median=a.median(),rest_median=b.median(),effect=a.median()-b.median(),stable_coverage=len(a)>=20 and len(b)>=20))
 return pd.DataFrame(rows)
def savefig(fig,name,folder=F):
 folder.mkdir(parents=True,exist_ok=True);fig.tight_layout();fig.savefig(folder/(name+'.png'),dpi=150,transparent=True);fig.savefig(folder/(name+'.pdf'),transparent=True);plt.close(fig)
def main():
 d=pd.read_parquet(P/'anchor_MP10_all_layers.parquet');a=d[d.accepted_label.ge(0)].copy();matched=a[a.epi_score_valid&a.spot_MP10_model_universe.notna()].copy();sets={'ALL_ACCEPTED':a,'FULL20_SPOT_PRIMARY':a[a.n_signature_available_raw.eq(20)],'MATCHED_FULL20_EPI':matched,'PARTIAL_SIGNATURE_SENSITIVITY':a[a.n_signature_available_raw.lt(20)|a.n_signature_available_model.lt(20)],'PATHOLOGY_SUPPORTED_TUMOR':a[a.pathology_label.eq('tumor_epithelium')],'PATHOLOGY_SUPPORTED_NON_TUMOR':a[a.pathology_label.eq('non_tumor_epithelium')],'COMPLETE_RINGS':matched[matched.complete_rings],'EPI_GE_0.5':matched[matched.epi_fraction_center_mean.ge(.5)],'EPI_GE_0.7':matched[matched.epi_fraction_center_mean.ge(.7)],'ALL_CANDIDATE':d[d.candidate_label.ge(0)],'LEAVE_LARGEST_DATASET_OUT':matched[matched.dataset_id.ne(a.dataset_id.value_counts().index[0])]}
 cov=[];summ=[]
 for name,s in sets.items():
  label='candidate_label' if name=='ALL_CANDIDATE' else 'accepted_label'
  for n in range(9):
   g=s[s[label].eq(n)];cov.append(dict(analysis_set=name,niche_name=f'niche{n}',n=len(g),sections=g.section_id.nunique(),known_patients=g.loc[g.source_patient_identity_known,'patient_uid'].nunique(),analysis_units=g.analysis_unit_uid.nunique(),n_matched=int(g.epi_score_valid.sum())))
   for score in SCORES+(['spot_MP10_model_partial','epi_MP10_model_partial'] if name=='PARTIAL_SIGNATURE_SENSITIVITY' else []):
    v=g[score].dropna();summ.append(dict(analysis_set=name,niche_name=f'niche{n}',score=score,n=len(v),mean=v.mean(),median=v.median(),q25=v.quantile(.25),q75=v.quantile(.75),zero_fraction=v.eq(0).mean() if len(v) else np.nan,n_sections=g.loc[v.index,'section_id'].nunique(),n_known_patients=g.loc[v.index].loc[lambda z:z.source_patient_identity_known,'patient_uid'].nunique(),n_units=g.loc[v.index,'analysis_unit_uid'].nunique(),median_epi_fraction=g.epi_fraction_center_mean.median(),median_nCount=g.center_nCount.median()))
 tab(P/'coverage_sensitivity.tsv',cov);tab(P/'niche_MP10_summary_all9.tsv',summ)
 eff=[]
 for score in SCORES:
  s=sets['FULL20_SPOT_PRIMARY'] if score==SCORES[0] else matched
  e=effects(s,score);e['analysis_set']='FULL20_SPOT_PRIMARY' if score==SCORES[0] else 'MATCHED_FULL20_EPI';eff.append(e)
 for name in ['COMPLETE_RINGS','EPI_GE_0.5','EPI_GE_0.7','ALL_CANDIDATE','LEAVE_LARGEST_DATASET_OUT']:
  for score in SCORES[1:]:
   e=effects(sets[name],score,'candidate_label' if name=='ALL_CANDIDATE' else 'accepted_label');e['analysis_set']=name;eff.append(e)
 ef=pd.concat(eff,ignore_index=True);tab(P/'section_niche_score_summary.tsv',ef)
 u=ef[ef.stable_coverage].groupby(['analysis_set','score','dataset_id','analysis_unit_uid','niche_name'],observed=True).agg(effect=('effect','mean'),n_sections=('section_id','nunique')).reset_index();tab(P/'independent_unit_effects.tsv',u)
 rows=[]
 for keys,g in u.groupby(['analysis_set','score','dataset_id','niche_name'],observed=True):
  rows.append(dict(zip(['analysis_set','score','dataset_id','niche_name'],keys),n_units=len(g),mean_effect=g.effect.mean(),median_effect=g.effect.median(),positive_fraction=g.effect.gt(0).mean()))
 tab(P/'dataset_niche_effects.tsv',rows)
 boot=[]
 for (aset,score,n),g in u.groupby(['analysis_set','score','niche_name'],observed=True):
  x=g.effect.to_numpy();ci=np.quantile(rng.choice(x,(1000,len(x)),replace=True).mean(1),[.025,.975]) if len(x)>=5 else [np.nan,np.nan];boot.append(dict(analysis_set=aset,score=score,niche_name=n,n_units=len(x),mean_effect=x.mean(),ci_low=ci[0],ci_high=ci[1],p_value=np.nan,interval_status='exploratory_unit_bootstrap' if len(x)>=5 else 'INSUFFICIENT_UNITS'))
 tab(P/'unit_bootstrap_effects.tsv',boot)
 # Predefined adjustment: coefficients only, no spot-based SE or P-values.
 adjust=[]
 for score in SCORES:
  s=(sets['FULL20_SPOT_PRIMARY'] if score==SCORES[0] else matched).dropna(subset=[score,'center_nCount','epi_fraction_center_mean']).copy()
  try:
   design=dmatrix('C(section_id)+C(accepted_label)+bs(epi_fraction_center_mean,df=3)+np.log1p(center_nCount)',s,return_type='dataframe');coef,resid,rank,_=np.linalg.lstsq(design.to_numpy(),s[score],rcond=None)
   for k,v in zip(design.columns,coef):
    if 'accepted_label' in k:adjust.append(dict(score=score,term=k,coefficient_vs_niche0=v,n=len(s),rank=rank,n_columns=design.shape[1],status='DESCRIPTIVE_COEFFICIENT_NO_IID_INFERENCE' if rank==design.shape[1] else 'RANK_DEFICIENT_NOT_ESTIMABLE'))
  except Exception as e:adjust.append(dict(score=score,status='NOT_ESTIMABLE',reason=repr(e)))
 tab(P/'technical_mixture_adjustment.tsv',adjust)
 shares=['signature_epi_count_share_of_observed','signature_nonepi_count_share_of_observed','signature_ambient_count_share_of_observed','signature_observed_counts','allgene_epi_count_share_of_observed']
 sh=a.groupby('niche_name')[shares].median().reindex([f'niche{i}' for i in range(9)]);tab(P/'signature_counts_attribution.tsv',sh.reset_index())
 diag=[]
 for sec,g in a.groupby('section_id'):
  x=g.epi_MP10_model_partial.dropna();diag.append(dict(section_id=sec,n=len(x),unique_scores=x.nunique(),std=x.std(),template_score=g.template_only_MP10.iloc[0],zero_fraction=x.eq(0).mean(),near_template_fraction=np.isclose(x,g.template_only_MP10.iloc[0],atol=1e-6).mean() if len(x) else np.nan,low_detected_gene_fraction=g.epi_detected_below_rank_cutoff.mean()))
 tab(OUT/'03_epi_attribution/template_and_information_diagnostics.tsv',diag)
 for dataset,s in [('ALL_DATASETS',a)]+list(a.groupby('dataset_id')):
  folder=F/dataset;folder.mkdir(exist_ok=True)
  for score in SCORES:
   q=s[s.n_signature_available_raw.eq(20)] if score==SCORES[0] else s[s.epi_score_valid];fig,ax=plt.subplots(figsize=(11,4));vals=[q.loc[q.accepted_label.eq(n),score].dropna().to_numpy() for n in range(9)];ax.boxplot([x if len(x) else [np.nan] for x in vals],positions=np.arange(9),showfliers=False,whis=(10,90));med=q.groupby(['section_id','accepted_label'])[score].median().reset_index();ax.scatter(med.accepted_label,med[score],s=7,alpha=.35,color='tab:blue');ax.set_xticks(range(9),[f'niche{n}' for n in range(9)]);ax.set_ylabel('MP10 AUCell');ax.set_title(f'{dataset}: {score}\nFull summaries; dots=section medians; boxes=anchors');savefig(fig,score,folder)
  sm=s[s.epi_score_valid].groupby('accepted_label')[SCORES[1:]].median().reindex(range(9));fig,ax=plt.subplots(figsize=(10,4));ax.plot(range(9),sm.iloc[:,0],'o-',label='observed, model genes');ax.plot(range(9),sm.iloc[:,1],'o-',label='model-attributed Epi, same centers/genes');ax.set_xticks(range(9),[f'niche{i}' for i in range(9)]);ax.legend();ax.set_title(dataset+' matched-universe medians (not source fractions)');savefig(fig,'matched_layers',folder)
  subef=ef[ef.analysis_set.eq('MATCHED_FULL20_EPI')&(ef.dataset_id.eq(dataset) if dataset!='ALL_DATASETS' else True)]
  fig,axes=plt.subplots(1,2,figsize=(13,5),sharey=True)
  for ax,sc in zip(axes,SCORES[1:]):
   z=subef[subef.score.eq(sc)&subef.stable_coverage];ax.scatter(z.niche_name.str.replace('niche','').astype(int),z.effect,s=12,alpha=.4);ax.axhline(0,c='gray');ax.set_xticks(range(9),[f'n{i}' for i in range(9)]);ax.set_title(sc);ax.set_ylabel('within-section niche minus rest median')
  savefig(fig,'section_effects',folder)
  pp=P/('all_datasets' if dataset=='ALL_DATASETS' else 'by_dataset/'+dataset);pp.mkdir(parents=True,exist_ok=True);tab(pp/'niche_summary.tsv',pd.DataFrame(summ) if dataset=='ALL_DATASETS' else s.groupby('niche_name')[SCORES].agg(['count','median','mean']).reset_index());tab(pp/'section_effects.tsv',subef)
 for score in SCORES:
  z=ef[ef.score.eq(score)&ef.analysis_set.isin(['FULL20_SPOT_PRIMARY','MATCHED_FULL20_EPI'])]
  for val in ['median','effect']:
   h=z.pivot(index='section_id',columns='niche_name',values=val).reindex(columns=[f'niche{i}' for i in range(9)]);fig,ax=plt.subplots(figsize=(10,13));im=ax.imshow(h,aspect='auto',cmap='plasma' if val=='median' else 'coolwarm');ax.set_xticks(range(9),h.columns);ax.set_yticks(range(len(h)),h.index,fontsize=6);fig.colorbar(im,ax=ax,label=score+' '+val);savefig(fig,score+'_'+val+'_heatmap')
 fig,ax=plt.subplots(figsize=(11,4));sh.iloc[:,:3].plot.bar(stacked=True,ax=ax);ax.set_ylabel('Median estimated signature RNA count shares');ax.set_title('Component medians need not sum exactly to one; zero-count rows excluded');savefig(fig,'signature_count_shares')
 log('Nine-niche comparisons and plots completed')
if __name__=='__main__':main()
