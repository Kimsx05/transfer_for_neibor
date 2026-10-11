from common import *
E=R/'05_evaluation'
if not all((outmap(x)/'MAPPING_DONE.json').exists() for x in ['DEV','DEV_OLD']):raise SystemExit('Paired baseline not yet available; no historical subtraction performed')
a=pd.read_csv(E/'DEV_spot_results.tsv.gz',sep='\t');b=pd.read_csv(E/'DEV_OLD_spot_results.tsv.gz',sep='\t');assert set(a.pseudospot_id)==set(b.pseudospot_id);fa=json.loads((outmap('DEV')/'run_fingerprint.json').read_text());fb=json.loads((outmap('DEV_OLD')/'run_fingerprint.json').read_text());assert fa['counts_sha256']==fb['counts_sha256'] and fa['spot_hash']==fb['spot_hash'] and fa['fit_seed']==fb['fit_seed'] and fa['post_seed']==fb['post_seed'];pan=pd.read_csv(R/'02_pseudospots/panel_references.tsv',sep='\t');ids=set(pan.loc[pan.panel.isin(['A','C'])&pan.mapped,'pseudospot_id']);keys=['pseudospot_id','sample','dataset','source_group','n_K32_true','n_epi_true','n_total_true'];p=a.merge(b,on=keys,suffixes=('_new','_old'),validate='one_to_one');rows=[]
def w(q):
 x=1/q.groupby(['source_group','sample'])['sample'].transform('size').to_numpy(float);g=q[['source_group','sample']].drop_duplicates().groupby('source_group').size();x/=q.source_group.map(g).to_numpy(float);return x/x.sum()
for target,t,pred in [('K32','n_K32_true','K32_mean'),('epithelial','n_epi_true','Epi_total_mean'),('total','n_total_true','total_cells_mean')]:
 p[target+'_absolute_error_delta_new_minus_old']=abs(p[pred+'_new']-p[t])-abs(p[pred+'_old']-p[t]);p[target+'_prediction_delta_new_minus_old']=p[pred+'_new']-p[pred+'_old']
 for scope,q in [('A_C_unique',p[p.pseudospot_id.isin(ids)]),('ALL_unique',p)]:
  for level in ['ALL','dataset','source_group','sample']:
   gs=[('ALL',q)] if level=='ALL' else q.groupby(level)
   for unit,z in gs:
    ww=w(z);newerr=abs(z[pred+'_new']-z[t]);olderr=abs(z[pred+'_old']-z[t]);rows.append({'target':target,'scope':scope,'level':level,'unit':unit,'n_bags':len(z),'n_samples':z['sample'].nunique(),'n_groups':z.source_group.nunique(),'new_MAE':float(np.sum(ww*newerr)),'old_MAE':float(np.sum(ww*olderr)),'paired_MAE_delta_new_minus_old':float(np.sum(ww*(newerr-olderr))),'new_mean_prediction':float(np.sum(ww*z[pred+'_new'])),'old_mean_prediction':float(np.sum(ww*z[pred+'_old']))})
wt(p,E/'paired_baseline_spots.tsv.gz');wt(pd.DataFrame(rows),E/'paired_baseline_summary.tsv');js(E/'PAIRED_BASELINE_STATUS.json',{'counts_identical':True,'seeds_identical':True,'n_bags':len(p),'new_reference':'new lower80 with withheld TRAIN nonK32 high excluded','old_reference':'existing TRAIN-only S31 with all frozen nonK32 Epi_other','not_historical_v6_MAE_subtraction':True});say('PAIRED BASELINE EVALUATED')
