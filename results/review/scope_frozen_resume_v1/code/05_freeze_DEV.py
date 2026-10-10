from common import *
import time
while (R/"05_deconvolution/DEV_convergence_probe_5000/PROBE_CONFIG.json").exists() and not (R/"06_evaluation/DEV_CONVERGENCE_REVIEW.json").exists():time.sleep(10)
while not all((R/f'05_deconvolution/{b}_DEV/MAPPING_DONE.json').exists() for b in ['MAIN','NO_P07T']):time.sleep(10)
cfg=json.loads((R/'CONFIG.yaml').read_text())
if (R/'06_evaluation/DEV_CONVERGENCE_REVIEW.json').exists():cfg['mapping'].update(json.loads((R/'06_evaluation/DEV_CONVERGENCE_REVIEW.json').read_text()).get('mapping_overrides',{}))
d=pd.read_csv(R/'04_pseudospots/DEV/manifest.tsv.gz',sep='\t');p=pd.read_csv(R/'05_deconvolution/MAIN_DEV/abundance_mean.tsv.gz',sep='\t');d=d.merge(p[['pseudospot_id','Epi_K32']],on='pseudospot_id',validate='one_to_one');neg=d[d.N_K32==0].copy();qs=[]
# Equal source weighting within each scenario; protect each difficult negative stratum from dilution.
for scene,z in neg.groupby('scenario'):
 z=z.sort_values('Epi_K32');w=1/z.groupby('source_group').source_group.transform('size');w=w/w.sum();q=float(z.Epi_K32.iloc[np.searchsorted(w.cumsum(),.95)]);qs.append({'scenario':scene,'n_spots':len(z),'n_sources':z.source_group.nunique(),'source_equal_95_percentile':q})
th=max(x['source_equal_95_percentile'] for x in qs);neg['detected']=neg.Epi_K32>th;per=neg.groupby(['scenario','source_group']).detected.mean().reset_index(name='FPR');wt(pd.DataFrame(qs),R/'06_evaluation/DEV_threshold_strata.tsv');wt(per,R/'06_evaluation/DEV_threshold_achieved_by_source.tsv');actual=per.groupby('scenario').FPR.mean().to_dict();states={b:json.loads((R/f'05_deconvolution/{b}_DEV/MAPPING_DONE.json').read_text()) for b in ['MAIN','NO_P07T']};assert all(v['converged'] for v in states.values()),'Do not freeze a biological decision rule on unconverged DEV'
f={'threshold_on_posterior_mean':th,'decision':'posterior_mean > threshold','FPR_target':.05,'DEV_negative_composition':'seven scenario strata, equal source weights within each; max scenario-specific 95% weighted quantile','achieved_source_equal_FPR_by_scenario':actual,'mapping_config':cfg['mapping'],'stopping_rule':'same DEV-amended loss plus posterior mean-abundance stability rule applied to TEST without truth','DEV_model_states':states,'MAIN_DEV_mean_sha256':sha(R/'05_deconvolution/MAIN_DEV/abundance_mean.tsv.gz'),'time':datetime.datetime.now().isoformat(),'TEST_predictions_seen':False}
if (R/'06_evaluation/DEV_FROZEN.json').exists():
 old=json.loads((R/'06_evaluation/DEV_FROZEN.json').read_text());assert old['MAIN_DEV_mean_sha256']==f['MAIN_DEV_mean_sha256'] and old['threshold_on_posterior_mean']==th
else:js(R/'06_evaluation/DEV_FROZEN.json',f)
say('DEV_FROZEN',th,actual)
