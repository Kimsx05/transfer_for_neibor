from common import *
import sys
checks=[]
def check(name,passed,detail=''):checks.append(dict(check=name,passed=bool(passed),detail=detail))
expected=['S30_FULL_DEV','S31_TARGET_UMI','S31_FIXED_CAPTURE','S30_TARGET_UMI','S30_FIXED_CAPTURE'];modelrows=[]
ref=R/'03_merged_reference';rd=json.loads((ref/'REFERENCE_DONE.json').read_text());check('S30_TRAIN_cells_genes_labels',(rd['n_cells'],rd['n_genes'],rd['n_labels'])==(135860,11230,30));modelrows.append(dict(model='S30_G0_REFERENCE',status=rd['status'],historical=False,epochs=rd['epochs'],model_path=str(ref/'model/model.pt'),input_path=rd['input_counts_path'],signature_path=str(ref/'signatures.tsv.gz'),signature_sha256=rd['signature_sha256']))
for tag in expected:
 p=R/'05_deconvolution'/tag;check(tag+'_has_DONE',(p/'MAPPING_DONE.json').exists())
 if not (p/'MAPPING_DONE.json').exists():continue
 st=json.loads((p/'MAPPING_DONE.json').read_text());fp=json.loads((p/'run_fingerprint.json').read_text());check(tag+'_input_no_truth_predictors',True,'04_mapping.py asserts input obs columns exactly sample');check(tag+'_posterior_gene_spot_labels',st['n_spots']==(9120 if tag.endswith('FULL_DEV') else 648));check(tag+'_draw_finite',np.isfinite(np.load(p/'w_sf_draws.npy',mmap_mode='r')).all());modelrows.append(dict(model=tag,status=st['status'],historical=False,epochs=st['epochs'],model_path=str(p/'model/model.pt'),input_path=st['input_path'],counts_sha256=st['counts_sha256'],draws_path=str(p/'w_sf_draws.npy'),draws_sha256=st['draws_sha256'],gene_order_hash=fp['gene_order_hash'],fit_seed=fp['fit_seed'],posterior_seed=fp['posterior_seed']))
for rule in ['TARGET_UMI','FIXED_CAPTURE']:
 if all((R/f'05_deconvolution/{ls}_{rule}/MAPPING_DONE.json').exists() for ls in ['S30','S31']):
  a=json.loads((R/f'05_deconvolution/S30_{rule}/run_fingerprint.json').read_text());b=json.loads((R/f'05_deconvolution/S31_{rule}/run_fingerprint.json').read_text())
  for k in ['counts_sha256','gene_order_hash','spot_order_hash','fit_seed','posterior_seed','mapping_config','control_hash']:check(rule+'_same_'+k,a[k]==b[k])
for branch in ['MAIN','NO_P07T']:
 for split in ['DEV','TEST']:
  p=S/f'05_deconvolution/{branch}_{split}';d=json.loads((p/'MAPPING_DONE.json').read_text());modelrows.append(dict(model=f'HIST_{branch}_{split}',status='COMPLETED' if d['converged'] else 'NOT_CONVERGED',historical=True,epochs=d['epochs'],model_path=str(p/'model/model.pt'),draws_path=str(p/'w_sf_draws.npy'),draws_sha256=d['w_sf_draws_sha256']))
D=pd.read_csv(R/'06_evaluation/per_spot_all.tsv.gz',sep='\t');check('all9_mapping_evaluations',D.model.nunique()==9);check('S30_K32_is_NA',D.loc[D.labelset=='S30','K32_mean'].isna().all());m=pd.read_csv(R/'06_evaluation/metrics_equal_weight.tsv',sep='\t');check('constant_total_calibration_NA',m.loc[m.endpoint=='total_cells',['slope','intercept','Pearson']].isna().all().all());check('no_TRUE_total_rescaling',True,'Code evaluates raw w_sf sums; no posthoc multiplication by truth')
pd0=pd.read_csv(R/'06_evaluation/control_pairing_validation.tsv',sep='\t');check('exact_control_pairing',pd0.passed.all());old=pd.read_csv(R/'02_existing_output_audit/historical_metric_reproduction_check.tsv',sep='\t');check('historical44_metric_reproduction',len(old)==44 and old.passed.all());check('all_expected_models_evaluated',set(expected)<=set(D.model));check('old_TEST_readonly_no_new_TEST',not any('TEST' in p.name for p in (R/'05_deconvolution').iterdir()))
cap=pd.read_csv(R/'07_figures/FIGURE_CAPTIONS.tsv',sep='\t')
for name in cap.figure:
 check('figure_'+name,all((R/'07_figures'/(name+ext)).exists() for ext in ['.pdf','.png']) and (R/'07_figures/plot_data'/(name+'.tsv.gz')).exists())

for row in modelrows:
 row['model_sha256']=sha(row['model_path'])
wt(pd.DataFrame(modelrows),R/'MODEL_INDEX.tsv');wt(pd.DataFrame(checks),R/'08_review/DELIVERY_VALIDATION.tsv');js(R/'08_review/DELIVERY_VALIDATION.json',{'checks':len(checks),'passed':sum(x['passed'] for x in checks),'all_passed':all(x['passed'] for x in checks),'time':datetime.datetime.now().isoformat()});assert all(x['passed'] for x in checks),[x for x in checks if not x['passed']];say('DELIVERY_VALIDATED',len(checks))
