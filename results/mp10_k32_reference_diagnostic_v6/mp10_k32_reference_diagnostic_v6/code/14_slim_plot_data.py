from common import *
O=R/'07_figures/plot_data';audit=[]
meta=['pseudospot_id','model','labelset','split','historical','dataset','sample','source_group','scenario','depth_condition','N_K32','N_epithelial','N_total','bag_id']
selected={'03_full_DEV_epithelial_calibration':['epithelial_mean'],'04_total_cells_and_depth':['post_full_UMI','total_cells_mean'],'06_K32_quantity_ratio_negatives':['K32_mean','K32_fraction_epithelial_true','K32_fraction_epithelial_mean'],'07_detection_background':['post_shared_UMI','detection_y_s_mean','background_fraction_from_means','epithelial_error']}
for p in sorted(O.rglob('*.tsv.gz')):
 name=p.name.removesuffix('.tsv.gz')
 if name not in selected and name not in ['05_control_source_pairing','08_historical_by_source']:continue
 d=pd.read_csv(p,sep='\t',low_memory=False);before=p.stat().st_size
 if name in ['05_control_source_pairing','08_historical_by_source']:
  cols=['epithelial_abs_error','total_cells_abs_error'] if name.startswith('05') else ['epithelial_abs_error','K32_abs_error']
  if 'sample' in d:d=d.groupby(['model','source_group','sample'])[cols].mean().groupby(['model','source_group']).mean().reset_index()
  else:assert not d.duplicated(['model','source_group']).any(), 'Already compacted group points must be unique'
 else:
  if name=='07_detection_background':d=d[d.model.isin(['HIST_MAIN_DEV','S30_FULL_DEV','S31_TARGET_UMI','S31_FIXED_CAPTURE','S30_TARGET_UMI','S30_FIXED_CAPTURE'])]
  d=d[[x for x in meta+selected[name] if x in d]]
 wt(d,p);audit.append(dict(plot_data=str(p.relative_to(R)),bytes_before=before,bytes_after=p.stat().st_size,rule='Only plotted variables retained; source-point plots store exact group means; full per-spot evaluation remains separately available'))
wt(pd.DataFrame(audit),R/'08_review/PLOT_DATA_COMPACTION.tsv');say('PLOT_DATA_COMPACTED',len(audit))
