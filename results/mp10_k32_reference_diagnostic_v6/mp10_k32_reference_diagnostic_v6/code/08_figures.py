from common import *
import matplotlib.pyplot as plt
O=R/'07_figures';E=R/'06_evaluation';D=pd.read_csv(E/'per_spot_all.tsv.gz',sep='\t');captions=[];plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False});depths=['low','medium','high']
def save(fig,name,data,caption,folder=O):
 folder.mkdir(exist_ok=True,parents=True);fig.tight_layout();fig.savefig(folder/(name+'.pdf'),transparent=True);fig.savefig(folder/(name+'.png'),dpi=300,transparent=True);plt.close(fig);rel=folder.relative_to(O);pdest=O/'plot_data'/rel;pdest.mkdir(parents=True,exist_ok=True);wt(data,pdest/(name+'.tsv.gz'));captions.append(dict(figure=str(rel/name),caption=caption))
def empty(ax,why):ax.text(.5,.5,why,ha='center',va='center',transform=ax.transAxes);ax.set_axis_off()
f=pd.read_csv(R/'01_feature_audit/gene_provenance_audit.tsv',sep='\t');f['assay_group']=f.dataset.str.strip()+' / '+f.assay.str.replace('Visium ','',regex=False);t=f.groupby(['assay_group','gene']).existing_h5ad_retained.mean().unstack();fig,ax=plt.subplots(figsize=(13,6));im=ax.imshow(t.values,vmin=0,vmax=1,cmap='Blues',aspect='auto');ax.set_xticks(range(len(t.columns)),t.columns,rotation=60,ha='right');ax.set_yticks(range(len(t)),t.index);fig.colorbar(im,ax=ax,label='Fraction of sections retaining feature');save(fig,'01_key_gene_matrix_availability',f,'Feature-list presence, not tissue expression positivity. Original assay and processing evidence is distinguished in the provenance table. HSD17B7, HMGCS1 and KRT18 are excluded from G0 by the63-section intersection.')
c=pd.read_csv(R/'01_feature_audit/frozen_MP_gene_coverage_by_assay.tsv',sep='\t');c=c[c.definition=='Top20'];c['assay_group']=c.dataset.str.strip()+' / '+c.assay.str.replace('Visium ','',regex=False);t=c.groupby(['assay_group','MP']).fraction_frozen_exactly_mapped_and_present.mean().unstack();fig,ax=plt.subplots(figsize=(14,6));im=ax.imshow(t.values,vmin=0,vmax=1,cmap='Blues',aspect='auto');ax.set_xticks(range(len(t.columns)),t.columns,rotation=60,ha='right');ax.set_yticks(range(len(t)),t.index);fig.colorbar(im,ax=ax,label='Frozen Top20 feature coverage');save(fig,'01_frozen_MP_Top20_coverage',c,'Frozen MP Top20 coverage from exact canonical mapping. Assay-group averages describe derived matrices; assay and processing limitations are supplied in gene_provenance_audit.tsv.')
sig=pd.read_csv(R/'02_existing_output_audit/epithelial_signature_comparison.tsv',sep='\t');fig,axes=plt.subplots(1,2,figsize=(10,4));axes[0].scatter(np.log1p(sig.S31_other),np.log1p(sig.S31_K32),s=3,alpha=.4,c='#173B75',rasterized=True,label='G0 genes');axes[0].set(xlabel='log1p S31 Epi_other',ylabel='log1p S31 Epi_K32');axes[0].legend();rc=pd.read_csv(R/'02_existing_output_audit/remaining29_signature_comparison.tsv',sep='\t');axes[1].scatter(rc.Pearson_log1p,rc.signature_sum_ratio,c='gray',label='Other29 labels');axes[1].axhline(1,c='black',ls='--',lw=.8);axes[1].set(xlabel='S31 / S30 log1p Pearson',ylabel='S30 / S31 signature sum');axes[1].legend();save(fig,'02_reference_signature_diagnostics',sig,'Left: frozen S31 epithelial signatures. Right: unchanged29 label identities after full S30 refitting. Correlation is descriptive. Other29 source values are in02_reference_other29.tsv.');wt(rc,O/'plot_data/02_reference_other29.tsv')
cor=pd.read_csv(R/'03_merged_reference/signature_Pearson_log1p.tsv',sep='\t',index_col=0);fig,ax=plt.subplots(figsize=(12,10));im=ax.imshow(cor,vmin=0,vmax=1,cmap='viridis');ax.set_xticks(range(len(cor)),cor.columns,rotation=90,fontsize=7);ax.set_yticks(range(len(cor)),cor.index,fontsize=7);fig.colorbar(im,ax=ax,label='Pearson after log1p');save(fig,'02_S30_reference_correlation',cor.reset_index(),'All30-label S30 signature correlations after log1p. No correlation threshold defines successful abundance recovery.')
full=D[D.model.isin(['HIST_MAIN_DEV','S30_FULL_DEV'])];maxepi=max(20,float(full.epithelial_mean.max())*1.05);maxtotal=max(25,float(full.total_cells_mean.max())*1.05)
def plots_for(d,folder):
 full=d[d.model.isin(['HIST_MAIN_DEV','S30_FULL_DEV'])];fig,axes=plt.subplots(2,3,figsize=(12,7),squeeze=False)
 for i,tag in enumerate(['HIST_MAIN_DEV','S30_FULL_DEV']):
  for j,dep in enumerate(depths):
   z=full[(full.model==tag)&(full.depth_condition==dep)];ax=axes[i,j]
   if z.empty:empty(ax,tag+' / '+dep+'\nUnavailable');continue
   for sample,q in z.groupby('sample'):ax.scatter(q.N_epithelial,q.epithelial_mean,s=6,alpha=.35,label=sample,rasterized=True)
   ax.plot([0,maxepi],[0,maxepi],'k--',lw=.8,label='Identity');ax.set(xlim=(0,21),ylim=(0,maxepi),xlabel='True epithelial cells',ylabel=('S31' if i==0 else 'S30')+' predicted / '+dep)
 if axes[0,0].has_data():axes[0,0].legend(fontsize=5,loc='upper left',bbox_to_anchor=(0,1.05),ncol=2)
 save(fig,'03_full_DEV_epithelial_calibration',full,'Identical fullDEV counts in S31/S30; each point is a technical pseudospot, colored by sample. All prediction extremes retained; identity is not fitted.',folder)
 fig,axes=plt.subplots(1,2,figsize=(10,4))
 for ax,tag in zip(axes,['HIST_MAIN_DEV','S30_FULL_DEV']):
  z=full[full.model==tag]
  if z.empty:empty(ax,'Unavailable');continue
  for dep in depths:
   q=z[z.depth_condition==dep];ax.scatter(q.post_full_UMI,q.total_cells_mean,s=5,alpha=.3,label=dep,rasterized=True)
  ax.axhline(20,c='black',ls='--',label='True20');ax.set(xlabel='Retained full RNA UMI',ylabel=tag+' total abundance',ylim=(0,maxtotal));ax.legend(fontsize=7)
 save(fig,'04_total_cells_and_depth',full,'Unnormalized total abundance; all truths equal20. No slope or correlation against constant total-cell truth is reported.',folder)
 ctl=d[~d.historical & d.model.str.contains('TARGET_UMI|FIXED_CAPTURE')];fig,axes=plt.subplots(1,2,figsize=(11,4))
 for ax,ep in zip(axes,['epithelial','total_cells']):
  if ctl.empty:empty(ax,'Controls unavailable');continue
  q=ctl.groupby(['model','source_group','sample'])[ep+'_abs_error'].mean().groupby(['model','source_group']).mean().unstack(0);xx=np.arange(len(q.columns))
  for gr,row in q.iterrows():ax.plot(xx,row.values,'o-',alpha=.7,lw=.7,label=gr)
  ax.set_xticks(xx,[x.replace('_','\n',1) for x in q.columns],fontsize=7);ax.set_ylabel(ep+' MAE / source group');ax.legend(fontsize=5,ncol=2)
 save(fig,'05_control_source_pairing',ctl,'Exactly paired20-cell bags across rules and references; source-group points average samples first. Unnormalized quantities; DEV technical diagnostic, not independent spatial validation.',folder)
 s31=d[(d.labelset=='S31') & d.model.isin(['HIST_MAIN_DEV','S31_TARGET_UMI','S31_FIXED_CAPTURE'])];fig,axes=plt.subplots(1,3,figsize=(13,4))
 if s31.empty:
  for ax in axes:empty(ax,'S31 unavailable')
 else:
  for tag,z in s31.groupby('model'):
   axes[0].scatter(z.N_K32,z.K32_mean,s=5,alpha=.25,c='#173B75',marker={'HIST_MAIN_DEV':'o','S31_TARGET_UMI':'x','S31_FIXED_CAPTURE':'+'}[tag],label=tag,rasterized=True);axes[1].scatter(z.K32_fraction_epithelial_true,z.K32_fraction_epithelial_mean,s=5,alpha=.25,c='#173B75',marker={'HIST_MAIN_DEV':'o','S31_TARGET_UMI':'x','S31_FIXED_CAPTURE':'+'}[tag],label=tag,rasterized=True)
  axes[0].plot([0,20],[0,20],'k--',lw=.7);axes[1].plot([0,1],[0,1],'k--',lw=.7);axes[0].set(xlabel='True frozen K32 cells',ylabel='S31 K32 abundance');axes[1].set(xlabel='True K32 / epithelial',ylabel='Posterior K32 / epithelial');axes[0].legend(fontsize=6)
  neg=s31[s31.N_K32==0];scenarios=sorted(neg.scenario.unique());axes[2].boxplot([neg.loc[neg.scenario==s,'K32_mean'] for s in scenarios],tick_labels=[s.replace('S2_','').replace('_non_K32','') for s in scenarios],showfliers=True,flierprops={'markersize':2});axes[2].axhline(8.367393,c='#173B75',ls='--',label='Historical DEV threshold');axes[2].tick_params(axis='x',rotation=80,labelsize=6);axes[2].set_ylabel('K32 abundance in true-zero spots');axes[2].legend(fontsize=6)
 save(fig,'06_K32_quantity_ratio_negatives',s31,'S31 K32 counts and drawwise ratios; true-zero tails retained. S30 has no K32 prediction. Boxplots pool displayed DEV collections descriptively; exact collections/source/scenario/depth tables remain separate.',folder)
 fig,axes=plt.subplots(1,2,figsize=(10,4))
 for tag,z in d[d.model.isin(['HIST_MAIN_DEV','S30_FULL_DEV','S31_TARGET_UMI','S31_FIXED_CAPTURE','S30_TARGET_UMI','S30_FIXED_CAPTURE'])].groupby('model'):
  axes[0].scatter(z.post_shared_UMI,z.detection_y_s_mean,s=4,alpha=.2,label=tag,rasterized=True);axes[1].scatter(z.background_fraction_from_means,z.epithelial_error,s=4,alpha=.2,label=tag,rasterized=True)
 axes[0].set(xlabel='Observed G0 UMI',ylabel='Per-location detection mean');axes[1].set(xlabel='Plug-in background RNA fraction',ylabel='Epithelial count bias');axes[0].legend(fontsize=5);save(fig,'07_detection_background',d,'Per-location detection is not one sample normalization. Background fractions are posterior-mean plug-in reconstruction, not posterior intervals. Descriptive associations do not establish a unique cause.',folder)
 # Historical evidence is also shown within every observed scRNA dataset.
 h=d[d.historical];fig,axes=plt.subplots(1,2,figsize=(10,4))
 for ax,col in zip(axes,['epithelial_abs_error','K32_abs_error']):
  q=h.groupby(['model','source_group','sample'])[col].mean().groupby(['model','source_group']).mean().reset_index()
  for tag,z in q.groupby('model'):ax.scatter([tag]*len(z),z[col],s=20,label=tag)
  ax.set_ylabel(col+' / group');ax.tick_params(axis='x',rotation=65,labelsize=7)
 axes[0].legend(fontsize=5);save(fig,'08_historical_by_source',h,'Previously inspected historical DEV/TEST only. Group-equal source points, no new fitting or untouched-test claim.',folder)
plots_for(D,O)
for dataset,d in D.groupby('dataset'):
 folder=O/'by_dataset'/str(dataset).replace('/','_');plots_for(d,folder);(folder/'STATUS.txt').write_text('Observed dataset only; source counts and unsupported conditions in tables. Technical bags are not independent patients.\n')
h=D[D.historical];q=h.groupby(['model','source_group','sample'])[['epithelial_abs_error','K32_abs_error','total_cells_abs_error']].mean().groupby(['model','source_group']).mean().reset_index();fig,axes=plt.subplots(1,3,figsize=(13,4))
for ax,col in zip(axes,['epithelial_abs_error','K32_abs_error','total_cells_abs_error']):
 for tag,z in q.groupby('model'):ax.scatter([tag]*len(z),z[col],s=18,label=tag)
 ax.tick_params(axis='x',rotation=75,labelsize=7);ax.set_ylabel(col+' / source group')
axes[0].legend(fontsize=5);save(fig,'08_historical_MAIN_NO_P07T_DEV_TEST',q,'Historical v5 models and previously inspected TEST, read-only re-summary. Conservative groups are not verified patients; no new untouched test set.')
a=pd.read_csv(E/'non_epithelial_group_equal.tsv',sep='\t');sel=a[a.model.isin(['HIST_MAIN_DEV','S30_FULL_DEV','S31_TARGET_UMI','S31_FIXED_CAPTURE','S30_TARGET_UMI','S30_FIXED_CAPTURE'])];t=sel.pivot(index='label',columns='model',values='error');v=max(abs(t.min().min()),abs(t.max().max()),.01);fig,ax=plt.subplots(figsize=(10,10));im=ax.imshow(t,aspect='auto',cmap='RdBu_r',vmin=-v,vmax=v);ax.set_yticks(range(len(t)),t.index,fontsize=8);ax.set_xticks(range(len(t.columns)),t.columns,rotation=60,ha='right');fig.colorbar(im,ax=ax,label='Predicted - true cells / group equal');save(fig,'09_non_epithelial_absorption',sel,'Signed errors for all29 unchanged non-epithelial labels; positive means excess assignment. Group-equal errors do not prove direct cell transitions.')
wt(pd.DataFrame(captions),O/'FIGURE_CAPTIONS.tsv');js(O/'FIGURES_DONE.json',{'figures':len(captions),'format':'PDF and transparent300dpiPNG','datasets':sorted(D.dataset.unique()),'all_extremes_retained':True});say('FIGURES_DONE',len(captions))
