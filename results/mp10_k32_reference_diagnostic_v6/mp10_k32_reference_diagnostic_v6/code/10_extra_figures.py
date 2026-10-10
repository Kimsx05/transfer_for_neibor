from common import *
import matplotlib.pyplot as plt
E=R/'06_evaluation';O=R/'07_figures';caps=pd.read_csv(O/'FIGURE_CAPTIONS.tsv',sep='\t').to_dict('records');plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
def save(fig,name,data,caption):
 fig.tight_layout();fig.savefig(O/(name+'.pdf'),transparent=True);fig.savefig(O/(name+'.png'),dpi=300,transparent=True);plt.close(fig);wt(data,O/'plot_data'/(name+'.tsv.gz'));caps.append(dict(figure=name,caption=caption))
c=pd.read_csv(E/'detection_curves.tsv.gz',sep='\t');fig,axes=plt.subplots(1,2,figsize=(11,4))
for ax,kind in zip(axes,['PR','ROC']):
 for tag,z in c[c.kind==kind].groupby('model'):ax.plot(z.x,z.y,label=tag,lw=1)
 ax.set(xlabel='Recall' if kind=='PR' else 'False positive rate',ylabel='Precision' if kind=='PR' else 'True positive rate',xlim=(0,1),ylim=(0,1));ax.legend(fontsize=6)
save(fig,'10_S31_PR_ROC',c,'Pooled simulated PR/ROC curves, with prevalence in plot data and per-source AUC/undefined denominators in detection tables. Precision at constructed prevalence is not real-Visium precision. No threshold was selected from these curves.')
d=pd.read_csv(E/'detection_group_equal_strata.tsv',sep='\t');d=d[d.stratification=='N_K32|depth_condition|N_epithelial'];m=pd.read_csv(E/'metrics_group_equal_all_strata.tsv',sep='\t');m=m[(m.stratification=='N_K32|depth_condition|N_epithelial')&(m.endpoint=='K32')];tags=sorted(d.model.unique());fig,axes=plt.subplots(2,len(tags),figsize=(max(12,len(tags)*2.5),7),squeeze=False)
for i,tag in enumerate(tags):
 for ax,src,value in [(axes[0,i],m,'MAE'),(axes[1,i],d,'recall')]:
  z=src[(src.model==tag)&(src.N_epithelial==10)];t=z.pivot(index='N_K32',columns='depth_condition',values=value).reindex(columns=['low','medium','high']);im=ax.imshow(t,aspect='auto',cmap='Blues',vmin=0,vmax=1 if value=='recall' else None);ax.set_xticks(range(3),['low','medium','high'],rotation=45);ax.set_yticks(range(len(t)),[int(x) for x in t.index]);ax.set_ylabel('True K32 cells');ax.set_title(tag,fontsize=8);fig.colorbar(im,ax=ax,label=value)
save(fig,'11_K32_quantity_depth',pd.concat([m.assign(measure='MAE'),d.assign(measure='recall')],ignore_index=True),'K32 MAE and frozen-threshold detection by quantity and depth at ten epithelial cells. Equal samples within conservative source groups then equal groups; unsupported cells remain blank. Negative spots have undefined recall.')
p=E/'intensity_by_sample.tsv'
if p.exists():
 d=pd.read_csv(p,sep='\t');fig,axes=plt.subplots(1,2,figsize=(11,4))
 for ax,col in zip(axes,['delta_high_minus_low__K32_mean','delta_high_minus_low__epithelial_mean']):
  for tag,z in d.groupby('model'):
   if col in z and z[col].notna().any():ax.scatter([tag]*len(z),z[col],label=tag,s=22)
  ax.axhline(0,c='black',lw=.8);ax.set_ylabel(col.replace('delta_high_minus_low__','High - low '));ax.tick_params(axis='x',rotation=70,labelsize=7)
 axes[0].legend(fontsize=5);save(fig,'12_existing_state_intensity',d,'Existing frozen S4 cells only, at fixed K32 counts and background. Points average paired pseudospots within source samples. RNA-depth differences remain confounders; no new high-core definition and no per-spot significance tests.')
rows=[];fig,axes=plt.subplots(1,2,figsize=(11,4))
for p in [R/'03_merged_reference/training_history.tsv']+sorted((R/'05_deconvolution').glob('*/training_history.tsv')):
 d=pd.read_csv(p,sep='\t');tag=p.parent.name;d['model']=tag;rows.append(d);ax=axes[0] if '03_merged' in str(p) else axes[1];ax.plot(d.epoch,d.loss/d.loss.iloc[0],label=tag)
for ax in axes:ax.set(xlabel='Training epochs',ylabel='Loss / first epoch loss');ax.legend(fontsize=6)
save(fig,'13_training_curves',pd.concat(rows),'Actual completed training histories. Final stopping additionally requires posterior-abundance stability for mapping; a flat loss alone is insufficient. NB and spatial mapping have separate epoch rules.')
sig=pd.read_csv(R/'02_existing_output_audit/epithelial_signature_comparison.tsv',sep='\t').set_index('gene');panel=['SREBF2','SREBF1','FASN','ELOVL5','ADIPOR2','ABHD3','INSIG1','ACSL3','ACSS2'];z=sig.loc[[x for x in panel if x in sig.index]].copy();fig,ax=plt.subplots(figsize=(11,4));xx=np.arange(len(z))
for k,(col,color,label) in enumerate([('S31_K32','#173B75','S31 Epi_K32'),('S31_other','#999999','S31 Epi_other'),('S30_Epithelial','#B45B2A','S30 Epithelial')]):ax.bar(xx+(k-1)*.24,np.log1p(z[col]),width=.24,color=color,label=label)
ax.set_xticks(xx,z.index,rotation=45,ha='right');ax.set_ylabel('log1p posterior mean signature');ax.legend();save(fig,'14_fixed_support_genes',z.reset_index(),'Fixed support genes retained in G0. HSD17B7/HMGCS1/KRT18 are absent from G0 and are not shown as zero expression. Signature contrasts are descriptive, not sample-level differential expression tests.')
neg=pd.read_csv(E/'negative_distributions_by_source.tsv',sep='\t');neg=neg[neg.model.isin(['HIST_MAIN_DEV','HIST_MAIN_TEST','S31_TARGET_UMI','S31_FIXED_CAPTURE'])];ng=neg.groupby(['model','source_group','scenario','depth_condition'])[['mean','FPR']].mean().groupby(['model','scenario','depth_condition']).mean().reset_index();tags=sorted(ng.model.unique());fig,axes=plt.subplots(2,len(tags),figsize=(len(tags)*3.8,8),squeeze=False)
for i,tag in enumerate(tags):
 for j,value in enumerate(['mean','FPR']):
  z=ng[ng.model==tag];t=z.pivot(index='scenario',columns='depth_condition',values=value).reindex(columns=['low','medium','high']);ax=axes[j,i];im=ax.imshow(t,aspect='auto',cmap='Blues',vmin=0,vmax=1 if value=='FPR' else ng['mean'].max());ax.set_xticks(range(3),['low','medium','high'],rotation=45);ax.set_yticks(range(len(t)),[x.replace('S2_','').replace('_non_K32','').replace('S1_regular','Regular').replace('S3_pure_epithelial','Pure other') for x in t.index],fontsize=7);ax.set_title(tag,fontsize=8);fig.colorbar(im,ax=ax,label='Negative K32 mean' if value=='mean' else 'False positive rate')
save(fig,'15_negative_background_depth',ng,'True-zero K32 means and frozen-threshold FPR, separate by historical/full or paired control collection and depth. Equal samples within groups then groups. Unavailable conditions blank. Source-specific tails and maxima are retained separately and must not be inferred from group averages.')
wt(pd.DataFrame(caps).drop_duplicates('figure',keep='last'),O/'FIGURE_CAPTIONS.tsv');js(O/'FIGURES_DONE.json',{'figures':len(pd.DataFrame(caps).drop_duplicates('figure')),'format':'PDF and transparent300dpiPNG','all_extremes_retained':True});say('EXTRA_FIGURES_DONE')
