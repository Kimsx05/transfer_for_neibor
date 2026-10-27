from common import *
import matplotlib.pyplot as plt
import re
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
def safe(x):return re.sub(r'[^a-zA-Z0-9_-]+','_',str(x))
def put(fig,ds,name):
 dest=O/'figures'/('overall' if ds=='ALL' else 'by_dataset/'+ds);dest.mkdir(parents=True,exist_ok=True)
 for ext in ['png','pdf']:fig.savefig(dest/(name+'.'+ext),dpi=300,transparent=True,bbox_inches='tight')
 plt.close(fig)
s=read('samples_v2.tsv');edge=read('network_node_edges.tsv.gz');io=read('network_incoming_outgoing.tsv');axes=read('selected_review_axes.tsv');context=read('selected_axes_all_network_context.tsv');expr=read('candidate_sample_expression.tsv.gz');evidence=read('candidate_dataset_evidence.tsv');colors=dict(read('cell_colors.tsv')[['cell_type','color']].values)
# Actual independently inferred major networks; main and secondary never mix their scales.
for ds in ['ALL']+sorted(s.dataset.unique()):
 for co,groups in [('epi_state_within_specimen',['shared']),('sample_MP10_group_allEpi',['low','high']),('Normal_Tumor',['Normal','Tumor']),('T1_T2',['T1','T2'])]:
  z=edge[edge.scope.eq(ds)&edge.comparison_type.eq(co)&edge.network_level.eq('major')];nodes=sorted(set(z.source)|set(z.target));fig,axs=plt.subplots(1,len(groups),figsize=(8*len(groups),7),squeeze=False);vmax=z.supported_sum.max() if len(z) else 1
  for j,g in enumerate(groups):
   a=z[z.condition.eq(g)];ax=axs[0,j]
   if a.empty:ax.text(.5,.5,'NE: no eligible condition object',ha='center');ax.axis('off');continue
   mat=a.pivot(index='source',columns='target',values='supported_sum').reindex(index=nodes,columns=nodes);im=ax.imshow(np.ma.masked_invalid(mat),cmap='viridis',vmin=0,vmax=max(vmax,1e-12));ax.set(xticks=range(len(nodes)),yticks=range(len(nodes)),xticklabels=nodes,yticklabels=nodes,title=ds+' / '+g);ax.tick_params(axis='x',rotation=70);fig.colorbar(im,ax=ax,shrink=.6,label='Sum of internally supported LR probabilities')
  fig.suptitle(co+'; actual major inference\n'+('Shared object, both epi states' if co=='epi_state_within_specimen' else 'Separate condition objects; object-specific expression maxima'));fig.tight_layout();put(fig,ds,'major_overview_'+co)
  a=io[io.scope.eq(ds)&io.comparison_type.eq(co)&io.network_level.eq('major')];fig,axs=plt.subplots(1,2,figsize=(14,5));nodes=sorted(a.node.unique())
  for j,field in enumerate(['outgoing_supported','incoming_supported']):
   for k,g in enumerate(groups):
    v=a[a.condition.eq(g)].set_index('node')[field].reindex(nodes);axs[j].bar(np.arange(len(nodes))+(k-(len(groups)-1)/2)*.36,v,width=.35,label=g)
   axs[j].set(xticks=range(len(nodes)),xticklabels=nodes,title=field,ylabel='Supported probability sum');axs[j].tick_params(axis='x',rotation=70);axs[j].legend()
   if a.empty:axs[j].text(.5,.5,'NE: no eligible condition object',ha='center',transform=axs[j].transAxes)
  fig.suptitle(ds+' '+co+'; missing sums stay NA');fig.tight_layout();put(fig,ds,'major_incoming_outgoing_'+co)
 # All selected B/Plasma LR raw contrasts across the four questions; NE explicit.
 a=context[context.scope.eq(ds)];cs=['epi_state_within_specimen','sample_MP10_group_allEpi','Normal_Tumor','T1_T2'];fig,ax=plt.subplots(figsize=(12,max(5,len(axes)*.55)));mat=a.pivot(index='axis_id',columns='comparison_type',values='raw_delta').reindex(index=axes.axis_id,columns=cs);v=mat.to_numpy();lim=np.nanmax(abs(v)) if np.isfinite(v).any() else 1
 im=ax.imshow(np.ma.masked_invalid(v),cmap='RdBu_r',vmin=-max(lim,1e-12),vmax=max(lim,1e-12),aspect='auto');ax.set(xticks=range(4),xticklabels=['Within-specimen epi state','Specimen MP10 high-low','Tumor-Normal','T2-T1'],yticks=range(len(axes)),yticklabels=[f'{x.axis_id} {x.TME_type} {x.direction} {x.LR_id}' for x in axes.itertuples()]);ax.tick_params(axis='x',rotation=20)
 for i in range(len(axes)):
  for j in range(4):ax.text(j,i,'NE' if not np.isfinite(v[i,j]) else f'{v[i,j]:.2g}',ha='center',va='center',fontsize=8)
 ax.set_title(ds+': raw differences only when both probabilities were computed\nSecondary objects have separate expression scaling; contrasts answer different questions');fig.colorbar(im,ax=ax,label='Raw probability difference');fig.tight_layout();put(fig,ds,'B_Plasma_LR_contrast_heatmap')
 # Sample paired expression: all evaluable samples, shared target >=30. Zero in one epi state remains visible.
 fig,axs=plt.subplots(len(axes),2,figsize=(12,max(5,len(axes)*2.3)),squeeze=False)
 for i,ar in enumerate(axes.itertuples()):
  a=expr[expr.network_level.eq(ar.network_level)&expr.TME_type.eq(ar.TME_type)&expr.direction.eq(ar.direction)&expr.LR_id.eq(ar.LR_id)&expr.coverage_eligible];a=a if ds=='ALL' else a[a.dataset.eq(ds)];ax=axs[i,0];bx=axs[i,1]
  if a.empty:ax.text(.5,.5,'NE: no matched specimen coverage',ha='center',transform=ax.transAxes);bx.text(.5,.5,'NE',ha='center',transform=bx.transAxes)
  for d,zz in a.groupby('dataset'):
   col=plt.get_cmap('tab10')(sorted(s.dataset.unique()).index(d))
   for row in zz.itertuples():ax.plot([0,1],[row.epi_low_mean,row.epi_high_mean],'-o',color=col,alpha=.45,ms=3,lw=.7)
  for k,row in enumerate(a.sort_values(['dataset','sample_uid']).itertuples()):bx.scatter(k,row.aux_product_delta,color=plt.get_cmap('tab10')(sorted(s.dataset.unique()).index(row.dataset)),s=20,marker='o' if row.high_supported or row.low_supported else 'x')
  ax.set(xticks=[0,1],xticklabels=['epi low','epi high'],ylabel='Geometric mean RNA:data',title=f'{ar.axis_id}: {ar.TME_type}; {ar.LR_id}; {ar.direction}');bx.axhline(0,color='gray',lw=.6);bx.set(title=f'{a.sample_uid.nunique()} specimens / {a.patient_uid_recorded.nunique()} recorded IDs',xlabel='Matched specimens (dataset, sample order)',ylabel='Auxiliary product high - low')
 fig.suptitle(ds+': paired expression, not CellChat probability; no patient-level CellChat P values');fig.tight_layout();put(fig,ds,'B_Plasma_paired_expression')
# Cross-dataset patient-direction summary, no invented significance.
fig,axs=plt.subplots(1,2,figsize=(17,max(5,len(axes)*.55)));ds=sorted(s.dataset.unique());vals=np.full((len(axes),len(ds)),np.nan)
for i,ar in enumerate(axes.itertuples()):
 a=evidence[evidence.network_level.eq(ar.network_level)&evidence.TME_type.eq(ar.TME_type)&evidence.direction.eq(ar.direction)&evidence.LR_id.eq(ar.LR_id)].set_index('dataset')
 for j,d in enumerate(ds):
  if d not in a.index:continue
  q=a.loc[d];vals[i,j]=q.n_positive-q.n_negative;axs[0].text(j,i,f'+{int(q.n_positive)}/−{int(q.n_negative)}/0:{int(q.n_zero)}',ha='center',va='center',fontsize=7);axs[1].text(j,i,f'n={int(q.n_patients)}; support={int(q.n_supported_patients)}\n'+('stable LOO' if q.leave_one_patient_direction_stable else 'unstable/NE LOO'),ha='center',va='center',fontsize=7)
lim=np.nanmax(abs(vals)) if np.isfinite(vals).any() else 1
for ax in axs:ax.imshow(np.ma.masked_invalid(vals),cmap='RdBu_r',vmin=-max(lim,1),vmax=max(lim,1),aspect='auto');ax.set(xticks=range(len(ds)),xticklabels=ds,yticks=range(len(axes)),yticklabels=axes.axis_id);ax.tick_params(axis='x',rotation=30)
axs[0].set_title('Recorded-patient direction counts (positive / negative / zero)');axs[1].set_title('Coverage and leave-one-patient description; not significance');fig.tight_layout();put(fig,'ALL','B_Plasma_patient_direction_summary')
print('Network figures done')
