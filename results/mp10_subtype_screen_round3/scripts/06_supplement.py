from common import *
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,FancyArrowPatch
plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
ss=read('supplementary_NECTIN2_TIGIT.tsv');ss=ss[ss.LR_id.eq('NECTIN2_TIGIT')];ex=read('supplementary_NECTIN2_TIGIT_expression.tsv');pat=read('supplementary_NECTIN2_TIGIT_patient_differences.tsv');lr=read('treg_LR_contribution_decomposition.tsv.gz');nd=read('node_order.tsv');colors=dict(nd[['node','color']].drop_duplicates().values);cov=read('node_coverage.tsv');folder=O/'figures/supplementary/NECTIN2_TIGIT';figs=[];plotedges=[];row=[]
def put(fig,name,scope,question):
 for ext in ['png','pdf']:fig.savefig(folder/(name+'.'+ext),dpi=300,transparent=True,bbox_inches='tight')
 figs.append(dict(path=str((folder/(name+'.pdf')).relative_to(O)),scope=scope,network_level='fine',question=question,measure='single_LR_only',self_loops=False,data_source='supplementary_NECTIN2_TIGIT.tsv; supplementary_NECTIN2_TIGIT_expression.tsv'));plt.close(fig)
for ds in ss.scope.unique():
 a=ss[ss.scope.eq(ds)].set_index('epi_state');coords={EPI[0]:(-.866,.5),EPI[1]:(.866,.5),'Treg':(0,-1)};fig,ax=plt.subplots(figsize=(8,7));ax.set(aspect='equal',xlim=(-1.6,1.6),ylim=(-1.5,1.2));ax.axis('off');patches={}
 for node,(x,y) in coords.items():
  circle=Circle((x,y),.08,facecolor=colors[node],edgecolor='white',lw=.7);ax.add_patch(circle);patches[node]=circle;cc=cov[cov.scope.eq(ds)&cov.network_level.eq('fine')&cov.node.eq(node)];label=node+(' [NE]' if cc.empty or not cc.computed.iloc[0] else '');ax.text(x,y+(.17 if node!='Treg' else -.17),label,ha='center',va='center',fontsize=12)
 vmax=a.display_strength.max();vmax=max(vmax,1e-12) if np.isfinite(vmax) else 1e-12
 for ep in EPI:
  z=a.loc[ep];v=z.display_strength
  if np.isfinite(v) and v>0:ax.add_patch(FancyArrowPatch(coords[ep],coords['Treg'],patchA=patches[ep],patchB=patches['Treg'],arrowstyle='-|>',mutation_scale=17,linewidth=.5+6*v/vmax,color=colors[ep],alpha=.65,connectionstyle='arc3,rad='+('-.16' if ep==EPI[0] else '.16')))
  ax.text(coords[ep][0],.88,('NE' if not np.isfinite(v) else f'raw p={z.raw_probability:.3g}\ninternal P={z.internal_P:.3g}\n{z.status}'),ha='center',fontsize=9);plotedges.append(dict(scope=ds,source=ep,target='Treg',LR_id='NECTIN2_TIGIT',raw_probability=z.raw_probability,internal_P=z.internal_P,status=z.status,display_strength=v,edge_weight_max=vmax))
 fig.suptitle(ds+': NECTIN2 → TIGIT (supplement only)\nBoth actual epi states in the same object',fontsize=13);put(fig,'circle_'+ds,ds,'supp_single_LR_circle')
 ev=ex[ex.coverage_eligible&ex.genes_observed_in_dataset];ev=ev if ds=='ALL' else ev[ev.dataset.eq(ds)];pp=pat if ds=='ALL' else pat[pat.dataset.eq(ds)];fig,axs=plt.subplots(1,3,figsize=(15,4.8));datasets=sorted(ex.dataset.unique())
 for dr,z in ev.groupby('dataset'):
  co=plt.get_cmap('tab10')(datasets.index(dr))
  for r in z.itertuples():axs[0].plot([0,1],[r.NECTIN2_low_mean,r.NECTIN2_high_mean],'-o',c=co,alpha=.55,ms=4,lw=.8)
 for i,r in enumerate(ev.sort_values(['dataset','sample_uid']).itertuples()):
  axs[1].scatter(i,r.NECTIN2_low_detection,c=colors[EPI[0]],s=20);axs[1].scatter(i,r.NECTIN2_high_detection,c=colors[EPI[1]],s=20);axs[1].scatter(i,r.TIGIT_Treg_detection,c=colors['Treg'],marker='x',s=25)
 for i,r in enumerate(pp.itertuples()):axs[2].bar(i,r.expression_delta,color=colors[EPI[1]] if r.expression_delta>=0 else colors[EPI[0]])
 axs[0].set(xticks=[0,1],xticklabels=['epi low','epi high'],ylabel='NECTIN2 mean RNA:data',title='Matched specimen expression');axs[1].axhline(.1,c='black',ls='--',lw=.6);axs[1].set(xlabel='Matched specimens',ylabel='Detection fraction',title='Blue/red: epi NECTIN2; ×: Treg TIGIT',ylim=(-.02,1.02));axs[2].axhline(0,c='black',lw=.6);axs[2].set(xlabel='Recorded IDs',ylabel='Mean paired NECTIN2 high − low',title='Repeated specimens averaged within ID');fig.suptitle(f'{ds}: {len(ev)} matched specimens / {len(pp)} recorded IDs; not validation of total Treg W');fig.tight_layout();put(fig,'paired_expression_'+ds,ds,'supp_single_LR_expression')
 # Contribution context is based on W contributions, not unfiltered raw differences.
 for dr in ['epi_to_TME','TME_to_epi']:
  aa=lr[lr.scope.eq(ds)&lr.direction.eq(dr)];aa=aa.copy();aa['positive_rank']=aa.contribution.where(aa.contribution.gt(0)).rank(ascending=False,method='min');den=aa.loc[aa.contribution.gt(0),'contribution'].sum();zz=aa[aa.LR_id.isin(['NECTIN2_TIGIT','TIGIT_NECTIN2'])]
  for rr in zz.itertuples():row.append(dict(scope=ds,direction=dr,LR_id=rr.LR_id,contribution=rr.contribution,positive_contribution_rank=rr.positive_rank,n_positive_LR=int(aa.contribution.gt(0).sum()),share_of_positive_contributions=rr.contribution/den if rr.contribution>0 and den>0 else np.nan,total_delta_W=aa.contribution.sum(),prob_low=rr.prob_low,prob_high=rr.prob_high,P_low=rr.P_low,P_high=rr.P_high))
fig,axs=plt.subplots(1,2,figsize=(13,5));dslist=[x for x in ss.scope.unique() if x!='ALL']+['ALL']
for i,ds in enumerate(dslist):
 a=ss[ss.scope.eq(ds)].set_index('epi_state')
 for j,metric in enumerate(['raw_probability','display_strength']):axs[j].plot([i-.14,i+.14],[a.loc[EPI[0],metric],a.loc[EPI[1],metric]],c='gray');axs[j].scatter(i-.14,a.loc[EPI[0],metric],c=colors[EPI[0]],s=40);axs[j].scatter(i+.14,a.loc[EPI[1],metric],c=colors[EPI[1]],s=40)
for j,metric in enumerate(['Raw probability','Internally supported display strength']):axs[j].set(xticks=range(len(dslist)),xticklabels=dslist,ylabel=metric,title=metric);axs[j].tick_params(axis='x',rotation=35);axs[j].axvline(len(dslist)-1.5,c='black',ls='--',lw=.5)
fig.suptitle('Supplement: NECTIN2 → TIGIT; compare state direction within objects, not absolute cohort values');fig.tight_layout();put(fig,'cohort_comparison','ALL','supp_cohort_comparison')
save(plotedges,'supplementary_circle_edges.tsv');save(row,'supplementary_NECTIN2_TIGIT_contribution_context.tsv');save(figs,'figure_index_supplement.tsv')
# A concise receiver comparison retains CD8 Exhausted and all other actually TIGIT-positive T/NK labels.
tg=read('supplementary_TIGIT_receiver_summary.tsv');net=read('supplementary_TIGIT_receivers_network.tsv');net=net[net.LR_id.eq('NECTIN2_TIGIT')&net.scope.ne('ALL')];tg=tg[tg.specimen_equal_TIGIT_mean.gt(0)];comparison=tg.merge(net,left_on=['dataset','node'],right_on=['scope','TME_type'],how='left');save(comparison,'supplementary_TIGIT_positive_receiver_comparison.tsv');print('SUPPLEMENT DONE')
