from common import *
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
tr=read('treg_total_summary.tsv');g=read('global_incoming_outgoing.tsv');c=read('epi_allTME_total_contrasts.tsv');cov=read('node_coverage.tsv');ord=read('node_order.tsv');lr=read('treg_LR_contribution_decomposition.tsv.gz');pc=read('treg_partner_categories.tsv');partners=read('treg_partner_breakdown.tsv');colors=dict(ord[['node','color']].drop_duplicates().values);blue=colors[EPI[0]];red=colors[EPI[1]];figs=[];plotlr=[]
def put(fig,ds,name,question,source,preview=False):
 folder=O/'figures'/('overall' if ds=='ALL' else 'by_dataset/'+ds)
 for ext in ['pdf','png']:fig.savefig(folder/(name+'.'+ext),dpi=300,transparent=True,bbox_inches='tight')
 if preview:fig.savefig(O/'logs'/('preview_'+ds+'_'+name+'.png'),dpi=140,transparent=False,facecolor='white',bbox_inches='tight')
 figs.append(dict(path=str((folder/(name+'.pdf')).relative_to(O)),scope=ds,network_level='fine',question=question,measure='W / C / share',self_loops='explicit in panels',data_source=source));plt.close(fig)
for r in tr.itertuples():
 ds=r.scope;cv=cov[cov.object_id.eq(r.object_id)].set_index('node');a=c[c.object_id.eq(r.object_id)&c.node.eq('Treg')].set_index('direction');cats=pc[pc.scope.eq(ds)];fig,axs=plt.subplots(2,3,figsize=(17,10))
 x=np.arange(2);axs[0,0].bar(x-.18,[r.incoming_all_W,r.outgoing_all_W],.36,label='Including self',color='#444444');axs[0,0].bar(x+.18,[r.incoming_other_W,r.outgoing_other_W],.36,label='Excluding self',color='#aaaaaa');axs[0,0].set(xticks=x,xticklabels=['Incoming','Outgoing'],ylabel='Total W',title='One shared Treg node: global totals');axs[0,0].legend(fontsize=8)
 catorder=[EPI[0],EPI[1],'other_TME','Treg_self'];catcolor=[blue,red,'#999999',colors['Treg']];bottom=np.zeros(2)
 for name,col in zip(catorder,catcolor):
  vals=np.array([cats.loc[cats.direction.eq(dr)&cats.category.eq(name),'W'].sum(min_count=1) for dr in ['incoming','outgoing']]);axs[0,1].bar(x,vals,bottom=bottom,color=col,label=name);bottom+=vals
 axs[0,1].set(xticks=x,xticklabels=['Incoming','Outgoing'],ylabel='W',title='Where global Treg totals come from / go');axs[0,1].legend(fontsize=7)
 dirs=['epi_to_TME','TME_to_epi'];lab=['Epi → Treg','Treg → Epi']
 for j,dr in enumerate(dirs):
  z=a.loc[dr];axs[0,2].plot([j-.16,j+.16],[z.low_W,z.high_W],c='gray',lw=1);axs[0,2].scatter(j-.16,z.low_W,c=blue,s=50);axs[0,2].scatter(j+.16,z.high_W,c=red,s=50);axs[0,2].text(j,max(z.low_W,z.high_W),f'\nC {int(z.low_C)} → {int(z.high_C)}',ha='center',va='bottom',fontsize=9);axs[1,0].plot([j-.16,j+.16],[z.low_share*100,z.high_share*100],c='gray',lw=1);axs[1,0].scatter(j-.16,z.low_share*100,c=blue,s=50);axs[1,0].scatter(j+.16,z.high_share*100,c=red,s=50)
 axs[0,2].set(xticks=x,xticklabels=lab,ylabel='W',title='Only epithelial aggregate edges vary by state');axs[0,2].margins(y=.2);axs[1,0].set(xticks=x,xticklabels=['Share of epi output','Share of epi input'],ylabel='Communication share (%)',title='Fraction among the same computable TME nodes')
 ranks=[r.rank_incoming_all_W,r.rank_outgoing_all_W,r.rank_incoming_other_W,r.rank_outgoing_other_W];axs[1,1].scatter(ranks,range(4),c='black',s=45);axs[1,1].set(yticks=range(4),yticklabels=['Incoming incl. self','Outgoing incl. self','Incoming excl. self','Outgoing excl. self'],xlim=(.5,r.n_TME_ranked+.5),xlabel=f'Rank among {r.n_TME_ranked} computable TME nodes (1 = largest)',title='Global Treg rank');axs[1,1].invert_yaxis()
 axs[1,2].axis('off');lines=['Coverage (recorded IDs are unverified)']
 for node in EPI+['Treg']:
  z=cv.loc[node];lines += [node,f'  cells raw / sampled: {int(z.n_raw_cells)} / {int(z.n_sampled_cells)}',f'  specimens / recorded IDs: {int(z.n_contributing_specimens)} / {int(z.n_contributing_recorded_IDs)}']
 lines += [f'Matched epi-low + epi-high + Treg ≥30:',f'  {int(r.n_matched_specimens)} specimens / {int(r.n_matched_recorded_IDs)} IDs','Blue = epi low; red = epi high'];axs[1,2].text(0,1,'\n'.join(lines),va='top',fontsize=10)
 fig.suptitle(ds+': Treg in the complete fine network — global totals are not state-specific',fontsize=15);fig.tight_layout();put(fig,ds,'Treg_reading_panel','E_Treg_global_vs_epi_share','treg_total_summary.tsv; treg_partner_categories.tsv; node_coverage.tsv',ds in ['ALL','InhouseData'])
 nodes=ord[ord.network_level.eq('fine')].sort_values('node_order').node.tolist();fig,axs=plt.subplots(1,2,figsize=(16,len(nodes)*.28),sharey=True)
 for j,dr in enumerate(['incoming','outgoing']):
  z=partners[partners.scope.eq(ds)&partners.direction.eq(dr)].set_index('partner').reindex(nodes);axs[j].barh(range(len(nodes)),z.W,color=[colors[n] for n in nodes]);axs[j].set(yticks=range(len(nodes)),yticklabels=[n+(' [NE]' if not z.loc[n,'computed'] else '') for n in nodes],xlabel='W',title='All partners → Treg' if dr=='incoming' else 'Treg → all partners');axs[j].tick_params(axis='y',labelsize=8)
 axs[0].invert_yaxis();fig.suptitle(ds+': every Treg partner, including self; no topN filtering');fig.tight_layout();put(fig,ds,'Treg_all_partners','Treg_partner_breakdown','treg_partner_breakdown.tsv')
 fig,axs=plt.subplots(1,2,figsize=(16,7))
 for j,dr in enumerate(dirs):
  z=lr[lr.scope.eq(ds)&lr.direction.eq(dr)].copy();p=z[z.contribution.gt(0)].nlargest(5,'contribution');n=z[z.contribution.lt(0)].nsmallest(5,'contribution');keep=set(p.LR_id)|set(n.LR_id);other=z[~z.LR_id.isin(keep)];a2=pd.concat([p[['LR_id','contribution']],n[['LR_id','contribution']],pd.DataFrame({'LR_id':['Other positive LR (sum)','Other negative LR (sum)'],'contribution':[other.loc[other.contribution.gt(0),'contribution'].sum(),other.loc[other.contribution.lt(0),'contribution'].sum()]})]).sort_values('contribution');assert np.isclose(a2.contribution.sum(),z.contribution.sum());a2['scope']=ds;a2['direction']=dr;plotlr.append(a2)
  axs[j].barh(range(len(a2)),a2.contribution,color=np.where(a2.contribution.ge(0),red,blue));axs[j].axvline(0,c='black',lw=.6);axs[j].set(yticks=range(len(a2)),yticklabels=a2.LR_id,xlabel='Supported W contribution: high − low',title=lab[j]+f'  total ΔW={z.contribution.sum():.4g}');axs[j].tick_params(axis='y',labelsize=8)
 fig.suptitle(ds+': LR contributions explain the aggregate change; both signs and remaining sum retained');fig.tight_layout();put(fig,ds,'Treg_LR_contributions','Treg_all_LR_decomposition','treg_LR_contribution_decomposition.tsv.gz')
save(pd.concat(plotlr),'treg_LR_contribution_plot_rows.tsv');save(figs,'figure_index_treg.tsv');print('TREG FIGURES DONE')
