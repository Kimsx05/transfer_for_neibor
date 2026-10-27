from common import *
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,ListedColormap
from matplotlib.lines import Line2D
plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'savefig.transparent':True})
obj=read('object_index.tsv');order=read('node_order.tsv');edge=read('total_node_edges.tsv.gz');g=read('global_incoming_outgoing.tsv');c=read('epi_allTME_total_contrasts.tsv');cross=read('cross_dataset_direction_long.tsv');cov=read('node_coverage.tsv');figs=[]
colors=dict(order[['node','color']].drop_duplicates().values);blue=colors[EPI[0]];red=colors[EPI[1]]
def put(fig,ds,lev,name,question,source,loops=False,preview=False):
 folder=O/'figures'/('overall' if ds=='ALL' else 'by_dataset/'+ds);folder.mkdir(parents=True,exist_ok=True)
 for ext in ['pdf','png']:fig.savefig(folder/(name+'.'+ext),dpi=300,transparent=True,bbox_inches='tight')
 if preview:fig.savefig(O/'logs'/('preview_'+ds+'_'+name+'.png'),dpi=120,transparent=False,facecolor='white',bbox_inches='tight')
 figs.append(dict(path=str((folder/(name+'.pdf')).relative_to(O)),scope=ds,network_level=lev,question=question,measure='see panels',self_loops=loops,data_source=source));plt.close(fig)
def sep(ax,nodes,lev):
 dd=order[order.network_level.eq(lev)].set_index('node');last=None
 for i,n in enumerate(nodes):
  pa=dd.loc[n,'parent']
  if last is not None and pa!=last:ax.axhline(i-.5,c='#aaaaaa',lw=.4)
  last=pa
for o in obj.itertuples():
 nd=order[order.network_level.eq(o.network_level)].sort_values('node_order');nodes=nd.node.tolist();fine=o.network_level=='fine';ee=edge[edge.object_id.eq(o.object_id)];gg=g[g.object_id.eq(o.object_id)].set_index('node');cc=c[c.object_id.eq(o.object_id)];height=max(5,len(nodes)*.27)
 # Complete directed source x target, zero is white and NE is gray.
 for metric in ['W','C']:
  mat=ee.pivot(index='source',columns='target',values=metric).reindex(index=nodes,columns=nodes);arr=mat.to_numpy();cmap=plt.get_cmap('YlGnBu').copy();cmap.set_bad('#c8c8c8');fig,ax=plt.subplots(figsize=(18 if fine else 8,16 if fine else 7));im=ax.imshow(np.ma.masked_invalid(arr),cmap=cmap,vmin=0,vmax=np.nanmax(arr),aspect='equal');labels=[n+(' [NE]' if not gg.loc[n,'computed'] else '') for n in nodes];ax.set(xticks=range(len(nodes)),yticks=range(len(nodes)),xticklabels=labels,yticklabels=labels,xlabel='Target / receiver',ylabel='Source / sender',title=f'{o.scope} {o.network_level}: '+('supported total strength W' if metric=='W' else 'supported LR count C')+'\nSelf loops included; white = computed zero, gray = NE');ax.tick_params(axis='x',rotation=90,labelsize=7 if fine else 9);ax.tick_params(axis='y',labelsize=8 if fine else 9);fig.colorbar(im,ax=ax,shrink=.65,label=metric);fig.tight_layout();put(fig,o.scope,o.network_level,'heatmap_'+o.network_level+'_'+metric,'D_complete_source_target','total_node_edges.tsv.gz',True)
 # Ranked incoming and outgoing are independent, not a combined score; all nodes shown, TME rank excludes epi.
 for metric,kind in [('W','all'),('W','other'),('C','all'),('C','other')]:
  fig,axs=plt.subplots(1,2,figsize=(16,height),sharex=False)
  for j,dr in enumerate(['incoming','outgoing']):
   fld=dr+'_'+kind+'_'+metric;a=gg.reset_index().sort_values(fld,ascending=False,na_position='last');vals=a[fld].to_numpy();yy=np.arange(len(a));axs[j].barh(yy,vals,color=[colors[n] for n in a.node],edgecolor=['black' if n=='Treg' else 'none' for n in a.node],linewidth=1.8);labels=[]
   for i,row in enumerate(a.itertuples()):
    rank=getattr(row,'rank_'+fld) if 'rank_'+fld in a else np.nan;labels.append(row.node+(' [NE]' if not row.computed else f'  #{int(rank)}/{o.n_TME_nodes}' if row.node not in EPI and np.isfinite(rank) else ''))
    if not row.computed:axs[j].text(0,i,'NE',ha='left',va='center',fontsize=7)
   axs[j].set(yticks=yy,yticklabels=labels,xlabel=f'{metric}: {dr}, '+('including self' if kind=='all' else 'excluding self'),title=dr);axs[j].invert_yaxis();axs[j].tick_params(axis='y',labelsize=8 if fine else 10)
  fig.suptitle(f'{o.scope} {o.network_level}: all-node {metric} ranks; TME-only rank labels');fig.tight_layout();put(fig,o.scope,o.network_level,f'global_IO_{o.network_level}_{kind}_{metric}','E_global_IO','global_incoming_outgoing.tsv',kind=='all')
 # Role scatter with a complete numbered node key (no topN), Treg explicitly highlighted.
 for kind in ['all','other']:
  fig,(ax,keyax)=plt.subplots(1,2,figsize=(16,max(7,len(nodes)*.25)),gridspec_kw={'width_ratios':[1.5,1]});keyax.axis('off');actual=gg[gg.computed].reindex([n for n in nodes if gg.loc[n,'computed']]);xx=actual['outgoing_'+kind+'_W'];yy=actual['incoming_'+kind+'_W'];ax.scatter(xx,yy,c=[colors[n] for n in actual.index],s=60,edgecolors=['black' if n=='Treg' else 'white' for n in actual.index],linewidths=[2 if n=='Treg' else .6 for n in actual.index]);mx=max(xx.max(),1e-9);my=max(yy.max(),1e-9);pos=np.column_stack([xx.to_numpy()/mx,yy.to_numpy()/my]);labs=pos.copy()+[.018,.018]
  # Spread short numeric labels in normalized coordinates; names remain in the complete key.
  for _ in range(160):
   for i in range(len(labs)):
    for j in range(i):
     delta=labs[i]-labs[j]
     if abs(delta[0])<.045 and abs(delta[1])<.027:
      shift=.003*(1 if delta[1]>=0 else -1);labs[i,1]+=shift;labs[j,1]-=shift
  for i,n in enumerate(actual.index):ax.annotate('Treg' if n=='Treg' else str(nodes.index(n)+1),xy=(xx[n],yy[n]),xytext=(labs[i,0]*mx,labs[i,1]*my),fontsize=10 if n=='Treg' else 7,fontweight='bold' if n=='Treg' else 'normal',arrowprops=dict(arrowstyle='-',color='black',lw=.35))
  ax.set(xlabel='Outgoing W',ylabel='Incoming W',title='Self included' if kind=='all' else 'Self excluded');ax.margins(.12)
  for i,n in enumerate(nodes):keyax.text(0,1-(i+.5)/len(nodes),f'{i+1:02d}',color=colors[n],fontsize=9,fontweight='bold');keyax.text(.08,1-(i+.5)/len(nodes),n+(' [NE]' if not gg.loc[n,'computed'] else ''),fontsize=8 if fine else 11,fontweight='bold' if n=='Treg' else 'normal')
  fig.suptitle(f'{o.scope} {o.network_level}: outgoing–incoming roles, all nodes');fig.tight_layout();put(fig,o.scope,o.network_level,f'roles_{o.network_level}_{kind}','E_roles','global_incoming_outgoing.tsv',kind=='all')
 # Low/high W, count C, W difference, and communication share; every TME node in fixed order.
 tme=[n for n in nodes if n not in EPI]
 for dr in ['epi_to_TME','TME_to_epi']:
  a=cc[cc.direction.eq(dr)].set_index('node').reindex(tme);fig,axs=plt.subplots(1,4,figsize=(20,max(5,len(tme)*.29)),sharey=True);yy=np.arange(len(tme))
  for j,metric in [(0,'W'),(1,'C')]:
   axs[j].hlines(yy,a['low_'+metric],a['high_'+metric],color='#777777',lw=.7);axs[j].scatter(a['low_'+metric],yy,c=blue,s=23,label='epi low');axs[j].scatter(a['high_'+metric],yy,c=red,s=23,label='epi high');axs[j].set(title='Supported strength W' if metric=='W' else 'Supported LR count C');axs[j].legend()
  axs[2].barh(yy,a.delta_W,color=np.where(a.delta_W.ge(0),red,blue));axs[2].axvline(0,c='black',lw=.7);axs[2].set(title='W high − low')
  axs[3].hlines(yy,a.low_share*100,a.high_share*100,color='#777777',lw=.7);axs[3].scatter(a.low_share*100,yy,c=blue,s=23);axs[3].scatter(a.high_share*100,yy,c=red,s=23);axs[3].set(title='Share of all TME (%)')
  for ax in axs:
   sep(ax,tme,o.network_level);ax.set(yticks=yy,yticklabels=[n+(' [NE]' if not a.loc[n,'computed'] else '') for n in tme]);ax.tick_params(axis='y',labelsize=8 if fine else 11)
  axs[0].invert_yaxis();fig.suptitle(f'{o.scope} {o.network_level}: {dr}; shared object, all TME');fig.tight_layout();put(fig,o.scope,o.network_level,f'epi_compare_{o.network_level}_{dr}','D_epi_totals_difference_share','epi_allTME_total_contrasts.tsv')
 # Self contribution, explicitly separate from incoming/outgoing.
 a=gg.reindex(nodes);fig,axs=plt.subplots(1,2,figsize=(15,height),sharey=True);yy=np.arange(len(nodes));axs[0].barh(yy,a.incoming_self_fraction*100,color=[colors[n] for n in nodes]);axs[1].barh(yy,a.outgoing_self_fraction*100,color=[colors[n] for n in nodes]);axs[0].set(title='Self / total incoming W (%)',yticks=yy,yticklabels=nodes);axs[1].set(title='Self / total outgoing W (%)');axs[0].invert_yaxis();fig.suptitle(f'{o.scope} {o.network_level}: self-loop contribution (NE stays missing)');fig.tight_layout();put(fig,o.scope,o.network_level,'self_fraction_'+o.network_level,'E_self_contribution','global_incoming_outgoing.tsv',True)
 print('Figures',o.object_id,flush=True)
# Direction and share matrices never average strengths across objects. ALL is a separate last column.
for lev in ['major','fine']:
 nodes=order[order.network_level.eq(lev)&~order.is_epi].sort_values('node_order').node.tolist();scopes=sorted(obj[obj.scope.ne('ALL')].scope.unique())+['ALL'];a=cross[cross.network_level.eq(lev)]
 for metric in ['delta_W','delta_share','delta_C']:
  fig,axs=plt.subplots(1,2,figsize=(17,max(5,len(nodes)*.3)),sharey=True)
  for j,dr in enumerate(['epi_to_TME','TME_to_epi']):
   z=a[a.direction.eq(dr)].pivot(index='node',columns='scope',values=metric).reindex(index=nodes,columns=scopes);n=a[a.direction.eq(dr)].pivot(index='node',columns='scope',values='n_matched_recorded_IDs').reindex(index=nodes,columns=scopes);vals=z.to_numpy();sg=np.sign(vals);cm=ListedColormap([blue,'white',red]);cm.set_bad('#c8c8c8');axs[j].imshow(np.ma.masked_invalid(sg),cmap=cm,vmin=-1,vmax=1,aspect='auto');axs[j].axvline(len(scopes)-1.5,c='black',lw=1.4)
   for i in range(len(nodes)):
    for k in range(len(scopes)):
     text='NE' if not np.isfinite(vals[i,k]) else ('+' if vals[i,k]>0 else '−' if vals[i,k]<0 else '0')+f' / {int(n.iloc[i,k])}'
     axs[j].text(k,i,text,ha='center',va='center',fontsize=7 if lev=='fine' else 10,color='black')
   axs[j].set(xticks=range(len(scopes)),xticklabels=scopes,yticks=range(len(nodes)),yticklabels=nodes,title=dr);axs[j].tick_params(axis='x',rotation=35);axs[j].tick_params(axis='y',labelsize=8 if lev=='fine' else 10);sep(axs[j],nodes,lev)
  fig.suptitle(f'{lev}: {metric} direction / matched recorded IDs; ALL is pooled, not replication');fig.tight_layout();put(fig,'ALL',lev,'direction_matrix_'+lev+'_'+metric,'D_cross_dataset_direction','cross_dataset_direction_long.tsv',preview=metric=='delta_W')
save(figs,'figure_index_global.tsv');print('GLOBAL FIGURES DONE')
