exec(open('mp10_subtype_screen_round1/scripts/05_abundance_figures.py').read().split('comps=')[0])
sel=pd.read_csv(T/'selected_axes.tsv',sep='\t');e=pd.read_csv(T/'CellChat_detected_edges.tsv.gz',sep='\t');d=pd.read_csv(T/'candidate_axis_dataset_support.tsv',sep='\t')
q=e[e.source.isin(sel.candidate_subtype)&~e.target.str.contains('epithelial',case=False)].copy();q.to_csv(T/'candidate_TME_to_TME_detected_edges.tsv.gz',sep='\t',index=False)
q.groupby(['scope','comparison','condition','source','target','category','pathway']).agg(total_strength=('probability','sum'),n_detected_LR=('LR_ID','nunique')).reset_index().to_csv(T/'candidate_TME_to_TME_pathways.tsv',sep='\t',index=False)
z=sel[['candidate_subtype','LR_ID','source','target','LR_display']].merge(d,on=['candidate_subtype','LR_ID','source','target']);z.to_csv(T/'candidate_LR_bubble_data.tsv',sep='\t',index=False,na_rep='NA')
ds=sorted(z.dataset.unique());labs=sel.candidate_subtype+' | '+sel.LR_display;fig,ax=plt.subplots(figsize=(10,4.6));mx=max(abs(z.network_difference.dropna()).max(),1e-10)
for i,a in sel.iterrows():
 b=z[z.candidate_subtype.eq(a.candidate_subtype)]
 for j,dataset in enumerate(ds):
  c=b[b.dataset.eq(dataset)].iloc[0];v=c.network_difference
  if pd.isna(v):ax.text(j,i,'NE',ha='center',va='center',fontsize=9);continue
  n=min(c.low_n_supported_patient_ids,c.high_n_supported_patient_ids);ax.scatter(j,i,s=40+n*55,c=[v],cmap='RdBu_r',vmin=-mx,vmax=mx,edgecolors='black' if c.replicated_expression_support else '#aaaaaa',linewidths=1.2 if c.replicated_expression_support else .4)
  ax.text(j,i+.22,f'{int(c.low_n_supported_patient_ids)}/{int(c.high_n_supported_patient_ids)}',ha='center',fontsize=8)
ax.set_xticks(range(len(ds)),ds,rotation=20,ha='right');ax.set_yticks(range(len(sel)),labs);ax.set_ylim(len(sel)-.45,-.5);ax.set_title('Selected LR: dataset network direction and sample support\nColor: high − low inferred strength; labels: supported mapped IDs low/high')
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
fig.colorbar(ScalarMappable(norm=Normalize(-mx,mx),cmap='RdBu_r'),ax=ax,shrink=.7,label='Network difference (exploratory)');fig.tight_layout();save(fig,'ALL','candidate_LR_dataset_support_bubbles')
print('Candidate TME context edges:',len(q))
