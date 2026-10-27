from common import *
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
def put(fig,ds,name):
 dest=O/'figures'/('overall' if ds=='ALL' else 'by_dataset/'+ds);dest.mkdir(parents=True,exist_ok=True)
 for ext in ['png','pdf']:fig.savefig(dest/(name+'.'+ext),dpi=300,transparent=True,bbox_inches='tight')
 plt.close(fig)
s=read('samples_v2.tsv');b=read('B_marker_sample_support.tsv');markers=read('B_marker_gene_coverage.tsv').gene.tolist();lim=b.groupby(['dataset','subtype','gene'])['mean'].mean().max()
for ds,z in b.groupby('dataset'):
 a=z.groupby(['subtype','gene']).agg(mean=('mean','mean'),detection=('detection','mean')).reset_index();types=sorted(z.subtype.unique());fig,ax=plt.subplots(figsize=(12,4))
 for i,ct in enumerate(types):
  q=a[a.subtype.eq(ct)].set_index('gene').reindex(markers);ax.scatter(np.arange(len(markers)),np.full(len(markers),i),s=q.detection*140,c=q['mean'],cmap='viridis',norm=Normalize(0,lim),edgecolor='gray',lw=.2)
 labels=[f'{ct}\ncells={int(z.loc[z.subtype.eq(ct),"n_cells"].sum()/len(markers))}; specimens={z.loc[z.subtype.eq(ct),"sample_uid"].nunique()}' for ct in types]
 ax.set(xticks=range(len(markers)),xticklabels=markers,yticks=range(len(types)),yticklabels=labels,ylim=(-.6,len(types)-.4),title=ds+': frozen B markers; size = detection, color = specimen-equal mean');ax.tick_params(axis='x',rotation=60);fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,lim),cmap='viridis'),ax=ax,label='Mean RNA:data');fig.tight_layout();put(fig,ds,'B_marker_comparability')
mat=read('matched_B_Plasma_LR_sensitivity.tsv');a=read('selected_axes_all_network_context.tsv');a=a[a.network_level.eq('major')&a.comparison_type.eq('epi_state_within_specimen')];z=mat.merge(a[['axis_id','scope','raw_delta']],on=['axis_id','scope'],suffixes=('_matched','_primary'));save(z,'matched_B_Plasma_sensitivity_comparison.tsv')
for ds in ['ALL']+sorted(s.dataset.unique()):
 zz=z if ds=='ALL' else z[z.scope.eq(ds)];zz=zz.sort_values(['axis_id','scope']);fig,ax=plt.subplots(figsize=(max(10,len(zz)*.42),5));pos=np.arange(len(zz));ax.bar(pos-.18,zz.raw_delta_primary,.36,label='Primary shared-state object');ax.bar(pos+.18,zz.raw_delta_matched,.36,label='Common B + Plasma specimens');ax.axhline(0,c='black',lw=.6);ax.set(xticks=pos,xticklabels=[f'{x.axis_id} {x.scope}\nn={x.n_patients} IDs' for x in zz.itertuples()],ylabel='Raw high - low probability',title='B/Plasma coverage sensitivity; compare direction, objects rescale expression separately');ax.tick_params(axis='x',rotation=60);ax.legend();fig.tight_layout();put(fig,ds,'B_Plasma_common_specimen_sensitivity')
