from common import *
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import re
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
s=read('samples_v2.tsv');p=read('proportions_v2.tsv.gz').merge(s[['sample_uid','sample_MP10_group_allEpi','MP10_sample_allEpi','n_parent_B','T_category','tissue_status']],on='sample_uid');r=read('proportion_models_v2.tsv');colors=dict(read('cell_colors.tsv')[['cell_type','color']].values);colors.update(Epi_MP10_low='#4575b4',Epi_MP10_high='#d73027');pd.DataFrame({'cell_type':list(colors),'color':list(colors.values())}).to_csv(T/'cell_colors.tsv',sep='\t',index=False)
def safe(x):return re.sub(r'[^a-zA-Z0-9_-]+','_',str(x))
def figsave(fig,ds,name):
 dest=O/'figures'/('overall' if ds=='ALL' else 'by_dataset/'+ds);dest.mkdir(parents=True,exist_ok=True)
 for ext in ['png','pdf']:fig.savefig(dest/(name+'.'+ext),dpi=300,transparent=True,bbox_inches='tight')
 plt.close(fig)
# Limited B marker check: specimen-equal mean and detection; missing vs zero distinguished.
b=read('B_marker_sample_support.tsv');markers=read('B_marker_gene_coverage.tsv').gene.tolist();a=b.groupby(['dataset','subtype','gene']).agg(mean=('mean','mean'),detection=('detection','mean'),n_specimens=('sample_uid','nunique')).reset_index();save(a,'B_marker_dataset_support.tsv');groups=a[['dataset','subtype']].drop_duplicates().sort_values(['dataset','subtype']);fig,ax=plt.subplots(figsize=(15,9));mx=a['mean'].max()
for i,row in enumerate(groups.itertuples(index=False)):
 aa=a[a.dataset.eq(row.dataset)&a.subtype.eq(row.subtype)].set_index('gene').reindex(markers)
 ax.scatter(np.arange(len(markers)),np.full(len(markers),i),s=aa.detection*130,c=aa['mean'],cmap='viridis',norm=Normalize(0,mx),edgecolors='gray',linewidths=.2)
ax.set(yticks=range(len(groups)),yticklabels=[f'{x.dataset}: {x.subtype} (cells={int(b.loc[b.dataset.eq(x.dataset)&b.subtype.eq(x.subtype), "n_cells"].sum()/len(markers))}; specimens={b.loc[b.dataset.eq(x.dataset)&b.subtype.eq(x.subtype), "sample_uid"].nunique()})' for x in groups.itertuples()],xticks=range(len(markers)),xticklabels=markers,title='Frozen B marker panel; specimen-equal expression; size = detection; no relabeling');ax.tick_params(axis='x',rotation=60);fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,mx),cmap='viridis'),ax=ax,label='Mean RNA:data');fig.tight_layout();figsave(fig,'ALL','B_marker_comparability')
print('Focus figures done')
