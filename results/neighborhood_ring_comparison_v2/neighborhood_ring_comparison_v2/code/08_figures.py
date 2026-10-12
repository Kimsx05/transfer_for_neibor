from common import *
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm,to_hex
from matplotlib.lines import Line2D
import colorsys

def main():
 m=pd.read_parquet(OUT/'03_features/common_anchors.parquet');alllabs={};palette={};records=[]
 configs=['R1_KM','R2_KM','R12_KM','R1_LEIDEN','R2_LEIDEN','R12_LEIDEN']
 for conf in configs:
  n,a=conf.rsplit('_',1);d=OUT/'04_clustering'/n/a;ll=pd.read_parquet(d/'labels.parquet');alllabs[conf]=ll.cluster_label_native_application.to_numpy();k=int(max(alllabs[conf])+1)
  colors=['#2878B5','#D95319'] if k==2 else [to_hex(colorsys.hsv_to_rgb((i*.61803398875)%1,.55+.15*(i%3),.7+.13*(i%2))) for i in range(k)]
  cmap=ListedColormap(colors);norm=BoundaryNorm(np.arange(k+1)-.5,k);palette[conf]=(cmap,norm)
  for c,col in enumerate(colors):records.append(dict(config=conf,cluster=c,color=col,meaning='configuration-specific ID; no cross-configuration homology implied'))
  fig,ax=plt.subplots(figsize=(8,max(2,.28*np.ceil(k/4))));ax.axis('off');handles=[Line2D([],[],marker='o',ls='',color=col,label=str(c)) for c,col in enumerate(colors)];ax.legend(handles=handles,ncol=4,loc='center',title=conf+' cluster IDs',frameon=False);fig.tight_layout();fig.savefig(d/'cluster_legend.png',dpi=150);plt.close(fig)
 tab(OUT/'05_comparison/cluster_palette.tsv',records)
 for n in ['R1','R2','R12']:
  od=OUT/'04_clustering'/n;u=np.load(od/'umap_reference.npy');fig,axes=plt.subplots(1,2,figsize=(16,7))
  for ax,a in zip(axes,['KM','LEIDEN']):
   conf=n+'_'+a;labs=np.load(od/a/'reference_native_labels.npy');cmap,norm=palette[conf];ax.scatter(u[:,0],u[:,1],c=labs,s=2,cmap=cmap,norm=norm,rasterized=True)
   for c in np.unique(labs):
    center=np.median(u[labs==c],axis=0);ax.text(*center,str(c),fontsize=8,ha='center',va='center',bbox=dict(boxstyle='round,pad=.12',fc='white',alpha=.85,ec='gray',lw=.3))
   ax.set_title(conf+f' | K={len(set(labs))} | same reference UMAP');ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2')
  fig.suptitle('Native reference labels | numbers are categorical IDs | UMAP not used for clustering');fig.tight_layout();fig.savefig(od/'UMAP_comparison.png',dpi=170);plt.close(fig)
 for (ds,s),mm in m.groupby(['dataset_id','section_id']):
  od=OUT/'05_comparison'/ds;od.mkdir(exist_ok=True);ix=mm.index.to_numpy();fig,axes=plt.subplots(2,3,figsize=(15,10))
  for ax,conf in zip(axes.ravel(),configs):
   cmap,norm=palette[conf];ax.scatter(mm.array_col,mm.array_row*np.sqrt(3),c=alllabs[conf][ix],s=4,cmap=cmap,norm=norm);ax.set_aspect('equal');ax.invert_yaxis();ax.set_title(conf+f' | K={cmap.N}');ax.set_xlabel('array_col');ax.set_ylabel('sqrt(3) array_row')
  fig.suptitle(f'{ds} / {s} | anchors={len(mm)} | cnv_candidate_v1\nCategorical colors: see cluster_palette.tsv and per-configuration cluster_legend.png');fig.tight_layout();fig.savefig(od/(s+'__six_schemes.png'),dpi=140);plt.close(fig)
 log('Categorical palettes and shared UMAP/spatial plots finalized')
if __name__=='__main__':main()
