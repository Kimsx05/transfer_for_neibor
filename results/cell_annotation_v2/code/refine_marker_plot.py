from common import *
import numpy as np,pandas as pd,anndata as ad
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
def save(fig,path):
 fig.tight_layout();fig.savefig(path,dpi=160,transparent=True,bbox_inches='tight');plt.close(fig)
for name in ['Cervilla','D1_Johnson']:
 out=ROOT/name;a=ad.read_h5ad(out/'cells.h5ad');o=a.obs
 genes=['EPCAM','KRT8','KRT19','COL1A1','DCN','LUM','COL11A1','POSTN','MMP11','CTHRC1','FAP','PI16','DPT','CCL19','RGS5','MYH11','LYZ','LST1','C1QA','FOLR2','SPP1','GPNMB','MS4A1','CD79A','CD19','BCL6','AICDA','RGS13','MZB1','JCHAIN','CD3D','CD3E','TRAC','TRBC1','TRBC2','NKG7','GNLY','PECAM1','VWF','GJA5','EFNB2','HEY1','ACKR1','PROX1']
 canon=a.var.reference_symbol.to_numpy();values=np.column_stack([np.asarray(a.X[:,canon==g].sum(1)).ravel() for g in genes]);groups=[(z,o.final_label.eq(z).to_numpy()) for z in sorted(o.final_label.unique())];fig,ax=plt.subplots(figsize=(21,max(5,len(groups)*.29)))
 dd=[]
 for i,(label,h) in enumerate(groups):
  fraction=(values[h]>0).mean(0);mean=np.log1p(values[h]/np.maximum(o.n_counts.to_numpy()[h,None],1)*1e4).mean(0);pts=ax.scatter(np.arange(len(genes)),np.repeat(i,len(genes)),s=fraction*130+1,c=mean,cmap='viridis',vmin=0,vmax=3)
  for g,f,v in zip(genes,fraction,mean):dd.append(dict(label=label,n=int(h.sum()),gene=g,detection_fraction=f,mean_log1p_10k=v))
 ax.set_xticks(range(len(genes)),genes,rotation=90,fontsize=8);ax.set_yticks(range(len(groups)),[f'{z} (n={h.sum():,})' for z,h in groups],fontsize=8);ax.invert_yaxis();ax.set_title('Measured markers: area=detection fraction; color=mean log1p(10k)');fig.colorbar(pts,ax=ax,shrink=.5);save(fig,out/'figures/marker_dotplot.png');pd.DataFrame(dd).to_csv(out/'metadata/final_label_marker_evidence.tsv',sep='\t',index=False)
 print('MARKER_PLOT',name,flush=True)
