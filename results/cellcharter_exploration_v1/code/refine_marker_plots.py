from common import *
import pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
m=pd.read_parquet(ROOT/'metadata/measured_marker_evidence.parquet')
for sec in ['Cervilla','D1_Johnson']:
 z=m[m.section_id.eq(sec)];gs=list(dict.fromkeys(['EPCAM','KRT5','KRT13','KRT14','TRAC','CD3D','MS4A1','CD79A','JCHAIN','MZB1','C1QA','LYZ','DCN','LUM','COL1A1','CTHRC1','RGS5','MYH11','VWF','PLVAP','CPA3']+z[z.representative].sort_values('difference_log2',ascending=False).groupby('niche').head(2).gene.tolist()));gs=[g for g in gs if z[z.gene.eq(g)].niche.nunique()==z.niche.nunique()];v=z.pivot(index='niche',columns='gene',values='mean_log1p10k').reindex(columns=gs);d=z.pivot(index='niche',columns='gene',values='detection').reindex(columns=gs);assert not v.isna().any().any();fig,ax=plt.subplots(figsize=(max(10,.3*len(gs)),max(4,.4*len(v))));xx,yy=np.meshgrid(range(len(gs)),range(len(v)));im=ax.scatter(xx.ravel(),yy.ravel(),s=80*d.to_numpy().ravel(),c=v.to_numpy().ravel(),cmap='viridis');ax.set_xticks(range(len(gs)),gs,rotation=90);ax.set_yticks(range(len(v)),v.index);ax.invert_yaxis();ax.set_title(sec+' observed markers: size=detection fraction; color=mean log1p10k');fig.colorbar(im,ax=ax);fig.savefig(ROOT/f'figures/{sec}_marker_heatmap.png',dpi=180,bbox_inches='tight',transparent=True);plt.close(fig)
