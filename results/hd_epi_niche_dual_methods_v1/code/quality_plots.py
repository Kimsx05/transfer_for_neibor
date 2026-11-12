from pathlib import Path
import pandas as pd,numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
font_manager.fontManager.addfont('/home/data/t070721/.local/share/fonts/corefonts-arial/Arial.TTF')
R=Path(__file__).resolve().parents[1];F=R/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':12,'pdf.fonttype':42,'figure.facecolor':'white'})
for s in ['Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2']:
 for r in [50,200]:
  d=pd.read_parquet(R/s/f'r{r}.parquet');m=pd.read_parquet(R/s/f'neighbor_marker_r{r}.parquet');assert d.bin_id.equals(m.bin_id);d=pd.concat([d,m.drop(columns='bin_id')],axis=1)
  cols=['blca_known_fraction','hpca_known_fraction','input_coverage_component','low_information_fraction','reference_conflict_fraction','neighbor_marker_epithelial','neighbor_marker_collagen','neighbor_marker_Myeloid','neighbor_marker_T/NK','neighbor_marker_Endothelial'];fig,axs=plt.subplots(2,5,figsize=(18,8))
  for ax,c in zip(axs.ravel(),cols):
   ax.scatter(d.x_um,d.y_um,c=d[c],vmin=0,vmax=1,cmap='viridis',s=.5,rasterized=True);ax.set_aspect('equal');ax.invert_yaxis();ax.set_xlabel(c.replace('_','\n'),fontsize=9);ax.set_ylabel('relative μm' if s.startswith('MPBC') else 'μm')
  # Identical color range, separate vector colorbar legend.
  import matplotlib as mpl
  cb=fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0,1),cmap='viridis'),ax=axs.ravel().tolist(),fraction=.015);cb.set_label('Fraction (denominator in METHODS)')
  stem=F/f'{s}_r{r}_support_raw_neighbor_markers';fig.savefig(str(stem)+'.png',dpi=300,bbox_inches='tight');fig.savefig(str(stem)+'.pdf',bbox_inches='tight');cb.remove();fig.savefig(str(stem)+'_panel_only.pdf',bbox_inches='tight');plt.close(fig)
  lf,la=plt.subplots(figsize=(2,4));lf.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0,1),cmap='viridis'),cax=la);lf.savefig(str(stem)+'_legend_only.pdf',bbox_inches='tight');plt.close(lf)
 print('QUALITY_PLOTS',s,flush=True)
