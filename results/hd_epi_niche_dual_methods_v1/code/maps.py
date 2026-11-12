from pathlib import Path
import numpy as np,pandas as pd,json
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import font_manager
font_manager.fontManager.addfont('/home/data/t070721/.local/share/fonts/corefonts-arial/Arial.TTF')
R=Path(__file__).resolve().parents[1];F=R/'figures';sections=['Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2'];types=json.loads((R/'type_mapping.json').read_text());labels={(ref,r):pd.read_parquet(R/f'{ref}_r{r}/labels.parquet') for ref in types for r in [50,200]};comp=pd.read_csv(R/'tables/niche_composition.tsv',sep='\t',dtype={'niche':str}).to_dict('records')
plt.rcParams.update({'font.family':'Arial','font.size':12,'pdf.fonttype':42})
def export(fig,stem,handles=None):
 fig.savefig(F/(stem+'.png'),dpi=300,bbox_inches='tight');fig.savefig(F/(stem+'.pdf'),bbox_inches='tight')
 if handles:
  for ax in fig.axes:
   leg=ax.get_legend()
   if leg:leg.remove()
  for leg in list(fig.legends):leg.remove()
  fig.savefig(F/(stem+'_panel_only.pdf'),bbox_inches='tight');
  from matplotlib.backends.backend_pdf import PdfPages
  with PdfPages(F/(stem+'_legend_only.pdf')) as pdf:
   for start in range(0,len(handles),100):
    lf,la=plt.subplots(figsize=(12,10));la.axis('off');la.legend(handles=handles[start:start+100],loc='center',frameon=False,ncol=4,fontsize=9);pdf.savefig(lf,bbox_inches='tight');plt.close(lf)
 plt.close(fig)
for s in sections:
 nodes=pd.read_parquet(R/s/'nodes.parquet');a=nodes[nodes.epi_anchor]
 for ref in types:
  fig,axs=plt.subplots(2,2,figsize=(12,11));handles=[]
  for row,method in enumerate(['A','B']):
   for col,radius in enumerate([50,200]):
    ax=axs[row,col];z=labels[ref,radius];z=z[z.section==s].set_index('bin_id').loc[a.bin_id];labs=sorted(set(labels[ref,radius][method+'_niche'])-{'unresolved'});pal=plt.get_cmap('hsv');colors={v:pal(((i*.61803398875)+(0 if method=='A' else .23+radius/1000+(0.11 if ref=='hpca' else 0)))%1) for i,v in enumerate(labs)};
    if method=='A' and len(labs)==2:
     cc=[v for v in comp if v['reference']==ref and v['radius']==radius and v['method']=='A' and v['niche']!='unresolved'];high=max(cc,key=lambda v:v['Epithelial'])['niche'];colors={k:('#f5b264' if k==high else '#8752a3') for k in labs}
    colors['unresolved']='#444444'
    ax.scatter(nodes.x_um,nodes.y_um,s=.12,c='#e8e0d0',rasterized=True);ax.scatter(a.x_um,a.y_um,c=[colors[x] for x in z[method+'_niche']],s=.4,rasterized=True);ax.set_aspect('equal');ax.invert_yaxis();ax.set_xlabel(f'{method}, r={radius} μm: x (μm)');ax.set_ylabel('Relative y (μm)' if s.startswith('MPBC') else 'y (μm)')
    hh=[Line2D([],[],marker='o',linestyle='',color=colors[k],label=f'{method} r{radius}: {k}') for k in labs];handles+=hh
  handles+=[Line2D([],[],marker='o',linestyle='',color='#444444',label='Unresolved anchor'),Line2D([],[],marker='o',linestyle='',color='#e8e0d0',label='Supplied background (not niche)')];fig.legend(handles=[h for prefix in ['A r50:', 'A r200:', 'B r50:', 'B r200:'] for h in [v for v in handles if v.get_label().startswith(prefix)][:6]]+handles[-2:],loc='center left',bbox_to_anchor=(1,.5),fontsize=8,title='Up to6 keys/panel; complete legend in PDF');export(fig,f'{s}_{ref}_methods_scales',handles)
