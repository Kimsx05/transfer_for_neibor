from pathlib import Path
import pandas as pd,numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
font_manager.fontManager.addfont('/home/data/t070721/.local/share/fonts/corefonts-arial/Arial.TTF');plt.rcParams.update({'font.family':'Arial','font.size':12,'pdf.fonttype':42})
R=Path(__file__).resolve().parents[1];F=R/'figures';d=pd.read_csv(R/'tables/niche_composition.tsv',sep='\t',dtype={'niche':str})
import json
types=json.loads((R/'type_mapping.json').read_text())
for (ref,r,m),z in d.groupby(['reference','radius','method']):
 cats=types[ref];c=z.set_index('niche').sort_values('n',ascending=False)
 if len(c)>32:
  tail=c.iloc[31:];other={t:np.average(tail[t],weights=tail.n) for t in cats};other['n']=tail.n.sum();c=pd.concat([c.iloc[:31],pd.DataFrame([other],index=['Other IDs (display only)'])])
 fig,ax=plt.subplots(figsize=(max(6,len(cats)*.8),max(3,len(c)*.24)));im=ax.imshow(c[cats],vmin=0,vmax=1,cmap='viridis',aspect='auto');ax.set_xticks(range(len(cats)),cats,rotation=45,ha='right');ax.set_yticks(range(len(c)),[m+' '+x for x in c.index]);cb=fig.colorbar(im,ax=ax,label='Mean known-neighbor fraction');stem=F/f'composition_{ref}_{r}_{m}';fig.savefig(str(stem)+'.png',dpi=300,bbox_inches='tight');fig.savefig(str(stem)+'.pdf',bbox_inches='tight');cb.remove();fig.savefig(str(stem)+'_panel_only.pdf',bbox_inches='tight');plt.close(fig);f,a=plt.subplots(figsize=(2,4));f.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0,1),cmap='viridis'),cax=a,label='Mean known-neighbor fraction');f.savefig(str(stem)+'_legend_only.pdf',bbox_inches='tight');plt.close(f)
