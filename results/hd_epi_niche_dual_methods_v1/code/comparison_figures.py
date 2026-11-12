from pathlib import Path
import pandas as pd,numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch,Rectangle
font_manager.fontManager.addfont('/home/data/t070721/.local/share/fonts/corefonts-arial/Arial.TTF');plt.rcParams.update({'font.family':'Arial','font.size':12,'pdf.fonttype':42})
R=Path(__file__).resolve().parents[1];F=R/'figures'
for ref in ['blca','hpca']:
 for method in ['A','B']:
  a=pd.read_parquet(R/f'{ref}_r50/labels.parquet');b=pd.read_parquet(R/f'{ref}_r200/labels.parquet');ok=a.evaluable&b.evaluable;x=a.loc[ok,method+'_niche'];y=b.loc[ok,method+'_niche']
  # Display aggregation only; full label transition table retained.
  topx=x.value_counts().head(10).index;topy=y.value_counts().head(10).index;x=x.where(x.isin(topx),'Other IDs (display)');y=y.where(y.isin(topy),'Other IDs (display)');c=pd.crosstab(x,y)/len(x);lv=c.sum(1);rv=c.sum(0);lpos={};rpos={};off=0
  for k,v in lv.items():lpos[k]=off;off+=v+.015
  off=0
  for k,v in rv.items():rpos[k]=off;off+=v+.015
  fig,ax=plt.subplots(figsize=(9,7));lo=lpos.copy();ro=rpos.copy()
  for i,lk in enumerate(c.index):
   color=plt.get_cmap('tab20')(i%20)
   for rk in c.columns:
    v=c.loc[lk,rk]
    if v<=0:continue
    p,q=lo[lk],ro[rk];verts=[(.2,p),(.45,p),(.55,q),(.8,q),(.8,q+v),(.55,q+v),(.45,p+v),(.2,p+v),(.2,p)];codes=[MPath.MOVETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.LINETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.CLOSEPOLY];ax.add_patch(PathPatch(MPath(verts,codes),facecolor=color,alpha=.55,linewidth=0));lo[lk]+=v;ro[rk]+=v
  for k,v in lv.items():ax.text(.18,lpos[k]+v/2,str(k),ha='right',va='center',fontsize=9)
  for k,v in rv.items():ax.text(.82,rpos[k]+v/2,str(k),ha='left',va='center',fontsize=9)
  ax.text(.2,-.05,'50 μm',ha='center');ax.text(.8,-.05,'200 μm',ha='center');ax.set_xlim(0,1);ax.set_ylim(-.08,max(max(lo.values()),max(ro.values()))+.03);ax.axis('off');stem=F/f'scale_flow_{ref}_{method}';fig.savefig(str(stem)+'.png',dpi=300,bbox_inches='tight');fig.savefig(str(stem)+'.pdf',bbox_inches='tight');plt.close(fig)
# Candidate evidence: all candidates as summary, not one complete figure set each.
fig,axs=plt.subplots(2,4,figsize=(16,7))
for j,(ref,r) in enumerate([(f,r) for f in ['blca','hpca'] for r in [50,200]]):
 d=pd.read_csv(R/f'{ref}_r{r}/candidate_metrics.tsv',sep='\t')
 for i,m in enumerate(['A','B']):
  z=d[d.method==m];ax=axs[i,j]
  for col in ['seed_ARI','block_ARI','silhouette']:ax.plot(z.parameter,z[col],marker='o',label=col,lw=1.1)
  ax.set_xlabel(f'{ref} {r} μm, '+('K' if m=='A' else 'resolution'));ax.set_ylabel(m+' diagnostic')
handles,labels=axs[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',ncol=3);stem=F/'parameter_stability';fig.savefig(str(stem)+'.png',dpi=300,bbox_inches='tight');fig.savefig(str(stem)+'.pdf',bbox_inches='tight');fig.legends[0].remove();fig.savefig(str(stem)+'_panel_only.pdf',bbox_inches='tight');plt.close(fig);lf,ax=plt.subplots(figsize=(5,2));ax.axis('off');ax.legend(handles,labels,frameon=False);lf.savefig(str(stem)+'_legend_only.pdf',bbox_inches='tight');plt.close(lf)
