from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from PIL import Image
import subprocess
for style in ['Regular','Bold']:
 path=subprocess.check_output(['fc-match','-f','%{file}','Arial:style='+style],text=True)
 assert font_manager.FontProperties(fname=path).get_name()=='Arial'
 font_manager.fontManager.addfont(path)
plt.rcParams.update({'font.family':'Arial','font.size':12,'pdf.fonttype':42,'axes.linewidth':1.1,'lines.linewidth':1.2,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
index=[]
def export(fig,name,ax,bar=None,handles=None):
 p=directory/name
 fig.savefig(str(p)+'.pdf',bbox_inches='tight');fig.savefig(str(p)+'.png',dpi=300,bbox_inches='tight')
 if bar:
  bar.ax.set_visible(False);fig.savefig(str(p)+'_panel_only.pdf',bbox_inches='tight');ax.set_visible(False);bar.ax.set_visible(True)
  fig.savefig(str(p)+'_legend_only.pdf',bbox_inches='tight')
 if handles:
  ax.get_legend().remove();fig.savefig(str(p)+'_panel_only.pdf',bbox_inches='tight')
  lg=plt.figure(figsize=(3.8,1.3));lg.legend(handles=handles,loc='center',frameon=False)
  lg.savefig(str(p)+'_legend_only.pdf',bbox_inches='tight');plt.close(lg)
 plt.close(fig)
 index.append(dict(section_id=sec,plot=name,png=str(p.relative_to(O))+'.png',pdf=str(p.relative_to(O))+'.pdf'))
def smooth(x,y):
 grid=np.linspace(np.min(x),np.max(x),140);bw=.08*np.ptp(x)
 w=np.exp(-.5*((grid[:,None]-x[None,:])/bw)**2)
 return grid,w@y/w.sum(axis=1)
def scatter(d,col,label):
 fig,ax=plt.subplots(figsize=(5.4,4.7))
 b=d[d['query']].dropna(subset=['EpiMP10_Top20_AUCell',col])
 x=b.EpiMP10_Top20_AUCell.to_numpy();y=b[col].to_numpy()
 ax.scatter(x,y,s=7,alpha=.38,color='#22567a',edgecolors='none');ax.plot(*smooth(x,y),color='#b2182b',lw=1.8)
 ax.set(xlabel='MP10 score',ylabel=label);export(fig,'scatter_'+col,ax)
for sec,r in M.iterrows():
 directory=O/'figures'/sec;directory.mkdir(exist_ok=True)
 d=read(O/f'spots/{sec}.tsv.gz');image=np.asarray(Image.open(r.image_path).convert('RGB'))
 xlim=(max(-.5,d.x.min()-.04*d.x.max()),min(image.shape[1]-.5,d.x.max()+.04*d.x.max()))
 ylim=(min(image.shape[0]-.5,d.y.max()+.04*d.y.max()),max(-.5,d.y.min()-.04*d.y.max()))
 def spatial():
  fig,ax=plt.subplots(figsize=(6.2,5.1));ax.imshow(image,alpha=.45,interpolation='bilinear')
  ax.set(xlim=xlim,ylim=ylim,aspect='equal');ax.axis('off');return fig,ax
 fig,ax=spatial();ax.images[0].set_alpha(1);export(fig,'tissue',ax)
 for col,label in [('EpiMP10_Top20_AUCell','MP10 score'),('H05_support','Epi high mapping support'),('H05_conditional_support','Conditional high support'),('Epi_support','Epi mapping support'),('Epithelial_fraction','c2l Epi fraction'),('nCount_RNA','nCount (UMI)')]:
  fig,ax=spatial()
  data=d[d['query']] if col=='EpiMP10_Top20_AUCell' else d
  ax.scatter(d.x,d.y,s=4,linewidths=0,color='#cccccc')
  data=data[data[col].notna()]
  lim={'vmin':0,'vmax':1} if col in ['Epi_support','H05_conditional_support','H05_support','Epithelial_fraction'] else {}
  if col=='H05_support':lim={'vmin':0,'vmax':data[col].max()}
  im=ax.scatter(data.x,data.y,c=data[col],s=8 if len(d)<2500 else 4,cmap='viridis',linewidths=0,**lim)
  cb=fig.colorbar(im,ax=ax,shrink=.74,pad=.02);cb.set_label(label);cb.solids.set_rasterized(False)
  export(fig,col,ax,bar=cb)
 fig,ax=spatial()
 ax.scatter(d.x,d.y,color='#cccccc',s=5,linewidths=0)
 small=d[d.low_denominator]
 ax.scatter(small.x,small.y,color='#d95f02',s=6,linewidths=0)
 handles=[Line2D([],[],marker='o',linestyle='',color=c,label=l) for c,l in [('#d95f02','Epi support < 0.10'),('#cccccc','Epi support ≥ 0.10')]]
 ax.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,-.01),frameon=False)
 export(fig,'low_denominator',ax,handles=handles)
 fig,ax=spatial()
 handles=[]
 for name,color in [('MP10 high / support low','#d73027'),('MP10 low / support high','#2166ac')]:
  g=d[d.discordance==name];ax.scatter(g.x,g.y,s=14,linewidths=0,color=color)
  handles.append(Line2D([],[],marker='o',color=color,linestyle='',label=name))
 ax.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,-.01),frameon=False)
 export(fig,'discordance',ax,handles=handles)
 scatter(d,'H05_support','Epi high mapping support');scatter(d,'H05_conditional_support','Conditional high support')
save(index,'figure_index.tsv')

