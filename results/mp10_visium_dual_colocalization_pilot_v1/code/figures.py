from shared import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib import font_manager
from PIL import Image
import subprocess
for style in ['Regular','Bold']:
 path=subprocess.check_output(['fc-match','-f','%{file}','Arial:style='+style],text=True)
 assert font_manager.FontProperties(fname=path).get_name()=='Arial'
 font_manager.fontManager.addfont(path)
plt.rcParams.update({'font.family':'Arial','font.size':11,'pdf.fonttype':42,'axes.linewidth':1.1,'axes.spines.right':False,'axes.spines.top':False,'savefig.facecolor':'white'})
idx=[]
def export(fig,name,kind,bar=None):
 p=out/name;fig.savefig(str(p)+'.pdf',bbox_inches='tight');fig.savefig(str(p)+'.png',dpi=300,bbox_inches='tight')
 if bar:
  bar.ax.set_visible(False);fig.savefig(str(p)+'_panel_only.pdf',bbox_inches='tight')
  for ax in fig.axes:ax.set_visible(ax is bar.ax)
  fig.savefig(str(p)+'_legend_only.pdf',bbox_inches='tight')
 else:
  legends=[a for a in fig.axes if a.get_legend() is not None]
  if legends:
   ax=legends[0];handles,labels=ax.get_legend_handles_labels();ax.get_legend().remove()
   fig.savefig(str(p)+'_panel_only.pdf',bbox_inches='tight')
   panel=plt.figure(figsize=(3.6,max(1,.32*len(handles))));panel.legend(handles,labels,loc='center',frameon=False)
   panel.savefig(str(p)+'_legend_only.pdf',bbox_inches='tight');plt.close(panel)
 plt.close(fig);idx.append(dict(section_id=s,kind=kind,figure=name,pdf=str(p.relative_to(O))+'.pdf',png=str(p.relative_to(O))+'.png'))
def heat(df,name,label='Pearson r',cor=True,height=None,width=None,kind='association'):
 fig,ax=plt.subplots(figsize=(width or max(4.5,min(14,2+len(df.columns)*.41)),height or max(2.7,1.7+len(df)*.29)))
 cmap=plt.get_cmap('RdBu_r' if cor else 'viridis').copy();cmap.set_bad('#e6e6e6')
 im=ax.imshow(df.to_numpy(float),aspect='auto',interpolation='none',cmap=cmap,vmin=-1 if cor else 0,vmax=1)
 ax.set_yticks(range(len(df)),df.index,fontsize=9);ax.set_xticks(range(len(df.columns)),df.columns,rotation=90 if len(df.columns)>8 else 25,ha='center' if len(df.columns)>8 else 'right',fontsize=9)
 cb=fig.colorbar(im,ax=ax,shrink=.75);cb.set_label(label);cb.solids.set_rasterized(False);export(fig,name,kind,cb)
def scatter(x,y,xlabel,ylabel,name,kind):
 ok=np.isfinite(x)&np.isfinite(y);x=np.asarray(x)[ok];y=np.asarray(y)[ok]
 fig,ax=plt.subplots(figsize=(5.2,4.3));ax.scatter(x,y,color='#255678',alpha=.33,s=7,linewidths=0)
 if len(x)>2 and np.ptp(x)>0:
  grid=np.linspace(x.min(),x.max(),140);w=np.exp(-.5*((grid[:,None]-x[None,:])/(.08*np.ptp(x)))**2);ax.plot(grid,w@y/w.sum(axis=1),color='#b2182b',lw=1.6)
 ax.set(xlabel=xlabel,ylabel=ylabel);export(fig,name,kind)
names={'Synthesis':'Cholesterol synthesis','Homeostasis':'Cholesterol homeostasis','Fatty_acid':'Fatty acid metabolism','Efflux':'Cholesterol efflux'}
def label(x):
 if x in [H,MP,N]:return {H:'Epi high mapping score',MP:'MP10 score',N:'Other Epi mapping score'}[x]
 if x.startswith(('log1p_c2l_mean__','TME_score__')):return x.split('__',1)[1]
 if '__' in x:
  f,v=x.split('__');return names[f]+(' (minus Top20)' if v=='minus20' else '')
 return x
four=read(T/'four_way.tsv');mapping=read(O/'annotation_mapping.tsv').query('subtype!="Epithelial"')
fine=list(mapping.subtype);fc=read(T/'functional_correlations.tsv');ft=read(T/'functional_TME.tsv.gz')
pal=read(R/'neighborhood_niche_visualization_v1/data_tables/niche_palette.tsv').sort_values('accepted_label')
ab={k:read(T/f'{k}_pairwise.tsv.gz') for k in ['A','B']}
for s,row in M.iterrows():
 out=O/'figures'/s;out.mkdir(exist_ok=True)
 d=read(O/f'spots/{s}.tsv.gz');e=d[d['query']]
 img=Image.open(row.image_path).convert('RGB');W,L=img.size;img.thumbnail((800,800))
 xpad=.04*np.ptp(d.x);ypad=.04*np.ptp(d.y)
 def spatial(col,ylabel,name,query_only=False,cat=False):
  fig,ax=plt.subplots(figsize=(6.2,5));ax.imshow(img,extent=(-.5,W-.5,L-.5,-.5),alpha=.5)
  ax.set(xlim=(d.x.min()-xpad,d.x.max()+xpad),ylim=(d.y.max()+ypad,d.y.min()-ypad),aspect='equal');ax.axis('off')
  sub=e if query_only else d;ax.scatter(d.x,d.y,s=4,color='#cccccc',linewidths=0)
  sub=sub[sub[col].notna()];kw={}
  if cat:kw=dict(cmap=ListedColormap(pal.color),norm=BoundaryNorm(np.arange(.5,10.5),9))
  elif col==H or col.startswith('TME_score'):kw=dict(vmin=0,vmax=sub[col].max() if sub[col].max()>0 else 1,cmap='viridis')
  else:kw=dict(cmap='viridis',vmin=0)
  im=ax.scatter(sub.x,sub.y,c=sub[col],s=8 if len(d)<2500 else 4,linewidths=0,**kw)
  cb=fig.colorbar(im,ax=ax,shrink=.76,pad=.02);cb.set_label(ylabel);cb.solids.set_rasterized(False)
  if cat:cb.set_ticks(range(1,10))
  export(fig,name,'spatial',cb)
 spatial(MP,'MP10 score','space_MP10',True);spatial(H,'Epi high mapping score','space_high');spatial('niche','Original niche','space_niche',True,True)
 for k in names:
  for v in ['full','minus20']:
   col=k+'__'+v;spatial(col,label(col)+' score','space_'+col)
   for x in [H,MP]:scatter(e[x],e[col],label(x),label(col)+' score','scatter_'+('high' if x==H else 'MP10')+'_'+col,'function')
 for t in ['Treg','SPP1 TAM','C1QC TAM','Matrix_CAF']:
  spatial(scol(t),t+' mapping score','space_transfer_'+t.replace(' ','_'))
  spatial(lcol(t),t+' log1p c2l mean abundance','space_c2l_'+t.replace(' ','_'))
 spatial('Epithelial_fraction','c2l Epi / all30','space_Epi');spatial('nCount_RNA','nCount (UMI)','space_depth')
 for model,tab in ab.items():
  g=tab[tab.section_id==s]
  cols=([H,N]+[scol(t) for t in fine]) if model=='A' else [MP]+[lcol(t) for t in fine]
  for metric in ['pearson','resid_pearson']:
   mat=pd.DataFrame(np.nan,index=cols,columns=cols)
   for r in g.itertuples():mat.loc[r.left,r.right]=mat.loc[r.right,r.left]=getattr(r,metric)
   mat.columns=mat.index=[label(c) for c in cols]
   heat(mat,'matrix_'+model+'_'+metric,('Residual r' if metric=='resid_pearson' else 'Pearson r')+(' · log1p c2l mean' if model=='B' else ' · mapping scores'),height=10,width=12,kind='matrix')
 for level in ['subtype','major']:
  g=four[(four.section_id==s)&(four.level==level)]
  for metric in ['pearson','resid_pearson']:
   mat=g.pivot(index='cell_type',columns='combination',values=metric).loc[:,['A','H05_c2l','MP10_transfer','B']]
   if level=='subtype':mat=mat.loc[fine]
   mat.index=[x.replace('major_','') for x in mat.index];mat.columns=['high × transfer','high × c2l','MP10 × transfer','MP10 × c2l']
   heat(mat,'four_'+level+'_'+metric,'Pearson r' if metric=='pearson' else 'Residual Pearson r',width=6.5)
 for model in ['A','B']:
  g=four[(four.section_id==s)&(four.level=='subtype')&(four.combination==model)].set_index('cell_type').loc[fine]
  fig,ax=plt.subplots(figsize=(6.2,8));y=np.arange(len(g))
  ax.hlines(y,np.minimum(g.pearson,g.resid_pearson),np.maximum(g.pearson,g.resid_pearson),color='#a8a8a8')
  ax.scatter(g.pearson,y,label='Unadjusted',color='#2166ac',s=25);ax.scatter(g.resid_pearson,y,label='Adjusted',color='#b2182b',s=15)
  ax.axvline(0,lw=.7,color='#444');ax.set(xlim=(-1,1),xlabel=('Epi high × TME mapping scores' if model=='A' else 'MP10 × log1p c2l mean')+' (Pearson r)');ax.set_yticks(y,fine);ax.invert_yaxis();ax.legend(frameon=False);export(fig,'TME_effects_'+model,'association')
 for x in [H,MP]:
  a=fc[(fc.section_id==s)&(fc.exposure==x)].set_index('module_score')[['pearson','resid_pearson','spearman','resid_spearman']]
  a.index=[label(c) for c in a.index];a.columns=['Pearson','Residual Pearson','Spearman','Residual Spearman'];heat(a,'function_'+('high' if x==H else 'MP10'),'Correlation',width=6.5,kind='function')
 for method in ['log1p_mean','transfer']:
  g=ft[(ft.section_id==s)&(ft.level=='subtype')&(ft.TME_method==method)];mat=g.pivot(index='cell_type',columns='module_score',values='resid_pearson').loc[fine];mat.columns=[label(c) for c in mat.columns]
  heat(mat,'function_TME_'+method,'Residual Pearson r',width=7.5,kind='function')
 assigned=e.dropna(subset=['niche'])
 cols=[MP,H]+[k+'__'+v for k in names for v in ['full','minus20']]
 a=assigned[cols].rank(pct=True).groupby(assigned.niche).mean().T;a.columns=['N'+str(int(i)) for i in a.columns];a.index=[label(c) for c in a.index]
 heat(a,'niche_scores','Mean within-section percentile',False,width=7,kind='niche')
 scatter(np.log1p(d.c2l_mean_TME29),d[[scol(t) for t in fine]].max(axis=1),'log1p c2l mean TME abundance','Largest TME mapping score','TME_depth_diagnostic','diagnostic')
 scatter(e.H05_support,e[H],'Joint Epi high support','Separate Epi high score','joint_vs_separate','diagnostic')
 # Keep the location of both discordant tails inspectable.
 fig,ax=plt.subplots(figsize=(5.5,5));ax.imshow(img,extent=(-.5,W-.5,L-.5,-.5),alpha=.55)
 for name,c in [('H05 high / MP10 low','#2166ac'),('H05 low / MP10 high','#d73027')]:
  a=e[e.discordance_new==name];ax.scatter(a.x,a.y,s=12,color=c,label=name,linewidths=0)
 ax.set(xlim=(d.x.min()-xpad,d.x.max()+xpad),ylim=(d.y.max()+ypad,d.y.min()-ypad));ax.axis('off');ax.legend(loc='upper center',bbox_to_anchor=(.5,-.01),frameon=False);export(fig,'discordant_locations','spatial')
 print(s,len(idx),flush=True)
s='four_sections';out=O/'figures/overview';out.mkdir(exist_ok=True)
for x,short in [(H,'high'),(MP,'MP10')]:
 z=fc[fc.exposure==x].pivot(index='module_score',columns='section_id',values='resid_pearson').loc[:,list(S)]
 z=z.loc[[k+'__'+v for k in names for v in ['full','minus20']]];z.index=[label(c) for c in z.index];z.columns=[c.split('_')[-1] for c in z.columns]
 heat(z,'four_sections_function_'+short,'Residual Pearson r',width=6.5,kind='functional comparison')
for method in ['A','B']:
 z=four[(four.level=='subtype')&(four.combination==method)].pivot(index='cell_type',columns='section_id',values='resid_pearson').loc[fine,list(S)]
 z.columns=[c.split('_')[-1] for c in z.columns];heat(z,'four_sections_TME_'+method,'Residual Pearson r',width=6,kind='TME comparison')
save(idx,'figure_index_primary.tsv')

