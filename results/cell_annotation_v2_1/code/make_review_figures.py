from common import *
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from programs import MAJOR,STATE
from expression_evidence import SUBTYPE
plt.rcParams.update({'figure.facecolor':'none','savefig.transparent':True,'axes.spines.top':False,'axes.spines.right':False,'font.size':8})
o=pd.read_parquet(ROOT/'cell_annotations.parquet');methods=json.loads((ROOT/'metadata/lineage_structure_methods.json').read_text());coords=np.load(ROOT/'metadata/existing_coordinates.npz');pc=np.load(ROOT/'metadata/lineage_pca.npz')['X_pca'];pal=json.loads((V2/'metadata/annotation_palette.json').read_text());pal.update(Unknown='#bdbdbd',Mural='#9c755f',Erythroid='#b22222',Neural_glial='#9467bd',major_unresolved='#bdbdbd',subtype_unresolved='#e6ab02',subtype_proposed='#1b9e77',mixed_or_artifact='#d95f02');labels=sorted(set(o.display_label)|set(o.proposed_major)|set(o.proposed_subtype))
for i,k in enumerate(labels):
 if k not in pal:pal[k]=pal.get(k.replace('_subtype_unknown',''),to_hex(plt.get_cmap('tab20')(i%20)))
js(ROOT/'metadata/annotation_palette.json',pal)
def save(fig,name):fig.savefig(ROOT/'figures'/name,dpi=180,bbox_inches='tight',transparent=True);plt.close(fig)
def tag(s):return s.replace('/','_')
def scatter(ax,xy,lab,title,spatial=False,background=None,size=.45):
 if background is not None:ax.scatter(background[:,0],background[:,1],s=.08,c='#e6e6e6',linewidths=0,rasterized=True)
 valid=np.isfinite(xy).all(1);ls=sorted(pd.unique(lab),key=lambda k:(k!='Unknown',k))
 for k in ls:
  h=valid&(lab==k);ax.scatter(xy[h,0],xy[h,1],s=size,c=pal.get(k,'#888888'),label=f'{k} ({h.sum():,})',linewidths=0,rasterized=True)
 ax.set_title(title,fontsize=10);ax.set_aspect('equal');ax.axis('off')
 if spatial:ax.invert_yaxis()
 ax.legend(loc='upper left',bbox_to_anchor=(1,1),markerscale=5,fontsize=7,frameon=False)
for sec in ['Cervilla','D1_Johnson']:
 ix=np.flatnonzero(o.section_id.eq(sec));z=o.iloc[ix];xy=coords['spatial_um'][ix];um=coords['X_umap'][ix];fig,ax=plt.subplots(figsize=(10,7));scatter(ax,xy,z.proposed_major.to_numpy(),f'{sec}: proposed major | all {len(ix):,} retained cells',True);save(fig,f'{sec}_major_spatial.png')
 fig,ax=plt.subplots(figsize=(10,7));scatter(ax,um,z.proposed_major.to_numpy(),f'{sec}: existing joint scVI UMAP | {np.isfinite(um[:,0]).sum():,} core only',size=1);save(fig,f'{sec}_major_core_UMAP.png')
 quality=np.where(z.proposed_major.eq('Unknown'),'major_unresolved',np.where(z.proposed_subtype.eq('Unknown'),'subtype_unresolved','subtype_proposed'));quality[z.mixed_observed_programs|z.artifact_review_flag]='mixed_or_artifact';fig,ax=plt.subplots(figsize=(10,7));scatter(ax,xy,quality,f'{sec}: identity resolution / quality | all retained cells',True);save(fig,f'{sec}_uncertainty_spatial.png')
 for lineage in methods:
  hit=z.proposed_major.eq(lineage).to_numpy();fig,axs=plt.subplots(1,2,figsize=(14,5));scatter(axs[0],xy[hit],z.display_label.to_numpy()[hit],f'{sec} {lineage}: {hit.sum():,} cells',True,xy,size=1);pcl=pc[ix[hit],:2];scatter(axs[1],pcl,z.display_label.to_numpy()[hit],f'Local expression PCA (informative {np.isfinite(pcl[:,0]).sum():,})',size=2);fig.subplots_adjust(wspace=.8);save(fig,f'{sec}_{tag(lineage)}_subtypes.png')
# States remain separate from identity, with per-section spatial views.
for sec in ['Cervilla','D1_Johnson']:
 ix=np.flatnonzero(o.section_id.eq(sec));z=o.iloc[ix];xy=coords['spatial_um'][ix]
 for lineage,states in {'Epithelial':['Epithelial_basal','Epithelial_luminal','Epithelial_squamous','Cycling','Stress','Interferon'],'Plasma':['Plasma_IgA','Plasma_IgG'],'Mast':['Mast_activation'],'Myeloid':['C1Q_program','SPP1_program','FOLR2_program','Myeloid_S100A8_A9_program']}.items():
  fig,axs=plt.subplots(1,len(states),figsize=(4*len(states),4),squeeze=False)
  for ax,state in zip(axs[0],states):
   lineagehit=z.proposed_major.eq(lineage).to_numpy();hit=z['state_'+state].to_numpy()&lineagehit;ax.scatter(xy[lineagehit,0],xy[lineagehit,1],s=.25,c='#cccccc',linewidths=0);ax.scatter(xy[hit,0],xy[hit,1],s=1.2,c='#b2182b',linewidths=0);ax.set_aspect('equal');ax.invert_yaxis();ax.axis('off');ax.set_title(f'{sec} {state}\n{hit.sum():,}/{lineagehit.sum():,} {lineage} cells')
  save(fig,f'{sec}_{tag(lineage)}_states.png')
# Marker dotplots: known identity/competition markers plus measured data-driven genes.
def dotplot(m,summary,gs,title,name):
 gs=list(dict.fromkeys(gs));gs=[g for g in gs if g in set(m.gene)];
 if not gs:return
 clusters=summary.cluster.astype(str).tolist();nmap=dict(zip(summary.cluster.astype(str),summary.n));mean=m.pivot(index='cluster',columns='gene',values='mean_log1p10k').reindex(index=clusters,columns=gs).fillna(0);det=m.pivot(index='cluster',columns='gene',values='detection_fraction').reindex(index=clusters,columns=gs).fillna(0);height=max(3,len(clusters)*.25+1.7);fig,ax=plt.subplots(figsize=(max(8,len(gs)*.28),height));xx,yy=np.meshgrid(np.arange(len(gs)),np.arange(len(clusters)));sc=ax.scatter(xx.ravel(),yy.ravel(),s=det.to_numpy().ravel()*65,c=mean.to_numpy().ravel(),cmap='viridis',vmin=0,vmax=max(.1,np.quantile(mean.to_numpy(),.98)),linewidths=0);ax.set_xticks(np.arange(len(gs)),gs,rotation=90,fontsize=7);ax.set_yticks(np.arange(len(clusters)),[f'{c} n={nmap[c]:,}' for c in clusters],fontsize=7);ax.invert_yaxis();ax.set_title(title+'\nDot size: observed fraction (0–100%); color: mean log1p(10k counts)',fontsize=10);fig.colorbar(sc,ax=ax,shrink=.45,pad=.02,label='Measured log expression');save(fig,name)
mapping=json.loads((ROOT/'metadata/reference.json').read_text())['major_mapping']
sets=[('global','global')]+[(l,'lineage_'+tag(l)) for l in methods]
for lineage,prefix in sets:
 m=pd.read_parquet(ROOT/f'metadata/{prefix}_markers.parquet');s=pd.read_csv(ROOT/f'metadata/{prefix}_summary.tsv',sep='\t');s.cluster=s.cluster.astype(str)
 known=sum([g[:3] for g in MAJOR.values()],[]) if lineage=='global' else MAJOR.get(lineage,[])+sum([g for label,g in SUBTYPE.items() if mapping.get(label)==lineage],[])
 if lineage=='Epithelial':known+=sum([g for k,g in STATE.items() if k.startswith('Epithelial_') or k in ['Cycling','Stress','Interferon']],[])
 if lineage=='Plasma':known+=sum([g for k,g in STATE.items() if k.startswith('Plasma_')],[])
 if lineage=='Mast':known+=STATE['Mast_activation']
 if lineage!='global':known+=['EPCAM','CD3D','MS4A1','MZB1','LST1','DCN','RGS5','PECAM1','TPSAB1']
 for scope in ['Cervilla','D1_Johnson']:
  ss=s[s.scope.eq(scope)].copy();mm=m[m.scope.eq(scope)].copy();
  if ss.empty:continue
  # Split long known panels into legible pages, never silently omit a lineage.
  kg=list(dict.fromkeys(known))
  for page,start in enumerate(range(0,len(kg),40),1):dotplot(mm,ss,kg[start:start+40],f'{scope} | {lineage} known / competitor markers | all assigned members',f'{scope}_{tag(lineage)}_known_markers_{page}.png')
  top=[]
  for value in ss.top_genes.fillna(''):top+=value.split(';')[:2]
  top=list(dict.fromkeys(top))
  for page,start in enumerate(range(0,len(top),40),1):dotplot(mm,ss,top[start:start+40],f'{scope} | {lineage} data-driven markers | descriptive',f'{scope}_{tag(lineage)}_data_markers_{page}.png')
 print('PLOTTED',lineage,flush=True)
print('ALL_FIGURES_COMPLETE',flush=True)
