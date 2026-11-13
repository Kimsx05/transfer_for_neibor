from common import *
import sys
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.font_manager import fontManager
from matplotlib.backends.backend_pdf import PdfPages
fontManager.addfont('/home/data/t070721/.local/share/fonts/corefonts-arial/Arial.TTF');plt.rcParams.update({'font.family':'Arial','font.size':12,'pdf.fonttype':42,'axes.linewidth':1.1})
F=R/'figures';F.mkdir(exist_ok=True);a=anchors();models=pd.read_csv(R/'tables/selected_models.tsv',sep='\t');comp=pd.read_csv(R/'tables/niche_composition.tsv',sep='\t');labels={v.tag:pd.read_parquet(R/'models'/v.tag/'assignments.parquet') for _,v in models.iterrows()}

def short(tag):
 return tag.replace('hellinger_kmeans','A').replace('adapted_leiden','B-adapted').replace('MENDER_adapted','MENDER').replace('CellCharter','CC').replace('ring_counts','rings').replace('scVI_without_Top20','scVI-exclusion')

def display_key(h,tag,section):
 if len(h)<=12:return h
 counts=labels[tag].loc[a.section==section,'label'].value_counts();keep=set(counts[counts.index>=0].head(4).index)
 return [v for v in h if int(v.get_label().rsplit(': ',1)[1]) in keep]

def export(fig,name,handles=None):
 if name.endswith('four_methods') or name.endswith('NONE_WEAK'):fig.text(.01,-.05,'All labels are plotted. Key: all groups if <=12, otherwise four largest per panel; complete key in legend-only PDF.',fontsize=10)
 if name.startswith('MPBC'):fig.text(.01,-.02,'Relative array coordinates; original tissue extent unverified',fontsize=10)
 fig.savefig(F/(name+'.png'),dpi=300,bbox_inches='tight',facecolor='white');fig.savefig(F/(name+'.pdf'),bbox_inches='tight',facecolor='white')
 if handles:
  for lg in list(fig.legends):lg.remove()
  for ax in fig.axes:
   if ax.get_legend():ax.get_legend().remove()
  fig.savefig(F/(name+'_panel_only.pdf'),bbox_inches='tight')
  with PdfPages(F/(name+'_legend_only.pdf')) as pdf:
   for st in range(0,len(handles),80):
    ff,aa=plt.subplots(figsize=(11,9));aa.axis('off');aa.legend(handles=handles[st:st+80],loc='center',frameon=False,ncol=3,fontsize=10);pdf.savefig(ff,bbox_inches='tight');plt.close(ff)
 plt.close(fig)

def palette(tag):
 labs=sorted(set(labels[tag].label)-{-1});offset=sum(tag.encode())/179.;pal={k:plt.get_cmap('hsv')((j*.61803398875+offset)%1) for j,k in enumerate(labs)};pal[-1]='#444444';return pal

def spatial(ax,tag,s,bg,region=None,reference='blca'):
 q=a.section==s
 if region is not None:
  x,y,X,Y=region;q&=a.x_um.between(x,X)&a.y_um.between(y,Y);bg=bg[bg.x_um.between(x,X)&bg.y_um.between(y,Y)]
 d=labels[tag];pal=palette(tag);ax.scatter(bg.x_um,bg.y_um,c=np.where(~bg.analysis_qc,'#f0f0f0',np.where(bg['type_'+reference].isna(),'#9fb8d0','#e8e0d0')),s=.12 if region is None else 4,rasterized=True);ax.scatter(a.loc[q,'x_um'],a.loc[q,'y_um'],c=[pal[k] for k in d.loc[q,'label']],s=.5 if region is None else 8,rasterized=True);ax.set_aspect('equal');ax.invert_yaxis();ax.set_ylabel('Relative y (μm)' if s.startswith('MPBC') else 'y (μm)')
 if region is not None:ax.set_xlim(x,X);ax.set_ylim(Y,y);ax.plot([x+10,x+74],[Y-10,Y-10],lw=1.2,c='black');ax.text(x+10,Y-18,'64 μm',fontsize=9)
 hs=[Line2D([],[],marker='o',ls='',color=pal[k],label=f'{tag.replace("hellinger_kmeans","A").replace("adapted_leiden","B-adapted").replace("MENDER_adapted","MENDER").replace("CellCharter","CC")}: {k}') for k in pal if k!=-1];return hs
for s in SECTIONS:
 bg=pd.read_parquet(V/s/'nodes.parquet',columns=['x_um','y_um','bin_id','analysis_qc','type_blca','type_hpca']);ri=pd.read_csv(V/'tables/ROI_review_index.tsv',sep='\t');ri=ri[ri.section==s]
 for ref in TYPES:
  fig,axs=plt.subplots(2,4,figsize=(20,11));handles=[];shown=[]
  for i,r in enumerate([50,200]):
   for j,(m,rep) in enumerate([('hellinger_kmeans','NONE'),('adapted_leiden','NONE'),('MENDER_adapted','ring_counts'),('CellCharter','scVI')]):
    tag=f'{m}_{"both" if m=="CellCharter" else ref}_r{r}_{rep}';ax=axs[i,j]
    if tag not in labels:ax.text(.5,.5,'Not completed',ha='center');continue
    h=spatial(ax,tag,s,bg,reference=ref);handles+=h;shown+=display_key(h,tag,s);ax.set_xlabel(f'{m}, {r} μm: x (μm)')
  handles+=[Line2D([],[],marker='o',ls='',c='#444444',label='Niche unresolved anchor'),Line2D([],[],marker='o',ls='',c='#e8e0d0',label='QC non-anchor background'),Line2D([],[],marker='o',ls='',c='#9fb8d0',label='Annotation Unknown background'),Line2D([],[],marker='o',ls='',c='#f0f0f0',label='Supplied outside analysis QC'),Line2D([],[],marker='s',ls='',markeredgecolor='#aaaaaa',c='white',label='White: no supplied input/background')];fig.legend(handles=shown+handles[-5:],loc='center left',bbox_to_anchor=(1,.5),fontsize=10);export(fig,f'{s}_{ref}_four_methods',handles)
  fig,axs=plt.subplots(2,4,figsize=(20,11));handles=[];shown=[]
  for i,rep in enumerate(['NONE','WEAK']):
   for j,(m,r) in enumerate([('hellinger_kmeans',50),('adapted_leiden',50),('hellinger_kmeans',200),('adapted_leiden',200)]):
    tag=f'{m}_{ref}_r{r}_{rep}'
    if tag not in labels:continue
    h=spatial(axs[i,j],tag,s,bg,reference=ref);handles+=h;shown+=display_key(h,tag,s);axs[i,j].set_xlabel(f'{m} {rep}, {r} μm')
  fig.legend(handles=shown,loc='center left',bbox_to_anchor=(1,.5),fontsize=10);export(fig,f'{s}_{ref}_NONE_WEAK',handles)
 # Existing fixed ROIs: four methods x radii, BLCA display; HPCA full sheets above.
 for _,rr in (ri.iloc[0:0] if '--maps-only' in sys.argv else ri).iterrows():
  region=[rr[k] for k in ['x','y','X','Y']]
  if not np.isfinite(region).all():continue
  fig,axs=plt.subplots(2,4,figsize=(17,9));handles=[]
  for i,r in enumerate([50,200]):
   for j,(m,rep) in enumerate([('hellinger_kmeans','NONE'),('adapted_leiden','NONE'),('MENDER_adapted','ring_counts'),('CellCharter','scVI')]):
    tag=f'{m}_{"both" if m=="CellCharter" else "blca"}_r{r}_{rep}'
    if tag not in labels:continue
    spatial(axs[i,j],tag,s,bg,region);axs[i,j].set_xlabel(f'{m}, {r} μm')
  export(fig,f'{s}_ROI_{rr.roi}')
if '--maps-only' in sys.argv:raise SystemExit(0)
for tag,g in comp.groupby('tag'):
 for ref,h in g.groupby('reference'):
  h=h[h.label>=0].sort_values('n',ascending=False);h=h.head(30);cats=TYPES[ref];fig,ax=plt.subplots(figsize=(7, max(3,len(h)*.25)));im=ax.imshow(h[cats].to_numpy(),aspect='auto',vmin=0,vmax=1,cmap='viridis');ax.set_xticks(range(len(cats)),cats,rotation=50,ha='right');ax.set_yticks(range(len(h)),[f'{k} (n={n})' for k,n in zip(h.label,h.n)]);ax.set_ylabel('Cluster (up to 30 largest; all in tables)');cb=fig.colorbar(im,ax=ax);cb.set_label('Mean original cumulative composition');fig.savefig(F/f'{tag}_{ref}_composition.png',dpi=300,bbox_inches='tight');fig.savefig(F/f'{tag}_{ref}_composition.pdf',bbox_inches='tight');cb.remove();fig.savefig(F/f'{tag}_{ref}_composition_panel_only.pdf',bbox_inches='tight');plt.close(fig);ff,aa=plt.subplots(figsize=(1.5,4));ff.colorbar(im,cax=aa);ff.savefig(F/f'{tag}_{ref}_composition_legend_only.pdf',bbox_inches='tight');plt.close(ff)
cm=pd.read_csv(R/'tables/candidate_metrics.tsv',sep='\t');fig,axs=plt.subplots(1,2,figsize=(12,5));handles=[]
for tag,g in cm.groupby('tag'):
 if 'seed_ARI' not in g or g.seed_ARI.isna().all():continue
 ln=axs[0].plot(g.parameter,g.seed_ARI,'o-',lw=1.1,label=short(tag))[0];axs[1].plot(g.parameter,g.silhouette,'o-',lw=1.1,color=ln.get_color());handles.append(ln)
axs[0].set(xlabel='K or resolution (branch-specific)',ylabel='Seed ARI');axs[1].set(xlabel='K or resolution (branch-specific)',ylabel='Within-representation silhouette');fig.legend(handles=handles,loc='center left',bbox_to_anchor=(1,.5),fontsize=10);export(fig,'candidate_stability',handles)
for name in ['technical_association','coverage']:
 d=pd.read_csv(R/'tables'/f'{name}.tsv',sep='\t');fig,ax=plt.subplots(figsize=(10,7))
 if name=='technical_association':
  cols=['section_NMI','eta2_nFeature_RNA','eta2_low_information_fraction','eta2_input_coverage_component','eta2_reference_conflict_fraction'];im=ax.imshow(d[cols],vmin=0,vmax=1,aspect='auto');ax.set_yticks(range(len(d)),d.tag.map(short),fontsize=10);ax.set_xticks(range(len(cols)),cols,rotation=45,ha='right');fig.colorbar(im,ax=ax)
 else:
  z=d.pivot_table(index='tag',columns='source',values='n',aggfunc='sum').fillna(0);z=z.div(z.sum(1),axis=0);z.index=z.index.map(short);z.plot.barh(stacked=True,ax=ax);ax.set_xlabel('Fraction of all shared anchors');ax.tick_params(axis='y',labelsize=10)
 if name=='coverage':export(fig,name,ax.get_legend_handles_labels()[0])
 else:
  fig.savefig(F/(name+'.png'),dpi=300,bbox_inches='tight');fig.savefig(F/(name+'.pdf'),bbox_inches='tight');fig.delaxes(fig.axes[-1]);fig.savefig(F/(name+'_panel_only.pdf'),bbox_inches='tight');plt.close(fig);ff,aa=plt.subplots(figsize=(1.5,4));ff.colorbar(im,cax=aa);ff.savefig(F/(name+'_legend_only.pdf'),bbox_inches='tight');plt.close(ff)
print('FIGURES_COMPLETE',flush=True)

def heat(mat,rows,cols,name,label,vmax=1):
 fig,ax=plt.subplots(figsize=(max(7,len(cols)*.25),max(3,len(rows)*.24)));im=ax.imshow(mat,aspect='auto',vmin=0,vmax=vmax,cmap='viridis');ax.set_xticks(range(len(cols)),cols,rotation=70,ha='right',fontsize=10);ax.set_yticks(range(len(rows)),rows,fontsize=10);cb=fig.colorbar(im,ax=ax);cb.set_label(label);fig.savefig(F/(name+'.png'),dpi=300,bbox_inches='tight');fig.savefig(F/(name+'.pdf'),bbox_inches='tight');cb.remove();fig.savefig(F/(name+'_panel_only.pdf'),bbox_inches='tight');plt.close(fig);ff,aa=plt.subplots(figsize=(1.7,4));ff.colorbar(im,cax=aa);ff.savefig(F/(name+'_legend_only.pdf'),bbox_inches='tight');plt.close(ff)
q=pd.read_csv(R/'tables/quality_marker_support.tsv',sep='\t');cols=[c for c in q if c.startswith('neighbor_marker_')]
for tag,d in q.groupby('tag'):
 d=d[d.label>=0].sort_values('n',ascending=False).head(30);heat(d[cols].to_numpy(),d.label.astype(str).tolist(),[c.replace('neighbor_marker_','') for c in cols],tag+'_raw_neighbor_markers','Raw neighbor panel detection fraction')
rf=pd.read_csv(R/'tables/MENDER_ring_composition.tsv',sep='\t')
for tag,g in rf.groupby('tag'):
 ref='blca' if '_blca_' in tag else 'hpca';g=g[g.label>=0];labs=comp[(comp.tag==tag)&(comp.label>=0)].sort_values('n',ascending=False).head(30).label.astype(int).tolist();rings=sorted(g.ring.unique());arr=np.array([g[g.label==k].sort_values('ring')[TYPES[ref]].to_numpy().ravel() for k in labs]);columns=[f'{25*j}-{25*(j+1)}:{t}' for j in rings for t in TYPES[ref]];heat(arr,list(map(str,labs)),columns,tag+'_radial_composition','Per-ring known-type composition')
pairs=pd.read_csv(R/'tables/method_comparisons.tsv',sep='\t');d=pairs[(pairs.kind.isin(['Harmony','feature_sensitivity']))&(pairs.scope=='common_native_T')];fig,ax=plt.subplots(figsize=(9,5));ax.barh(np.arange(len(d)),d.ARI);ax.set_yticks(range(len(d)),d.left.map(short),fontsize=10);ax.set_xlabel('Native T partition ARI (not accuracy)');export(fig,'Harmony_Top20_partition_sensitivity')
nd=pd.read_csv(R/'tables/native_graph_diagnostics.tsv',sep='\t');gd=pd.read_csv(R/'tables/identical_input_diagnostics.tsv',sep='\t').set_index('tag');corrected=[int(gd.loc[f'adapted_leiden_{v.reference}_r{v.radius}_NONE','identical_full_inputs_split']) for v in nd.itertuples()];fig,ax=plt.subplots(figsize=(8,4));ax.bar(np.arange(len(nd)),nd.identical_profiles_split,color='#a6611a',label='Native v1 raw P');ax.bar(np.arange(len(nd))+.3,corrected,width=.3,color='#018571',label='Adapted full-input splits (see all branches table)');ax.set_xticks(range(len(nd)),nd.reference+' '+nd.radius.astype(str));ax.set_ylabel('Identical profiles split across labels');ax.legend(frameon=False);export(fig,'native_duplicate_diagnostic',ax.get_legend_handles_labels()[0])

manifest=[]
for (tag,ref),g in comp.groupby(['tag','reference']):
 h=g[g.label>=0].sort_values('n',ascending=False);shown=h.head(30);manifest.append(dict(tag=tag,reference=ref,total_clusters=len(h),heatmap_clusters=len(shown),heatmap_anchor_fraction=float(shown.n.sum()/h.n.sum()),all_labels_retained_in_spatial_maps=True,complete_table='tables/niche_composition.tsv',legend_note='Full categorical key in legend-only PDFs; ROI colors share matching full-sheet legend'))
pd.DataFrame(manifest).to_csv(R/'tables/figure_display_scope.tsv',sep='\t',index=False)
