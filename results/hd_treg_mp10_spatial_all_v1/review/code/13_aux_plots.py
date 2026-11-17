from common import *
import importlib.util
spec=importlib.util.spec_from_file_location('plot',R/'code/10_plot.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
plt=m.plt
fig,axs=plt.subplots(2,5,figsize=(17,6.5))
for j,s in enumerate(SECTIONS):
 d=pd.read_parquet(R/s/'T32_neighborhoods.parquet')
 for i,rad in enumerate([50,200]):
  v=d[d.radius_um.eq(rad)&d.source_rule.eq('exclude_shared')&d.association_eligible];ax=axs[i,j];ax.scatter(v.Treg_main,v.Epi_MP10,s=5,alpha=.4,c='#4674a1',rasterized=True,lw=0);ax.set_xlabel(s+f' ({rad} µm)\nTreg main UCell',fontsize=10);ax.set_ylabel('Neighbor Epi MP10'if j==0 else '');ax.tick_params(labelsize=9)
fig.subplots_adjust(wspace=.4,hspace=.55);m.save(fig,'continuous_source_disjoint')
d=pd.read_csv(R/'tables/cell_bin_common_effects.tsv',sep='\t');v=d[d.radius_um.eq(200)];fig,ax=plt.subplots(figsize=(6,4));handles=[]
for j,unit in enumerate(['bin_enriched','cell_limited']):
 for k,block in enumerate([400,800]):
  q=v[v.unit.eq(unit)&v.block_um.eq(block)].iloc[0];y=j+(k-.5)*.16;col=['#4674a1','#cc8348'][k];ax.plot([q.low,q.high],[y,y],c=col);h=ax.plot(q.effect,y,'o',c=col,label=f'{block} µm blocks')[0]
  if j==0:handles.append(h)
ax.axvline(0,c='gray',ls='--');ax.set_yticks([0,1],['Bin-enriched MP10','Cell-limited MP10']);ax.set_xlabel('Cervilla 200 µm: same 247 T regions\nAdjusted rank effect (95% block CI)');ax.legend(handles=handles,frameon=False);m.save(fig,'Cervilla_cell_bin_common',legends=handles)
d=pd.read_csv(R/'tables/spottedpy_hotspot_distances.tsv.gz',sep='\t');fig,ax=plt.subplots(figsize=(8,4.5));handles=[]
for k,g in enumerate(['Epi_L50','Epi_H50','Epi_H10']):
 col=['#6b8db4','#d89549','#b95464'][k]
 for j,s in enumerate(SECTIONS):
  z=d[(d.section==s)&d.comparison.eq(g)&d.observable].distance_um
  if len(z):
   x=j+(k-1)*.2;ax.scatter(np.repeat(x,len(z)),z,s=9,c=col,alpha=.5,rasterized=True);ax.plot([x-.07,x+.07],[z.median()]*2,color=col,lw=2)
 handles.append(m.Line2D([],[],marker='o',c=col,ls='',label=g))
ax.set_xticks(range(5),SECTIONS,rotation=20);ax.set_ylabel('Source-disjoint nearest distance (µm)');ax.legend(handles=handles,frameon=False);m.save(fig,'SpottedPy_auxiliary_distance',legends=handles)
d=pd.read_csv(R/'tables/candidate_expression_support.tsv',sep='\t');d=d[d.radius_um.eq(200)&d.Epi_group.eq('all')];fig,axs=plt.subplots(1,2,figsize=(11,4));handles=[]
for k,axis in enumerate(d.axis.unique()):
 v=d[d.axis.eq(axis)].set_index('section').loc[SECTIONS];col=['#4674a1','#cc8348'][k]
 for ax,field in zip(axs,['T_regions_with_neighbor_expression','Treg_like_supported_T_regions']):ax.bar(np.arange(5)+(k-.5)*.35,v[field],width=.35,color=col,label=axis)
 handles.append(m.Patch(color=col,label=axis))
for ax,label in zip(axs,['T-lineage regions with both sides detected','Treg-like subset with both sides detected']):ax.set_xticks(range(5),SECTIONS,rotation=25,ha='right');ax.set_ylabel(label,fontsize=10)
axs[-1].legend(handles=handles,frameon=False,fontsize=9);m.save(fig,'candidate_axis_expression_support',legends=handles)
# Direct ROI single-gene evidence: no normalized or imputed counts used.
d=pd.read_csv(R/'tables/roi_manifest.tsv',sep='\t');genes=['CD3D','TRAC','FOXP3','IL2RA','CTLA4','IKZF2','CCR8','ICOS','CXCL13','PDCD1','TIGIT','CD96'];values=[];labels=[]
for s in SECTIONS:
 t=pd.read_parquet(R/s/'T32_regions.parquet').set_index('region_id')
 for q in d[d.section==s].itertuples():values.append(t.loc[q.region_id,['raw_'+g for g in genes]].to_numpy(float));labels.append(s+' '+q.roi)
values=np.array(values);fig,ax=plt.subplots(figsize=(9,8));im=ax.imshow(np.log1p(values),aspect='auto',cmap='viridis',interpolation='nearest');ax.set_xticks(range(len(genes)),genes,rotation=55,ha='right');ax.set_yticks(range(len(labels)),labels,fontsize=8);fig.colorbar(im,ax=ax,label='Raw regional log1p UMI');m.save(fig,'ROI_raw_marker_evidence',cbar_spec=[('Raw regional log1p UMI','viridis',0,np.log1p(values).max())]);put(pd.DataFrame(values,index=labels,columns=genes).reset_index(names='ROI'),'ROI_raw_marker_plot_data.tsv')
print('AUX PLOTS OK')
