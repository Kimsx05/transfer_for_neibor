from fu_common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.colors import Normalize
from matplotlib.patches import Circle,Patch
from matplotlib.lines import Line2D
fm.fontManager.addfont('/home/data/t070721/.local/share/fonts/corefonts-arial/Arial.TTF')
plt.rcParams.update({'font.family':'Arial','font.size':12,'axes.linewidth':1,'lines.linewidth':1.1,'pdf.fonttype':42,'figure.facecolor':'white','savefig.facecolor':'white'})
spec=importlib.util.spec_from_file_location('oldplot',V/'code/10_plot.py');vp=importlib.util.module_from_spec(spec);spec.loader.exec_module(vp);vp.FIG=R/'figures'
def save_roi(fig,name,cbar_spec):
 handles=[Patch(color='#e5e5e5',label='Supplied background (not MP10-low)'),Patch(facecolor='white',edgecolor='gray',label='No supplied input'),Line2D([],[],c='#af799f',marker='x',ls='',label='Unmeasurable score'),Line2D([],[],c=plt.cm.viridis(0),marker='s',ls='',label='Measured zero (color scale 0)'),Line2D([],[],c='#e58e28',marker='s',mfc='none',ls='',label='Center T source bins'),Line2D([],[],c='#3582a5',marker='s',mfc='none',ls='',label='Bin epithelial source'),Line2D([],[],c='#4a9d62',marker='s',mfc='none',ls='',label='Cell epithelial source'),Line2D([],[],c='#c65555',marker='x',ls='',label='Shared source excluded'),Line2D([],[],c='#c0444c',label='50 µm radius'),Line2D([],[],c='#4d4d4d',ls='--',label='200 µm radius'),Line2D([],[],c='black',label='Scale bar = 64 µm')]
 leg=fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,-.075),ncol=4,frameon=False,fontsize=11)
 fig.savefig(R/'figures'/(name+'.png'),dpi=300,bbox_inches='tight');fig.savefig(R/'figures'/(name+'.pdf'),bbox_inches='tight');leg.remove()
 for ax in fig.axes:
  if ax.get_label()=='<colorbar>':ax.set_visible(False)
 fig.savefig(R/'figures'/(name+'_panel_only.pdf'),bbox_inches='tight');plt.close(fig)
 lf=plt.figure(figsize=(11,12));a=lf.add_axes([0,.72,1,.27]);a.axis('off');a.legend(handles=handles,loc='center',ncol=2,frameon=False,fontsize=12)
 for j,(label,cmap,lo,hi)in enumerate(cbar_spec):
  a=lf.add_axes([.3,.66-j*.052,.6,.013]);cb=lf.colorbar(plt.cm.ScalarMappable(norm=Normalize(lo,hi),cmap=cmap),cax=a,orientation='horizontal',extend='max');cb.set_label(label,fontsize=10);cb.ax.tick_params(labelsize=9)
 lf.savefig(R/'figures'/(name+'_legend_only.pdf'),bbox_inches='tight');plt.close(lf)
ROI=pd.read_csv(R/'tables/roi_manifest_display512.tsv',sep='\t');descriptions=[];neigh=[];summary=[]
for s in SECTIONS:
 b=pd.read_parquet(V/s/'bin_evidence.parquet');bc=b.set_index('bin_id');t=pd.read_parquet(V/s/'T32_regions.parquet');ttab=t.set_index('region_id');e=pd.read_parquet(V/s/'Epi64_regions.parquet');em=pd.read_parquet(V/s/'Epi64_members.parquet');tm=pd.read_parquet(V/s/'T32_members.parquet');ev=pd.read_parquet(R/'tables'/f'{s}_evidence.parquet').set_index('region_id');pn=pd.read_parquet(V/s/'region_pairs.parquet');pn=pn[pn.T_grid_um.eq(32)&pn.radius_um.eq(200)];n=pd.read_parquet(V/s/'T32_neighborhoods.parquet');n=n[n.source_rule.eq('exclude_shared')]
 if s=='Cervilla':
  ce=pd.read_parquet(R/'tables/Cervilla_cell_regions.parquet');cs=pd.read_parquet(R/'tables/Cervilla_cell_region_source8.parquet');cp=pd.read_parquet(R/'tables/Cervilla_cell_region_pairs.parquet');cp=cp[cp.radius_um.eq(200)]
 for rr in ROI[ROI.section.eq(s)].itertuples():
  bounds=(rr.display_x0,rr.display_x1,rr.display_y0,rr.display_y1);inside=lambda d:d.x_um.between(bounds[0],bounds[1])&d.y_um.between(bounds[2],bounds[3]);bb=b[inside(b)];tt=t[inside(t)];ee=e[inside(e)];center=ttab.loc[rr.region_id];bi=pn[pn.T_region_id.eq(rr.region_id)];binsrc=bi[~bi.shared_original_bins];escore=e[e.region_id.isin(binsrc.Epi_region_id)];csrc=cp[cp.T_region_id.eq(rr.region_id)]if s=='Cervilla'else pd.DataFrame();cvalid=csrc[~csrc.shared_original_bins]if len(csrc)else pd.DataFrame();cells=ce[ce.region_id.isin(cvalid.cell_region_id)]if len(cvalid)else pd.DataFrame()
  fig,ax=plt.subplots(3,5,figsize=(22,13));axs=ax.ravel();specs=[]
  def base(a):
   a.scatter(bb.x_um,bb.y_um,c='#e5e5e5',s=4,marker='s',lw=0,rasterized=True);a.set_xlim(bounds[:2]);a.set_ylim(bounds[3],bounds[2]);a.set_aspect('equal');a.tick_params(labelsize=8)
  def colored(a,d,values,label,vmax,eligible=None):
   base(a);ok=np.ones(len(d),bool)if eligible is None else np.asarray(eligible);v=np.asarray(values);a.scatter(d.loc[ok,'x_um'],d.loc[ok,'y_um'],c=v[ok],s=12 if len(d)<1000 else 5,marker='s',cmap='viridis',vmin=0,vmax=vmax,lw=0,rasterized=True)
   if (~ok).any():a.scatter(d.loc[~ok,'x_um'],d.loc[~ok,'y_um'],s=10,marker='x',c='#af799f',lw=.7,rasterized=True)
   a.set_xlabel(label+'\nx (µm)',fontsize=10);fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,vmax),cmap='viridis'),ax=a,fraction=.03,extend='max');specs.append((label+' (cap '+str(vmax)+')','viridis',0,vmax))
  im=vp.image_crop(s,bounds)
  if im is not None:axs[0].imshow(im,extent=[bounds[0],bounds[1],bounds[3],bounds[2]]);axs[0].set_xlabel('H&E (hires only)'if s=='UC_SCC'else'H&E (registered)')
  else:base(axs[0]);axs[0].set_xlabel('Expression coverage\nRelative µm; tissue unknown',fontsize=10)
  colored(axs[1],bb,bb[['raw_'+g for g in TL]].sum(1),'Raw T-lineage UMI',10)
  for a,g in zip(axs[2:6],CORE):colored(a,bb,bb['raw_'+g],g+' raw UMI',3)
  colored(axs[6],tt,tt.Treg_main_detected,'Signature detected genes',11)
  for a,field,cap in [(axs[7],'Treg_main',.4),(axs[8],'Tfh',.5),(axs[9],'Activated_T',.3)]:colored(a,tt,tt[field],field+' UCell',cap,tt.measurement_eligible)
  colored(axs[10],ee,ee.MP10,'Bin-enriched MP10',.25,ee.measurement_eligible)
  if s=='Cervilla':
   z=ce[inside(ce)];colored(axs[11],z,z.MP10,'Cell-limited MP10',.25,z.coordinate_supported)
  else:base(axs[11]);axs[11].text(.5,.5,'No existing cell-limited\nMP10 region result',transform=axs[11].transAxes,ha='center',fontsize=11);axs[11].set_xlabel('Not measured (not zero)')
  base(axs[12]);mm=bc.loc[tm.loc[tm.region_id.eq(rr.region_id),'bin_id']];axs[12].scatter(mm.x_um,mm.y_um,s=20,c='#e58e28',marker='s',label='T members')
  eb=bc.loc[em.loc[em.region_id.isin(escore.region_id),'bin_id']];axs[12].scatter(eb.x_um,eb.y_um,s=9,c='#3582a5',marker='s',alpha=.7,rasterized=True);axs[12].scatter(escore.x_um,escore.y_um,s=24,facecolors='none',edgecolors='#124d6e',lw=1.1)
  excluded=e[e.region_id.isin(bi.loc[bi.shared_original_bins,'Epi_region_id'])];axs[12].scatter(excluded.x_um,excluded.y_um,s=24,c='#c65555',marker='x',lw=1.1)
  if s=='Cervilla'and len(cells):
   cb=bc.loc[cs.loc[cs.region_id.isin(cells.region_id),'bin_id'].unique()];axs[12].scatter(cb.x_um,cb.y_um,s=12,facecolors='none',edgecolors='#4a9d62',marker='s',lw=.5,rasterized=True);axs[12].scatter(cells.x_um,cells.y_um,s=20,facecolors='none',edgecolors='#236c38',marker='^',lw=1.1)
  axs[12].set_xlabel('Exact source members / neighbor centers\nShared-bin exclusions marked x',fontsize=10)
  # Overlay exact source members in score panels as noninterpolated support outlines.
  for a in [axs[7],axs[8],axs[9]]:a.scatter(mm.x_um,mm.y_um,s=18,facecolors='none',edgecolors='#e58e28',marker='s',lw=.7)
  ebdisp=bc.loc[em.loc[em.region_id.isin(ee.loc[ee.measurement_eligible,'region_id']),'bin_id']];axs[10].scatter(ebdisp.x_um,ebdisp.y_um,s=9,facecolors='none',edgecolors='#3582a5',marker='s',lw=.5,rasterized=True)
  if s=='Cervilla':
   cbdisp=bc.loc[cs.loc[cs.region_id.isin(ce.loc[inside(ce),'region_id']),'bin_id'].unique()];axs[11].scatter(cbdisp.x_um,cbdisp.y_um,s=9,facecolors='none',edgecolors='#4a9d62',marker='s',lw=.5,rasterized=True)
  colored(axs[13],tt,np.log10(1+tt.nCount),'Regional log10(UMI+1)',4.5)
  textrows=[]
  for rad in [50,200]:
   q=n[n.T_region_id.eq(rr.region_id)&n.radius_um.eq(rad)];bn=binsrc[binsrc.distance_um<=rad+1e-8]if 'distance_um'in binsrc else binsrc
   # Pair distance is persisted; field may be distance_um or distance.
   dc='centroid_distance_um'
   bn=binsrc[binsrc[dc]<=rad+1e-8];cn=cvalid[cvalid.distance_um<=rad+1e-8]if len(cvalid)else cvalid
   z=dict(section=s,roi=rr.roi,T_region_id=rr.region_id,radius_um=rad,bin_score_regions=len(bn),cell_score_regions=len(cn)if s=='Cervilla'else np.nan,bin_mean_depth=float(e.set_index('region_id').loc[bn.Epi_region_id,'nCount'].mean())if len(bn)else np.nan,cell_mean_depth=float(ce.set_index('region_id').loc[cn.cell_region_id,'nCount_RNA'].mean())if len(cn)else np.nan,input_observed_fraction=float(q.observed_circle_fraction.iloc[0])if len(q)else np.nan,bin_shared_excluded=int((bi.shared_original_bins&(bi[dc]<=rad+1e-8)).sum()),cell_shared_excluded=int((csrc.shared_original_bins&(csrc.distance_um<=rad+1e-8)).sum())if len(csrc)else 0);summary.append(z);textrows.append(f'{rad} µm: bin n={len(bn)}, cell n='+ (str(len(cn))if s=='Cervilla'else'NA')+f'\nInput fraction={z["input_observed_fraction"]:.2f}')
  for typ,pair,idc in [('bin',bi,'Epi_region_id'),('cell',csrc,'cell_region_id')]:
   if len(pair):
    for qq in pair.to_dict('records'):neigh.append(dict(section=s,roi=rr.roi,source=typ,**qq))
  axs[14].axis('off');old=ev.loc[rr.region_id,'old_evidence'];new=ev.loc[rr.region_id,'followup_evidence'];txt='\n\n'.join(textrows)+f'\n\nT depth={int(center.nCount)}\nT genes={int(center.T_detected)}\nMain genes={int(center.Treg_main_detected)}\nCore genes={int(center.Treg_core_detected)}\n\nEvidence unchanged'
  if s=='Cervilla'and rr.roi.startswith('candidate'):txt+='\nAdjacent measurements\nof one hotspot'
  axs[14].text(0,1,txt,va='top',fontsize=11)
  for a in axs[:14]:
   a.set_xlim(bounds[:2]);a.set_ylim(bounds[3],bounds[2]);a.set_aspect('equal');a.tick_params(labelsize=8)
   for rad,col,ls in [(50,'#c0444c','-'),(200,'#4d4d4d','--')]:a.add_patch(Circle((center.x_um,center.y_um),rad,fill=False,ec=col,lw=1.1,ls=ls))
   a.plot([bounds[0]+16,bounds[0]+80],[bounds[3]-16]*2,c='black',lw=1.4)
  fig.subplots_adjust(wspace=.4,hspace=.4);save_roi(fig,s+'_'+rr.roi+'_512atlas',cbar_spec=specs)
  descriptions.append(dict(section=s,roi=rr.roi,center_id=rr.region_id,old_evidence=old,new_evidence=new,raw_main_genes=center.Treg_main_detected_genes,raw_core_genes=center.Treg_core_detected_genes,Treg_main=center.Treg_main,Tfh=center.Tfh,Activated_T=center.Activated_T,bin_neighbors200=len(binsrc),cell_neighbors200=len(cvalid)if s=='Cervilla'else np.nan,shared_bin_excluded200=int(bi.shared_original_bins.sum()),comment='unchanged center; full200um circle shown; raw multigene signal does not establish cell identity; no morphology claim without HE'))
  print('ROI',s,rr.roi,flush=True)
put(pd.DataFrame(summary),'ROI_neighborhood_coverage.tsv');put(pd.DataFrame(neigh),'ROI_neighbor_sources.tsv.gz');put(pd.DataFrame(descriptions),'ROI_evidence_notes.tsv')
