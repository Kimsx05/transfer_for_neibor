from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
fm.fontManager.addfont('/home/data/t070721/.local/share/fonts/corefonts-arial/Arial.TTF')
from matplotlib.lines import Line2D
from matplotlib.patches import Circle,Patch
from matplotlib.colors import Normalize
from PIL import Image
Image.MAX_IMAGE_PIXELS=None
plt.rcParams.update({'font.family':'Arial','font.size':12,'axes.linewidth':1,'lines.linewidth':1.2,'pdf.fonttype':42,'savefig.facecolor':'white','figure.facecolor':'white'})
FIG=R/'figures';FIG.mkdir(exist_ok=True)
def save(fig,name,legends=None,cbar_spec=None):
 fig.savefig(FIG/(name+'.png'),dpi=300,bbox_inches='tight');fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight')
 if legends or cbar_spec:
  for ax in fig.axes:
   leg=ax.get_legend()
   if leg:leg.set_visible(False)
  for ax in list(fig.axes):
   if ax.get_label()=='<colorbar>':ax.set_visible(False)
  fig.savefig(FIG/(name+'_panel_only.pdf'),bbox_inches='tight')
  n=len(cbar_spec or []);lf=plt.figure(figsize=(5,max(1.2,.65*n+(.8 if legends else 0))))
  if legends:
   ax=lf.add_axes([0,.02,1,.9]);ax.axis('off');ax.legend(handles=legends,loc='center',frameon=False,ncol=2)
  for j,(label,cmap,vmin,vmax)in enumerate(cbar_spec or []):
   ax=lf.add_axes([.2,.88-(j+1)*.75/max(n,1),.65,.18/max(n,1)]);cb=lf.colorbar(plt.cm.ScalarMappable(norm=Normalize(vmin,vmax),cmap=cmap),cax=ax,orientation='horizontal');cb.set_label(label)
  lf.savefig(FIG/(name+'_legend_only.pdf'),bbox_inches='tight');plt.close(lf)
 plt.close(fig)
def spatial(ax,d,c=None,cmap='viridis',vmax=None,point=4,alpha=1,bounds=None):
 if c is None:ax.scatter(d.x_um,d.y_um,c='#dddddd',s=point,marker='s',lw=0,rasterized=True,alpha=alpha)
 else:ax.scatter(d.x_um,d.y_um,c=c,s=point,cmap=cmap,vmin=0,vmax=vmax,marker='s',lw=0,rasterized=True,alpha=alpha)
 ax.set_aspect('equal');ax.set_xlabel('x (µm)');ax.set_ylabel('y (µm)')
 if bounds:ax.set_xlim(bounds[0],bounds[1]);ax.set_ylim(bounds[3],bounds[2])
 else:ax.invert_yaxis()
def image_crop(s,bounds):
 if s in SECTIONS[:2]:
  tr=json.loads((A/'segmentation_assignment_v2'/s/'metadata/transform.json').read_text());p=Path(tr['source_image']);mpp=tr['native_mpp'];ox,oy=tr['fullres_crop_xyxy'][:2];scale=1
 elif s=='UC_SCC':
  tr=json.loads((A/'hd_bin8_banksy_extension_v1/UC_SCC/transform.json').read_text());p=H/'data/GSE299573/GSM9040999/square_008um/spatial/tissue_hires_image.png';mpp=tr['microns_per_pixel'];scale=tr['scalefactors']['tissue_hires_scalef'];ox=oy=0
 else:return None
 box=[int((bounds[0]/mpp-ox)*scale),int((bounds[2]/mpp-oy)*scale),int((bounds[1]/mpp-ox)*scale),int((bounds[3]/mpp-oy)*scale)]
 with Image.open(p)as im:return im.crop(box).convert('RGB')
def forest():
 d=pd.read_csv(R/'tables/association_primary.tsv',sep='\t');fig,ax=plt.subplots(figsize=(7,5.5));colors={50:'#3577a8',200:'#ce7544'}
 for i,s in enumerate(SECTIONS):
  for r,off in [(50,-.12),(200,.12)]:
   v=d[(d.section==s)&(d.radius_um==r)].iloc[0];y=i+off
   if np.isfinite(v.effect):
    ax.plot(v.effect,y,'o',color=colors[r],markerfacecolor=colors[r] if np.isfinite(v.low) else 'white');
    if np.isfinite(v.low):ax.plot([v.low,v.high],[y,y],color=colors[r])
   else:ax.text(0,y,'NE',ha='center',va='center',fontsize=9)
 ax.axvline(0,color='gray',ls='--');ax.set_yticks(range(5),[s+(' (mask unknown)'if s.startswith('MPBC')else '')for s in SECTIONS]);ax.invert_yaxis();ax.set_xlabel('Adjusted standardized rank effect (95% block CI)');handles=[Line2D([],[],marker='o',color=c,label=f'{r} µm')for r,c in colors.items()];handles.append(Line2D([],[],marker='o',color='gray',markerfacecolor='white',ls='',label='Descriptive only (no CI/P)'));ax.legend(handles=handles,frameon=False,loc='upper left',bbox_to_anchor=(1,1),fontsize=10);save(fig,'primary_forest',legends=handles)
 # Main vs common-source and unadjusted effect, no selective significance display.
 z=pd.read_csv(R/'tables/association_all_results.tsv',sep='\t');families=['main_unadjusted','main','shared_source_inclusion','larger_blocks','competing_state_adjustment'];fig,axs=plt.subplots(1,2,figsize=(12,5),sharey=True);hh=[]
 for k,f in enumerate(families):
  col=plt.cm.tab10(k);hh.append(Line2D([],[],marker='o',color=col,label=f.replace('_',' ')))
  for ax,rad in zip(axs,[50,200]):
   for i,s in enumerate(SECTIONS):
    v=z[(z.family==f)&(z.section==s)&(z.radius_um==rad)&z.T_grid_um.eq(32)]
    if len(v) and np.isfinite(v.effect.iloc[0]):ax.plot(v.effect.iloc[0],i+(k-2)*.12,'o',color=col,ms=4)
 for ax,rad in zip(axs,[50,200]):ax.axvline(0,c='gray',ls='--');ax.set_xlabel(f'{rad} µm rank effect');ax.set_yticks(range(5),SECTIONS);ax.invert_yaxis()
 axs[1].legend(handles=hh,frameon=False,fontsize=8,loc='upper left',bbox_to_anchor=(1,1));save(fig,'source_and_quality_sensitivity',legends=hh)
def validation():
 d=pd.read_csv(R/'tables/signature_validation.tsv',sep='\t');v=d[d.signature.eq('Treg_main')&d.comparison.isin(['all_nonTreg','T_competitors','Tfh/CXCL13-like'])];fig,ax=plt.subplots(figsize=(6,4));hh=[]
 for i,g in enumerate(['all_nonTreg','T_competitors','Tfh/CXCL13-like']):
  z=v[v.comparison==g].set_index('depth');vals=[z.loc[str(dep),'AUROC']for dep in [500,2000,5000,'full']];line=ax.plot(range(4),vals,'o-',label=g)[0];hh.append(line)
 ax.axhline(.5,c='gray',ls='--');ax.set_ylim(.45,1);ax.set_xticks(range(4),['500','2000','5000','Original']);ax.set_xlabel('Expected UMI thinning (held-out samples)');ax.set_ylabel('Treg discrimination AUROC');ax.legend(frameon=False,fontsize=9);save(fig,'reference_discrimination',legends=hh)
 d=pd.read_csv(R/'tables/reference_mixture_validation.tsv',sep='\t');fig,axs=plt.subplots(1,3,figsize=(12,3.6),sharey=True);hh=[]
 for ax,g in zip(axs,d.competitor.unique()):
  for dep,col in zip([500,2000,5000],['#5b8ab2','#e09f4c','#7aa873']):
   z=d[(d.competitor==g)&(d.depth==dep)].sort_values('Treg_RNA_fraction');ln=ax.plot(z.Treg_RNA_fraction,z.score_median,'o-',color=col,label=str(dep)+' UMI')[0];ax.fill_between(z.Treg_RNA_fraction,z.score_q10,z.score_q90,color=col,alpha=.15)
   if len(hh)<3:hh.append(ln)
  ax.set_xlabel(g+'\nTreg RNA mixture fraction',fontsize=10)
 axs[0].set_ylabel('Treg main UCell (median, 10–90%)');axs[-1].legend(handles=hh,frameon=False,fontsize=9);save(fig,'reference_negative_mixtures',legends=hh)
 d=pd.read_csv(R/'tables/targeted_T_depth_validation.tsv',sep='\t');fig,ax=plt.subplots(figsize=(6,4));hh=[]
 for j,k in enumerate(['T32','T16']):
  v=d[d.branch==k].set_index('section').loc[SECTIONS];h=ax.bar(np.arange(5)+(j-.5)*.35,v.Treg_half_rho,width=.35,label=k);hh.append(Patch(color=h[0].get_facecolor(),label=k))
 ax.set_xticks(range(5),SECTIONS,rotation=25,ha='right');ax.set_ylabel('Half-depth score Spearman correlation');ax.set_ylim(0,1);ax.legend(frameon=False);save(fig,'Treg_depth_stability',legends=hh)
def whole_and_roi():
 roi=pd.read_csv(R/'tables/roi_manifest.tsv',sep='\t');cfg=yaml.safe_load((R/'frozen_config.yaml').read_text())
 for s in SECTIONS:
  b=pd.read_parquet(R/s/'bin_evidence.parquet');t=pd.read_parquet(R/s/'T32_regions.parquet');e=pd.read_parquet(R/s/'Epi64_regions.parquet');bounds=(b.x_um.min()-8,b.x_um.max()+8,b.y_um.min()-8,b.y_um.max()+8);limitT=max(.1,float(t.loc[t.measurement_eligible,'Treg_main'].quantile(.99)));limitE=max(.02,float(e.loc[e.measurement_eligible,'MP10'].quantile(.99)))
  fig,axs=plt.subplots(1,3,figsize=(15,5));spec=[]
  fields=[('T_lineage',t,t.measurement_eligible,'T-lineage UCell'),('Treg_main',t,t.measurement_eligible,'Treg main UCell'),('MP10',e,e.measurement_eligible,'Epi-enriched MP10')]
  for ax,(field,d,good,label)in zip(axs,fields):
   spatial(ax,b.iloc[::4],point=.15,bounds=bounds);v=d[good];lim=limitE if field=='MP10'else limitT if field=='Treg_main'else max(.1,float(v[field].quantile(.99)));spatial(ax,v,c=v[field],vmax=lim,point=3,bounds=bounds);ax.set_xlabel(label+'\nx (µm)');fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,lim),cmap='viridis'),ax=ax,fraction=.03);spec.append((label,'viridis',0,lim))
  save(fig,s+'_whole_scores',cbar_spec=spec)
  # Coverage and unscorable regions are distinct from genuine score zero.
  fig,axs=plt.subplots(1,2,figsize=(10,5));handles=[Patch(color='#dddddd',label='Supplied input'),Patch(color='#4674a1',label='T measurement eligible'),Patch(color='#dd9448',label='Epi MP10 eligible'),Line2D([],[],marker='x',ls='',color='#a44b69',label='Unscorable region')]
  for ax,d,col in zip(axs,[t,e],['#4674a1','#dd9448']):
   spatial(ax,b.iloc[::4],point=.2,bounds=bounds);v=d[d.measurement_eligible];ax.scatter(v.x_um,v.y_um,s=2,c=col,rasterized=True);v=d[~d.measurement_eligible];ax.scatter(v.x_um,v.y_um,s=2,c='#a44b69',marker='x',rasterized=True);ax.set_xlim(bounds[:2]);ax.set_ylim(bounds[3],bounds[2])
  axs[-1].legend(handles=handles,frameon=False,fontsize=8,loc='upper left',bbox_to_anchor=(1,1));save(fig,s+'_measurement_coverage',legends=handles)
  tm=pd.read_parquet(R/s/'T32_members.parquet');em=pd.read_parquet(R/s/'Epi64_members.parquet');bc=b.set_index('bin_id')
  for rr in roi[roi.section==s].itertuples():
   bounds=(rr.x0,rr.x1,rr.y0,rr.y1);bb=b[b.x_um.between(rr.x0,rr.x1)&b.y_um.between(rr.y0,rr.y1)];tt=t[t.x_um.between(rr.x0,rr.x1)&t.y_um.between(rr.y0,rr.y1)];ee=e[e.x_um.between(rr.x0,rr.x1)&e.y_um.between(rr.y0,rr.y1)];fig,axs=plt.subplots(2,4,figsize=(15,8));axs=axs.ravel();spec=[];im=image_crop(s,bounds)
   if im is not None:axs[0].imshow(im,extent=[rr.x0,rr.x1,rr.y1,rr.y0]);axs[0].set_xlabel('Registered H&E\nx (µm)')
   else:spatial(axs[0],bb,point=10,bounds=bounds);axs[0].set_xlabel('Supplied grid (tissue unknown)\nx (µm)')
   rawT=bb[['raw_'+g for g in PANELS['T']]].sum(axis=1);rawcore=bb[['raw_'+g for g in ['FOXP3','IL2RA','CTLA4','IKZF2']]].sum(axis=1)
   items=[(bb,np.log1p(rawT),'Raw T-lineage log1p UMI',max(1,np.log1p(rawT).max())),(tt[tt.measurement_eligible],tt.loc[tt.measurement_eligible,'Treg_main'],'Treg main UCell',limitT),(tt[tt.measurement_eligible],tt.loc[tt.measurement_eligible,'Treg_literature'],'Literature UCell',max(.1,float(t.Treg_literature.quantile(.99)))),(bb,np.log1p(rawcore),'Raw Treg core log1p UMI',max(1,np.log1p(rawcore).max())),(tt[tt.measurement_eligible],tt.loc[tt.measurement_eligible,'Tfh'],'Tfh UCell',max(.1,float(t.Tfh.quantile(.99)))),(ee[ee.measurement_eligible],ee.loc[ee.measurement_eligible,'MP10'],'Epi-enriched MP10',limitE),(tt,np.log10(1+tt.nCount),'Region log10(UMI+1)',max(2,float(np.log10(1+tt.nCount).max())))]
   for ax,(d,col,label,lim)in zip(axs[1:],items):
    spatial(ax,bb,point=8,bounds=bounds);spatial(ax,d,c=col,vmax=lim,point=28 if len(d)<100 else 10,bounds=bounds);ax.set_xlabel(label+'\nx (µm)',fontsize=10);fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,lim),cmap='viridis'),ax=ax,fraction=.04);spec.append((label,'viridis',0,lim))
   # Exact contributing bin centers; no interpolated region or convex-hull support.
   members=bc.loc[tm.loc[tm.region_id.eq(rr.region_id),'bin_id']]
   for ax in [axs[2],axs[3],axs[5],axs[7]]:ax.scatter(members.x_um,members.y_um,s=18,marker='s',facecolors='none',edgecolors='#e68b30',linewidths=.6)
   eid=ee.loc[ee.measurement_eligible,'region_id'];ep=bc.loc[em.loc[em.region_id.isin(eid),'bin_id']]
   axs[6].scatter(ep.x_um,ep.y_um,s=12,marker='s',facecolors='none',edgecolors='#21a7a0',linewidths=.5)
   center=t[t.region_id==rr.region_id].iloc[0]
   for ax in axs:
    ax.set_xlim(rr.x0,rr.x1);ax.set_ylim(rr.y1,rr.y0);ax.set_aspect('equal');ax.tick_params(labelsize=8);ax.add_patch(Circle((center.x_um,center.y_um),50,fill=False,color='#c44e52',lw=1));ax.plot([rr.x0+12,rr.x0+76],[rr.y1-12]*2,color='black',lw=2)
   fig.subplots_adjust(wspace=.45,hspace=.3);save(fig,s+'_'+rr.roi+'_atlas',cbar_spec=spec)
  print('PLOTTED',s,flush=True)
if __name__=='__main__':
 validation();forest();whole_and_roi()
