from run_v3_3 import *
def plots(sec,d,poly,im,ed,un,co,out):
 xy=d[['x','y']].to_numpy();pad=co['median_neighbor_image_px'];bounds=[xy[:,0].min()-pad,xy[:,0].max()+pad,xy[:,1].max()+pad,xy[:,1].min()-pad]
 def setup(ax,title):
  ax.set_title(title,fontsize=9);ax.set_xlim(bounds[:2]);ax.set_ylim(bounds[2:]);ax.set_aspect('equal');ax.axis('off');x=bounds[0]+pad*.4;y=bounds[2]-pad*.3;ax.plot([x,x+pad*5],[y,y],c='black',lw=1);ax.text(x,y-pad*.15,'500 um (nominal)',fontsize=5)
 def cat(ax,v,pal,labels=None):
  ax.add_collection(PolyCollection(poly,facecolors=[pal.get(x,'#999999') for x in v],edgecolors='none'));ax.legend(handles=[Patch(color=c,label=(labels or {}).get(l,l)) for l,c in pal.items() if l in set(v)],fontsize=6,loc='upper center',bbox_to_anchor=(.5,-.01),ncol=1,frameon=False)
 def lines(ax):
  for typ,col,style,label in [('candidate_observed_outer','#151515','-','Observed candidate outer'),('unlocalized_contact','#825995',':','Unlocalized contact'),('sampling_or_tissue_cut_edge','#d0912b','--','Tissue / sampling cut'),('QC_or_missing_edge','#777777',':','QC / missing'),('internal_defect_or_unknown_gap','#777777',':','Defect / unknown')]:
   e=ed[ed.boundary_type.eq(typ)]
   if len(e):ax.add_collection(LineCollection(e[['x1_px','y1_px','x2_px','y2_px']].to_numpy().reshape(-1,2,2),colors=col,linestyles=style,linewidths=.5))
 def cont(ax,v,title,limits,cmap='viridis'):
  setup(ax,title);a=np.asarray(v,float);p=PolyCollection(poly,array=np.ma.masked_invalid(a),edgecolors='none',cmap=cmap);p.set_clim(*limits);p.cmap=p.cmap.copy();p.cmap.set_bad('#bdbdbd');ax.add_collection(p)
  if not np.isfinite(a).any():ax.text(.05,.05,'UNAVAILABLE',transform=ax.transAxes,fontsize=8)
  plt.colorbar(p,ax=ax,fraction=.028,pad=.02);lines(ax)
 def save(fig,name):
  fig.savefig(out/(name+'.png'),dpi=300,bbox_inches='tight');fig.savefig(out/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
 fig,axs=plt.subplots(1,3,figsize=(16,6))
 for ax,c,title in zip(axs,['v3_1_lesion_location','v3_2_lesion_location','lesion_location'],['v3.1 candidate location','v3.2: position/support were conflated','v3.3 candidate location']):
  setup(ax,title);cat(ax,d[c],LP,LL)
 lines(axs[2]);fig.suptitle(sec+' | identical coordinates; red/blue are BOTH computational estimates');fig.subplots_adjust(bottom=.18,wspace=.1);save(fig,'01_location_versions')
 fig,axs=plt.subplots(1,3,figsize=(16,6));setup(axs[0],'v3.3 candidate position (not confirmed core)');cat(axs[0],d.lesion_location,LP,LL);lines(axs[0]);setup(axs[1],'H&E + current candidate / missing edges')
 if im is not None:axs[1].imshow(im)
 else:axs[1].text(.1,.8,'H&E unavailable',transform=axs[1].transAxes)
 axs[1].add_collection(PolyCollection(poly[d.lesion_location.eq('inside')],facecolors='#bd3345',alpha=.12,edgecolors='none'));lines(axs[1]);axs[1].legend(handles=[plt.Line2D([0],[0],color=c,ls=l,label=t) for c,l,t in [('#151515','-','Observed candidate outer'),('#825995',':','Unlocalized contact'),('#d0912b','--','Tissue / sampling cut'),('#777777',':','QC / missing / defect')]],fontsize=6,loc='upper center',bbox_to_anchor=(.5,-.01),frameon=False)
 setup(axs[2],'Position support (separate from red/blue location)');sp={'baseline_candidate_estimate':'#c9bd9f','observed_local_context':'#3b9e90','author_coarse_calibrated':'#5b8dc5','unstable_small_candidate':'#e59d47','reference_edge_conflict_candidate':'#c05b9c','original_unlocalized':'#999999','original_quality_excluded':'#333333'};cat(axs[2],d.position_support,sp);fig.suptitle(sec+' | author information is coarse calibration, not independent truth');fig.subplots_adjust(bottom=.2,wspace=.12);save(fig,'02_location_HE_support')
 fig,axs=plt.subplots(2,3,figsize=(16,11))
 for ax,w in zip(axs[0],[100,200,300]):
  setup(ax,f'Half-width {w} um | same new geometry');v=d[f'distance_band_{w}um'];cat(ax,v,BP);ix=v.eq('Band unresolved');ax.add_collection(PolyCollection(poly[ix],facecolors='none',edgecolors='#555555',linewidths=0,hatch='////'));lines(ax)
 cont(axs[1,0],d.signed_distance_um,'Exact signed distance; red inside / blue outside',(-2000,2000),'RdBu');axs[1,0].text(.02,-.02,'um; clipped at ±2000; gray = unavailable',transform=axs[1,0].transAxes,fontsize=6)
 statuses={'estimable':'#329b78','outer_edge_not_observed':'#e5d7b5','unlocalized':'#999999','path_occluded':'#614d30','outer_sampling_or_tissue_truncation':'#df8a15','unresolved_internal_gap':'#9d75b3','image_supported_cavity_or_defect':'#533264','unresolved_location_or_gap':'#aaaaaa','quality_excluded':'#333333'}
 setup(axs[1,1],'Exact-distance reliability');cat(axs[1,1],d.exact_distance_status,statuses)
 setup(axs[1,2],'Position retained when no outer edge is observed');cat(axs[1,2],d.lesion_location,LP,LL);lines(axs[1,2]);fig.suptitle(sec+' | distances CONDITIONAL ON CANDIDATE ROI; no edge => no margin/core');fig.subplots_adjust(hspace=.45,bottom=.12);save(fig,'03_bands_distance_reliability')
 fig,axs=plt.subplots(2,4,figsize=(19,10));sc=json.loads((O/'color_scales.json').read_text());features=['marker_epithelial']+['gene_'+g for g in GENES]+['malignant_fraction','native_spot_MSQ']
 for ax,f in zip(axs.ravel(),features):
  title=f.replace('gene_','').replace('marker_epithelial','EPI: mean available epithelial genes')
  if f=='native_spot_MSQ':title='CopyKAT spot mean(CNA²); baseline caution'
  if f=='malignant_fraction':title='SpaCET estimated malignant fraction'
  cont(ax,d[f],title,tuple(sc[f]['limits']))
 fig.suptitle(sec+' | log1p(CP10k) markers; common cohort color limits; gray = unavailable\nLimits fixed across sections, above-limit values clipped; EPI is not independent validation',fontsize=11);fig.subplots_adjust(hspace=.15,wspace=.12);save(fig,'04_epithelial_molecular_evidence')
 fig,axs=plt.subplots(1,3,figsize=(16,6));setup(axs[0],'CNV support: original confidence retained');cnpal={'strong':'#8c1d36','moderate':'#d67863','moderate_baseline_caution':'#c9a053','coherent_pattern_native_diploid_conflict':'#966cb4','weak_or_no_detectable_pattern':'#b8d1d8','unavailable':'#999999'};cat(axs[0],d.CNV_support,cnpal);lines(axs[0])
 setup(axs[1],'Buffered author coarse reference');rp={'author_coarse_inside':'#bd3345','author_coarse_outside':'#4389b6','reference_edge_or_unassigned':'#e4e4e4','Others_not_assigned':'#aaaaaa','unavailable':'#cccccc'};cat(axs[1],d.author_coarse_reference,rp)
 setup(axs[2],'Author correspondence / conflict');rp={'no_matched_coarse_reference':'#dddddd','coarse_reference_agrees_with_v3.1':'#3c9c87','coarse_reference_conflicts_with_v3.1_corrected':'#7654ad','author_edge_or_nonbinary_reference_candidate_retained':'#d8c9ac','possible_author_conflict_in_registration_boundary_zone':'#e58c3f'};cat(axs[2],d.reference_status,rp);fig.suptitle(sec+' | both author tumour and normal compartments checked');fig.subplots_adjust(bottom=.25,wspace=.1);save(fig,'05_reference_CNV_support')
 rp=O/'registration'/sec/'registration_check.png'
 if rp.exists():
  # Keep the real registration diagnostic accessible alongside every same-section comparison.
  shutil.copyfile(rp,out/'06_author_registration.png');shutil.copyfile(rp.with_suffix('.pdf'),out/'06_author_registration.pdf')
  meta=json.loads((rp.parent/'registration.json').read_text());whole=np.array(Image.open(PRE/'sources/image1.png'));x,y,w,h=meta['annotation_crop'];xx,yy,ww,hh=meta['HE_crop'];au=whole[y:y+h,x:x+w];he=whole[yy:yy+hh,xx:xx+ww]
  fig,axs=plt.subplots(1,3,figsize=(16,6));setup(axs[0],'v3.1 candidate');cat(axs[0],d.v3_1_lesion_location,LP,LL);setup(axs[1],'v3.3 candidate + separate reference support');cat(axs[1],d.lesion_location,LP,LL);lines(axs[1]);axs[2].imshow(np.concatenate([au,he],axis=1));axs[2].axis('off');axs[2].set_title('Published S1 classes + H&E (own scale)')
  la=fig.add_axes([.68,.12,.29,.10]);la.imshow(Image.open(PRE/'sources/S1_legend.png'));la.axis('off');fig.suptitle(sec+' | '+meta['status']);fig.text(.05,.035,'Buffered coarse references may calibrate BOTH inside and outside. Uncertain image edges/Others are not author barcode labels.',fontsize=8);save(fig,'07_author_comparison')
