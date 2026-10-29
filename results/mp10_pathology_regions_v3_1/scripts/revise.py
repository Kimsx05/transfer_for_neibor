import sys,os,json,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3];OLD=ROOT/'results/mp10_pathology_regions_v3';OUT=ROOT/'results/mp10_pathology_regions_v3_1'
sys.path.insert(0,str(OLD/'scripts'))
from run_v3 import np,pd,M,FULL,OFF,VERT,graph,components,exterior_empty,foreground,nearest,Image,plt,PolyCollection,LineCollection,Patch,cKDTree
CONFIG={'widths_um':[100,200,300],'foreground_min':.65,'foreground_blank_max':.15,'pilot':['GSE238145_B16','Inhousedata_T4413900','GSE319536_BC01','GSE285715_N1'],'ROI_revision':'none_without_demonstrated_artifact','unknown_model':'Any unobserved space can contain an extension; observed location labels fixed; qualified enclosed technical holes inherit surrounding macro-location for geometry only.'}
(OUT/'config.json').write_text(json.dumps(CONFIG,indent=2))
LP={'inside':'#bd3345','outside':'#4389b6','unlocalized':'#999999','quality_excluded':'#333333','no_lesion':'#4389b6'}
BP={'Interior (>w confirmed)':'#a92742','Inner margin':'#f0956b','Outer margin':'#e9bb51','Outside (>w confirmed)':'#4389b6','Band unresolved':'#eeeeee','No supported lesion':'#4389b6','Location unresolved':'#999999','QC excluded':'#333333'}
def segs(e):return e[['x1_um','y1_um','x2_um','y2_um']].to_numpy().reshape(-1,2,2)
def process(sec):
 r=M.loc[sec];src=OLD/'by_dataset'/r.dataset_id/sec;out=OUT/'by_dataset'/r.dataset_id/sec;out.mkdir(parents=True,exist_ok=True)
 if (out/'complete.json').exists() and '--force' not in sys.argv:return
 d=pd.read_parquet(src/'spot_regions_v3.parquet');n=len(d);ed=pd.read_csv(src/'lesion_boundary_edges.tsv.gz',sep='\t');co=json.loads((src/'coordinate_scale.json').read_text());aff=np.array(co['array_to_image_affine']);im=np.asarray(Image.open(r.image_path).convert('RGB'));adj,lookup=graph(d);ext=exterior_empty(lookup);uv=d[['x_um','y_um']].to_numpy();loc=d.lesion_location.to_numpy();tc,tcs=components(adj,np.ones(n,bool));keys=np.c_[3*d.array_col,3*d.array_row].astype(int)[:,None,:]+VERT;poly=np.concatenate([keys[:,:,:1]/6,keys[:,:,1:]/(2*np.sqrt(3)),np.ones((n,6,1))],axis=2)@aff
 # Enumerate every internal absent component, not only missing cells next to the ROI.
 rr=d.array_row;cc=d.array_col;missing={(a,b) for a in range(int(rr.min()),int(rr.max())+1) for b in range(int(cc.min()),int(cc.max())+1) if (a+b)%2==int(rr.iloc[0]+cc.iloc[0])%2 and (a,b) not in lookup and (a,b) not in ext};seen=set();tech={};holes=[]
 for p in sorted(missing):
  if p in seen:continue
  q=[p];seen.add(p)
  for a,b in q:
   for da,db in OFF:
    v=(a+da,b+db)
    if v in missing and v not in seen:seen.add(v);q.append(v)
  nb=sorted({lookup[v] for a,b in q for da,db in OFF if (v:=(a+da,b+db)) in lookup});px=np.array([[b/2,a*np.sqrt(3)/2,1] for a,b in q])@aff;fg=foreground(pd.DataFrame(px,columns=['x','y']),im,co['median_neighbor_image_px']);full=all(v in FULL.get(sec,set()) for v in q);same=bool(nb) and len(set(tc[nb]))==1 and len(set(loc[nb]))==1 and loc[nb[0]] in ['inside','outside'];same=same and (loc[nb[0]]!='inside' or len(set(d.lesion_roi_id.iloc[nb]))==1)
  supported=full and same and np.isfinite(fg).all() and min(fg)>=CONFIG['foreground_min'] and co['affine_residual_max_pitch']<=.15
  typ='supported_internal_technical_missing' if supported else ('image_supported_cavity_or_defect' if np.nanmax(fg)<CONFIG['foreground_blank_max'] else 'unresolved_internal_gap')
  hid=f'{sec}:H{len(holes)+1}'
  for k,v in enumerate(q):
   holes.append({'section_id':sec,'gap_component_id':hid,'array_row':v[0],'array_col':v[1],'x':px[k,0],'y':px[k,1],'HE_foreground':fg[k],'in_original_tissue_grid':v in FULL.get(sec,set()),'gap_type':typ,'same_surrounding_location':same,'molecular_values_imputed':False,'QC_reason':'not_in_v3_retained_set; upstream individual exclusion reason unavailable'})
   if supported:tech[v]=(tc[nb[0]],loc[nb[0]])
 pd.DataFrame(holes,columns=['section_id','gap_component_id','array_row','array_col','x','y','HE_foreground','in_original_tissue_grid','gap_type','same_surrounding_location','molecular_values_imputed','QC_reason']).to_csv(out/'internal_gap_assessment.tsv',sep='\t',index=False)
 # Boundary of ALL observed known tissue, augmented only by qualified geometry-only technical cells.
 # This closes the outside-bound loophole: unknown territory on the outside side also counts.
 augmented={p:(tc[i],loc[i]) for p,i in lookup.items() if loc[i] in ['inside','outside']};augmented.update(tech);unknown=[]
 def shared(p,q):
  a=np.array([3*p[1],3*p[0]])+VERT;b=np.array([3*q[1],3*q[0]])+VERT;k=np.array(sorted(set(map(tuple,a))&set(map(tuple,b))));return np.c_[k[:,0]/6,k[:,1]/(2*np.sqrt(3))]*100
 htype={(h['array_row'],h['array_col']):h['gap_type'] for h in holes}
 for p,(comp,_) in augmented.items():
  for da,db in OFF:
   q=(p[0]+da,p[1]+db)
   if q in augmented:continue
   typ=htype.get(q,'outer_sampling_or_tissue_truncation' if q in ext else 'unresolved_location_or_gap');v=shared(p,q);unknown.append([comp,typ,*v.ravel()])
 un=pd.DataFrame(unknown,columns=['component','unknown_type','x1_um','y1_um','x2_um','y2_um']);un.to_csv(out/'unknown_region_boundary.tsv.gz',sep='\t',index=False)
 olddist=d.signed_distance_um.to_numpy().copy();d['v3_signed_distance_um']=olddist;d['v3_distance_status']=d.distance_status;d['v3_main_region']=d.main_region
 for w in CONFIG['widths_um']:d[f'v3_distance_band_{w}um']=d[f'distance_band_{w}um']
 lo=np.full(n,np.nan);hi=lo.copy();exact=lo.copy();status=np.full(n,'no_valid_outer_boundary',object);reason=np.full(n,'',object)
 groups=[(ix,ed[ed.lesion_roi_id==rid]) for rid,ix in d[d.lesion_location=='inside'].groupby('lesion_roi_id').groups.items()]
 groups += [(np.array(ix)[loc[ix]=='outside'],ed[ed.tissue_component==k]) for k,ix in enumerate(tcs)]
 for ix,e in groups:
  ix=np.array(ix,dtype=int);e=e[e.boundary_type=='biological_lesion_outer']
  if not len(ix) or not len(e):continue
  dd,jj,foot=nearest(uv[ix],segs(e));u=un[un.component==tc[ix[0]]];bad,bj,_=nearest(uv[ix],segs(u));lo[ix]=np.minimum(dd,bad)
  pts=np.array([[p[1]*50,p[0]*np.sqrt(3)*50] for p,(comp,_) in augmented.items() if comp==tc[ix[0]]]);tree=cKDTree(pts)
  for k,i in enumerate(ix):
   path=np.linspace(uv[i],foot[k],max(2,int(np.ceil(dd[k]/25))+1));visible=tree.query(path)[0].max()<=58
   if not visible:status[i]='path_occluded';continue
   hi[i]=dd[k]
   if bad[k]>dd[k]+1e-6:exact[i]=dd[k];status[i]='estimable'
   else:status[i]=u.unknown_type.iloc[bj[k]]
 # No lesion / unresolved position are not assigned numerical boundary bounds.
 for state in ['unlocalized','quality_excluded']:status[loc==state]=state;lo[loc==state]=np.nan;hi[loc==state]=np.nan
 d['signed_distance_um']=np.where(loc=='inside',-exact,exact);d['exact_distance_status']=status;d['distance_status']=status;d['abs_distance_lower_bound_um']=lo;d['abs_distance_upper_bound_um']=hi
 recovered=np.isnan(olddist)&np.isfinite(exact);d['exact_distance_change_reason']=np.where(recovered,'technical_gap_geometry_support_and_no_closer_unknown_front',np.where(np.isfinite(olddist)&np.isnan(exact),'more_conservative_all_unknown_region_model','unchanged_estimability'))
 for w in CONFIG['widths_um']:
  band=np.full(n,'Band unresolved',object);bs=np.full(n,'interval_crosses_threshold_or_insufficient_bounds',object);near=np.isfinite(hi)&(hi<=w);far=np.isfinite(lo)&(lo>w)
  band[near&(loc=='inside')]='Inner margin';band[near&(loc=='outside')]='Outer margin';band[far&(loc=='inside')]='Interior (>w confirmed)';band[far&(loc=='outside')]='Outside (>w confirmed)';bs[near]='upper_bound_confirms_margin';bs[far]='lower_bound_confirms_non_margin';bs[np.isfinite(exact)]='exact_distance'
  if not (loc=='inside').any():band[loc=='outside']='No supported lesion';bs[loc=='outside']='no_supported_lesion_not_a_distance_band'
  band[loc=='unlocalized']='Location unresolved';bs[loc=='unlocalized']='location_unresolved';band[loc=='quality_excluded']='QC excluded';bs[loc=='quality_excluded']='QC_excluded'
  d[f'distance_band_{w}um']=band;d[f'band_status_{w}um']=bs
 d['band_status']=d.band_status_200um;d['distance_band']=d.distance_band_200um;d['main_region']=d.distance_band;d['region_version']='v3.1';d['distance_definition']='Euclidean_to_unchanged_valid_outer_with_all_unknown_region_bounds';d.to_parquet(out/'spot_regions_v3_1.parquet',index=False)
 cols=['dataset_id','section_id','barcode','x','y','lesion_roi_id','lesion_location','signed_distance_um','exact_distance_status','abs_distance_lower_bound_um','abs_distance_upper_bound_um','distance_band_100um','distance_band_200um','distance_band_300um','band_status_100um','band_status_200um','band_status_300um','exact_distance_change_reason','v3_signed_distance_um','v3_distance_status','v3_main_region'];d[cols].to_csv(out/'spot_labels.tsv.gz',sep='\t',index=False)
 for f in ['lesion_rois_um.geojson','lesion_rois_image_px.geojson','lesion_outer_boundaries_um.geojson','lesion_boundary_edges.tsv.gz','coordinate_scale.json']:shutil.copyfile(src/f,out/f)
 # Small ROI evidence, no size-based edits.
 small=[]
 for rid,z in d[d.lesion_location=='inside'].groupby('lesion_roi_id'):
  if len(z)>=10 and not (sec=='GSE238145_B16' and len(z)<d.lesion_roi_id.value_counts().drop('',errors='ignore').max()):continue
  small.append({'section_id':sec,'lesion_roi_id':rid,'n_spots':len(z),'barcodes':';'.join(z.barcode),'x':z.x.mean(),'y':z.y.mean(),'malignant_median':z.malignant_fraction.median(),'CNV_support':';'.join(sorted(set(z.CNV_support))),'domain_components':';'.join(sorted(set(z.domain_component_id))),'n_valid_edges':int(((ed.lesion_roi_id==rid)&ed.boundary_type.eq('biological_lesion_outer')).sum()),'decision':'retain_candidate_pending_review','basis':'no demonstrated generation artifact; size alone insufficient','changed_spots':0,'changed_edges':0})
 pd.DataFrame(small,columns=['section_id','lesion_roi_id','n_spots','barcodes','x','y','malignant_median','CNV_support','domain_components','n_valid_edges','decision','basis','changed_spots','changed_edges']).to_csv(out/'small_ROI_review.tsv',sep='\t',index=False)
 plot(sec,d,poly,im,ed,un,holes,aff,out)
 for w in CONFIG['widths_um']:
  assert not ((d[f'band_status_{w}um']=='lower_bound_confirms_non_margin') & ~(lo>w)).any()
 assert np.array_equal(d.lesion_location,pd.read_parquet(src/'spot_regions_v3.parquet',columns=['lesion_location']).lesion_location)
 info={'section_id':sec,'dataset_id':r.dataset_id,'n_spots':n,'technical_geometry_cells':len(tech),'exact_v3':int(np.isfinite(olddist).sum()),'exact_v3_1':int(np.isfinite(exact).sum()),'exact_restored':int(recovered.sum()),'exact_newly_censored':int((np.isfinite(olddist)&np.isnan(exact)).sum()),'exact_unestimated':int(np.isnan(exact).sum()),'band_200_determined':int(np.isin(d.distance_band_200um,list(BP)[:4]).sum()),'band_200_unresolved':int(d.distance_band_200um.eq('Band unresolved').sum()),'band_only_restored':int((np.isnan(exact)&np.isin(d.distance_band_200um,list(BP)[:4])&~d.v3_distance_band_200um.isin(['Tumor interior','Inner margin','Outer margin','Outside'])).sum()),'ROI_changed_spots':0};(out/'complete.json').write_text(json.dumps(info,indent=2));print(info,flush=True)
def plot(sec,d,poly,im,ed,un,holes,aff,out):
 xy=d[['x','y']].to_numpy();bounds=[xy[:,0].min()-10,xy[:,0].max()+10,xy[:,1].max()+10,xy[:,1].min()-10]
 def setup(ax,title):ax.set_title(title,fontsize=10);ax.set_xlim(bounds[:2]);ax.set_ylim(bounds[2:]);ax.set_aspect('equal');ax.axis('off');length=500/100*np.linalg.norm(aff[0]);x=bounds[0]+10;y=bounds[2]-10;ax.plot([x,x+length],[y,y],c='black',lw=1.4);ax.text(x,y-4,'500 um',fontsize=6)
 def cat(ax,vals,pal):ax.add_collection(PolyCollection(poly,facecolors=[pal.get(x,'#999999') for x in vals],edgecolors='none'));ax.legend(handles=[Patch(color=c,label=l) for l,c in pal.items() if l in set(vals)],fontsize=6,loc='upper center',bbox_to_anchor=(.5,-.01),ncol=2,frameon=False)
 def lines(ax):
  for typ,col,ls,label in [('biological_lesion_outer','#202020','-','Valid outer'),('QC_gap','#8b559c',':','Internal QC / defect'),('trunc','#df8a15','--','External truncation')]:
   e=ed[ed.boundary_type.eq(typ)] if typ!='trunc' else ed[ed.boundary_type.str.contains('truncation|surface_edge')]
   if typ=='QC_gap':e=ed[ed.boundary_type.isin(['QC_gap','internal_tissue_defect'])]
   if len(e):ax.add_collection(LineCollection(e[['x1_px','y1_px','x2_px','y2_px']].to_numpy().reshape(-1,2,2),colors=col,linestyles=ls,linewidths=.65))
 def location(ax):
  values=d.lesion_location.copy()
  if not values.eq('inside').any():values=values.replace({'outside':'no_lesion'})
  cat(ax,values,LP);lines(ax)
 def band(ax,w):
  vals=d[f'distance_band_{w}um'];cat(ax,vals,BP);ix=vals.eq('Band unresolved').to_numpy();ax.add_collection(PolyCollection(poly[ix],facecolors='none',edgecolors='#555555',linewidths=0,hatch='////'));lines(ax)
 def save(fig,name):fig.savefig(out/(name+'.png'),dpi=300,bbox_inches='tight');fig.savefig(out/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
 fig,axs=plt.subplots(1,2,figsize=(12,6));setup(axs[0],'H&E + unchanged ROI / edge types');axs[0].imshow(im);lines(axs[0]);axs[0].add_collection(PolyCollection(poly[d.lesion_location.eq('inside')],facecolors='#bd3345',alpha=.15,edgecolors='none'));setup(axs[1],'Lesion location only (inside is not deep core)');location(axs[1]);axs[0].legend(handles=[plt.Line2D([0],[0],color=c,ls=ls,label=l) for c,ls,l in [('#202020','-','Valid outer'),('#df8a15','--','External truncation'),('#8b559c',':','Internal QC / defect')]],fontsize=7,loc='upper center',bbox_to_anchor=(.5,0));fig.suptitle(sec);save(fig,'01_location_HE')
 fig,axs=plt.subplots(1,3,figsize=(16,6))
 for ax,w in zip(axs,CONFIG['widths_um']):setup(ax,f'w={w} um | hatch: band unresolved');band(ax,w)
 fig.suptitle(sec+' | same ROI and bounds at all widths');save(fig,'02_bands')
 fig,axs=plt.subplots(1,2,figsize=(12,6));setup(axs[0],'Exact signed distance: inside red / outside blue');v=d.signed_distance_um.to_numpy();pc=PolyCollection(poly,array=np.ma.masked_invalid(v),cmap='RdBu',edgecolors='none');pc.set_clim(-2000,2000);pc.cmap=pc.cmap.copy();pc.cmap.set_bad('#dddddd');axs[0].add_collection(pc);fig.colorbar(pc,ax=axs[0],fraction=.03,label='um; clipped at ±2000; gray = unavailable');lines(axs[0]);setup(axs[1],'Exact distance reliability (separate from location)');sp={'estimable':'#329b78','outer_sampling_or_tissue_truncation':'#df8a15','unresolved_internal_gap':'#9d75b3','image_supported_cavity_or_defect':'#533264','unresolved_location_or_gap':'#999999','path_occluded':'#614d30','no_valid_outer_boundary':'#dddddd','unlocalized':'#999999','quality_excluded':'#333333'};cat(axs[1],d.exact_distance_status,sp);fig.suptitle(sec);save(fig,'03_distance_reliability')
 if sec in CONFIG['pilot']:
  fig,axs=plt.subplots(1,3,figsize=(16,6));setup(axs[0],'v3: position and distance status mixed');from run_v3 import PAL;cat(axs[0],d.v3_main_region,PAL);setup(axs[1],'v3.1: position only');location(axs[1]);setup(axs[2],'v3.1: 200 um bands / unresolved hatch');band(axs[2],200);fig.suptitle(sec);save(fig,'04_v3_v3_1_comparison')
  sm=pd.read_csv(out/'small_ROI_review.tsv',sep='\t');fig,ax=plt.subplots(figsize=(8,8));setup(ax,'Small ROI and technical gap review');ax.imshow(im);lines(ax)
  for z in sm.itertuples():ax.plot(z.x,z.y,'o',mfc='none',mec='cyan',ms=8);ax.annotate(z.lesion_roi_id.split(':')[-1],(z.x,z.y),fontsize=7,color='blue')
  hh=[h for h in holes if h['gap_type']=='supported_internal_technical_missing']
  if hh:ax.scatter([h['x'] for h in hh],[h['y'] for h in hh],s=12,marker='+',c='lime',label='Geometry-only technical support');ax.legend(fontsize=7)
  fig.suptitle(sec);save(fig,'05_local_review')
if __name__=='__main__':
 for sec in (CONFIG['pilot'] if '--pilot' in sys.argv else list(M.index)):process(sec)
