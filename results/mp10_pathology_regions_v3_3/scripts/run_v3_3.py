"""Revision of v3.1 candidates, v3.2 geometry and plot templates.
Candidate positions, author calibration and geometric estimability are independent.
"""
import sys,os,json,shutil,traceback,time
from pathlib import Path
R=Path(__file__).resolve().parents[3];O=R/'results/mp10_pathology_regions_v3_3';OLD=R/'results/mp10_pathology_regions_v3_1';PRE=R/'results/mp10_pathology_regions_v3_2';V3=R/'results/mp10_pathology_regions_v3'
sys.path.insert(0,str(V3/'scripts'))
from run_v3 import np,pd,M,FULL,OFF,VERT,graph,components,exterior_empty,foreground,nearest,Image,plt,PolyCollection,LineCollection,Patch,cKDTree,Polygon,LineString,unary_union,feature,writegeo
CFG=json.loads((O/'config.json').read_text())
LP={'inside':'#bd3345','outside':'#4389b6','unlocalized':'#999999','quality_excluded':'#333333','extent_unresolved':'#9b919e'}
LL={'inside':'Candidate inside','outside':'Candidate outside','unlocalized':'Locally unlocalized','quality_excluded':'QC excluded','extent_unresolved':'v3.2 extent unresolved'}
BP={'Interior (>w confirmed)':'#a92742','Inner margin':'#f0956b','Outer margin':'#e9bb51','Outside (>w confirmed)':'#4389b6','Inside; outer edge unobserved':'#bd3345','Outside; no usable lesion edge':'#4389b6','Band unresolved':'#eeeeee','Location unresolved':'#999999','QC excluded':'#333333'}
GENES=['EPCAM','KRT8','KRT19','KRT7','KRT18']
def segments(e):return e[['x1_um','y1_um','x2_um','y2_um']].to_numpy(float).reshape(-1,2,2)

def revise_positions(d,prev,adj,lookup,sec,out):
 n=len(d);base=d.v3_1_lesion_location.to_numpy();loc=base.copy();support=np.full(n,'baseline_candidate_estimate',object);basis=np.full(n,'v3.1 computational candidate retained symmetrically; no pathology confirmation',object);change=np.full(n,'baseline_candidate_retained',object)
 restored=(prev.lesion_location.eq('extent_unresolved').to_numpy() & (base=='outside'));change[restored]='v3.2_restriction_only_restored';basis[restored]='restore v3.1 candidate outside; insufficient biological support is not geometric or positional absence'
 for c in ['domain_local_evidence_score','domain_shared_evidence_score','shared_only_threshold_crossing']:d[c]=prev[c].to_numpy()
 # Review each previous small ROI as an object, against immutable v3.1 neighbours.
 reviews=[];sm=pd.read_csv(PRE/'by_dataset'/M.loc[sec,'dataset_id']/sec/'small_roi_decisions.tsv',sep='\t')
 for row in sm.itertuples():
  ix=np.flatnonzero(d.v3_1_lesion_roi_id.eq(row.old_roi_id));S=set(ix);nb=sorted({j for i in ix for j in adj[i] if j not in S});complete=all(len(adj[i])==6 for i in ix);states=set(base[nb]);weak=row.decision=='weak_borrowing_artifact_corrected';ls=d.domain_local_evidence_score.iloc[ix]
  decision='retain_local_supported_small_candidate';reason='size alone is not a deletion criterion'
  if weak:
   if complete and states=={'outside'} and ls.max()<.42:
    loc[ix]='outside';support[ix]='observed_local_context';change[ix]='weak_borrowed_roi_to_surrounding_outside';reason='weak borrowed ROI; every external neighbour observed and unanimously v3.1 candidate outside';basis[ix]=reason;decision='remove_independent_roi_inherit_outside'
   elif complete and states=={'inside'}:
    support[ix]='observed_local_context';basis[ix]='internal observed lesion compartment; positional inheritance is not malignant identity';decision='retain_inside_internal_compartment';reason=basis[ix[0]]
   else:
    support[ix]='unstable_small_candidate';basis[ix]='weak borrowed ROI retained at v3.1 location: incomplete/mixed external neighbourhood; no justified relabel';decision='retain_unstable_v3.1_candidate';reason=basis[ix[0]]
  reviews.append({'section_id':sec,'old_roi_id':row.old_roi_id,'n_spots':len(ix),'v3_2_weak_roi':weak,'complete_external_neighbourhood':complete,'n_external_observed_neighbours':len(nb),'external_v3_1_locations':';'.join(sorted(states)),'local_score_min':ls.min(),'local_score_max':ls.max(),'decision':decision,'basis':reason,'changed_from_v3_1':int(np.sum(loc[ix]!=base[ix])),'barcodes':';'.join(d.barcode.iloc[ix])})
 # No novel single-spot ROI from reverting local/shared score shrinkage.
 # This baseline rule applies to every v3.2 inside gain, including the specified B16 point.
 gained=prev.lesion_location.eq('inside').to_numpy()&(base=='outside')
 for i in np.flatnonzero(gained):
  if sec=='GSE238145_B16' and d.barcode.iloc[i]=='AGATAACACGGTATAA-1_4':
   assert d.array_row.iloc[i]==89 and d.array_col.iloc[i]==175
   assert len(adj[i])==6 and set(base[adj[i]])=={'outside'}
   basis[i]='correct v3.2 new singleton: local score crossed 0.42 while shared score did not; all six v3.1 neighbours outside; no independent local lesion structure';change[i]='v3.2_unsupported_new_singleton_reverted';support[i]='observed_local_context'
 # Author coarse reference: derived exclusively from reliable H&E registration.
 # Both tumour/stroma and normal stroma/urothelium constrain candidates identically.
 refpath=O/'registration'/sec/'coarse_reference.tsv.gz';ref=np.full(n,'unavailable',object);classes=np.full(n,'unavailable',object);refstatus=np.full(n,'no_matched_coarse_reference',object)
 if refpath.exists():
  z=pd.read_csv(refpath,sep='\t').set_index('barcode').reindex(d.barcode);assert z.section_id.eq(sec).all();ref=z.author_coarse_reference.to_numpy();classes=z.author_panel_class_at_center.to_numpy()
  for label,target in [('author_coarse_inside','inside'),('author_coarse_outside','outside')]:
   ix=np.flatnonzero((ref==label)&~np.isin(base,['quality_excluded','unlocalized']))
   flip=base[ix]!=target;was_local_relabel=change[ix]=='weak_borrowed_roi_to_surrounding_outside';loc[ix]=target;support[ix]='author_coarse_calibrated';refstatus[ix]=np.where(base[ix]==target,'coarse_reference_agrees_with_v3.1','coarse_reference_conflicts_with_v3.1_corrected');change[ix[flip]]='author_coarse_region_corrected';change[ix[(~flip)&was_local_relabel]]='baseline_candidate_retained';basis[ix]='candidate constrained by reliably registered buffered author coarse '+target+' region; not original barcode truth'
  edge=~np.isin(ref,['author_coarse_inside','author_coarse_outside','unavailable']);refstatus[edge]='author_edge_or_nonbinary_reference_candidate_retained'
  rawtarget=np.where(np.isin(classes,['Tumor Cell','Tumor Stroma']),'inside',np.where(np.isin(classes,['Normal Stroma','Urothelium']),'outside',''))
  conflict=edge&(rawtarget!='')&(loc!=rawtarget);refstatus[conflict]='possible_author_conflict_in_registration_boundary_zone';support[conflict]='reference_edge_conflict_candidate'
  # Others and non-binary compartments remain non-binary references, never automatic tumour/normal.
 d['author_coarse_reference']=ref;d['author_panel_class_at_center']=classes;d['reference_status']=refstatus
 for state in ['unlocalized','quality_excluded']:
  sel=base==state;loc[sel]=state;support[sel]='original_'+state;change[sel]='original_'+state+'_preserved';basis[sel]='original v3.1 location/QC limitation retained; no automatic recovery'
 d['lesion_location']=loc;d['location_confidence']=support;d['position_support']=support;d['assignment_basis']=basis;d['change_class']=change
 d['boundary_interpretation']='computed_candidate_contour_not_histological_outer_confirmation'
 if sec=='GSE285715_T3':d['boundary_interpretation']='candidate_interface_may_reflect_intralesional_composition_or_state; macro_outer_unconfirmed'
 d['domain_annotation']=np.array([{'inside':'Candidate lesion-associated domain','outside':'Candidate outside domain','unlocalized':'Locally unlocalized','quality_excluded':'QC excluded'}[x] for x in loc])
 
 for row in reviews:
  ix=np.flatnonzero(d.v3_1_lesion_roi_id.eq(row['old_roi_id']));row['final_locations']=';'.join(sorted(set(loc[ix])));row['final_changed_from_v3_1']=int(np.sum(loc[ix]!=base[ix]));row['author_coarse_constrained_spots']=int(np.isin(ref[ix],['author_coarse_inside','author_coarse_outside']).sum())
 pd.DataFrame(reviews,columns=['section_id','old_roi_id','n_spots','v3_2_weak_roi','complete_external_neighbourhood','n_external_observed_neighbours','external_v3_1_locations','local_score_min','local_score_max','decision','basis','changed_from_v3_1','barcodes','final_locations','final_changed_from_v3_1','author_coarse_constrained_spots']).to_csv(out/'small_roi_decisions.tsv',sep='\t',index=False)
 return d

def new_geometry(d,adj,lookup,tc,tcs,aff,co,im,sec,out):
 n=len(d);loc=d.lesion_location.to_numpy();roi=loc=='inside';rid,regions=components(adj,roi);ids=np.array([f'{sec}:v3_3:L{x+1}' if x>=0 else '' for x in rid],object);d['lesion_roi_id']=ids
 keys=np.c_[3*d.array_col,3*d.array_row].astype(int)[:,None,:]+VERT;vp=np.concatenate([keys[:,:,:1]/6,keys[:,:,1:]/(2*np.sqrt(3)),np.ones((n,6,1))],axis=2);poly=vp@aff;polyu=vp[:,:,:2]*100;ext=exterior_empty(lookup)
 # Positive outside islands fully surrounded by lesion cannot automatically be outer fronts.
 _,outside_comps=components(adj,loc=='outside');outer=np.zeros(n,bool)
 for ix in outside_comps:
  outer[ix]=any((q:=(int(d.array_row.iloc[i])+a,int(d.array_col.iloc[i])+b)) in ext and q not in FULL.get(sec,set()) for i in ix for a,b in OFF)
 def shared_edge(p,q):
  aa=np.array([3*p[1],3*p[0]])+VERT;bb=np.array([3*q[1],3*q[0]])+VERT;kk=np.array(sorted(set(map(tuple,aa))&set(map(tuple,bb))));v=np.c_[kk[:,0]/6,kk[:,1]/(2*np.sqrt(3)),np.ones(2)];return v[:,:2]*100,v@aff
 edges=[]
 for i in np.flatnonzero(roi):
  p=(int(d.array_row.iloc[i]),int(d.array_col.iloc[i]))
  for da,db in OFF:
   q=(p[0]+da,p[1]+db);j=lookup.get(q)
   if j is not None and roi[j]:continue
   if j is not None:
    typ='candidate_observed_outer' if loc[j]=='outside' and outer[j] else ('internal_nonlesion_compartment_contact' if loc[j]=='outside' else 'unlocalized_contact')
   elif q in FULL.get(sec,set()):typ='QC_or_missing_edge'
   elif q not in ext:typ='internal_defect_or_unknown_gap'
   else:typ='sampling_or_tissue_cut_edge'
   su,sp=shared_edge(p,q);edges.append([i,j if j is not None else -1,ids[i],int(tc[i]),typ,*su.ravel(),*sp.ravel()])
 ed=pd.DataFrame(edges,columns=['inside_index','outside_index','lesion_roi_id','tissue_component','boundary_type','x1_um','y1_um','x2_um','y2_um','x1_px','y1_px','x2_px','y2_px']);ed['segment_id']=[f'{sec}:v3_3:E{i+1}' for i in range(len(ed))];ed['boundary_scope']='lesion_outer';ed['evidence_status']='candidate_geometry_not_histological_confirmation'
 # Re-assess technical gaps against NEW location and ROI. No observed rows are added.
 oldgap=OLD/'by_dataset'/M.loc[sec,'dataset_id']/sec/'internal_gap_assessment.tsv';hh=pd.read_csv(oldgap,sep='\t');tech={};hrows=[]
 for hid,z in hh.groupby('gap_component_id'):
  q=list(zip(z.array_row.astype(int),z.array_col.astype(int)));nb=sorted({lookup[v] for a,b in q for da,db in OFF if (v:=(a+da,b+db)) in lookup});same=bool(nb) and len(set(tc[nb]))==1 and len(set(loc[nb]))==1 and loc[nb[0]] in ['inside','outside'];same=same and (loc[nb[0]]!='inside' or len(set(ids[nb]))==1) if nb else False
  supported=bool(same and z.in_original_tissue_grid.all() and z.HE_foreground.notna().all() and z.HE_foreground.min()>=.65 and co['affine_residual_max_pitch']<=.15)
  typ='supported_internal_technical_missing' if supported else ('image_supported_cavity_or_defect' if z.HE_foreground.max()<.15 else 'unresolved_internal_gap')
  for a,b in q:
   hrows.append({'section_id':sec,'gap_component_id':hid,'array_row':a,'array_col':b,'gap_type':typ,'geometry_only':supported,'molecular_values_imputed':False})
   if supported:tech[(a,b)]=(tc[nb[0]],loc[nb[0]])
 gaps=pd.DataFrame(hrows,columns=['section_id','gap_component_id','array_row','array_col','gap_type','geometry_only','molecular_values_imputed']);gaps.to_csv(out/'internal_gap_assessment.tsv',sep='\t',index=False)
 aug={p:(tc[i],loc[i]) for p,i in lookup.items() if loc[i] in ['inside','outside']};aug.update(tech);htype={(z['array_row'],z['array_col']):z['gap_type'] for z in hrows};unknown=[]
 for p,(k,_) in aug.items():
  for a,b in OFF:
   q=(p[0]+a,p[1]+b)
   if q in aug:continue
   typ=htype.get(q,'outer_sampling_or_tissue_truncation' if q in ext else 'unresolved_location_or_gap');su,sp=shared_edge(p,q);unknown.append([k,typ,*su.ravel(),*sp.ravel()])
 un=pd.DataFrame(unknown,columns=['component','unknown_type','x1_um','y1_um','x2_um','y2_um','x1_px','y1_px','x2_px','y2_px'])
 uv=d[['x_um','y_um']].to_numpy();exact=np.full(n,np.nan);lower=exact.copy();upper=exact.copy();candidate=exact.copy();qdist=exact.copy();ref=np.full(n,'',object);refseg=ref.copy();status=np.full(n,'outer_edge_not_observed',object)
 groups=[(ix,ed[ed.lesion_roi_id==ids[ix[0]]]) for ix in regions];groups += [(ix[loc[ix]=='outside'],ed[ed.tissue_component==k]) for k,ix in enumerate(tcs)]
 for ix,e in groups:
  e=e[e.boundary_type=='candidate_observed_outer'];ix=np.asarray(ix,int)
  if not len(ix) or not len(e):continue
  g,jj,foot=nearest(uv[ix],segments(e));u=un[un.component==tc[ix[0]]];q,qj,_=nearest(uv[ix],segments(u));lower[ix]=np.minimum(g,q);candidate[ix]=g;qdist[ix]=q;ref[ix]=e.lesion_roi_id.to_numpy()[jj];refseg[ix]=e.segment_id.to_numpy()[jj]
  tree=cKDTree(np.array([[p[1]*50,p[0]*np.sqrt(3)*50] for p,(k,_) in aug.items() if k==tc[ix[0]]]))
  for k,i in enumerate(ix):
   path=np.linspace(uv[i],foot[k],max(2,int(np.ceil(g[k]/25))+1))
   if tree.query(path)[0].max()>58:status[i]='path_occluded';continue
   upper[i]=g[k]
   if q[k]>g[k]+1e-6:exact[i]=g[k];status[i]='estimable'
   else:status[i]=u.unknown_type.iloc[qj[k]] if len(u) else 'unknown_region'
 status[loc=='unlocalized']='unlocalized';status[loc=='quality_excluded']='quality_excluded'
 d['signed_distance_um']=np.where(roi,-exact,exact);d['exact_distance_status']=status;d['distance_status']=status;d['abs_distance_lower_bound_um']=lower;d['abs_distance_upper_bound_um']=upper;d['reference_lesion_roi_id']=ref;d['reference_boundary_segment_id']=refseg;d['distance_to_observed_outer_candidate_um']=candidate;d['distance_to_unobserved_edge_um']=qdist
 d['distance_definition']='Euclidean_to_v3.3_observed_candidate_outer; v3.1 all-unknown-region bounds';d['distance_zero']='hex_edge_not_spot_center';d['distance_unit']='um';d['boundary_scope']='lesion_outer'
 for w in [100,200,300]:
  band=np.full(n,'Band unresolved',object);bs=np.full(n,'interval_crosses_threshold_or_insufficient_bounds',object);near=np.isfinite(upper)&(upper<=w);far=np.isfinite(lower)&(lower>w)
  band[near&roi]='Inner margin';band[near&(loc=='outside')]='Outer margin';band[far&roi]='Interior (>w confirmed)';band[far&(loc=='outside')]='Outside (>w confirmed)';bs[near]='upper_bound_confirms_margin';bs[far]='lower_bound_confirms_non_margin';bs[np.isfinite(exact)]='exact_distance'
  band[(status=='outer_edge_not_observed')&roi]='Inside; outer edge unobserved';band[(status=='outer_edge_not_observed')&(loc=='outside')]='Outside; no usable lesion edge';bs[status=='outer_edge_not_observed']='no_outer_edge_no_distance_band';band[loc=='unlocalized']='Location unresolved';bs[loc=='unlocalized']='unlocalized';band[loc=='quality_excluded']='QC excluded';bs[loc=='quality_excluded']='quality_excluded'
  d[f'distance_band_{w}um']=band;d[f'band_status_{w}um']=bs
 d['band_status']=d.band_status_200um;d['distance_band']=d.distance_band_200um;d['main_region']=d.distance_band;d['exact_distance_change_reason']='geometry rebuilt from v3.3 candidates; biological confidence does not censor geometry'
 # Replace ambiguous inherited geometry fields instead of leaving old v2/v3 results under current names.
 for name in ['tumor_region_id','boundary_side','boundary_side_two_rings','boundary_type','distance_to_boundary','signed_distance_dNN','boundary_graph_ring','display_region']:
  if name in d:
   d['legacy_'+name]=d[name];d=d.drop(columns=name)
 d['display_region']=d.lesion_location;d['boundary_type']=np.where(np.isfinite(candidate),'candidate_observed_outer','no_candidate_observed_outer');d['distance_to_boundary']=abs(d.signed_distance_um)
 # Assign contextual support AFTER geometry and distance computation, never as an admission filter.
 ed['boundary_support']='molecular_spatial_candidate_only'
 for comp in set(tc):
  ii=np.flatnonzero((tc==comp)&d.author_coarse_reference.eq('author_coarse_inside').to_numpy()&roi)
  oo=np.flatnonzero((tc==comp)&d.author_coarse_reference.eq('author_coarse_outside').to_numpy()&(loc=='outside'))
  ex=ed.index[(ed.tissue_component==comp)&ed.boundary_type.eq('candidate_observed_outer')]
  if not len(ex):continue
  mid=segments(ed.loc[ex]).mean(1);di=cKDTree(uv[ii]).query(mid)[0] if len(ii) else np.full(len(ex),np.inf);do=cKDTree(uv[oo]).query(mid)[0] if len(oo) else np.full(len(ex),np.inf)
  ed.loc[ex[(di<=350)|(do<=350)],'boundary_support']='author_coarse_one_side_context_contour_uncertain'
  ed.loc[ex[(di<=350)&(do<=350)],'boundary_support']='author_coarse_bilateral_neighborhood_not_exact_contour'
 if sec=='GSE285715_T3':ed.loc[ed.boundary_type.eq('candidate_observed_outer'),'boundary_support']='candidate_composition_or_state_interface_possible'
 d['reference_boundary_support']=d.reference_boundary_segment_id.map(ed.set_index('segment_id').boundary_support).fillna('no_candidate_reference')
 ed.to_csv(out/'lesion_boundary_edges.tsv.gz',sep='\t',index=False);un.to_csv(out/'unknown_region_boundary.tsv.gz',sep='\t',index=False)
 for px,polys,unit in [('um',polyu,'local_Visium_nominal_um'),('image_px',poly,'aligned_image_pixels_y_down')]:
  ff=[feature(unary_union([Polygon(polys[i]) for i in ix]),{'lesion_roi_id':ids[ix[0]],'section_id':sec,'n_spots':len(ix),'status':'computational_lesion_extent; open/uncertain edges separate'}) for ix in regions]
  writegeo(out/f'lesion_rois_{px}.geojson',ff,unit)
 writegeo(out/'lesion_outer_boundaries_um.geojson',[feature(LineString(v),{'segment_id':z.segment_id,'lesion_roi_id':z.lesion_roi_id,'boundary_type':z.boundary_type,'boundary_support':z.boundary_support}) for v,z in zip(segments(ed),ed.itertuples())],'local_Visium_nominal_um')
 return d,ed,un,poly,len(tech)

def process(sec):
 r=M.loc[sec];src=OLD/'by_dataset'/r.dataset_id/sec;out=O/'by_dataset'/r.dataset_id/sec;out.mkdir(parents=True,exist_ok=True)
 if (out/'complete.json').exists() and '--force' not in sys.argv:return
 print('START',sec,flush=True);old=pd.read_parquet(src/'spot_regions_v3_1.parquet');prev=pd.read_parquet(PRE/'by_dataset'/r.dataset_id/sec/'spot_regions_v3_2.parquet');d=old.copy();n=len(d);assert d.barcode.equals(prev.barcode)
 for c in ['lesion_location','lesion_roi_id','assignment_basis','signed_distance_um','exact_distance_status','distance_band_100um','distance_band_200um','distance_band_300um','main_region']:
  d['v3_1_'+c]=old[c];d['v3_2_'+c]=prev[c]
 adj,lookup=graph(d);tc,tcs=components(adj,np.ones(n,bool));co=json.loads((src/'coordinate_scale.json').read_text());aff=np.array(co['array_to_image_affine']);im=np.asarray(Image.open(r.image_path).convert('RGB')) if Path(r.image_path).exists() else None
 d=revise_positions(d,prev,adj,lookup,sec,out);d,ed,un,poly,nt=new_geometry(d,adj,lookup,tc,tcs,aff,co,im,sec,out)
 comp=[];d['interface_status']='no_observable_candidate_outer'
 for k,ix in enumerate(tcs):
  ls=set(d.lesion_location.iloc[ix]);e=ed[(ed.tissue_component==k)&ed.boundary_type.eq('candidate_observed_outer')]
  status='candidate_interface_observed' if len(e) else ('all_inside_outer_edge_unobserved' if ls=={'inside'} else ('all_outside_no_lesion_interface' if ls=={'outside'} else 'no_observable_candidate_outer'))
  d.loc[ix,'interface_status']=status;comp.append({'dataset_id':r.dataset_id,'section_id':sec,'tissue_component_id':f'{sec}:T{k+1}','n_spots':len(ix),'status':status,'inside':int(d.lesion_location.iloc[ix].eq('inside').sum()),'outside':int(d.lesion_location.iloc[ix].eq('outside').sum()),'unlocalized':int(d.lesion_location.iloc[ix].eq('unlocalized').sum()),'n_candidate_edges':len(e),'exact_n':int(d.signed_distance_um.iloc[ix].notna().sum())})
 pd.DataFrame(comp).to_csv(out/'tissue_component_summary.tsv',sep='\t',index=False)
 # Original numerical data and identity fields stay unchanged.
 original=[c for c in old if c.startswith(('gene_','marker_','detected_','c2l_','SpaCET_')) or c in ['MP10_original','malignant_fraction','CNV_support','baseline_confidence','patient_id','barcode','x','y','array_row','array_col']]
 pd.testing.assert_frame_equal(d[original],old[original],check_exact=True)
 assert not d[['section_id','barcode']].duplicated().any()
 assert d.loc[old.lesion_location.isin(['unlocalized','quality_excluded']),'lesion_location'].equals(old.loc[old.lesion_location.isin(['unlocalized','quality_excluded']),'lesion_location'])
 d['region_version']='v3.3';d['roi_review_status']='computational_candidate; author_reference_is_coarse_calibration_not_independent_truth';d.to_parquet(out/'spot_regions_v3_3.parquet',index=False)
 cols=['dataset_id','section_id','barcode','x','y','array_row','array_col','x_um','y_um','tissue_component_id','spatial_domain','domain_component_id','v3_1_lesion_location','v3_2_lesion_location','lesion_location','position_support','reference_status','author_coarse_reference','assignment_basis','change_class','lesion_roi_id','boundary_interpretation','interface_status','reference_lesion_roi_id','reference_boundary_segment_id','reference_boundary_support','signed_distance_um','exact_distance_status','distance_to_observed_outer_candidate_um','abs_distance_lower_bound_um','abs_distance_upper_bound_um','distance_band_100um','distance_band_200um','distance_band_300um','band_status_100um','band_status_200um','band_status_300um','domain_local_evidence_score','domain_shared_evidence_score','v3_1_signed_distance_um','v3_2_signed_distance_um']
 d[cols].to_csv(out/'spot_labels.tsv.gz',sep='\t',index=False);shutil.copyfile(src/'coordinate_scale.json',out/'coordinate_scale.json')
 # Save updated mapping on the unchanged secondary contact object; never use it as macro outer.
 ne=pd.read_csv(V3/'by_dataset'/r.dataset_id/sec/'tumor_nest_contact_edges.tsv.gz',sep='\t');ne['inside_lesion_roi']=[d.lesion_location.iloc[int(i)]=='inside' and d.lesion_location.iloc[int(j)]=='inside' for i,j in zip(ne.inside_index,ne.outside_index)];ne.to_csv(out/'tumor_nest_contact_edges.tsv.gz',sep='\t',index=False)
 evidence=[]
 for version,col in [('v3.1','v3_1_lesion_location'),('v3.3','lesion_location')]:
  for state,z in d.groupby(col):
   for f in ['marker_epithelial','malignant_fraction','native_spot_MSQ']+['gene_'+g for g in GENES]:
    evidence.append({'section_id':sec,'version':version,'location':state,'feature':f,'n':len(z),'n_available':int(z[f].notna().sum()),'median':z[f].median(),'q25':z[f].quantile(.25),'q75':z[f].quantile(.75)})
 pd.DataFrame(evidence).to_csv(out/'evidence_distributions.tsv',sep='\t',index=False)
 from plot_v3_3 import plots
 plots(sec,d,poly,im,ed,un,co,out)
 edge=ed[ed.boundary_type.eq('candidate_observed_outer')];inside=d.lesion_location.eq('inside');info={'section_id':sec,'dataset_id':r.dataset_id,'n_spots':n,'inside':int(inside.sum()),'outside':int(d.lesion_location.eq('outside').sum()),'unlocalized':int(d.lesion_location.eq('unlocalized').sum()),'changed_from_v3_1':int(d.lesion_location.ne(old.lesion_location).sum()),'restored_from_v3_2_only':int(d.change_class.eq('v3.2_restriction_only_restored').sum()),'local_small_roi_changed':int(d.change_class.eq('weak_borrowed_roi_to_surrounding_outside').sum()),'author_changed':int(d.change_class.eq('author_coarse_region_corrected').sum()),'author_inside_to_outside':int((old.lesion_location.eq('inside')&d.lesion_location.eq('outside')&d.change_class.eq('author_coarse_region_corrected')).sum()),'author_outside_to_inside':int((old.lesion_location.eq('outside')&inside&d.change_class.eq('author_coarse_region_corrected')).sum()),'candidate_edges':len(edge),'bilateral_author_context_edges':int(edge.boundary_support.eq('author_coarse_bilateral_neighborhood_not_exact_contour').sum()),'exact_estimable':int(d.signed_distance_um.notna().sum()),'band_200_determined':int(d.distance_band_200um.isin(list(BP)[:4]).sum()),'band_only_200':int((d.signed_distance_um.isna()&d.distance_band_200um.isin(list(BP)[:4])).sum()),'rois':int(d.loc[inside,'lesion_roi_id'].nunique()),'technical_geometry_cells':nt,'checks':'passed'}
 (out/'complete.json').write_text(json.dumps(info,indent=2));print('DONE',sec,json.dumps(info),flush=True)
if __name__=='__main__':
 secs=CFG['pilot'] if '--pilot' in sys.argv else ([s for s in sys.argv[1:] if not s.startswith('--')] or list(M.index));failed=[]
 for sec in secs:
  try:process(sec)
  except Exception:
   failed.append(sec);(O/'logs'/f'{sec}.error').write_text(traceback.format_exc());print(traceback.format_exc(),flush=True)
 if failed:raise SystemExit(1)
