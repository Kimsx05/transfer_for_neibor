"""Macro-location revision. No upstream fitting and no downstream outcome selection."""
import os,sys,json,shutil,traceback,time
from pathlib import Path
R=Path(__file__).resolve().parents[3];O=R/'results/mp10_pathology_regions_v3_2';OLD=R/'results/mp10_pathology_regions_v3_1';V3=R/'results/mp10_pathology_regions_v3'
sys.path.insert(0,str(V3/'scripts'))
from run_v3 import np,pd,M,FULL,OFF,VERT,graph,components,exterior_empty,foreground,nearest,Image,plt,PolyCollection,LineCollection,Patch,cKDTree,Polygon,LineString,unary_union,feature,writegeo
CFG=json.loads((O/'config.json').read_text())
LP={'inside':'#bd3345','outside':'#4389b6','extent_unresolved':'#9b919e','quality_excluded':'#333333'}
LL={'inside':'Inside lesion ROI (computational)','outside':'Outside: non-lesion tissue support','extent_unresolved':'Extent unresolved (not outside)','quality_excluded':'QC excluded'}
BP={'Interior (>w confirmed)':'#a92742','Inner margin':'#f0956b','Outer margin':'#e9bb51','Outside (>w confirmed)':'#4389b6','Inside; outer edge unobserved':'#bd3345','Outside; no usable lesion edge':'#4389b6','Band unresolved':'#eeeeee','Extent unresolved':'#9b919e','QC excluded':'#333333'}
CNMAP={'strong':1.,'moderate':.55,'moderate_baseline_caution':.35,'coherent_pattern_native_diploid_conflict':.3,'weak_or_no_detectable_pattern':0.,'unavailable':0.}
GENES=['EPCAM','KRT8','KRT19','KRT7','KRT18']
def segments(e):return e[['x1_um','y1_um','x2_um','y2_um']].to_numpy(float).reshape(-1,2,2)
def local_annotation(d,adj,lookup,tc,tcs,sec):
 n=len(d);oldloc=d.lesion_location.to_numpy();valid=d.algorithm_qc.eq('retained').to_numpy();dl,dc=components(adj,valid,d.spatial_domain.astype(str).to_numpy());ext=exterior_empty(lookup)
 cn=d.CNV_support.map(CNMAP).fillna(0).to_numpy();contract=(d.c2l_SMC.ge(.5)&d.detected_smooth_muscle.ge(3)&d.marker_epithelial.lt(.8)).to_numpy();score=.65*d.malignant_fraction.fillna(0).to_numpy()+.20*cn+.15*np.clip(d.marker_epithelial.fillna(0).to_numpy()/2.5,0,1)-.45*contract
 loc=np.full(n,'extent_unresolved',object);basis=np.full(n,'low malignant support does not establish macro outside',object);local=np.zeros(n);shared=np.zeros(n);cross=np.zeros(n,bool);domrows=[];normal=np.zeros(n,bool);seed=np.zeros(n,bool)
 for k,ix in enumerate(dc):
  z=d.iloc[ix];gix=np.flatnonzero((tc==tc[ix[0]])&d.spatial_domain.eq(z.spatial_domain.iloc[0]).to_numpy());ls=float(score[ix].mean());ss=.25*ls+.75*score[gix].mean() if len(ix)<5 else ls;local[ix]=ls;shared[ix]=ss;cross[ix]=(len(ix)<5 and ls<.42<=ss)
  # Local molecular support describes lesion-containing tissue, never a calibrated probability.
  epithelial=z.marker_epithelial.mean()>=.8
  cna=z.CNV_support.isin(['strong','moderate']).mean()>=.5
  malignant=z.malignant_fraction.mean()>=.5
  seqonly=z.SpaCET_signature_source.astype(str).str.contains('seq_depth',case=False).any()
  tumour=ls>=.42 and not seqonly
  # Positive non-lesion lineage evidence, not the complement of a tumour score.
  muscle=z.c2l_SMC.mean()>=.5 and z.detected_smooth_muscle.mean()>=3 and z.marker_epithelial.mean()<.8 and z.malignant_fraction.mean()<.2 and z.CNV_support.eq('strong').mean()<.1
  if tumour:seed[ix]=True;loc[ix]='inside';basis[ix]='local connected-domain composite tumour evidence; candidate inside, not calibrated pathology probability'
  elif muscle:normal[ix]=True;loc[ix]='outside';basis[ix]='organized smooth-muscle program + SMC composition + low malignant support; computational non-lesion candidate'
  domrows.append(dict(section_id=sec,dataset_id=d.dataset_id.iloc[0],domain_component_id=f'{sec}:D{k+1}',spatial_domain=z.spatial_domain.iloc[0],tissue_component=int(tc[ix[0]]),n_spots=len(ix),local_score=ls,shared_score=ss,shared_only_threshold_crossing=bool(cross[ix[0]]),local_tumour_support=tumour,nonlesion_muscle_support=muscle,malignant_mean=z.malignant_fraction.mean(),epi_mean=z.marker_epithelial.mean(),CNV_strong_fraction=z.CNV_support.eq('strong').mean(),SMC_mean=z.c2l_SMC.mean()))
 # Preserve established adjacent-sample adjudication; the raw malignant estimates stay untouched.
 if sec in ['GSE285715_N1','GSE285715_N2','GSE285715_N3']:
  loc[:]='extent_unresolved';seed[:]=False;normal[:]=False
  keep=(oldloc=='outside');loc[keep]='outside';normal[keep]=True;basis[keep]='retained v3.1 adjacent-source urothelial/muscle context adjudication; molecular conflict retained separately';basis[~keep]='retained adjacent-source malignant/context conflict'
 else:
  # Preserve a previously connected lesion body and its observed internal compartments
  # when the SAME observed ROI contains local (not borrowed) support. This is
  # positional inheritance along observed adjacency, not malignant identity transfer.
  _,old_rois=components(adj,valid&(oldloc=='inside'))
  for ix in old_rois:
   if seed[ix].any():
    inherit=ix[~seed[ix]];loc[inherit]='inside';basis[inherit]='retained observed connected lesion body/internal compartment with local-supported anchor in same ROI; composition remains independent'
  # Existing enclosed OBSERVED compartments are admitted only when all observed contact
  # is with local-supported inside and no boundary faces any missing/unresolved cell.
  for _ in range(2):
   _,ucs=components(adj,valid&(loc=='extent_unresolved'))
   add=[]
   for ix in ucs:
    nb={j for i in ix for j in adj[i] if loc[j]!='extent_unresolved'}
    closed=all(len(adj[i])==6 for i in ix)
    if closed and nb and all(loc[j]=='inside' for j in nb) and d.nFeature.iloc[ix].median()>=500 and d.HE_foreground_fraction.iloc[ix].median()>=.35:
     add.extend(ix);basis[ix]='observed viable internal compartment enclosed by locally supported lesion; no missing-cell route assumed'
   loc[add]='inside'
  # Same-sample author S1 supports Tumor Stroma, but gives no per-barcode truth.
  # A limited contextual extension is permitted only for an interleaved domain
  # contacting >=2 distinct locally supported domains, with >=60% of observed
  # cross-domain contact edges to those domains and nonzero malignant/EPI context.
  # This is a computational calibration, not image-derived annotation.
  if sec in CFG['author_tumour_stroma_context']:
   add=[]
   for k,ix in enumerate(dc):
    if loc[ix[0]]!='extent_unresolved' or cross[ix[0]]:continue
    nb=[j for i in ix for j in adj[i] if dl[j]!=k];nz=d.iloc[ix]
    frac=np.mean(seed[nb]) if nb else 0;nd=len(set(dl[j] for j in nb if seed[j]))
    if nd>=2 and frac>=.6 and nz.malignant_fraction.mean()>=.2 and nz.marker_epithelial.mean()>=.5 and nz.c2l_SMC.mean()<.5:
     add.extend(ix);basis[ix]='author S1 Tumor Stroma context + interleaved local domain contacts; calibrated candidate, not barcode truth'
   loc[add]='inside'
  # A weak small component can inherit a CLOSED observed location, not remote-domain identity.
  # No majority filtering, and missing neighbors prevent this operation.
  for k,ix in enumerate(dc):
   if not cross[ix[0]] or len(ix)>=5:continue
   nb={j for i in ix for j in adj[i] if dl[j]!=k}
   if all(len(adj[i])==6 for i in ix) and nb and len({loc[j] for j in nb})==1 and next(iter({loc[j] for j in nb})) in ['inside','outside']:
    loc[ix]=loc[next(iter(nb))];basis[ix]='weak shared-only tiny domain inherits unanimous fully observed immediate spatial context; no malignant upgrade'
   elif not seed[ix[0]] and not np.all(loc[ix]=='inside'):
    loc[ix]='extent_unresolved';basis[ix]='remote same-domain borrowing alone cannot create an independent lesion; local extent unresolved'
 loc[~valid]='quality_excluded';basis[~valid]='upstream QC exclusion retained'
 d['domain_local_evidence_score']=local;d['domain_shared_evidence_score']=shared;d['shared_only_threshold_crossing']=cross;d['local_lesion_support']=seed;d['positive_nonlesion_support']=normal;d['lesion_location']=loc;d['assignment_basis']=basis
 d['location_confidence']=np.where(loc=='extent_unresolved','unresolved_macro_extent','computational_context_not_pathology_confirmed')
 d['domain_component_id']=[f'{sec}:D{x+1}' if x>=0 else '' for x in dl]
 d['domain_annotation']=[{'inside':'Lesion-associated tissue domain','outside':'Non-lesion-supported tissue domain','extent_unresolved':'Macro extent unresolved','quality_excluded':'QC excluded'}[x] for x in loc]
 dt=pd.DataFrame(domrows);dt['v3_2_location']=dt.domain_component_id.map(d.groupby('domain_component_id').lesion_location.first());dt['assignment_basis']=dt.domain_component_id.map(d.groupby('domain_component_id').assignment_basis.first())
 # Check reconstruction of the previous shared score, independent of the new location decisions.
 assert np.allclose(shared[valid],d.domain_evidence_score.to_numpy()[valid],atol=1e-10),sec
 return dt

def new_geometry(d,adj,lookup,tc,tcs,aff,co,im,sec,out):
 n=len(d);loc=d.lesion_location.to_numpy();roi=loc=='inside';rid,regions=components(adj,roi);ids=np.array([f'{sec}:v3_2:L{x+1}' if x>=0 else '' for x in rid],object);d['lesion_roi_id']=ids
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
    typ='supported_macro_outer' if loc[j]=='outside' and outer[j] else ('internal_nonlesion_compartment_contact' if loc[j]=='outside' else 'extent_uncertain_contact')
   elif q in FULL.get(sec,set()):typ='QC_or_missing_edge'
   elif q not in ext:typ='internal_defect_or_unknown_gap'
   else:typ='sampling_or_tissue_cut_edge'
   su,sp=shared_edge(p,q);edges.append([i,j if j is not None else -1,ids[i],int(tc[i]),typ,*su.ravel(),*sp.ravel()])
 ed=pd.DataFrame(edges,columns=['inside_index','outside_index','lesion_roi_id','tissue_component','boundary_type','x1_um','y1_um','x2_um','y2_um','x1_px','y1_px','x2_px','y2_px']);ed['segment_id']=[f'{sec}:v3_2:E{i+1}' for i in range(len(ed))];ed['boundary_scope']='lesion_outer';ed['evidence_status']='computational_supported_not_histologically_confirmed'
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
  e=e[e.boundary_type=='supported_macro_outer'];ix=np.asarray(ix,int)
  if not len(ix) or not len(e):continue
  g,jj,foot=nearest(uv[ix],segments(e));u=un[un.component==tc[ix[0]]];q,qj,_=nearest(uv[ix],segments(u));lower[ix]=np.minimum(g,q);candidate[ix]=g;qdist[ix]=q;ref[ix]=e.lesion_roi_id.to_numpy()[jj];refseg[ix]=e.segment_id.to_numpy()[jj]
  tree=cKDTree(np.array([[p[1]*50,p[0]*np.sqrt(3)*50] for p,(k,_) in aug.items() if k==tc[ix[0]]]))
  for k,i in enumerate(ix):
   path=np.linspace(uv[i],foot[k],max(2,int(np.ceil(g[k]/25))+1))
   if tree.query(path)[0].max()>58:status[i]='path_occluded';continue
   upper[i]=g[k]
   if q[k]>g[k]+1e-6:exact[i]=g[k];status[i]='estimable'
   else:status[i]=u.unknown_type.iloc[qj[k]] if len(u) else 'unknown_region'
 status[loc=='extent_unresolved']='extent_unresolved';status[loc=='quality_excluded']='quality_excluded'
 d['signed_distance_um']=np.where(roi,-exact,exact);d['exact_distance_status']=status;d['distance_status']=status;d['abs_distance_lower_bound_um']=lower;d['abs_distance_upper_bound_um']=upper;d['reference_lesion_roi_id']=ref;d['reference_boundary_segment_id']=refseg;d['distance_to_observed_outer_candidate_um']=candidate;d['distance_to_unobserved_edge_um']=qdist
 d['distance_definition']='Euclidean_to_new_supported_macro_edge; v3.1 all-unknown-region bounds';d['distance_zero']='hex_edge_not_spot_center';d['distance_unit']='um';d['boundary_scope']='lesion_outer'
 for w in [100,200,300]:
  band=np.full(n,'Band unresolved',object);bs=np.full(n,'interval_crosses_threshold_or_insufficient_bounds',object);near=np.isfinite(upper)&(upper<=w);far=np.isfinite(lower)&(lower>w)
  band[near&roi]='Inner margin';band[near&(loc=='outside')]='Outer margin';band[far&roi]='Interior (>w confirmed)';band[far&(loc=='outside')]='Outside (>w confirmed)';bs[near]='upper_bound_confirms_margin';bs[far]='lower_bound_confirms_non_margin';bs[np.isfinite(exact)]='exact_distance'
  band[(status=='outer_edge_not_observed')&roi]='Inside; outer edge unobserved';band[(status=='outer_edge_not_observed')&(loc=='outside')]='Outside; no usable lesion edge';bs[status=='outer_edge_not_observed']='no_outer_edge_no_distance_band';band[loc=='extent_unresolved']='Extent unresolved';bs[loc=='extent_unresolved']='extent_unresolved';band[loc=='quality_excluded']='QC excluded';bs[loc=='quality_excluded']='quality_excluded'
  d[f'distance_band_{w}um']=band;d[f'band_status_{w}um']=bs
 d['band_status']=d.band_status_200um;d['distance_band']=d.distance_band_200um;d['main_region']=d.distance_band;d['exact_distance_change_reason']='geometry rebuilt from v3.2 location and supported macro interface'
 # Replace ambiguous inherited geometry fields instead of leaving old v2/v3 results under current names.
 for name in ['tumor_region_id','boundary_side','boundary_side_two_rings','boundary_type','distance_to_boundary','signed_distance_dNN','boundary_graph_ring','display_region']:
  if name in d:
   d['legacy_'+name]=d[name];d=d.drop(columns=name)
 d['display_region']=d.lesion_location;d['boundary_type']=np.where(np.isfinite(candidate),'supported_macro_outer','no_supported_macro_outer');d['distance_to_boundary']=abs(d.signed_distance_um)
 ed.to_csv(out/'lesion_boundary_edges.tsv.gz',sep='\t',index=False);un.to_csv(out/'unknown_region_boundary.tsv.gz',sep='\t',index=False)
 for px,polys,unit in [('um',polyu,'local_Visium_nominal_um'),('image_px',poly,'aligned_image_pixels_y_down')]:
  ff=[feature(unary_union([Polygon(polys[i]) for i in ix]),{'lesion_roi_id':ids[ix[0]],'section_id':sec,'n_spots':len(ix),'status':'computational_lesion_extent; open/uncertain edges separate'}) for ix in regions]
  writegeo(out/f'lesion_rois_{px}.geojson',ff,unit)
 writegeo(out/'lesion_outer_boundaries_um.geojson',[feature(LineString(v),{'segment_id':z.segment_id,'lesion_roi_id':z.lesion_roi_id,'boundary_type':z.boundary_type}) for v,z in zip(segments(ed),ed.itertuples())],'local_Visium_nominal_um')
 return d,ed,un,poly,len(tech)

def process(sec):
 r=M.loc[sec];src=OLD/'by_dataset'/r.dataset_id/sec;out=O/'by_dataset'/r.dataset_id/sec;out.mkdir(parents=True,exist_ok=True)
 if (out/'complete.json').exists() and '--force' not in sys.argv:return
 print('START',sec,flush=True);old=pd.read_parquet(src/'spot_regions_v3_1.parquet');d=old.copy();n=len(d)
 for c in ['lesion_location','lesion_roi_id','domain_annotation','assignment_basis','reference_lesion_roi_id','signed_distance_um','exact_distance_status','abs_distance_lower_bound_um','abs_distance_upper_bound_um','distance_band_100um','distance_band_200um','distance_band_300um','band_status','main_region']:
  d['v3_1_'+c]=d[c]
 adj,lookup=graph(d);tc,tcs=components(adj,np.ones(n,bool));co=json.loads((src/'coordinate_scale.json').read_text());aff=np.array(co['array_to_image_affine']);im=np.asarray(Image.open(r.image_path).convert('RGB')) if Path(r.image_path).exists() else None
 dt=local_annotation(d,adj,lookup,tc,tcs,sec);d,ed,un,poly,nt=new_geometry(d,adj,lookup,tc,tcs,aff,co,im,sec,out)
 comp=[];d['interface_status']='extent_or_interface_unresolved'
 for k,ix in enumerate(tcs):
  ls=set(d.lesion_location.iloc[ix]);good=ed[ed.tissue_component.eq(k)&ed.boundary_type.eq('supported_macro_outer')]
  if ls=={'inside'}:st='all_inside_outer_edge_unobserved'
  elif ls=={'outside'}:st='all_outside_positive_nonlesion_context'
  elif len(good):st='supported_local_macro_interface_with_extent_caution' if 'extent_unresolved' in ls else 'supported_macro_interface'
  else:st='extent_or_interface_unresolved'
  d.loc[ix,'interface_status']=st;comp.append({'section_id':sec,'dataset_id':r.dataset_id,'tissue_component_id':f'{sec}:T{k+1}','n_spots':len(ix),'status':st,'inside':int(d.lesion_location.iloc[ix].eq('inside').sum()),'outside':int(d.lesion_location.iloc[ix].eq('outside').sum()),'extent_unresolved':int(d.lesion_location.iloc[ix].eq('extent_unresolved').sum()),'supported_outer_segments':len(good),'basis':'; '.join(sorted(set(d.assignment_basis.iloc[ix])))})
 ct=pd.DataFrame(comp);ct.to_csv(out/'tissue_component_summary.tsv',sep='\t',index=False);dt.to_csv(out/'domain_annotations.tsv',sep='\t',index=False)
 small=[]
 for rid,z in old[old.lesion_location.eq('inside')].groupby('lesion_roi_id'):
  if len(z)>=10:continue
  ix=z.index;weak=d.shared_only_threshold_crossing.iloc[ix];changed=d.lesion_location.iloc[ix].ne('inside');small.append({'section_id':sec,'old_roi_id':rid,'n_spots':len(ix),'shared_only_crossing_spots':int(weak.sum()),'n_changed':int(changed.sum()),'n_weak_changed':int((weak&changed).sum()),'new_locations':';'.join(sorted(set(d.lesion_location.iloc[ix]))),'barcodes':';'.join(z.barcode),'decision':'weak_borrowing_artifact_corrected' if (weak&changed).any() else ('local_evidence_reassessment' if changed.any() else 'retain_local_supported_or_enclosed_candidate'),'basis':';'.join(sorted(set(d.assignment_basis.iloc[ix]))),'old_roi_total_edge_segments':int(pd.read_csv(src/'lesion_boundary_edges.tsv.gz',sep='\t').lesion_roi_id.eq(rid).sum())})
 pd.DataFrame(small,columns=['section_id','old_roi_id','n_spots','shared_only_crossing_spots','n_changed','n_weak_changed','new_locations','barcodes','decision','basis','old_roi_total_edge_segments']).to_csv(out/'small_roi_decisions.tsv',sep='\t',index=False)
 # Keep original numerical observations and original identity columns, not just selected markers.
 d['region_version']='v3.2';d['roi_review_status']='computational_revision_no_new_manual_pathology_truth';d.to_parquet(out/'spot_regions_v3_2.parquet',index=False)
 cols=['dataset_id','section_id','barcode','x','y','array_row','array_col','x_um','y_um','tissue_component_id','spatial_domain','domain_component_id','domain_local_evidence_score','domain_shared_evidence_score','shared_only_threshold_crossing','v3_1_lesion_location','lesion_location','lesion_roi_id','assignment_basis','location_confidence','interface_status','reference_lesion_roi_id','reference_boundary_segment_id','signed_distance_um','exact_distance_status','abs_distance_lower_bound_um','abs_distance_upper_bound_um','distance_band_100um','distance_band_200um','distance_band_300um','band_status_100um','band_status_200um','band_status_300um','v3_1_signed_distance_um','v3_1_distance_band_200um']
 d[cols].to_csv(out/'spot_labels.tsv.gz',sep='\t',index=False);shutil.copyfile(src/'coordinate_scale.json',out/'coordinate_scale.json')
 # Unchanged v2 nest-contact object is a secondary comparator, with membership freshly linked.
 ne=pd.read_csv(V3/'by_dataset'/r.dataset_id/sec/'tumor_nest_contact_edges.tsv.gz',sep='\t');ne['inside_lesion_roi']=[bool(d.lesion_location.iloc[int(i)]=='inside' and d.lesion_location.iloc[int(j)]=='inside') for i,j in zip(ne.inside_index,ne.outside_index)];ne['region_version']='v3.2_membership_on_unchanged_v2_contact';ne.to_csv(out/'tumor_nest_contact_edges.tsv.gz',sep='\t',index=False)
 evidence=[]
 for version,lcol in [('v3.1','v3_1_lesion_location'),('v3.2','lesion_location')]:
  for loc,z in d.groupby(lcol):
   for f in ['marker_epithelial','malignant_fraction','native_spot_MSQ','marker_urothelial','marker_smooth_muscle']+['gene_'+g for g in GENES]:
    v=z[f];evidence.append({'section_id':sec,'version':version,'location':loc,'feature':f,'n':len(z),'n_available':int(v.notna().sum()),'median':v.median(),'q25':v.quantile(.25),'q75':v.quantile(.75)})
 pd.DataFrame(evidence).to_csv(out/'evidence_distributions.tsv',sep='\t',index=False)
 side=[];good=ed[ed.boundary_type.eq('supported_macro_outer')]
 for f in ['marker_epithelial','malignant_fraction']+['gene_'+g for g in GENES]:
  for ver,ee in [('v3.1',pd.read_csv(src/'lesion_boundary_edges.tsv.gz',sep='\t').query("boundary_type == 'biological_lesion_outer'")),('v3.2',good)]:
   ii=ee.inside_index.astype(int).unique();jj=ee.outside_index.astype(int).unique();side.append({'section_id':sec,'boundary_version':ver,'feature':f,'inside_n':len(ii),'outside_n':len(jj),'inside_median':d[f].iloc[ii].median(),'outside_median':d[f].iloc[jj].median(),'interpretation':'descriptive shared-expression consistency, not independent accuracy'})
 pd.DataFrame(side).to_csv(out/'boundary_evidence.tsv',sep='\t',index=False)
 # Focused checks: identities, coordinates, every original molecular value including missingness.
 assert len(d)==len(old) and not d[['section_id','barcode']].duplicated().any()
 original=[c for c in old if c.startswith(('gene_','marker_','detected_','c2l_','SpaCET_')) or c in ['MP10_original','malignant_fraction','native_spot_MSQ','CNV_support','baseline_confidence','patient_id','barcode','x','y','array_row','array_col']]
 pd.testing.assert_frame_equal(d[original],old[original],check_exact=True)
 assert (good.inside_index.map(lambda i:tc[int(i)]).to_numpy()==good.outside_index.map(lambda j:tc[int(j)]).to_numpy()).all()
 assert set(d.reference_lesion_roi_id)-{''}<=set(ed.lesion_roi_id)
 no=d.exact_distance_status.eq('outer_edge_not_observed');assert d.loc[no,'signed_distance_um'].isna().all();assert not d.loc[no,'distance_band_200um'].isin(['Inner margin','Outer margin','Interior (>w confirmed)']).any()
 for w in [100,200,300]:
  assert not ((d[f'band_status_{w}um']=='lower_bound_confirms_non_margin')&~d.abs_distance_lower_bound_um.gt(w)).any()
 from plot_v3_2 import plots
 plots(sec,d,poly,im,ed,un,co,out)
 info={'section_id':sec,'dataset_id':r.dataset_id,'n_spots':n,'inside':int(d.lesion_location.eq('inside').sum()),'outside':int(d.lesion_location.eq('outside').sum()),'extent_unresolved':int(d.lesion_location.eq('extent_unresolved').sum()),'changed_positions':int(d.lesion_location.ne(old.lesion_location).sum()),'old_outside_to_inside':int((old.lesion_location.eq('outside')&d.lesion_location.eq('inside')).sum()),'old_outside_to_unresolved':int((old.lesion_location.eq('outside')&d.lesion_location.eq('extent_unresolved')).sum()),'n_rois':int(d.loc[d.lesion_location.eq('inside'),'lesion_roi_id'].nunique()),'supported_outer_segments':len(good),'exact_estimable':int(d.signed_distance_um.notna().sum()),'band_200_determined':int(d.distance_band_200um.isin(list(BP)[:4]).sum()),'shared_crossing_spots':int(d.shared_only_threshold_crossing.sum()),'weak_shared_inside_removed':int((d.shared_only_threshold_crossing&old.lesion_location.eq('inside')&d.lesion_location.ne('inside')).sum()),'technical_geometry_cells':nt,'all_inside_blocks':int(ct.status.eq('all_inside_outer_edge_unobserved').sum()),'all_outside_blocks':int(ct.status.eq('all_outside_positive_nonlesion_context').sum()),'checks':'passed'}
 (out/'complete.json').write_text(json.dumps(info,indent=2));print('DONE',sec,json.dumps(info),flush=True)
if __name__=='__main__':
 secs=CFG['pilot'] if '--pilot' in sys.argv else ([x for x in sys.argv[1:] if not x.startswith('--')] or list(M.index))
 failed=[]
 for sec in secs:
  try:process(sec)
  except Exception:
   failed.append(sec);(O/'logs'/f'{sec}.error').write_text(traceback.format_exc());print(traceback.format_exc(),flush=True)
 if failed:raise SystemExit(1)
