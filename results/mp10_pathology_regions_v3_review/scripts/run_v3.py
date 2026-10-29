import os,sys,json,time,warnings
from pathlib import Path
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
R=Path(__file__).resolve().parents[3];O=R/'results/mp10_pathology_regions_v3';V=R/'mp10_pathology_regions_v2';sys.path.insert(0,str(R/'spatial_ecology_top20_v2/runtime'))
import numpy as np,pandas as pd
from scipy.spatial import cKDTree
from shapely.geometry import Polygon,mapping,LineString
from shapely.ops import unary_union
from collections import deque,defaultdict
from PIL import Image
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection,LineCollection
from matplotlib.patches import Patch
warnings.filterwarnings('ignore',category=RuntimeWarning)
warnings.filterwarnings('ignore',category=pd.errors.PerformanceWarning)
C=json.loads((O/'config.json').read_text());M=pd.read_csv(V/'tables/section_manifest.tsv',sep='\t').set_index('section_id');FULL=pd.read_csv(R/'neighborhood_ring_comparison_v2/00_preflight/full_image_tissue_coordinates.csv.gz',usecols=['section','array_row','array_col']);FULL={s:set(zip(d.array_row.astype(int),d.array_col.astype(int))) for s,d in FULL.groupby('section')}
OFF=[(0,2),(1,1),(1,-1),(0,-2),(-1,-1),(-1,1)];VERT=np.array([(0,2),(3,1),(3,-1),(0,-2),(-3,-1),(-3,1)])
PAL={'Tumor interior':'#a92742','Margin':'#f5a14d','Outside':'#4389b6','Inside (depth unestimated)':'#d77b8c','Outside (distance unestimated)':'#a8c9dd','No supported lesion':'#4389b6','Locally unlocalized':'#8c8c8c','QC excluded':'#303030'}
DPAL={'Tumor-related domain':'#bd3345','Interface-related domain':'#f5a14d','Non-tumor-related domain':'#4389b6','Context conflict domain':'#999999'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'axes.titlesize':10})
def graph(d):
 lookup={(int(r),int(c)):i for i,(r,c) in enumerate(zip(d.array_row,d.array_col))};adj=[[] for _ in range(len(d))]
 for (r,c),i in lookup.items():
  for dr,dc in OFF:
   j=lookup.get((r+dr,c+dc))
   if j is not None:adj[i].append(j)
 return adj,lookup

def components(adj,allowed,groups=None):
 lab=np.full(len(adj),-1,dtype=int);out=[]
 for i in np.flatnonzero(allowed):
  if lab[i]>=0:continue
  idx=len(out);q=[i];lab[i]=idx
  for j in q:
   for k in adj[j]:
    if allowed[k] and lab[k]<0 and (groups is None or groups[k]==groups[j]):lab[k]=idx;q.append(k)
  out.append(np.array(q,dtype=int))
 return lab,out

def exterior_empty(lookup):
 rr=[p[0] for p in lookup];cc=[p[1] for p in lookup];r0,r1=min(rr)-2,max(rr)+2;c0,c1=min(cc)-4,max(cc)+4;empty={(r,c) for r in range(r0,r1+1) for c in range(c0,c1+1) if (r+c)%2==sum(next(iter(lookup)))%2 and (r,c) not in lookup};seeds=[p for p in empty if p[0] in [r0,r1] or p[1]<=c0+1 or p[1]>=c1-1];ext=set(seeds);q=deque(seeds)
 while q:
  r,c=q.popleft()
  for dr,dc in OFF:
   p=(r+dr,c+dc)
   if p in empty and p not in ext:ext.add(p);q.append(p)
 return ext

def nearest(xy,segs):
 dist=np.full(len(xy),np.inf);ids=np.full(len(xy),-1,dtype=int);feet=np.full((len(xy),2),np.nan)
 if not len(segs):return dist,ids,feet
 a=segs[:,0];v=segs[:,1]-a;v2=(v*v).sum(1)
 for s in range(0,len(xy),64):
  z=xy[s:s+64,None,:]-a;t=np.clip((z*v).sum(2)/v2,0,1);dd=((z-t[:,:,None]*v)**2).sum(2);ii=dd.argmin(1);dist[s:s+64]=np.sqrt(dd[np.arange(len(ii)),ii]);ids[s:s+64]=ii;feet[s:s+64]=a[ii]+t[np.arange(len(ii)),ii,None]*v[ii]
 return dist,ids,feet

def foreground(d,im,dnn):
 h,w=im.shape[:2];rad=max(1,int(round(.275*dnn)));out=[]
 for x,y in zip(d.x,d.y):
  x,y=int(round(x)),int(round(y));p=im[max(0,y-rad):min(h,y+rad+1),max(0,x-rad):min(w,x+rad+1),:3].astype(float)/255
  if not p.size:out.append(np.nan);continue
  sat=p.max(2)-p.min(2);od=-np.log(np.clip(p.mean(2),.001,1));out.append(float(((sat>.07)|(od>.15)).mean()))
 return np.asarray(out)

def feature(geom,props):return {'type':'Feature','geometry':mapping(geom),'properties':props}
def writegeo(p,features,unit):p.write_text(json.dumps({'type':'FeatureCollection','coordinate_system':unit,'features':features},separators=(',',':')))

def process(sec,force=False):
 r=M.loc[sec];src=Path(r.output_dir);out=O/'by_dataset'/r.dataset_id/sec;out.mkdir(parents=True,exist_ok=True)
 if (out/'complete.json').exists() and not force:return
 if force and (out/'complete.json').exists():(out/'complete.json').unlink()
 print(time.strftime('%H:%M:%S'),sec,flush=True);d=pd.read_parquet(src/'spot_region_labels.parquet');n=len(d);assert n==r.n_spots and not d.barcode.duplicated().any();d=d.rename(columns={'tissue_region':'v2_tissue_region','boundary_side':'v2_boundary_side','signed_distance':'v2_signed_distance_image_px','main_tumor_mask':'v2_tumor_mask','BANKSY_domain':'spatial_domain'});d.spatial_domain=d.spatial_domain.astype(str)
 adj,lookup=graph(d);valid=d.algorithm_qc.eq('retained').to_numpy();tc,tcomps=components(adj,valid);ext=exterior_empty(lookup);xy=d[['x','y']].to_numpy();uv=np.c_[d.array_col/2,d.array_row*np.sqrt(3)/2];um=uv*100;aff=np.linalg.lstsq(np.c_[uv,np.ones(n)],xy,rcond=None)[0];keys=np.c_[3*d.array_col,3*d.array_row][:,None,:].astype(int)+VERT;vp=np.concatenate([keys[:,:,:1]/6,keys[:,:,1:]/(2*np.sqrt(3)),np.ones((n,6,1))],axis=2);polys_px=vp@aff;polys_um=vp[:,:,:2]*100;dnn=float(np.median([np.linalg.norm(xy[i]-xy[j]) for i,a in enumerate(adj) for j in a if i<j]));resid=np.linalg.norm(np.c_[uv,np.ones(n)]@aff-xy,axis=1)/dnn
 imagepath=Path(r.image_path);im=np.asarray(Image.open(imagepath).convert('RGB')) if imagepath.exists() else None;fg=foreground(d,im,dnn) if im is not None else np.full(n,np.nan);d['HE_foreground_fraction']=fg
 # The score is region-level evidence, never a probability or target-outcome optimization.
 cn=d.CNV_support.map({'strong':1.,'moderate':.55,'moderate_baseline_caution':.35,'coherent_pattern_native_diploid_conflict':.3,'weak_or_no_detectable_pattern':0.,'unavailable':0.}).fillna(0).to_numpy();epi=np.clip(d.marker_epithelial.fillna(0).to_numpy()/2.5,0,1);contract=(d.c2l_SMC.ge(.5)&d.detected_smooth_muscle.ge(3)&d.marker_epithelial.lt(.8)).to_numpy();weights=C['domain_score_weights'];score=weights['SpaCET_mean']*d.malignant_fraction.fillna(0).to_numpy()+weights['CNA_pattern_mean']*cn+weights['epithelial_program']*epi+weights['contractile_penalty']*contract
 dl,dc=components(adj,valid,d.spatial_domain.to_numpy());dom=[];dscore=np.zeros(n);da=np.full(n,'Non-tumor-related domain',dtype=object)
 for k,ix in enumerate(dc):
  z=d.iloc[ix];global_ix=np.flatnonzero(d.spatial_domain.eq(z.spatial_domain.iloc[0]).to_numpy()&(tc==tc[ix[0]]));local=float(score[ix].mean());ds=(.25*local+.75*float(score[global_ix].mean())) if len(ix)<C['tiny_domain_component_n'] else local;dscore[ix]=ds;label='Tumor-related domain' if ds>=C['tumor_related_domain_score'] else ('Interface-related domain' if ds>=C['interface_related_domain_score'] else 'Non-tumor-related domain');da[ix]=label
  dom.append(dict(section_id=sec,dataset_id=r.dataset_id,spatial_domain=z.spatial_domain.iloc[0],domain_component_id=f'{sec}:D{k+1}',n_spots=len(ix),tissue_component=int(tc[ix[0]]),domain_evidence_score=ds,domain_annotation=label,malignant_mean=z.malignant_fraction.mean(),malignant_q25=z.malignant_fraction.quantile(.25),malignant_median=z.malignant_fraction.median(),malignant_q75=z.malignant_fraction.quantile(.75),CNV_strong_fraction=z.CNV_support.eq('strong').mean(),epithelial_fraction=z.c2l_Epithelial.median(),immune_fraction=z.c2l_immune.median(),stromal_fraction=z.c2l_stroma.median(),SMC_fraction=z.c2l_SMC.median(),urothelial_marker=z.marker_urothelial.median(),HE_foreground_median=float(np.nanmedian(fg[ix])),v2_Tumor_fraction=z.v2_tumor_mask.mean(),v2_tissue_explanation=z.v2_tissue_region.value_counts().index[0]))
 roi=valid&(da=='Tumor-related domain');reason=np.full(n,'domain molecular and spatial context',dtype=object);location_uncertain=np.zeros(n,dtype=bool);presence=C['normal_source_context'].get(sec,'candidate_lesions')
 # Known adjacent contexts are adjudicated using the already inspected morphology/markers, not artificial malignant zeros.
 if sec in ['GSE285715_N1','GSE285715_N2']:
  roi[:]=False;da[:]='Non-tumor-related domain';reason[:]='adjacent-source morphology: urothelial lining / muscle; malignant evidence remains separate'
  if sec.endswith('N2'):location_uncertain=d.v2_tumor_mask.to_numpy();da[location_uncertain]='Context conflict domain'
 if sec=='GSE285715_N3':
  roi[:]=False;location_uncertain=d.v2_tissue_region.eq('Uncertain/conflicting').to_numpy();da[location_uncertain]='Context conflict domain';reason[:]='adjacent source and epithelial malignant-model conflict; candidate location locally unresolved'
 # Transition-domain components surrounded by lesion domains can be internal; no smoothing of individual spots.
 promotions=[]
 if roi.any():
  for k,ix in enumerate(dc):
   if da[ix[0]]!='Interface-related domain':continue
   outside_neighbors=[j for i in ix for j in adj[i] if dl[j]!=k];ratio=np.mean(roi[outside_neighbors]) if outside_neighbors else 0;touch_exterior=any((int(d.array_row.iloc[i])+dr,int(d.array_col.iloc[i])+co) in ext for i in ix for dr,co in OFF)
   if ratio>=.6 and not touch_exterior:promotions.extend(ix.tolist());reason[ix]='interface-related domain enclosed by lesion-related domains'
  roi[promotions]=True
  _,noncomps=components(adj,valid&~roi)
  for ix in noncomps:
   touches=any((int(d.array_row.iloc[i])+dr,int(d.array_col.iloc[i])+co) in ext for i in ix for dr,co in OFF);neighbors=[j for i in ix for j in adj[i] if roi[j]];fgok=im is None or float(np.nanmedian(fg[ix]))>=.35;viable=float(d.nFeature.iloc[ix].median())>=500
   if not touches and neighbors and fgok and viable:
    roi[ix]=True;reason[ix]='enclosed observed viable tissue compartment; composition does not define external location'
 # No spatial joining across missing spots, cracks or separate fragments; no convex hull.
 rid,regions=components(adj,roi);roi_ids=np.array([f'{sec}:L{x+1}' if x>=0 else '' for x in rid],dtype=object);d['domain_component_id']=[f'{sec}:D{x+1}' if x>=0 else '' for x in dl];d['domain_annotation']=da;d['domain_evidence_score']=dscore;d['lesion_roi_id']=roi_ids;d['lesion_location']=np.where(roi,'inside','outside');d.loc[location_uncertain,'lesion_location']='unlocalized';d.loc[~valid,'lesion_location']='quality_excluded';d['location_confidence']=np.where(location_uncertain,'unresolved',np.where(roi,'computational_ROI_context','computational_outside_context'));d['assignment_basis']=reason
 if resid.max()>.15:d['location_confidence']=d.location_confidence+'_image_alignment_caution'
 d['tissue_component_id']=[f'{sec}:T{x+1}' for x in tc];d['malignant_evidence_confidence']=d.evidence_level;d['quality_status']=d.algorithm_qc;d['boundary_scope']='lesion_outer';d['lesion_presence_status']=presence;d['roi_review_status']=C['roi_review_status'];d['x_um']=um[:,0];d['y_um']=um[:,1];d['distance_scale_source']='ordinary_Visium_nominal_100um_pitch_with_verified_array_mapping'
 # Classify all ROI edges: observed exterior, internal compartment or unobserved/censored edge.
 _,outs=components(adj,valid&~roi);outside_outer=np.zeros(n,dtype=bool)
 for ix in outs:
  touch=any((int(d.array_row.iloc[i])+dr,int(d.array_col.iloc[i])+co) in ext for i in ix for dr,co in OFF)
  outside_outer[ix]=touch
 edges=[]
 def edge_seg(i,j=None,missing=None):
  other=keys[j] if j is not None else np.array([3*missing[1],3*missing[0]])+VERT;kk=np.array(sorted(set(map(tuple,keys[i]))&set(map(tuple,other))));assert len(kk)==2;v=np.c_[kk[:,0]/6,kk[:,1]/(2*np.sqrt(3)),np.ones(2)];return v[:,:2]*100,v@aff
 for i in np.flatnonzero(roi):
  rr,cc=int(d.array_row.iloc[i]),int(d.array_col.iloc[i])
  for dr,co in OFF:
   p=(rr+dr,cc+co);j=lookup.get(p)
   if j is not None and roi[j]:continue
   if j is not None:
    typ='biological_lesion_outer' if outside_outer[j] and valid[j] and not location_uncertain[j] else ('local_contour_uncertain' if location_uncertain[j] or not valid[j] else 'internal_compartment_boundary')
   else:
    typ='QC_gap' if p in FULL.get(sec,set()) else ('unresolved_capture_or_tissue_edge' if p in ext else 'internal_tissue_defect')
    if typ=='unresolved_capture_or_tissue_edge' and im is not None:
     predicted=np.array([p[1]/2,p[0]*np.sqrt(3)/2,1])@aff
     xx,yy=predicted
     if xx<0 or yy<0 or xx>=im.shape[1] or yy>=im.shape[0]:typ='image_or_capture_truncation'
     else:
      missing_fg=foreground(pd.DataFrame({'x':[xx],'y':[yy]}),im,dnn)[0]
      typ='tissue_cut_or_surface_edge' if missing_fg<.35 else 'capture_or_sampling_truncation'
    # These are image-supported candidate edge types, not verified instrument boundaries.
   su,sp=edge_seg(i,j,missing=p);edges.append(dict(inside_index=int(i),outside_index=int(j) if j is not None else -1,lesion_roi_id=roi_ids[i],roi_number=int(rid[i]),tissue_component=int(tc[i]),boundary_scope='lesion_outer',boundary_type=typ,x1_um=su[0,0],y1_um=su[0,1],x2_um=su[1,0],y2_um=su[1,1],x1_px=sp[0,0],y1_px=sp[0,1],x2_px=sp[1,0],y2_px=sp[1,1]))
 ed=pd.DataFrame(edges,columns=['inside_index','outside_index','lesion_roi_id','roi_number','tissue_component','boundary_scope','boundary_type','x1_um','y1_um','x2_um','y2_um','x1_px','y1_px','x2_px','y2_px']);dist=np.full(n,np.nan);rawdist=np.full(n,np.nan);state=np.full(n,'no_local_lesion_boundary',dtype=object);refroi=np.full(n,'',dtype=object);invalid_dist=np.full(n,np.nan)
 def evaluate(ix,good,bad):
  if not len(ix):return
  if not len(good):state[ix]='boundary_not_observed';return
  seg=good[['x1_um','y1_um','x2_um','y2_um']].to_numpy().reshape(-1,2,2);dd,jj,foot=nearest(um[ix],seg);refroi[ix]=good.lesion_roi_id.to_numpy()[jj];rawdist[ix]=dd;state[ix]='estimable';allowed=np.ones(len(ix),bool)
  if len(bad):
   bd=nearest(um[ix],bad[['x1_um','y1_um','x2_um','y2_um']].to_numpy().reshape(-1,2,2))[0];invalid_dist[ix]=bd;badclose=bd<=dd+1e-6;allowed[badclose]=False;state[ix[badclose]]='censored_by_unobserved_or_defect_edge'
  # A straight distance path must stay in this retained-tissue component. Do not silently substitute geodesic distances.
  tree=cKDTree(um[tc==tc[ix[0]]])
  for k in np.flatnonzero(allowed):
   steps=max(2,int(np.ceil(dd[k]/25))+1);line=np.linspace(um[ix[k]],foot[k],steps);near=tree.query(line)[0]
   if near.max()>58.0:allowed[k]=False;state[ix[k]]='occluded_by_tissue_gap'
  dist[ix[allowed]]=dd[allowed]
 for k,ix in enumerate(regions):
  e=ed[ed.roi_number==k];evaluate(ix,e[e.boundary_type=='biological_lesion_outer'],e[~e.boundary_type.isin(['biological_lesion_outer','internal_compartment_boundary'])])
 for k,ix in enumerate(tcomps):
  ix=ix[~roi[ix]&~location_uncertain[ix]];e=ed[ed.tissue_component==k];evaluate(ix,e[e.boundary_type=='biological_lesion_outer'],e[~e.boundary_type.isin(['biological_lesion_outer','internal_compartment_boundary'])])
 state[location_uncertain]='location_unresolved';state[~valid]='quality_excluded';dist[location_uncertain|~valid]=np.nan;d['signed_distance_um']=np.where(roi,-dist,dist);d['distance_to_observed_outer_candidate_um']=rawdist;d['distance_to_unobserved_edge_um']=invalid_dist;d['distance_status']=state;d['reference_lesion_roi_id']=refroi;d['distance_definition']='Euclidean_to_valid_hex_ROI_edge_visible_in_same_tissue_component';d['distance_zero']='ROI_edge_not_spot_center'
 for w in C['sensitivity_half_widths_um']:
  band=np.full(n,'Outside (distance unestimated)',dtype=object);band[roi]='Inside (depth unestimated)';known=np.isfinite(dist);band[known&roi&(dist>w)]='Tumor interior';band[known&roi&(dist<=w)]='Inner margin';band[known&~roi&(dist<=w)]='Outer margin';band[known&~roi&(dist>w)]='Outside';band[location_uncertain]='Locally unlocalized';band[~valid]='QC excluded'
  if not roi.any():band[~location_uncertain&valid]='No supported lesion'
  d[f'distance_band_{w}um']=band
 d['distance_band']=d[f"distance_band_{C['main_half_width_um']}um"];d['main_region']=d.distance_band.replace({'Inner margin':'Margin','Outer margin':'Margin'})
 # Independent local tumor-nest contact object: retains v2 malignant identity only as an explicit secondary comparator.
 nest=d.v2_tumor_mask.to_numpy()&valid;nl,nc=components(adj,nest);ne=[]
 for i in np.flatnonzero(nest):
  for j in adj[i]:
   if nest[j] or not valid[j] or d.v2_tissue_region.iloc[j]=='Uncertain/conflicting':continue
   su,sp=edge_seg(i,j);ne.append(dict(inside_index=int(i),outside_index=int(j),tumor_nest_id=f'{sec}:N{nl[i]+1}',tissue_component=int(tc[i]),boundary_scope='tumor_nest_contact',inside_lesion_roi=bool(roi[i] and roi[j]),x1_um=su[0,0],y1_um=su[0,1],x2_um=su[1,0],y2_um=su[1,1],x1_px=sp[0,0],y1_px=sp[0,1],x2_px=sp[1,0],y2_px=sp[1,1]))
 ne=pd.DataFrame(ne,columns=['inside_index','outside_index','tumor_nest_id','tissue_component','boundary_scope','inside_lesion_roi','x1_um','y1_um','x2_um','y2_um','x1_px','y1_px','x2_px','y2_px']);nd=np.full(n,np.nan);nref=np.full(n,'',object);nst=np.full(n,'no_observed_contact',object)
 def nest_evaluate(ix,ee):
  if not len(ix) or not len(ee):return
  dd,jj,foot=nearest(um[ix],ee[['x1_um','y1_um','x2_um','y2_um']].to_numpy().reshape(-1,2,2));nref[ix]=ee.tumor_nest_id.to_numpy()[jj];tree=cKDTree(um[tc==tc[ix[0]]])
  for k,i in enumerate(ix):
   path=np.linspace(um[i],foot[k],max(2,int(np.ceil(dd[k]/25))+1))
   if tree.query(path)[0].max()<=58:nd[i]=dd[k];nst[i]='observed_contact_distance'
   else:nst[i]='occluded_by_tissue_gap'
 for kk,ix in enumerate(nc):nest_evaluate(ix,ne[ne.tumor_nest_id==f'{sec}:N{kk+1}'])
 for kk,ix in enumerate(tcomps):nest_evaluate(ix[~nest[ix]],ne[ne.tissue_component==kk])
 contact=pd.DataFrame({'dataset_id':r.dataset_id,'section_id':sec,'barcode':d.barcode,'boundary_scope':'tumor_nest_contact','tumor_nest_id':nref,'inside_tumor_nest':nest,'signed_distance_um':np.where(nest,-nd,nd),'distance_status':nst,'distance_definition':'Euclidean_to_observed_v2_contact_own_nest_inside_same_tissue_visible_path_outside','inside_lesion_roi':roi,'lesion_roi_id':roi_ids,'biological_status':np.where(location_uncertain,'unresolved_malignant_candidate','v2_computational_malignant_identity')});contact.to_csv(out/'tumor_nest_contact_spots.tsv.gz',sep='\t',index=False);ne.to_csv(out/'tumor_nest_contact_edges.tsv.gz',sep='\t',index=False);ed.to_csv(out/'lesion_boundary_edges.tsv.gz',sep='\t',index=False)
 # Full composition is retained, and no target-cell subtype is used above to make the ROI.
 ab=pd.read_csv(r.c2l_mean_path,index_col=0).reindex(d.barcode);assert ab.notna().all().all();fr=ab.div(ab.sum(1),axis=0)
 for col in ab:d['c2l_mean_'+col]=ab[col].to_numpy();d['c2l_fraction_'+col]=fr[col].to_numpy()
 sp=pd.read_csv(src/'spacet_major_fractions.tsv',sep='\t').set_index('barcode').reindex(d.barcode)
 for col in sp:d['SpaCET_major_fraction_'+col]=sp[col].to_numpy()
 # Original MP10 is deliberately joined only after the region rules and geometry have finished.
 fm=pd.read_csv(R/'spatial_ecology_top20_v2/00_inputs/section_manifest.tsv',sep='\t').set_index('section_id');mp=pd.read_csv(fm.loc[sec,'score_path'],index_col=0);mc='EpiMP10_Top20_AUCell';assert mc in mp and set(d.barcode)<=set(mp.index);assert mp.reindex(d.barcode)['section'].eq(sec).all();d['MP10_original']=mp.reindex(d.barcode)[mc].to_numpy();d['MP10_source_column']=mc
 d.to_parquet(out/'spot_regions_v3.parquet',index=False);d.to_csv(out/'spot_regions_v3.tsv.gz',sep='\t',index=False)
 dt=pd.DataFrame(dom);dt['domain_annotation']=dt.domain_component_id.map(d.groupby('domain_component_id').domain_annotation.first());dt['inside_ROI_fraction']=dt.domain_component_id.map(d.groupby('domain_component_id').lesion_location.apply(lambda x:x.eq('inside').mean()));dt.to_csv(out/'domain_annotations.tsv',sep='\t',index=False)
 daedges=defaultdict(int)
 for i,a in enumerate(adj):
  for j in a:
   if i<j and dl[i]!=dl[j]:daedges[tuple(sorted((int(dl[i]),int(dl[j]))))]+=1
 pd.DataFrame([{'domain_component_a':f'{sec}:D{k[0]+1}','domain_component_b':f'{sec}:D{k[1]+1}','shared_edges':v} for k,v in daedges.items()]).to_csv(out/'domain_adjacency.tsv',sep='\t',index=False)
 geou=[];geop=[]
 for k,ix in enumerate(regions):
  props={'lesion_roi_id':f'{sec}:L{k+1}','section_id':sec,'n_spots':len(ix),'review_status':C['roi_review_status']};geou.append(feature(unary_union([Polygon(polys_um[i]) for i in ix]),props));geop.append(feature(unary_union([Polygon(polys_px[i]) for i in ix]),props))
 writegeo(out/'lesion_rois_um.geojson',geou,'local_Visium_array_micrometers_not_geographic');writegeo(out/'lesion_rois_image_px.geojson',geop,'original_aligned_image_pixels_y_down')
 for ef,scope in [(ed,'lesion_outer'),(ne,'tumor_nest_contact')]:
  ff=[feature(LineString([[x.x1_um,x.y1_um],[x.x2_um,x.y2_um]]),{'boundary_scope':scope,'boundary_type':getattr(x,'boundary_type','tumor_nest_contact'),'roi_id':getattr(x,'lesion_roi_id',getattr(x,'tumor_nest_id',''))}) for x in ef.itertuples()];writegeo(out/(scope+'_boundaries_um.geojson'),ff,'local_Visium_array_micrometers_not_geographic')
 coord={'section_id':sec,'pitch_um':100,'scale_basis':'manifest ordinary_Visium and 10x manufacturer 100um center spacing','array_to_image_affine':aff.tolist(),'median_neighbor_image_px':dnn,'nominal_um_per_image_pixel':100/dnn,'affine_residual_max_pitch':float(resid.max()),'affine_residual_median_pitch':float(np.median(resid)),'coordinate_direction':'image x right/y down; physical array u=col*50, v=row*sqrt(3)*50','image_available':im is not None};(out/'coordinate_scale.json').write_text(json.dumps(coord,indent=2))
 summary={'section_id':sec,'dataset_id':r.dataset_id,'n_spots':n,'n_rois':len(regions),'n_inside':int(roi.sum()),'n_outer_segments':int(ed.boundary_type.eq('biological_lesion_outer').sum()),'n_nest_contact_segments':len(ne),'n_internal_contact_segments':int(ne.inside_lesion_roi.sum()) if len(ne) else 0,'n_distance_estimable':int(np.isfinite(dist).sum()),'n_location_unresolved':int(location_uncertain.sum()),'n_inside_depth_unestimated':int((roi&~np.isfinite(dist)).sum()),'v2_gray_resolved_inside':int((d.v2_tissue_region.isin(['Other tissue','Uncertain/conflicting'])&roi).sum()),'v2_gray_resolved_outside':int((d.v2_tissue_region.isin(['Other tissue','Uncertain/conflicting'])&~roi&~location_uncertain).sum()),'presence_status':presence,'image_alignment_caution':bool(resid.max()>.15),'n_promoted_enclosed':int(np.isin(reason,['enclosed observed viable tissue compartment; composition does not define external location','interface-related domain enclosed by lesion-related domains']).sum())}
 summary.update({'count_'+x:int(d.main_region.eq(x).sum()) for x in PAL});plot(sec,d,polys_px,im,ed,ne,out,dnn);(out/'complete.json').write_text(json.dumps(summary,indent=2));print('DONE',sec,summary['n_rois'],summary['n_inside'],summary['n_distance_estimable'],flush=True)

def plot(sec,d,polys,im,ed,ne,out,dnn):
 xy=d[['x','y']].to_numpy();pad=dnn*1.5;bounds=[xy[:,0].min()-pad,xy[:,0].max()+pad,xy[:,1].max()+pad,xy[:,1].min()-pad]
 def setup(ax,title):ax.set_title(title,loc='left');ax.set_xlim(bounds[:2]);ax.set_ylim(bounds[2:]);ax.set_aspect('equal');ax.axis('off')
 def cats(ax,val,pal,legend=True,alpha=1):
  ax.add_collection(PolyCollection(polys,facecolors=[pal.get(x,'#999999') for x in val],edgecolors='none',alpha=alpha))
  if legend:ax.legend(handles=[Patch(color=v,label=k) for k,v in pal.items() if k in set(val)],loc='upper left',bbox_to_anchor=(1,1),frameon=False,fontsize=6)
 def cont(ax,val,label,lim=None,cmap='viridis'):
  a=np.asarray(val,float);p=PolyCollection(polys,array=np.ma.masked_invalid(a),edgecolors='none',cmap=cmap);p.cmap=p.cmap.copy();p.cmap.set_bad('#bcbcbc');ax.add_collection(p)
  if lim:p.set_clim(*lim)
  plt.colorbar(p,ax=ax,fraction=.035,pad=.01,label=label)
 def lines(ax,e,color,lw=.5):
  if len(e):ax.add_collection(LineCollection(e[['x1_px','y1_px','x2_px','y2_px']].to_numpy().reshape(-1,2,2),colors=color,linewidths=lw))
 fig,axs=plt.subplots(2,3,figsize=(18,10));ax=axs.ravel();setup(ax[0],'A  H&E + lesion ROI / valid outer boundary')
 if im is not None:ax[0].imshow(im)
 else:ax[0].text(.1,.9,'H&E missing',transform=ax[0].transAxes)
 inside=d.lesion_location.eq('inside');ax[0].add_collection(PolyCollection(polys[inside],facecolors='#bd3345',alpha=.22,edgecolors='none'));lines(ax[0],ed[ed.boundary_type=='biological_lesion_outer'],'#c02040',.8);lines(ax[0],ed[ed.boundary_type!='biological_lesion_outer'],'#5b5b5b',.35)
 # A manufacturer-calibrated scale bar in the aligned image coordinate frame.
 length_px=500/100*dnn;xx=bounds[0]+pad;yy=bounds[2]-pad*.5;ax[0].plot([xx,xx+length_px],[yy,yy],c='black',lw=2);ax[0].text(xx,yy-pad*.2,'500 um (nominal Visium pitch)',fontsize=7)
 setup(ax[1],'B  BANKSY-connected domain annotation');cats(ax[1],d.domain_annotation,DPAL)
 setup(ax[2],'C  Lesion interior / margin / outside');cats(ax[2],d.main_region,PAL)
 setup(ax[3],'D  Signed distance to valid LESION OUTER edge');v=d.signed_distance_um;
 if not v.notna().any():ax[3].text(.02,.02,'No estimable lesion outer boundary',transform=ax[3].transAxes,fontsize=9)
 lim=max(300,float(np.nanquantile(abs(v),.95))) if v.notna().any() else 300;cont(ax[3],v,'um; inside < 0; gray = not estimable',(-lim,lim),'RdBu_r')
 setup(ax[4],'E  SpaCET malignant component (independent field)');cont(ax[4],d.malignant_fraction,'estimated fraction',(0,1))
 setup(ax[5],'F  Secondary nest contacts vs lesion outer edge');cats(ax[5],d.main_region,PAL,False,.25);lines(ax[5],ne,'#7752a3',.45);lines(ax[5],ed[ed.boundary_type=='biological_lesion_outer'],'#cc203a',.85);ax[5].legend(handles=[Patch(color='#cc203a',label='Lesion outer'),Patch(color='#7752a3',label='Tumor-nest contact')],loc='upper left',bbox_to_anchor=(1,1),fontsize=7,frameon=False)
 fig.suptitle(sec+' | v3 computational lesion ROI; w = 200 um',fontsize=13);fig.subplots_adjust(wspace=.55,hspace=.25,right=.87,top=.92);fig.savefig(out/'01_v3_core.png',dpi=300);fig.savefig(out/'01_v3_core.pdf');plt.close(fig)
 fig,axs=plt.subplots(1,3,figsize=(17,5));v2pal={'Tumor':'#bd3345','Normal urothelium':'#e2b637','Non-tumor stromal region':'#479b8d','Immune-rich region':'#477dc1','Smooth-muscle-rich region':'#8b67ad','Other tissue':'#d3c6b6','Uncertain/conflicting':'#999999'}
 for a,title in zip(axs,['v2 cell/tissue identity','SpaCET native 0.5 (one-ring outside)','v3 lesion-space partition; 200 um']):setup(a,title)
 cats(axs[0],d.v2_tissue_region,v2pal);cats(axs[1],d.SpaCET_native_interface,{'Tumor':'#bd3345','Interface':'#f5a14d','Stroma':'#4389b6'});cats(axs[2],d.main_region,PAL);fig.suptitle(sec);fig.subplots_adjust(right=.86,wspace=.65,top=.88);fig.savefig(out/'02_v2_SpaCET_v3_comparison.png',dpi=300);fig.savefig(out/'02_v2_SpaCET_v3_comparison.pdf');plt.close(fig)
 fig,axs=plt.subplots(1,3,figsize=(15,5))
 for ax,w in zip(axs,C['sensitivity_half_widths_um']):setup(ax,f'w = {w} um');cats(ax,d[f'distance_band_{w}um'].replace({'Inner margin':'Margin','Outer margin':'Margin'}),PAL,False)
 fig.legend(handles=[Patch(color=v,label=k) for k,v in PAL.items() if k in set(d.main_region)],loc='lower center',ncol=4,fontsize=7,frameon=False);fig.suptitle(sec+' | fixed ROI, bandwidth sensitivity');fig.subplots_adjust(bottom=.18,top=.85);fig.savefig(out/'03_bandwidth_sensitivity.png',dpi=300);fig.savefig(out/'03_bandwidth_sensitivity.pdf');plt.close(fig)
 fig,axs=plt.subplots(1,2,figsize=(12,6));setup(axs[0],'BANKSY domains (original IDs)');ids=sorted(d.spatial_domain.unique());pal={s:matplotlib.colors.to_hex(plt.get_cmap('tab20')(i%20)) for i,s in enumerate(ids)};cats(axs[0],d.spatial_domain,pal,True);setup(axs[1],'Connected domain components → tissue-region interpretation');cats(axs[1],d.domain_annotation,DPAL);fig.suptitle(sec);fig.subplots_adjust(right=.8,wspace=.5);fig.savefig(out/'04_domains.png',dpi=300);fig.savefig(out/'04_domains.pdf');plt.close(fig)

if __name__=='__main__':
 secs=C['pilot_sections'] if '--pilot' in sys.argv else ([s for s in sys.argv[1:] if not s.startswith('--')] or list(M.index))
 failures=0
 for sec in secs:
  try:process(sec,'--force' in sys.argv)
  except Exception:
   failures+=1
   import traceback;p=O/'logs'/f'{sec}.error';p.write_text(traceback.format_exc());print(p.read_text(),flush=True)

 sys.exit(1 if failures else 0)
