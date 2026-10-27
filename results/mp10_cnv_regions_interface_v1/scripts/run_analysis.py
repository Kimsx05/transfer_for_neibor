#!/usr/bin/python3.12
import os,sys,json,gzip,time,hashlib,warnings
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/'mp10_cnv_regions_interface_v1'; ECO=ROOT/'spatial_ecology_top20_v2'
sys.path.insert(0,str(ECO/'runtime'))
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
import numpy as np,pandas as pd
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from shapely.geometry import Polygon,LineString,mapping
from shapely.ops import unary_union
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from PIL import Image
warnings.filterwarnings('ignore',category=RuntimeWarning)
plt.rcParams.update({'font.size':8,'axes.titlesize':10,'axes.titlecolor':'black','pdf.fonttype':42})
def tab(d,path):
 pd.DataFrame(d).to_csv(OUT/path,sep='\t',index=False,na_rep='NA')
def read(p):return pd.read_csv(p,sep='\t')
def log(s):print(time.strftime('%H:%M:%S'),s,flush=True)
paths={'frozen_positions':ECO/'00_inputs/frozen_spot_inputs.parquet','cnv_positions':ECO/'03_copykat_geography/spot_geography.parquet','cnv_sections':ECO/'03_copykat_geography/summary_table1_sections.tsv','section_manifest':ECO/'00_inputs/section_manifest.tsv','composition_manifest':ECO/'00_inputs/c2l_manifest.tsv','images':ECO/'00_inputs/histology_manifest.tsv','identity':ECO/'patient_identity_audit_20261002/revised_by_user_v2/section_patient_and_analysis_unit_crosswalk.tsv','major_map':ECO/'00_inputs/major_map.tsv','gene_coverage':ECO/'00_inputs/gene_coverage.tsv','local_meta_niche':ECO/'04B_meta_niche/v1/spot_meta_niche_candidates.parquet','R12_niche':ROOT/'neighborhood_ring_harmony_v3/03_main/SAMPLE_R12_LEIDEN/labels.parquet','full_tissue_grid':ROOT/'neighborhood_ring_comparison_v2/00_preflight/full_image_tissue_coordinates.csv.gz','legacy_regions':ROOT/'neighborhood_ring_comparison_v2/03_features/anchor_master.parquet','R12_scores':ROOT/'neighborhood_mp10_source_v4/04_niche_comparison/anchor_MP10_all_layers.parquet'}
# Frozen BEFORE inspecting response contrasts; all qualifying sections processed.
params={'version':'1.0','primary_region':'Epi-rich region','epi_threshold':.25,'epi_source':'frozen c2l mean Epithelial / sum all 30 modeled types','cnv_role':'diagnostic annotation only; no externally validated baseline or transferable threshold','connectivity':'six true lattice neighbors AND actual distance <=1.35*d_NN','mask_geometry':'affine fitted hexagonal lattice Voronoi cells; no hull','max_fit_error_tolerance_dNN':.35,'small_components':'retain including singletons','holes':'retain all; no smoothing or filling','bandwidth':'per-section median nearest-neighbor distance in recorded image pixels','signed_convention':'inside negative; outside positive','primary_exclusion':'distance to any unsupported tissue/grid edge <=w; signed distance retained but not used in primary contrasts','distance':'exact Euclidean point-to-segment distance, within tissue connected component; inner uses own region','far_tissue':'only when distance to support edge >w','high_low':'new section-relative Epi-anchor median; >median high, <median low, ==median tied_unclassified; not frozen original state','statistics':'descriptive section paired contrasts; average repeated sections within verified patient; no spot-level inferential P','representatives':'one per dataset, lexicographically first section with image, coordinates, composition; selected before response analysis','sensitivity':'2w for representatives only','physical_area':'not available; image-pixel grid area only','threads':1}
(OUT/'regions/frozen_rules.json').write_text(json.dumps(params,indent=2))
base=pd.read_parquet(paths['frozen_positions']); cnv=pd.read_parquet(paths['cnv_positions']); sm=read(paths['section_manifest']); cm=read(paths['composition_manifest']).set_index('section_id'); hm=read(paths['images']).set_index('section_id'); ids=read(paths['identity']).set_index('section_id'); cs=read(paths['cnv_sections']).set_index('section_id'); mm=read(paths['major_map']).rename(columns={'subtype':'spatial_label'}); gc=read(paths['gene_coverage']);gc=gc[gc.MP=='EpiMP10'].set_index('section_id')
keys=['section_id','barcode'];assert not base.duplicated(keys).any() and not cnv.duplicated(keys).any()
ct=base[keys+['x','y']].merge(cnv[keys+['x','y']],on=keys,validate='one_to_one',suffixes=('_base','_cnv'));assert len(ct)==len(base) and np.allclose(ct[['x_base','y_base']],ct[['x_cnv','y_cnv']])
cnvcols=[c for c in cnv if c not in base]+keys;base=base.merge(cnv[cnvcols],on=keys,validate='one_to_one')
niche=pd.read_parquet(paths['local_meta_niche']);nc=['local_niche_id','candidate_meta_niche_id','mapping_status'];base=base.merge(niche[keys+nc],on=keys,validate='one_to_one')
r12=pd.read_parquet(paths['R12_niche']);r12['barcode']=r12.spot_uid.str.split('||',regex=False).str[-1];assert not r12.duplicated(keys).any();base=base.merge(r12[keys+['accepted_label','assignment_status']],on=keys,how='left',validate='one_to_one');base['R12_niche']=base.accepted_label.map(lambda x:'R12_niche'+str(int(x)+1) if pd.notna(x) and x>=0 else ('Unassigned' if pd.notna(x) else 'Not_in_anchor_universe'))
fullgrid=pd.read_csv(paths['full_tissue_grid']);legacy=pd.read_parquet(paths['legacy_regions'],columns=['section_id','barcode','roi_id','region_auto','region_final','region_review_status']);base=base.merge(legacy,on=keys,how='left',validate='one_to_one')
# Component-resolved score is explicitly distinguished from raw mixed AUCell; reused only on existing anchors.
scores=pd.read_parquet(paths['R12_scores']);scores['barcode']=scores.spot_uid.str.split('||',regex=False).str[-1]; scols=[c for c in ['spot_MP10_model_universe','epi_MP10_model_universe','epi_score_valid','n_signature_available_model','model_definition_id'] if c in scores];base=base.merge(scores[keys+scols],on=keys,how='left',validate='one_to_one');base['epi_MP10_valid_analysis']=base.epi_MP10_model_universe.where(base.epi_score_valid.eq(True))
reps=sm.sort_values('section_id').groupby('dataset').first().section_id.tolist();(OUT/'input_review/representative_selection.json').write_text(json.dumps({'rule':params['representatives'],'sections':reps},indent=2))
index=[dict(role=k,path=str(p),exists=p.exists(),bytes=p.stat().st_size,source_mode='read_only') for k,p in paths.items()]
for s,r in cm.iterrows():index.append(dict(role='c2l_mean',section_id=s,path=r.mean_path,exists=Path(r.mean_path).exists(),bytes=Path(r.mean_path).stat().st_size,source_mode='read_only'))
for s,r in hm.iterrows():index.append(dict(role='H&E',section_id=s,path=r.image_path,exists=Path(r.image_path).exists(),bytes=Path(r.image_path).stat().st_size,source_mode='read_only'))
tab(index,'input_review/input_index.tsv');cross=mm.rename(columns={'subtype':'spatial_label'}).copy();cross['observation_unit']='spot';cross['quantity']='c2l posterior mean expected abundance, not observed cells';cross['scRNA_mapping']='project 30-type merged spatial reference; no recovery of absent finer labels';tab(cross,'input_review/spatial_label_crosswalk.tsv')
allpos=[];alldist=[];bound=[];reg=[];coverage=[];coord=[];paramrows=[];mp=[];target=[];profiles=[];nichereg=[];nichetarget=[];nicheoverlap=[];effect=[];sens=[];cnvev=[];figindex=[];featuretables=[];checks=[]
offsets=[(0,2),(1,1),(1,-1),(0,-2),(-1,-1),(-1,1)]
vertex_offsets=np.array([[0,2],[3,1],[3,-1],[0,-2],[-3,-1],[-3,1]])
# Exact vectorized distances in small chunks. All segments retained; no nearest-midpoint approximation.
def distances(points,A,B):
 if len(A)==0:return np.full(len(points),np.nan),np.full(len(points),-1,dtype=int)
 ds=[]; js=[];v=B-A;vv=(v*v).sum(1)
 for p in np.array_split(points,max(1,int(np.ceil(len(points)/96)))):
  delta=p[:,None,:]-A[None,:,:];t=np.clip((delta*v).sum(2)/vv,0,1);sq=((delta-t[:,:,None]*v)**2).sum(2);j=sq.argmin(1);ds.extend(np.sqrt(sq[np.arange(len(p)),j]));js.extend(j)
 return np.array(ds),np.array(js)
def savefig(fig,name,desc):
 fig.savefig(OUT/'figures'/f'{name}.pdf',transparent=True,bbox_inches='tight');fig.savefig(OUT/'figures'/f'{name}.png',dpi=300,transparent=True,bbox_inches='tight');plt.close(fig);figindex.append(dict(figure=name,pdf=f'figures/{name}.pdf',png=f'figures/{name}.png',description=desc))

def process(s):
 log('PROCESS '+s);g=base[base.section_id==s].copy().reset_index(drop=True);n=len(g);ds=g.dataset.iloc[0];ident=ids.loc[s];ab=pd.read_csv(cm.loc[s,'mean_path'],index_col=0);assert ab.index.is_unique and set(g.barcode)==set(ab.index);ab=ab.loc[g.barcode].reset_index(drop=True);fine=list(ab);tot=ab.sum(axis=1);props=ab.div(tot,axis=0);assert np.allclose(props.Epithelial,g.epi_fraction_gate,atol=1e-8)
 majors={m:mm[mm.major==m].spatial_label.tolist() for m in mm.major.unique()};abmaj=pd.DataFrame({m:ab[cols].sum(axis=1) for m,cols in majors.items()});pm=abmaj.div(tot,axis=0)
 g['position_id']=g.barcode;g['dataset_id']=ds;g['patient_id']=ident.patient_uid;g['sample_id']=s;g['record_id']=ident.analysis_unit_uid;g['identity_status']=ident.patient_id_status;g['observation_unit']='spot';g['coordinate_unit']='lowres_image_pixel';g['coordinate_reference']=hm.loc[s,'image_path'];g['original_region_label']=g.region_auto;g['region_evidence']='c2l Epithelial proportion >=0.25; no aligned pathology; CNV diagnostic only';g['region_confidence']='derived_composition_only'
 xy=g[['x','y']].to_numpy();w=float(np.median(cKDTree(xy).query(xy,k=2)[0][:,1]));X=np.c_[g.array_col/2,g.array_row*np.sqrt(3)/2,np.ones(n)];coef=np.linalg.lstsq(X,xy,rcond=None)[0];err=np.linalg.norm(X@coef-xy,axis=1);assert max(err)/w<=params['max_fit_error_tolerance_dNN']
 mask=g.epi_fraction_gate.ge(.25).to_numpy();g['candidate_mask']=mask;g['final_mask']=mask;g['epi_anchor']=g.analysis_eligible_legacy;assert np.array_equal(mask,g.epi_anchor)
 lookup={(int(r.array_row),int(r.array_col)):i for i,r in g.iterrows()};pairs=[]
 for (r,c),i in lookup.items():
  for dr,dc in offsets:
   j=lookup.get((r+dr,c+dc))
   if j is not None and j>i and np.linalg.norm(xy[i]-xy[j])<=1.35*w:pairs.append((i,j))
 ij=np.array(pairs,dtype=int);row=np.r_[ij[:,0],ij[:,1]];col=np.r_[ij[:,1],ij[:,0]];adj=coo_matrix((np.ones(len(row)),(row,col)),shape=(n,n)).tocsr();_,tissuecomp=connected_components(adj);kept=mask[row]&mask[col];mat=coo_matrix((np.ones(kept.sum()),(row[kept],col[kept])),shape=(n,n)).tocsr();_,comp=connected_components(mat);insidecomps=sorted(set(comp[mask]));rmap={c:f'{s}::EPI{j+1:04d}' for j,c in enumerate(insidecomps)};rids=np.array([rmap.get(c,'Outside') if m else 'Outside' for c,m in zip(comp,mask)],dtype=object);g['region_id']=rids;g['region_label']=np.where(mask,'Epi-rich region','Outside Epi-rich regions')
 # Hex cells share identical vertex keys; no convex hull and no imputation across missing spots.
 edges={};polygons=[];edge_neighbors={};full=fullgrid[fullgrid.section==s];fullset=set(zip(full.array_row,full.array_col));edge_offsets=[(1,1),(0,2),(-1,1),(-1,-1),(0,-2),(1,-1)]
 def vertex(k):return np.array([k[0]/6,k[1]/(2*np.sqrt(3)),1.])@coef
 for i,r in g.iterrows():
  vertices=vertex_offsets+np.array([3*int(r.array_col),3*int(r.array_row)]);ks=[tuple(v) for v in vertices];polygons.append(Polygon([vertex(k) for k in ks]))
  for j in range(6):
   ek=tuple(sorted([ks[j],ks[(j+1)%6]]));edges.setdefault(ek,[]).append(i);dr,dc=edge_offsets[j];edge_neighbors[ek]=(int(r.array_row)+dr,int(r.array_col)+dc)
 validpairs=set(pairs);segments=[];support=[];allboundary=[]
 for (ka,kb),owners in edges.items():
  A=vertex(ka);B=vertex(kb)
  if len(owners)==2 and tuple(sorted(owners)) in validpairs:
   i,j=owners
   if mask[i]==mask[j]:continue
   ins=i if mask[i] else j;outs=j if mask[i] else i;btype='Epi-rich_region_interface';rid=rids[ins];segments.append((A,B,rid,tissuecomp[ins],ins,outs))
  else:
   # Pathology cannot distinguish true tissue edge, QC hole, crack and cut edge.
   i=owners[0];btype=('missing_QC_observation_edge' if edge_neighbors[(ka,kb)] in fullset else 'tissue_or_capture_edge_unresolved');rid=rids[i];support.append((A,B,tissuecomp[i]))
  bid=f'{s}::B{len(allboundary)+1:06d}';allboundary.append(dict(dataset_id=ds,section_id=s,boundary_id=bid,region_id=rid,boundary_type=btype,x1=A[0],y1=A[1],x2=B[0],y2=B[1],length_image_pixels=float(np.linalg.norm(B-A)),coordinate_unit='lowres_image_pixel',coordinate_reference=hm.loc[s,'image_path'],eligible_interface=btype=='Epi-rich_region_interface',edge_subtype_status=('full_tissue_node_without_retained_observation' if btype=='missing_QC_observation_edge' else 'tissue_edge_vs_capture_truncation_unresolved') if btype!='Epi-rich_region_interface' else 'between_two_valid_observed_spots'))
  if btype=='Epi-rich_region_interface':segments[-1]=segments[-1]+(bid,)
 bound.extend(allboundary);segA=np.array([v[0] for v in segments]);segB=np.array([v[1] for v in segments]);segR=np.array([v[2] for v in segments]);segT=np.array([v[3] for v in segments]);segID=np.array([v[6] for v in segments]);supA=np.array([v[0] for v in support]);supB=np.array([v[1] for v in support]);supT=np.array([v[2] for v in support]);dist=np.full(n,np.nan);near=np.full(n,-1,dtype=int);dsupport=np.full(n,np.nan)
 for tc in np.unique(tissuecomp):
  ix=np.flatnonzero(tissuecomp==tc);ss=np.flatnonzero(supT==tc)
  if len(ss):dsupport[ix]=distances(xy[ix],supA[ss],supB[ss])[0]
  ix=ix[~mask[ix]];ss=np.flatnonzero(segT==tc)
  if len(ix) and len(ss):dd,jj=distances(xy[ix],segA[ss],segB[ss]);dist[ix]=dd;near[ix]=ss[jj]
 for rid in rmap.values():
  ix=np.flatnonzero(rids==rid);ss=np.flatnonzero(segR==rid)
  if len(ss):dd,jj=distances(xy[ix],segA[ss],segB[ss]);dist[ix]=-dd;near[ix]=ss[jj]
 available=np.isfinite(dist);clean=available&(dsupport>w);band=np.full(n,'not_evaluable',dtype=object);band[available&(dist<-w)]='inner_core';band[available&(dist>=-w)&(dist<=0)]='inner_edge';band[available&(dist>0)&(dist<=w)]='outer_neighborhood';band[available&(dist>w)&(dsupport>w)]='far_tissue';status=np.where(~available,'no_valid_interface_in_component',np.where(dsupport<=w,'censored_near_support_edge','evaluable'))
 g['distance_band']=band;g['signed_distance']=dist;g['distance_dNN']=dist/w;g['distance_to_support_edge']=dsupport;g['evaluation_status']=status;g['tissue_component']=tissuecomp;g['nearest_interface_region_id']=[segR[j] if j>=0 else '' for j in near];g['boundary_id']=[segID[j] if j>=0 else '' for j in near];g['boundary_type']=np.where(available,'Epi-rich_region_interface','unavailable');g['signed_distance_available']=available;g['distance_to_boundary']=np.abs(dist);g['distance_unit']='lowres_image_pixel';g['missing_reason']=np.where(~available,'no Epi/non-Epi interface within connected retained tissue',np.where(dsupport<=w,'within one bandwidth of unsupported tissue/QC/cut edge',''))
 med=g.loc[mask,'EpiMP10_Top20_AUCell'].median();g['MP10_relative_anchor_state']='not_epi_anchor';g.loc[mask,'MP10_relative_anchor_state']='tied_unclassified';g.loc[mask&g.EpiMP10_Top20_AUCell.gt(med),'MP10_relative_anchor_state']='relative_high';g.loc[mask&g.EpiMP10_Top20_AUCell.lt(med),'MP10_relative_anchor_state']='relative_low';g['MP10_relative_threshold']=med
 g['MP10_linear_epi_residual']=np.nan;fit=clean&np.isfinite(g.EpiMP10_Top20_AUCell)&np.isfinite(g.epi_fraction_gate)
 if fit.sum()>=10:
  design=np.c_[np.ones(fit.sum()),g.loc[fit,'epi_fraction_gate']];beta=np.linalg.lstsq(design,g.loc[fit,'EpiMP10_Top20_AUCell'],rcond=None)[0];g.loc[fit,'MP10_linear_epi_residual']=g.loc[fit,'EpiMP10_Top20_AUCell']-design@beta
 meta=dict(dataset_id=ds,section_id=s,patient_id=ident.patient_uid,record_id=ident.analysis_unit_uid,identity_status=ident.patient_id_status)
 cellarea=abs(np.linalg.det(coef[:2]))*np.sqrt(3)/2
 for rid in rmap.values():
  ix=rids==rid;cnvev.append(dict(**meta,region_id=rid,n_positions=int(ix.sum()),profile_fraction=g.loc[ix,'profile_available'].mean(),native_aneuploid_fraction=g.loc[ix,'native_copykat_label'].eq('aneuploid').mean(),median_spot_MSQ=g.loc[ix,'spot_MSQ'].median(),baseline_confidence=cs.loc[s,'baseline_confidence'],cnv_supported_status='Unknown_unvalidated_baseline',pathology_concordance='not_evaluable_no_pathology',evidence='Epi-only boundary; native CNV overlap descriptive'))
  reg.append(dict(**meta,region_id=rid,region_label='Epi-rich region',n_positions=int(ix.sum()),tissue_position_fraction=float(ix.mean()),grid_area_image_pixel_squared=cellarea*ix.sum(),physical_area_um2=np.nan,tissue_component=int(tissuecomp[np.flatnonzero(ix)[0]]),n_valid_interface_segments=int((segR==rid).sum()),n_evaluable_core=int((ix&clean&(band=='inner_core')).sum()),n_evaluable_inner_edge=int((ix&clean&(band=='inner_edge')).sum()),source='frozen Epi gate; six-neighbor components; raw mask unchanged',core_status='evaluable' if (ix&clean&(band=='inner_core')).any() else 'not_evaluable'))
 # Region masks and valid/support boundary geometries use a local image coordinate reference, not geographic WGS84.
 features=[]
 for rid in rmap.values():
  geom=unary_union([polygons[i] for i in np.flatnonzero(rids==rid)]);features.append(dict(type='Feature',properties={'region_id':rid,'region_label':'Epi-rich region','coordinate_unit':'lowres_image_pixel'},geometry=mapping(geom)))
 tissue=unary_union(polygons);features.append(dict(type='Feature',properties={'region_id':'valid_observation_support','region_label':'QC retained tissue lattice support'},geometry=mapping(tissue)))
 with gzip.open(OUT/'regions/masks'/f'{s}.geojson.gz','wt') as f:json.dump(dict(type='FeatureCollection',coordinate_reference=hm.loc[s,'image_path'],coordinate_unit='lowres_image_pixel',features=features),f)
 with gzip.open(OUT/'regions/masks'/f'{s}_outlines.geojson.gz','wt') as f:json.dump(dict(type='FeatureCollection',coordinate_reference=hm.loc[s,'image_path'],features=[dict(type='Feature',properties={k:v for k,v in b.items() if k not in ['x1','x2','y1','y2']},geometry=mapping(LineString([(b['x1'],b['y1']),(b['x2'],b['y2'])]))) for b in allboundary]),f)
 coord.append(dict(**meta,coordinate_reference=hm.loc[s,'image_path'],coordinate_unit='lowres_image_pixel',transform_grid_to_image=json.dumps(coef.tolist()),max_fit_error_px=max(err),max_fit_error_dNN=max(err)/w,median_fit_error_px=np.median(err),determinant=np.linalg.det(coef[:2]),image_transform='none; x right y down; origin upper; mirror/rotation retained in affine',geometry_precision='inferred_lattice_pixel_geometry',image_in_bounds=bool(((xy[:,0]>=0)&(xy[:,0]<hm.loc[s,'image_width'])&(xy[:,1]>=0)&(xy[:,1]<hm.loc[s,'image_height'])).all())))
 paramrows.append(dict(**meta,epi_threshold=.25,bandwidth_w=w,adjacency_cap_px=1.35*w,connectivity='six lattice neighbors',min_component_size=1,holes_filled=0,smoothing='none',relative_MP10_median=med,coordinate_unit='lowres_image_pixel',physical_scale_status='not_verified',version='1.0'))
 coverage.append(dict(**meta,n_positions=n,n_cnv_profiles=int(g.profile_available.sum()),n_regions=len(rmap),n_epi_anchors=int(mask.sum()),n_valid_interfaces=len(segments),n_signed_distance=int(available.sum()),n_primary_evaluable=int(clean.sum()),n_censored=int((available&~clean).sum()),n_R12_assigned=int(g.accepted_label.ge(0).sum()),n_local_niche=int(g.local_niche_id.notna().sum()),n_MP10_available=int(g.EpiMP10_Top20_AUCell.notna().sum()),MP10_genes_available=int(gc.loc[s,'n_available']),observation_unit='spot',platform='ordinary_Visium',resolution='original spot grid; physical scale not newly calibrated',status='completed_with_limitations',missing_reason='aligned pathology absent; CNV baseline not externally validated; support edges unresolved'))
 # Parallel dimensions of abundance: mean estimate, mean local fraction, pooled abundance fraction and lineage fraction.
 types=[('major',m,m,abmaj[m],pm[m],np.ones(n)) for m in majors]
 types += [('fine',t,mm.set_index('spatial_label').loc[t,'major'],ab[t],props[t],ab[t]/abmaj[mm.set_index('spatial_label').loc[t,'major']]) for t in fine]
 groups=[('distance_band',b,(band==b)&clean) for b in ['inner_core','inner_edge','outer_neighborhood','far_tissue']]+[('region',rid,rids==rid) for rid in rmap.values()]+[('whole_tissue','all',np.ones(n,dtype=bool))]
 def mp_record(scope,label,ix):
  vals=g.loc[ix,'EpiMP10_Top20_AUCell'];a=ix&mask;return dict(**meta,scope=scope,group=label,n_positions=int(ix.sum()),n_epi_anchors=int(a.sum()),n_MP10=int(vals.notna().sum()),MP10_mean=vals.mean(),MP10_median=vals.median(),MP10_anchor_mean=g.loc[a,'EpiMP10_Top20_AUCell'].mean(),epi_fraction_mean=g.loc[ix,'epi_fraction_gate'].mean(),n_high_anchors=int((a&g.MP10_relative_anchor_state.eq('relative_high')).sum()),n_low_anchors=int((a&g.MP10_relative_anchor_state.eq('relative_low')).sum()),n_tied_anchors=int((a&g.MP10_relative_anchor_state.eq('tied_unclassified')).sum()),high_fraction_among_anchors=((a&g.MP10_relative_anchor_state.eq('relative_high')).sum()/a.sum() if a.sum() else np.nan),score_meaning='mixed spot raw Top20 AUCell',epi_component_MP10_mean=g.loc[ix,'epi_MP10_valid_analysis'].mean() if 'epi_MP10_model_universe' in g else np.nan,n_epi_component_scores=int(g.loc[ix,'epi_MP10_valid_analysis'].notna().sum()) if 'epi_MP10_model_universe' in g else 0,evaluation_status='evaluable' if ix.sum() else 'not_evaluable')
 for scope,label,ix in groups:
  mp.append(mp_record(scope,label,ix))
  for level,t,line,a,p,lp in types:
   target.append(dict(**meta,scope=scope,group=label,level=level,target=t,lineage=line,n_positions=int(ix.sum()),n_epi_anchors=int((ix&mask).sum()),mean_model_abundance=a[ix].mean(),mean_local_fraction=p[ix].mean(),pooled_modeled_fraction=(a[ix].sum()/tot[ix].sum() if ix.sum() else np.nan),mean_lineage_fraction=np.asarray(lp)[ix].mean() if ix.sum() else np.nan,units='c2l modeled expected abundance per spot; proportions use all 30 types or named lineage',evaluation_status='evaluable' if ix.sum() else 'not_evaluable'))
 # Continuous profiles in fixed dNN intervals; no smoothing or tail extrapolation.
 for left,right in zip([-np.inf,-4,-3,-2,-1,0,1,2,3,4],[-4,-3,-2,-1,0,1,2,3,4,np.inf]):
  ix=clean&(dist/w>left)&(dist/w<=right)
  for level,t,line,a,p,lp in types:
   profiles.append(dict(**meta,distance_left_dNN=left,distance_right_dNN=right,target=t,level=level,n_positions=int(ix.sum()),n_epi_anchors=int((ix&mask).sum()),mean_model_abundance=a[ix].mean(),mean_local_fraction=p[ix].mean(),mean_lineage_fraction=np.asarray(lp)[ix].mean() if ix.sum() else np.nan,MP10_mean=g.loc[ix,'EpiMP10_Top20_AUCell'].mean(),epi_fraction_mean=g.loc[ix,'epi_fraction_gate'].mean(),evaluation_status='evaluable' if ix.sum() else 'not_evaluable'))
 metrics={'MP10_mixed_AUCell':g.EpiMP10_Top20_AUCell,'Epi_fraction':g.epi_fraction_gate,'MP10_linear_epi_residual':g.MP10_linear_epi_residual}
 if 'epi_MP10_model_universe' in g:metrics['MP10_epi_component_existing']=g.epi_MP10_valid_analysis
 for level,t,line,a,p,lp in types:metrics[f'{level}:{t}:all_modeled_fraction']=p
 metrics['fine:Treg:T_NK_fraction']=ab.Treg/abmaj.T_NK
 for lhs,rhs in [('inner_edge','inner_core'),('outer_neighborhood','inner_edge')]:
  for metric,v in metrics.items():
   a=v[(band==lhs)&clean].dropna();b=v[(band==rhs)&clean].dropna();effect.append(dict(**meta,level='section',contrast=lhs+' minus '+rhs,metric=metric,n_left=len(a),n_right=len(b),mean_left=a.mean(),mean_right=b.mean(),effect=a.mean()-b.mean() if len(a) and len(b) else np.nan,status='evaluable' if len(a) and len(b) else 'not_evaluable'))
 for coln in ['local_niche_id','candidate_meta_niche_id','R12_niche']:
  for lab,ix0 in g.groupby(coln,dropna=False).groups.items():
   ix=np.zeros(n,bool);ix[list(ix0)]=True
   for b in ['inner_core','inner_edge','outer_neighborhood','far_tissue','not_evaluable']:
    for regionlabel in ['Epi-rich region','Outside Epi-rich regions']:
     z=ix&(band==b)&g.region_label.eq(regionlabel).to_numpy();nichereg.append(dict(**meta,niche_system=coln,niche_label=lab,region_label=regionlabel,distance_band=b,n_positions=int(z.sum()),n_primary_evaluable=int((z&clean).sum()),niche_status='candidate_not_confirmatory' if coln=='candidate_meta_niche_id' else 'original_frozen_assignment'))
   for level,t,line,a,p,lp in types:nichetarget.append(dict(**meta,niche_system=coln,niche_label=lab,level=level,target=t,n_positions=int(ix.sum()),MP10_mean=g.loc[ix,'EpiMP10_Top20_AUCell'].mean(),epi_fraction_mean=g.loc[ix,'epi_fraction_gate'].mean(),mean_model_abundance=a[ix].mean(),mean_local_fraction=p[ix].mean(),mean_lineage_fraction=np.asarray(lp)[ix].mean()))
 for (l,r),z in g.groupby(['local_niche_id','R12_niche']):nicheoverlap.append(dict(**meta,local_niche_id=l,R12_niche=r,n_positions=len(z)))
 common=['dataset_id','patient_id','sample_id','record_id','identity_status','section_id','position_id','observation_unit','x','y','coordinate_unit','coordinate_reference','region_id','region_label','region_evidence','region_confidence']
 allpos.append(g[common+['candidate_mask','final_mask','original_region_label','roi_id','region_final','region_review_status','evaluation_status','tissue_component']]);dd=g[common+['nearest_interface_region_id','boundary_id','boundary_type','distance_to_boundary','signed_distance','distance_dNN','distance_unit','signed_distance_available','distance_to_support_edge','distance_band','evaluation_status','missing_reason']].copy();dd['position_region_id']=dd.region_id;dd['region_id']=dd.nearest_interface_region_id.where(dd.signed_distance_available,dd.region_id);dd['region_label']=np.where(dd.signed_distance_available,'Epi-rich region',dd.region_label);alldist.append(dd)
 featcols=['dataset_id','section_id','position_id','epi_fraction_gate','analysis_eligible_legacy','rctd_anchor_legacy','EpiMP10_Top20_AUCell','MP10_relative_anchor_state','MP10_relative_threshold','native_copykat_label','profile_available','spot_MSQ','cluster_consensus_MSQ','baseline_confidence','mixture_nonconvergence','native_consensus_correlation','local_niche_id','candidate_meta_niche_id','mapping_status','R12_niche']+scols
 featuretables.append(g[featcols]);checks.append(dict(section_id=s,n_positions=n,unique_keys=not g.position_id.duplicated().any(),composition_matches_frozen_gate=True,mask_preserved=bool(np.array_equal(mask,g.final_mask)),signed_inside_negative=bool((dist[mask&available]<0).all()),signed_outside_positive=bool((dist[~mask&available]>0).all()),n_inner_without_interface=int((mask&~available).sum()),max_affine_error_dNN=max(err)/w))
 if s in reps:
  for scale in [1,2]:
   ww=w*scale;ok=available&(dsupport>ww);bb=np.where(dist<-ww,'inner_core',np.where(dist<=0,'inner_edge',np.where(dist<=ww,'outer_neighborhood','far_tissue')))
   for b in ['inner_core','inner_edge','outer_neighborhood','far_tissue']:
    ix=ok&(bb==b);sens.append(dict(**meta,bandwidth_multiplier=scale,distance_band=b,n_positions=int(ix.sum()),MP10_mean=g.loc[ix,'EpiMP10_Top20_AUCell'].mean(),Treg_all_modeled_fraction=props.loc[ix,'Treg'].mean(),Treg_T_NK_fraction=(ab.Treg/abmaj.T_NK)[ix].mean(),Epi_fraction=g.loc[ix,'epi_fraction_gate'].mean()))
  fig,axs=plt.subplots(3,4,figsize=(16,12));axs=axs.ravel();img=Image.open(hm.loc[s,'image_path']);ext=[0,img.width,img.height,0];allseg=np.array([[v[0],v[1]] for v in segments]);supseg=np.array([[v[0],v[1]] for v in support]);nichevals=pd.Categorical(g.R12_niche).codes
  vals=[None,g.spot_MSQ, g.epi_fraction_gate,mask.astype(int),g.EpiMP10_Top20_AUCell,g.MP10_relative_anchor_state.map({'not_epi_anchor':0,'relative_low':1,'tied_unclassified':2,'relative_high':3}),pm.T_NK,props.Treg,pm.Endothelial,pm.Fibroblast,pm.Mural,nichevals]
  titles=['H&E + retained tissue support','CopyKAT spot MSQ (diagnostic)','Epithelial modeled fraction','Derived Epi-rich regions','MP10 mixed-spot Top20 AUCell','MP10 relative anchor state','T/NK modeled fraction','Treg / all 30 modeled types','Endothelial modeled fraction','Fibroblast modeled fraction','Pericyte + SMC modeled fraction','Existing R1/R2 niche codes']
  for j,(ax,v,title) in enumerate(zip(axs,vals,titles)):
   ax.imshow(img,origin='upper',alpha=1 if j==0 else .12)
   if v is not None:
    sc=ax.scatter(xy[:,0],xy[:,1],c=v,s=max(2,(240*w/np.ptp(xy[:,0]))**2*.55),cmap='tab10' if j in [3,5,11] else 'viridis',linewidths=0);fig.colorbar(sc,ax=ax,shrink=.65,pad=.01)
   else:ax.scatter(xy[:,0],xy[:,1],s=.3,c='cyan',alpha=.6)
   if len(allseg):ax.add_collection(LineCollection(allseg,colors='black',linewidths=.3))
   if j in [0,3] and len(supseg):ax.add_collection(LineCollection(supseg,colors='#d95f02',linewidths=.35))
   ax.set_xlim(max(0,xy[:,0].min()-w*2),min(img.width,xy[:,0].max()+w*2));ax.set_ylim(min(img.height,xy[:,1].max()+w*2),max(0,xy[:,1].min()-w*2));ax.set_aspect('equal');ax.set_title(title,loc='center',color='black');ax.set_xlabel('image x (pixel)');ax.set_ylabel('image y (pixel)')
  fig.suptitle(f'{s} | {n:,} Visium spots | Epi-rich edges, no pathology confirmation',color='black');fig.text(.5,.01,'Black: valid Epi-rich interface. Orange: unsupported tissue/QC/cut edge. CNV baseline unvalidated. Fractions are modeled composition.\nMP10 state codes: 0 non-anchor; 1 relative low; 2 tied; 3 relative high. Niche codes are display-only; mapping in table.',ha='center');fig.tight_layout(rect=[0,.035,1,.965]);savefig(fig,f'{s}_spatial_overview','Aligned source H&E, inferred masks, native CNV, mixed MP10 and modeled target composition; 12 panels')
  tab([dict(display_code=i,R12_niche=l) for i,l in enumerate(pd.Categorical(g.R12_niche).categories)],f'niche_review/{s}_display_codes.tsv')
 return g

if __name__=='__main__':
 order=reps+[s for s in sorted(sm.section_id) if s not in reps]
 # --pilot performs only the frozen representatives so visual validation precedes batch.
 if '--pilot' in sys.argv:order=reps
 for s in order:process(s)
 tab(pd.concat(allpos,ignore_index=True),'regions/position_regions.tsv.gz');tab(pd.concat(alldist,ignore_index=True),'interfaces/position_interface_distances.tsv.gz');tab(pd.concat(featuretables,ignore_index=True),'input_review/position_features.tsv.gz');tab(reg,'regions/region_summary.tsv');tab(bound,'interfaces/boundary_summary.tsv.gz');tab(paramrows,'regions/region_parameters.tsv');tab(coverage,'input_review/section_coverage.tsv');tab(coord,'input_review/coordinate_alignment.tsv');tab(cnvev,'cnv_review/cnv_region_evidence.tsv');tab(mp,'target_profiles/region_mp10_summary.tsv');tab(target,'target_profiles/target_region_summary.tsv.gz');tab(profiles,'target_profiles/target_distance_profiles.tsv');tab(nichereg,'niche_review/niche_region_crosswalk.tsv');tab(nichetarget,'niche_review/niche_target_summary.tsv.gz');tab(nicheoverlap,'niche_review/local_R12_overlap.tsv');tab(effect,'target_profiles/patient_section_effects.tsv');tab(sens,'interfaces/representative_2w_sensitivity.tsv');tab(figindex,'figures/figure_index.tsv');tab(checks,'logs/section_validation.tsv')
 c=cs.reset_index();c['analysis_unit']='observed mixed-expression Visium spot';c['source']=str(paths['cnv_positions']);c['space_mapping']='exact section+full barcode; x/y identical to frozen positions';c['reference']='section automatic tissue baseline';c['score_definition']='mean(native final CNA profile squared), no rescale';c['original_threshold']='native CopyKAT labels retained; no transferable MSQ threshold';c['derived_CNV_threshold']='none';c['spatial_diagnostic_available']=True;c['sufficient_for_primary_CNV_region']=False;c['reason']='41 automatic but externally unvalidated baselines; 22 additional uncertain/nonconvergent; no validated pathology or scale threshold';c['scale_policy']='do not pool raw MSQ across sections';tab(c,'cnv_review/cnv_availability.tsv')
 log('ANALYSIS DONE')
