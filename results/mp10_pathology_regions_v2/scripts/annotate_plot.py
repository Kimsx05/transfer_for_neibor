import os,sys,json,time,warnings
from pathlib import Path
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
R=Path(__file__).resolve().parents[2];O=R/'mp10_pathology_regions_v2';sys.path.insert(0,str(R/'spatial_ecology_top20_v2/runtime'))
import numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection,LineCollection
from matplotlib.patches import Patch
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from collections import deque,defaultdict
from PIL import Image
warnings.filterwarnings('ignore',category=RuntimeWarning)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'axes.titlesize':10,'savefig.facecolor':'white'})
PAL={'Tumor':'#bd3345','Normal urothelium':'#e2b637','Non-tumor stromal region':'#479b8d','Immune-rich region':'#477dc1','Smooth-muscle-rich region':'#8b67ad','Necrotic/acellular region':'#bba98a','Other tissue':'#d3c6b6','Uncertain/conflicting':'#999999','Algorithm QC excluded':'#333333'}
BPAL={'tumor_interior':'#9f1734','tumor_inner_margin':'#ed7566','peritumoral_outer_interface':'#2ca7bd','distal_non_tumor':'#d0d9ce','unassessable':'#999999'}
CPAL={'strong':'#8d1b2d','moderate':'#e9904f','moderate_baseline_caution':'#d8b54c','coherent_pattern_native_diploid_conflict':'#966aaf','weak_or_no_detectable_pattern':'#e0e1d8','unavailable':'#999999'}
OFF=[(0,2),(1,1),(1,-1),(0,-2),(-1,-1),(-1,1)]
VERT=np.array([(0,2),(3,1),(3,-1),(0,-2),(-3,-1),(-3,1)])
man=pd.read_csv(O/'tables/section_manifest.tsv',sep='\t').set_index('section_id');mapping=pd.read_csv(R/'spatial_ecology_top20_v2/00_inputs/major_map.tsv',sep='\t')
def graph(df):
 lookup={(int(r),int(c)):i for i,(r,c) in enumerate(zip(df.array_row,df.array_col))};adj=[[] for _ in range(len(df))]
 for (r,c),i in lookup.items():
  for dr,dc in OFF:
   j=lookup.get((r+dr,c+dc))
   if j is not None:adj[i].append(j)
 return adj,lookup

def hops(adj,seeds,allowed):
 d=np.full(len(adj),np.inf);q=deque()
 for i in seeds:d[i]=1;q.append(i)
 while q:
  i=q.popleft()
  for j in adj[i]:
   if allowed[j] and d[j]>d[i]+1:d[j]=d[i]+1;q.append(j)
 return d

def segment_distance(xy,seg):
 out=np.full(len(xy),np.inf)
 if not len(seg):return out
 a=seg[:,0];v=seg[:,1]-a;v2=(v*v).sum(1)
 for s in range(0,len(xy),128):
  z=xy[s:s+128,None,:]-a[None,:,:];t=np.clip((z*v).sum(2)/v2,0,1);dd=np.sum((z-t[:,:,None]*v)**2,2);out[s:s+128]=np.sqrt(dd.min(1))
 return out

def annotate(sec,force=False):
 r=man.loc[sec];dest=Path(r.output_dir)
 if (dest/'annotation.done').exists() and (dest/'annotation.done').stat().st_mtime>(dest/'cnv_v3.done').stat().st_mtime and not force:return
 if not (dest/'marker_evidence.tsv').exists():return
 print(time.strftime('%H:%M:%S'),'annotate',sec,flush=True)
 d=pd.read_csv(dest/'algorithm_qc.tsv',sep='\t');n=len(d)
 for fn in ['spacet_spots.tsv','cnv_spot_support.tsv','marker_evidence.tsv','banksy_domains.tsv']:
  if (dest/fn).exists():
   t=pd.read_csv(dest/fn,sep='\t');d=d.merge(t.drop(columns=[c for c in t if c in d and c!='barcode']),on='barcode',how='left',validate='one_to_one')
 for c,default in [('malignant_fraction',np.nan),('SpaCET_signature_source','unavailable'),('SpaCET_major_label','unavailable'),('SpaCET_native_interface','unavailable'),('CNV_support','unavailable')]:
  if c not in d:d[c]=default
  if isinstance(default,str):d[c]=d[c].fillna(default)
 domaincols=[c for c in d if c.startswith('clust_')];main=[c for c in domaincols if 'res0.6' in c or 'res_0.6' in c];main=main or domaincols[-1:];d['BANKSY_domain']=d[main[0]].astype(str) if main else 'unavailable'
 ab=pd.read_csv(r.c2l_mean_path,index_col=0).reindex(d.barcode);fr=ab.div(ab.sum(axis=1),axis=0)
 for g in mapping.major.unique():d['c2l_'+g]=fr[mapping.loc[mapping.major==g,'subtype']].sum(axis=1).to_numpy()
 d['c2l_SMC']=fr.SMC.to_numpy();d['c2l_immune']=d[['c2l_'+g for g in ['T_NK','B','Plasma','Myeloid','Mast']]].sum(axis=1)
 d['c2l_stroma']=d[['c2l_'+g for g in ['Fibroblast','Endothelial','Mural']]].sum(axis=1)
 d['domain_malignant_median']=d.groupby('BANKSY_domain').malignant_fraction.transform('median')
 d['domain_native_aneuploid_fraction']=d.native_copykat_label.eq('aneuploid').groupby(d.BANKSY_domain).transform('mean')
 mal=d.malignant_fraction.fillna(-1);strong=d.CNV_support.eq('strong');cnvpattern=d.CNV_support.isin(['strong','moderate','moderate_baseline_caution','coherent_pattern_native_diploid_conflict']);depth=d.SpaCET_signature_source.str.contains('seq_depth',case=False)
 epiexpr=(d.detected_epithelial>=3)&(d.marker_epithelial>=.8)
 imm=(d.c2l_immune>=.5)&(d.detected_immune>=2);sm=(d.c2l_SMC>=.25)&(d.detected_smooth_muscle>=3)&(d.gene_MYH11>.5)&((d.gene_ACTG2>.5)|(d.gene_CNN1>.5))
 strom=(d.c2l_stroma>=.45)&((d.detected_fibroblast>=3)|(d.detected_endothelial>=2)|(d.detected_smooth_muscle>=3))
 d['tissue_region']='Other tissue';d['annotation_reason']='mixed composition; insufficient evidence for a specific tissue identity'
 for mask,label,reason in [(strom,'Non-tumor stromal region','stromal composition with matching multi-gene expression'),(imm,'Immune-rich region','immune major-lineage fraction >=0.5 with multi-gene immune expression'),(sm,'Smooth-muscle-rich region','SMC fraction >=0.25 with MYH11 and ACTG2/CNN1 plus >=3 muscle markers')]:d.loc[mask,['tissue_region','annotation_reason']]=[label,reason]
 initial=(mal>=.5)&~depth
 # Major-lineage maximum supports mixed tumor; CNA can rescue low-epithelial spots without an epithelial fraction gate.
 mixed=(mal>=.25)&strong&(d.SpaCET_major_label.eq('Malignant')|d.domain_malignant_median.ge(.35))
 cnvrescue=strong&d.native_copykat_label.eq('aneuploid')&epiexpr&d.domain_malignant_median.ge(.15)
 tumor=initial|mixed|cnvrescue
 d.loc[tumor,['tissue_region','annotation_reason']]=['Tumor','SpaCET BLCA malignant estimate >=0.5; computational tumor inference']
 d.loc[mixed&~initial,'annotation_reason']='mixed malignant estimate with coherent native-aneuploid CNA support and domain context'
 d.loc[cnvrescue&~initial&~mixed,'annotation_reason']='coherent native-aneuploid CNA with epithelial expression and malignant domain context; no epithelial fraction gate'
 # Composition inconsistency is localized, not propagated to an entire BANKSY domain.
 conflict=tumor&((d.c2l_immune>=.8)&~epiexpr&~strong)
 conflict|=tumor&sm&(d.c2l_SMC>=.5)&~epiexpr
 conflict|=depth&(mal>=.5)&~strong
 normalcontext=d.source_context.eq('adjacent_normal_source')
 conflict|=normalcontext&tumor&~strong
 uro=(d.detected_urothelial>=2)&(d.marker_urothelial>=.5)&epiexpr
 normal=normalcontext&uro&(mal<.35)&~strong&(d.marker_proliferation<1.5)
 d.loc[normal,['tissue_region','annotation_reason']]=['Normal urothelium','inferred from adjacent-normal source plus differentiated urothelial expression and low malignant estimate; not histology-confirmed']
 # Urothelial differentiation alone is not evidence of benignity in a tumor specimen.
 epithelial_unresolved=~tumor&~normal&epiexpr&(d.c2l_Epithelial>=.35)&(d.tissue_region=='Other tissue')
 conflict|=epithelial_unresolved
 d.loc[conflict,['tissue_region','annotation_reason']]=['Uncertain/conflicting','malignancy or urothelial identity conflicts with composition/source; local review required']
 d['annotation_conflict']=np.where(conflict,'local_malignancy_or_identity_conflict','none')
 morph=pd.read_csv(O/'tables/morphology_context.tsv',sep='\t',dtype={'BANKSY_domain':str});mdom=set(morph.loc[morph.section_id.eq(sec),'BANKSY_domain'])
 morphnormal=normalcontext&d.BANKSY_domain.isin(mdom)&uro&d.native_copykat_label.eq('diploid')&(d.marker_proliferation<1.5)&~strong
 d.loc[morphnormal,['tissue_region','annotation_reason']]=['Normal urothelium','inferred normal urothelial lining: adjacent-normal source, surface lining on H&E comparison, differentiated urothelial markers and low proliferation; SpaCET conflict retained']
 d.loc[morphnormal&(mal>=.5),'annotation_conflict']='SpaCET_high_malignant_fraction_in_context_supported_normal_urothelium'
 normal=normal|morphnormal
 morphmuscle=normalcontext&('all_muscle_supported_spots' in mdom)&sm&~epiexpr
 d.loc[morphmuscle,['tissue_region','annotation_reason']]=['Smooth-muscle-rich region','H&E muscle-bundle context with SMC composition and multi-gene contractile program; low epithelial expression; discordant malignant/CNA estimates retained']
 d.loc[morphmuscle&((mal>=.5)|strong),'annotation_conflict']='SpaCET_or_CNA_malignant_support_in_context_supported_smooth_muscle'
 conflict=conflict&~morphnormal&~morphmuscle
 excluded=d.algorithm_qc.ne('retained');d.loc[excluded,['tissue_region','annotation_reason']]=['Uncertain/conflicting','algorithm QC excluded; original position retained']
 d['evidence_level']='computational_composition_marker';d.loc[d.tissue_region.eq('Tumor'),'evidence_level']='computational_moderate';d.loc[d.tissue_region.eq('Tumor')&strong&(mal>=.25),'evidence_level']='computational_high_concordance';d.loc[normal,'evidence_level']='computational_context_supported';d.loc[conflict|excluded,'evidence_level']='unresolved_local';d['review_status']='no_documented_human_pathology_review';d['main_tumor_mask']=d.tissue_region.eq('Tumor');d['initial_spacet_cutoff_mask']=mal>=.5;d['mask_adjustment_reason']=np.where(d.main_tumor_mask==d.initial_spacet_cutoff_mask,'unchanged',d.annotation_reason)
 adj,lookup=graph(d);t=d.main_tumor_mask.to_numpy();valid=~d.tissue_region.eq('Uncertain/conflicting').to_numpy();xy=d[['x','y']].to_numpy();ideal=np.c_[d.array_col/2,d.array_row*np.sqrt(3)/2,np.ones(n)];aff=np.linalg.lstsq(ideal,xy,rcond=None)[0];fit=ideal@aff;edge_lengths=[np.linalg.norm(xy[i]-xy[j]) for i,js in enumerate(adj) for j in js if i<j];dnn=float(np.median(edge_lengths)) if edge_lengths else 1.;keys=(np.c_[3*d.array_col,3*d.array_row][:,None,:]+VERT).astype(int);polys=np.concatenate([keys[:,:,:1]/6,keys[:,:,1:]/(2*np.sqrt(3)),np.ones((n,6,1))],axis=2)@aff
 ti=np.flatnonzero(t);tid=np.full(n,-1);reg=0
 for i in ti:
  if tid[i]>=0:continue
  q=[i];tid[i]=reg
  for a in q:
   for j in adj[a]:
    if t[j] and tid[j]<0:tid[j]=reg;q.append(j)
  reg+=1
 edges=[];segments=[];inner=set();outer=set();types=defaultdict(set);tumoredges=0;gapcount=0
 short={'Non-tumor stromal region':'stroma','Immune-rich region':'immune-rich tissue','Normal urothelium':'normal epithelium','Smooth-muscle-rich region':'smooth-muscle-rich tissue','Other tissue':'other'}
 for i in ti:
  for j in adj[i]:
   if t[j]:continue
   tumoredges+=1
   if not valid[j]:gapcount+=1;continue
   common=set(map(tuple,keys[i]))&set(map(tuple,keys[j]));assert len(common)==2
   kk=np.array(sorted(common));seg=np.c_[kk[:,0]/6,kk[:,1]/(2*np.sqrt(3)),np.ones(2)]@aff;segments.append(seg);bt='tumor–'+short.get(d.tissue_region.iloc[j],'other');edges.append((i,j,int(tid[i]),int(d.tissue_component.iloc[i]),bt));inner.add(i);outer.add(j);types[i].add(bt);types[j].add(bt)
 segments=np.asarray(segments).reshape(-1,2,2);dist=np.full(n,np.nan);outside_region=np.full(n,-1)
 for rid in range(reg):
  ei=[k for k,e in enumerate(edges) if e[2]==rid];ix=np.flatnonzero(tid==rid)
  if ei:dist[ix]=segment_distance(xy[ix],segments[ei])
 for tc in d.tissue_component.unique():
  ei=[k for k,e in enumerate(edges) if e[3]==tc];ix=np.flatnonzero((d.tissue_component.to_numpy()==tc)&~t&valid)
  if ei and len(ix):dist[ix]=segment_distance(xy[ix],segments[ei])
 ih=hops(adj,inner,t);oh=hops(adj,outer,~t&valid);band=np.full(n,'unassessable',dtype=object);band[t&np.isfinite(ih)]='tumor_interior';band[t&(ih==1)]='tumor_inner_margin';band[~t&np.isfinite(oh)]='distal_non_tumor';band[~t&(oh==1)]='peritumoral_outer_interface';b2=band.copy();b2[t&(ih<=2)]='tumor_inner_margin';b2[~t&(oh<=2)]='peritumoral_outer_interface'
 d['tumor_region_id']=[f'{sec}:T{x+1}' if x>=0 else '' for x in tid];d['boundary_side']=band;d['boundary_side_two_rings']=b2;d['boundary_type']=[';'.join(sorted(types[i])) if i in types else '' for i in range(n)];d['distance_to_boundary']=dist;d['signed_distance']=np.where(t,-dist,dist);d['distance_unit']='aligned_image_pixel';d['signed_distance_dNN']=d.signed_distance/dnn;d['boundary_graph_ring']=np.where(t,ih,oh);d.loc[~np.isfinite(d.boundary_graph_ring),'boundary_graph_ring']=np.nan;d['display_region']=d.tissue_region+' | '+d.boundary_side
 d.to_csv(dest/'spot_region_labels.tsv.gz',sep='\t',index=False);d.to_parquet(dest/'spot_region_labels.parquet',index=False)
 ds=d.groupby('BANKSY_domain').agg(n_spots=('barcode','size'),malignant_median=('malignant_fraction','median'),tumor_fraction=('main_tumor_mask','mean'),dominant_region=('tissue_region',lambda z:z.value_counts().index[0]),immune_fraction=('c2l_immune','median'),stromal_fraction=('c2l_stroma','median'),epithelial_fraction=('c2l_Epithelial','median'),urothelial_marker=('marker_urothelial','median'),muscle_marker=('marker_smooth_muscle','median'));ds.to_csv(dest/'domain_annotation_summary.tsv',sep='\t')
 pd.DataFrame([dict(tumor_barcode=d.barcode.iloc[e[0]],outer_barcode=d.barcode.iloc[e[1]],tumor_region_id=f'{sec}:T{e[2]+1}',boundary_type=e[4],x1=seg[0,0],y1=seg[0,1],x2=seg[1,0],y2=seg[1,1]) for e,seg in zip(edges,segments)]).to_csv(dest/'assessable_boundary_edges.tsv.gz',sep='\t',index=False)
 old=d.epi_fraction_gate.ge(.25).to_numpy();summ={'section_id':sec,'dataset_id':r.dataset_id,'n_spots':n,'n_algorithm_retained':int((~excluded).sum()),'n_tumor':int(t.sum()),'n_tumor_components':reg,'n_inner':len(inner),'n_outer':len(outer),'n_interior':int((band=='tumor_interior').sum()),'n_distal':int((band=='distal_non_tumor').sum()),'n_assessable_edges':len(edges),'n_uncertain':int((~valid).sum()),'uncertain_fraction':float((~valid).mean()),'tumor_below_old_epi_gate':int((t&~old).sum()),'old_epi_not_tumor':int((old&~t).sum()),'mask_jaccard':float((old&t).sum()/max(1,(old|t).sum())),'signature_source':str(d.SpaCET_signature_source.iloc[0]),'n_normal_urothelium':int(normal.sum()),'affine_max_residual_dNN':float(np.linalg.norm(fit-xy,axis=1).max()/dnn),'dNN_image_pixels':dnn,'boundary_coverage_among_observed_non_tumor_neighbors':len(edges)/max(1,tumoredges),'n_unknown_adjacent_edges':gapcount,'n_unobserved_neighbor_sides_at_tumor':int(sum(6-len(adj[i]) for i in ti)),'boundary_analysis_status':'eligible_computational' if len(inner)>=20 and len(outer)>=20 and (band=='tumor_interior').sum()>=20 else ('limited_no_core_or_small_interface' if len(edges) else 'no_assessable_tumor_interface')}
 summ.update({'count_'+k:int((d.tissue_region==k).sum()) for k in PAL if k!='Algorithm QC excluded'});(dest/'section_summary.json').write_text(json.dumps(summ,indent=2));plot(sec,d,polys,segments,r,summ);(dest/'annotation.done').write_text('complete')

def plot(sec,d,polys,segs,r,summ):
 dest=Path(r.output_dir);im=np.asarray(Image.open(r.image_path));h,w=im.shape[:2];xy=d[['x','y']].to_numpy();pad=.035*max(np.ptp(xy[:,0]),np.ptp(xy[:,1]));limits=(max(0,xy[:,0].min()-pad),min(w,xy[:,0].max()+pad),min(h,xy[:,1].max()+pad),max(0,xy[:,1].min()-pad))
 def setup(ax,title):ax.set_title(title,loc='left');ax.set_aspect('equal');ax.set_xlim(limits[:2]);ax.set_ylim(limits[2:]);ax.set_axis_off()
 def cats(ax,values,pal,alpha=1,legend=True):
  vals=pd.Series(values).fillna('unavailable');ax.add_collection(PolyCollection(polys,facecolors=[pal.get(v,'#999999') for v in vals],edgecolors='none',alpha=alpha))
  if legend:ax.legend(handles=[Patch(color=pal[k],label=k) for k in pal if k in set(vals)],loc='upper left',bbox_to_anchor=(1.0,1),fontsize=6,frameon=False,handlelength=.8)
 def cont(ax,values,label,vmin=None,vmax=None,cmap='viridis'):
  a=np.asarray(values,float);c=PolyCollection(polys,array=np.ma.masked_invalid(a),cmap=cmap,edgecolors='none');c.cmap=c.cmap.copy();c.cmap.set_bad('#a6a6a6');c.set_clim(vmin,vmax);ax.add_collection(c);cb=plt.colorbar(c,ax=ax,fraction=.035,pad=.01);cb.set_label(label,fontsize=7)
 fig,axs=plt.subplots(2,3,figsize=(18,10));axs=axs.ravel()
 setup(axs[0],'A  H&E + retained tissue locations');axs[0].imshow(im);axs[0].scatter(xy[:,0],xy[:,1],s=.3,c='cyan',alpha=.35)
 setup(axs[1],'B  H&E + computational tissue regions');axs[1].imshow(im);cats(axs[1],d.tissue_region,PAL,.45,False)
 if len(segs):axs[1].add_collection(LineCollection(segs,colors='#a31330',linewidths=.45))
 setup(axs[2],'C  Tissue identity (not pathology-confirmed)');cats(axs[2],d.tissue_region,PAL)
 setup(axs[3],'D  SpaCET BLCA malignant fraction');cont(axs[3],d.malignant_fraction,'estimated component proportion',0,1);axs[3].text(0,-.03,'Signature: '+str(d.SpaCET_signature_source.iloc[0]),transform=axs[3].transAxes,fontsize=8)
 setup(axs[4],'E  CNA pattern support / label discordance');cats(axs[4],d.CNV_support,CPAL)
 setup(axs[5],'F  Main mask: first-ring interface');cats(axs[5],d.boundary_side,BPAL)
 fig.suptitle(f'{sec}  |  tumor {summ["n_tumor"]:,}/{len(d):,}; uncertain {summ["uncertain_fraction"]:.1%}',fontsize=13);fig.subplots_adjust(wspace=.62,hspace=.2,right=.88,top=.92);fig.savefig(dest/'01_core_regions.png',dpi=300);fig.savefig(dest/'01_core_regions.pdf');plt.close(fig)
 fig,axs=plt.subplots(1,3,figsize=(17,5));setup(axs[0],'Signed distance to evaluable main-mask interface');lim=np.nanquantile(np.abs(d.signed_distance_dNN),.95) if d.signed_distance_dNN.notna().any() else 1;cont(axs[0],d.signed_distance_dNN,'distance / median neighbor spacing; inside < 0',-max(lim,1),max(lim,1),'RdBu_r');setup(axs[1],'SpaCET native 0.5 classification');cats(axs[1],d.SpaCET_native_interface,{'Tumor':'#bd3345','Interface':'#2ca7bd','Stroma':'#d0d9ce','unavailable':'#999999'});setup(axs[2],'BANKSY spatial domains (unannotated IDs)');ids=sorted(d.BANKSY_domain.unique());pal={v:matplotlib.colors.to_hex(plt.get_cmap('tab20')(i%20)) for i,v in enumerate(ids)};cats(axs[2],d.BANKSY_domain,pal,legend=len(ids)<20);fig.suptitle(sec);fig.subplots_adjust(wspace=.5,right=.9);fig.savefig(dest/'02_distance_native_domains.png',dpi=300);fig.savefig(dest/'02_distance_native_domains.pdf');plt.close(fig)
 z=np.load(dest/'cnv_plot_profiles.npz');a=z['profiles'];ch=z['chrom'];lim=max(.05,float(np.quantile(np.abs(a),.98)));fig,axs=plt.subplots(2,1,figsize=(13,6),gridspec_kw={'height_ratios':[3,1]});v=axs[0].imshow(a.T,aspect='auto',cmap='RdBu_r',vmin=-lim,vmax=lim,rasterized=True);axs[0].set_ylabel('Spots ordered by CNA cluster\n(deterministic sample up to 600)');axs[1].imshow(z['consensus'].T,aspect='auto',cmap='RdBu_r',vmin=-lim,vmax=lim,rasterized=True);axs[1].set_ylabel('CNA cluster median');axs[1].set_yticks(np.arange(len(z['cluster_ids'])),z['cluster_ids']);ticks=[np.flatnonzero(ch==c).mean() for c in sorted(set(ch))];labs=[str(c) if c<23 else 'X' for c in sorted(set(ch))]
 for ax in axs:
  ax.set_xticks(ticks,labs);ax.set_xlabel('Chromosome (native ordered CNA profile)')
  for p in np.flatnonzero(np.diff(ch)!=0):ax.axvline(p+.5,color='black',lw=.3)
 fig.colorbar(v,ax=axs,label='native CopyKAT CNA value; no artificial rescaling',fraction=.02);fig.suptitle(sec+' | CNA pattern evidence; inferred CNAs are not malignant fractions');fig.savefig(dest/'03_CNA_chromosome_profiles.png',dpi=300);fig.savefig(dest/'03_CNA_chromosome_profiles.pdf');plt.close(fig)
 fig,ax=plt.subplots(figsize=(8,6));setup(ax,sec);cats(ax,d.tissue_region,PAL);fig.subplots_adjust(right=.66);fig.savefig(dest/'04_final_tissue_regions.png',dpi=300);fig.savefig(dest/'04_final_tissue_regions.pdf');plt.close(fig)
if __name__=='__main__':
 secs=[a for a in sys.argv[1:] if not a.startswith('--')] or list(man.index)
 for sec in secs:
  try:annotate(sec,force='--force' in sys.argv)
  except Exception:
   import traceback;print(sec,traceback.format_exc(),flush=True)
