import os,sys,json,zipfile
from pathlib import Path
R=Path(__file__).resolve().parents[2];O=R/'mp10_pathology_regions_v2';sys.path.insert(0,str(R/'spatial_ecology_top20_v2/runtime'))
import pandas as pd,numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from PIL import Image
from annotate_plot import PAL,VERT
from matplotlib.collections import PolyCollection
def adjusted_rand_score(x,y):
 t=pd.crosstab(x,y).to_numpy(dtype=float);n=t.sum();choose=lambda z:(z*(z-1)/2).sum();nij=choose(t);ai=choose(t.sum(1));bj=choose(t.sum(0));tot=n*(n-1)/2;expected=ai*bj/tot if tot else 0;den=(ai+bj)/2-expected;return float((nij-expected)/den) if den else 1.0
m=pd.read_csv(O/'tables/section_manifest.tsv',sep='\t');rows=[];allspots=[];checks=[];sens=[];status=[]
for r in m.itertuples():
 p=Path(r.output_dir);status.append(dict(section_id=r.section_id,SpaCET=(p/'spacet.done').exists(),BANKSY=(p/'banksy.done').exists(),CNA=(p/'cnv_v3.done').exists(),annotation=(p/'annotation.done').exists()))
 if not (p/'annotation.done').exists():continue
 s=json.loads((p/'section_summary.json').read_text());d=pd.read_parquet(p/'spot_region_labels.parquet');assert len(d)==r.n_spots and not d.barcode.duplicated().any();assert d.main_tumor_mask.equals(d.tissue_region.eq('Tumor'));assert (d.loc[d.signed_distance.notna()&d.main_tumor_mask,'signed_distance']<=0).all();assert (d.loc[d.signed_distance.notna()&~d.main_tumor_mask,'signed_distance']>=0).all()
 if r.source_context=='adjacent_normal_source' and s['n_tumor']>0:s['boundary_analysis_status']='review_required_adjacent_source_conflict'
 if s['uncertain_fraction']>.4 and s['n_assessable_edges']>0:s['boundary_analysis_status']='review_required_high_uncertainty'
 edges=pd.read_csv(p/'assessable_boundary_edges.tsv.gz',sep='\t') if s['n_assessable_edges'] else pd.DataFrame()
 if len(edges):
  x=d.set_index('barcode');a=x.loc[edges.tumor_barcode];b=x.loc[edges.outer_barcode];assert a.main_tumor_mask.all() and not b.main_tumor_mask.any();assert not b.tissue_region.eq('Uncertain/conflicting').any();dr=np.abs(a.array_row.to_numpy()-b.array_row.to_numpy());dc=np.abs(a.array_col.to_numpy()-b.array_col.to_numpy());assert (((dr==0)&(dc==2))|((dr==1)&(dc==1))).all()
 checks.append(dict(section_id=r.section_id,n_spots=len(d),true_six_neighbor_edges=len(edges),checks='PASS: keys, mask identity, distance signs, known opposite-side neighbors'))
 cols=[c for c in d if c.startswith('clust_')]
 if len(cols)>1:sens.append(dict(section_id=r.section_id,n_domains_primary=d[cols[-1]].nunique(),n_domains_res04=d[cols[0]].nunique(),adjusted_rand_index=adjusted_rand_score(d[cols[0]],d[cols[-1]])))
 if s['affine_max_residual_dNN']>.15:s['boundary_analysis_status']='coordinate_precision_caution'
 rows.append(s);allspots.append(d)
s=pd.DataFrame(rows);s.to_csv(O/'tables/section_region_summary.tsv',sep='\t',index=False);pd.DataFrame(status).to_csv(O/'tables/completion_status.tsv',sep='\t',index=False);pd.DataFrame(checks).to_csv(O/'tables/definition_checks.tsv',sep='\t',index=False);pd.DataFrame(sens).to_csv(O/'tables/pilot_BANKSY_sensitivity.tsv',sep='\t',index=False)
if not allspots:sys.exit()
a=pd.concat(allspots,ignore_index=True);a.to_parquet(O/'tables/all_spot_region_labels.parquet',index=False)
core=['dataset_id','section_id','barcode','x','y','tissue_region','display_region','tumor_region_id','boundary_side','boundary_type','distance_to_boundary','signed_distance','distance_unit','signed_distance_dNN','malignant_fraction','SpaCET_signature_source','CNV_support','evidence_level','review_status','annotation_conflict','annotation_reason','main_tumor_mask','boundary_side_two_rings','algorithm_qc'];a[core].to_csv(O/'tables/all_spot_region_labels.tsv.gz',sep='\t',index=False)
a.groupby(['dataset_id','section_id','tissue_region','boundary_side']).size().rename('n_spots').to_csv(O/'tables/region_band_counts.tsv',sep='\t');a.groupby(['dataset_id','section_id','native_copykat_label']).size().rename('n_spots').to_csv(O/'tables/native_CopyKAT_label_counts.tsv',sep='\t');a.groupby(['dataset_id','section_id','boundary_side','boundary_side_two_rings']).size().rename('n_spots').to_csv(O/'tables/one_vs_two_ring_counts.tsv',sep='\t')
for ds,ss in s.groupby('dataset_id'):
 n=len(ss);nc=4;nr=int(np.ceil(n/nc));fig,axs=plt.subplots(nr,nc,figsize=(13,3.1*nr),squeeze=False)
 for ax,(_,q) in zip(axs.ravel(),ss.iterrows()):
  z=a[a.section_id.eq(q.section_id)];xy=z[['x','y']].to_numpy();ideal=np.c_[z.array_col/2,z.array_row*np.sqrt(3)/2,np.ones(len(z))];aff=np.linalg.lstsq(ideal,xy,rcond=None)[0];keys=np.c_[3*z.array_col,3*z.array_row][:,None,:]+VERT;poly=np.concatenate([keys[:,:,:1]/6,keys[:,:,1:]/(2*np.sqrt(3)),np.ones((len(z),6,1))],axis=2)@aff;ax.add_collection(PolyCollection(poly,facecolors=z.tissue_region.map(PAL),edgecolors='none'));ax.autoscale_view();ax.set_aspect('equal');ax.invert_yaxis();ax.axis('off');ax.set_title(q.section_id+'\nTumor '+str(q.n_tumor)+'; uncertain '+f'{q.uncertain_fraction:.0%}',fontsize=8)
 for ax in axs.ravel()[n:]:ax.axis('off')
 fig.legend(handles=[Patch(color=v,label=k) for k,v in PAL.items() if k in set(a.tissue_region)],loc='lower center',ncol=4,fontsize=8,frameon=False);fig.suptitle(ds+' | computational tissue regions',fontsize=13);fig.subplots_adjust(bottom=.11 if nr<3 else .07,top=.90,hspace=.35);fig.savefig(O/'figures'/f'{ds}_overview.png',dpi=300);fig.savefig(O/'figures'/f'{ds}_overview.pdf');plt.close(fig)
import subprocess
subprocess.run([sys.executable,str(O/'scripts/write_report.py')],check=True)
# Review bundle intentionally excludes expression, full CNA arrays, package libraries and repeated objects.
with zipfile.ZipFile(O/'review.zip','w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in [O/'METHODS.md',O/'RESULTS.md']:
  if p.exists():z.write(p,p.name)
 for p in (O/'tables').glob('*'):
  if p.suffix in ['.tsv','.json'] or p.name=='all_spot_region_labels.tsv.gz':z.write(p,str(p.relative_to(O)))
 for p in (O/'figures').glob('*.pdf'):z.write(p,str(p.relative_to(O)))
 for r in m.itertuples():
  p=Path(r.output_dir)
  for name in ['01_core_regions.pdf','02_distance_native_domains.pdf','03_CNA_chromosome_profiles.pdf','cnv_cluster_support.tsv','domain_annotation_summary.tsv']:
   f=p/name
   if f.exists():z.write(f,str(f.relative_to(O)))
print(json.dumps({'complete':len(s),'spots':len(a),'boundary_any':int(s.n_assessable_edges.gt(0).sum()),'boundary_eligible':int(s.boundary_analysis_status.eq('eligible_computational').sum()),'tumor_spots':int(a.main_tumor_mask.sum()),'normal_urothelium_spots':int(a.tissue_region.eq('Normal urothelium').sum()),'review_zip_MB':round((O/'review.zip').stat().st_size/1e6,1)},indent=2))
