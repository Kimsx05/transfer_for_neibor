from pathlib import Path
import pandas as pd,numpy as np,json,sys
r=Path(__file__).resolve().parents[1];s=sys.argv[1];o=r/s;b=r.parent/'hd_bin8_comparison_v1'/s
d=pd.read_csv(b/'analysis_metadata.tsv.gz',sep='\t');t=pd.read_csv(o/'type_metadata.tsv.gz',sep='\t');a=pd.read_csv(o/'spatial_metadata.tsv.gz',sep='\t');assert d.bin_id.equals(t.bin_id) and d.bin_id.equals(a.bin_id)
for c in t:
 if c not in d:d[c]=t[c]
d['native_edge_distance_um']=a.native_edge_distance_um
tr=json.loads((b/'transform.json').read_text());sc=tr['local_to_model_scale_xy'];cr=tr['fullres_crop_xyxy'];mp=tr['native_mpp'];mx=(d.x_um/mp-cr[0]+.5)*sc[0]-.5;my=(d.y_um/mp-cr[1]+.5)*sc[1]-.5
rows=[]
for _,roi in pd.read_csv(o/'fixed_rois.tsv',sep='\t').iterrows():
 mask=(mx>=roi.x)&(mx<roi.X)&(my>=roi.y)&(my<roi.Y)
 for group,z in d[mask].groupby('reference_conflict_group'):
  rr=dict(roi=roi.roi,group=group,n=len(z),UMI=int(z.nCount_RNA.sum()),ROI_n=int(mask.sum()),nFeature_median=z.nFeature_RNA.median(),epi_ge2=(z.detected_Epithelial>=2).mean(),fib_ge2=(z.detected_Fibroblast>=2).mean(),both_ge2=z.epi_fib_ge2.mean(),epi_strength=z.strength_Epithelial.mean(),fib_strength=z.strength_Fibroblast.mean());rows.append(rr)
pd.DataFrame(rows).to_csv(o/'conflict_ROI_summary.tsv',sep='\t',index=False)
rows=[]
for (group,band),z in d.groupby(['reference_conflict_group','banksy_information_band']):
 rows.append(dict(group=group,information_band=band,n=len(z),UMI=int(z.nCount_RNA.sum()),nFeature_median=z.nFeature_RNA.median(),epi_ge2=(z.detected_Epithelial>=2).mean(),fib_ge2=(z.detected_Fibroblast>=2).mean(),both_ge2=z.epi_fib_ge2.mean(),edge16_fraction=(z.native_edge_distance_um<=16).mean(),epi_strength=z.strength_Epithelial.mean(),fib_strength=z.strength_Fibroblast.mean()))
pd.DataFrame(rows).to_csv(o/'early_conflict_summary.tsv',sep='\t',index=False)
print(s);print(pd.DataFrame(rows).query("group in ['BLCA_Fibroblast_HPCA_Epithelial','BLCA_Epithelial_HPCA_Fibroblast']").to_string(index=False));print(pd.DataFrame(rows).query("group=='Agree_Epithelial'").to_string(index=False))
