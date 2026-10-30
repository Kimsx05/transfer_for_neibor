import sys
from pathlib import Path
R=Path(__file__).resolve().parents[3];sys.path.insert(0,str(R/'results/mp10_pathology_regions_v3/scripts'))
from run_v3 import *
O=R/'results/mp10_pathology_regions_v3_2'
rows=[]
for sec in M.index:
 d=pd.read_parquet(R/'results/mp10_pathology_regions_v3_1/by_dataset'/M.loc[sec,'dataset_id']/sec/'spot_regions_v3_1.parquet')
 for dom,z in d.groupby('domain_component_id'):
  rows.append(dict(section=sec,dom=dom,n=len(z),old_inside=z.lesion_location.eq('inside').mean(),mal=z.malignant_fraction.mean(),epi=z.marker_epithelial.mean(),epidet=z.detected_epithelial.mean(),CNV=z.CNV_support.eq('strong').mean(),SMC=z.c2l_SMC.mean(),stroma=z.c2l_stroma.mean(),immune=z.c2l_immune.mean(),uro=z.marker_urothelial.mean(),smooth=z.marker_smooth_muscle.mean(),musdet=z.detected_smooth_muscle.mean(),score=z.domain_evidence_score.mean(),tissue=z.v2_tissue_region.mode().iloc[0],source=';'.join(z.SpaCET_signature_source.unique().astype(str))))
q=pd.DataFrame(rows);q.to_csv(O/'tables/input_domain_summary.tsv',sep='\t',index=False)
print(q[q.section.isin(['GSE319536_BC18','GSE285715_T3','GSE319536_BC01','GSE285715_N1']) & q.n.ge(10)].round(3).to_string(index=False))
print('columns',list(d.columns))
