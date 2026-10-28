import sys,json
from pathlib import Path
R=Path(__file__).resolve().parents[2];O=R/'mp10_pathology_regions_v2';sys.path.insert(0,str(R/'spatial_ecology_top20_v2/runtime'))
import pandas as pd,numpy as np
rows=[]
for p in (O/'by_dataset').glob('*/*/spot_region_labels.parquet'):
 d=pd.read_parquet(p);cols=[c for c in d if c.startswith('clust_')]
 if len(cols)<2:continue
 mal=d.malignant_fraction.fillna(-1);strong=d.CNV_support.eq('strong');depth=d.SpaCET_signature_source.str.contains('seq_depth',case=False);epi=d.detected_epithelial.ge(3)&d.marker_epithelial.ge(.8);med=d.groupby(cols[0]).malignant_fraction.transform('median');alt=(mal.ge(.5)&~depth)|(mal.ge(.25)&strong&(d.SpaCET_major_label.eq('Malignant')|med.ge(.35)))|(strong&d.native_copykat_label.eq('aneuploid')&epi&med.ge(.15));conflict=alt&d.c2l_immune.ge(.8)&~epi&~strong;conflict|=alt&d.c2l_SMC.ge(.5)&d.detected_smooth_muscle.ge(3)&d.gene_MYH11.gt(.5)&(d.gene_ACTG2.gt(.5)|d.gene_CNN1.gt(.5))&~epi;conflict|=depth&mal.ge(.5)&~strong;conflict|=d.source_context.eq('adjacent_normal_source')&alt&~strong;alt&=~conflict;alt&=d.algorithm_qc.eq('retained');alt&=~d.tissue_region.eq('Normal urothelium');main=d.main_tumor_mask
 rows.append(dict(section_id=d.section_id.iloc[0],n_spots=len(d),main_tumor=int(main.sum()),res04_tumor=int(alt.sum()),n_changed=int((alt!=main).sum()),changed_fraction=float((alt!=main).mean()),tumor_mask_jaccard=float((main&alt).sum()/(main|alt).sum()) if (main|alt).any() else 1.0))
pd.DataFrame(rows).to_csv(O/'tables/pilot_mask_domain_sensitivity.tsv',sep='\t',index=False);print(pd.DataFrame(rows).to_string(index=False))
