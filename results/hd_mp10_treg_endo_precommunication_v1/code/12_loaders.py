from pathlib import Path
import pandas as pd,json
ROOT=Path(__file__).resolve().parents[1]
def load_precommunication(root=ROOT,include_expression=False):
 root=Path(root);d=pd.read_parquet(root/'tables/Cervilla_cell_formal_metadata.parquet').set_index('cell_id');assert d.index.is_unique
 for r in [50,200]:
  n=pd.read_parquet(root/'tables'/f'Cervilla_cell_neighborhood_r{r}.parquet').set_index('cell_id');keep=['niche','niche_K2','niche_assignment','neighborhood_eligible','neighborhood_reason','known_fraction','n_known','observed_area_fraction','distance_Endothelial_observed_um'];d=d.join(n[keep].add_prefix(f'r{r}_'))
 if not include_expression:return d
 import anndata as ad
 source=root.parent/'segmentation_assignment_v2/Cervilla/cells.h5ad';a=ad.read_h5ad(source,backed='r');z=d.add_prefix('pc_').reindex(a.obs_names);a.obs=a.obs.join(z);a.obs['pc_in_analysis']=a.obs_names.isin(d.index);return a
if __name__=='__main__':
 import sys
 load_precommunication().reset_index().to_csv(sys.argv[1],sep='\t',index=False)
