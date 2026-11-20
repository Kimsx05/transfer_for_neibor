from pathlib import Path
import pandas as pd
def load_section(section,branch='proseg_observed',include_counts=False):
 from scipy import sparse
 if section not in ('Cervilla','D1_Johnson') or branch not in ('proseg_observed','conservative'):raise ValueError('Unknown section or branch')
 base=Path(__file__).resolve().parents[1]
 if branch=='proseg_observed':base=base/branch
 folder=base/section
 cells=pd.read_parquet(folder/'cells_scored.parquet').set_index('id',verify_integrity=True)
 data={'cells':cells,'matches':pd.read_parquet(folder/'matches.parquet'),'neighborhoods':pd.read_parquet(folder/'neighborhoods.parquet'),'sources':pd.read_parquet(folder/'sources.parquet'),'review':pd.read_parquet(folder/'review_index.parquet')}
 if include_counts:
  X=sparse.load_npz(folder/'counts.npz');assert X.shape[0]==len(cells);data['counts']=X
  data['genes']=pd.read_csv(Path(__file__).resolve().parents[2]/'hd_treg_roi_refinement_v5'/section/'genes.tsv',sep='\t')
 return data
