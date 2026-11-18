"""Load lightweight region results without modifying upstream expression objects."""
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
def load_enrichment(section,root=ROOT):
 root=Path(root);allowed=['Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2']
 if section not in allowed:raise ValueError(section)
 out={k:pd.read_parquet(root/section/f)for k,f in [('T','T_regions.parquet'),('Epi','Epi_regions.parquet'),('neighbors','neighborhoods.parquet')]}
 for key in ['T','Epi']:
  if not out[key].region_id.is_unique:raise ValueError('Duplicated region IDs')
 out['effects']=pd.read_csv(root/'tables/enrichment_effects.tsv',sep='\t').query('section==@section')
 out['dependencies']={'T_members':root.parent/'hd_treg_mp10_spatial_all_v1'/section/'T32_members.parquet','Epi_members':root.parent/'hd_treg_mp10_spatial_all_v1'/section/'Epi64_members.parquet','pairs':root.parent/'hd_treg_mp10_spatial_all_v1'/section/'region_pairs.parquet'}
 return out

def attach_to_regions(metadata,section,kind='T',id_col='region_id',root=ROOT):
 """Only region-level join. Never interpret inherited labels as per-bin/cell identities."""
 d=load_enrichment(section,root)[kind].set_index('region_id');z=metadata.copy();ids=z[id_col]if id_col in z else z.index
 if pd.Index(ids).has_duplicates:raise ValueError('Caller region IDs must be unique')
 for col in d.columns:
  dest='v3_'+col
  if dest in z:raise ValueError('Refusing to overwrite '+dest)
  z[dest]=d[col].reindex(ids).to_numpy()
 return z
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--section',required=True);p.add_argument('--kind',choices=['T','Epi','neighbors'],default='T');p.add_argument('--tsv',required=True);a=p.parse_args();load_enrichment(a.section)[a.kind].to_csv(a.tsv,sep='\t',index=False,na_rep='NA')
