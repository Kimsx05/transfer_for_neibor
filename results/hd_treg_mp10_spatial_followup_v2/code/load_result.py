"""Light regional results. No expression copy or regional-score replication to bins."""
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
SECTIONS=['Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2']
def load_followup(section, kind='regions',root=ROOT):
 root=Path(root)
 if section not in SECTIONS:raise ValueError(section)
 if kind=='regions':return pd.read_parquet(root/'tables'/f'{section}_evidence.parquet')
 if kind=='membership':
  d=pd.read_parquet(root.parent/'hd_treg_mp10_spatial_all_v1'/section/'T32_members.parquet')
  return d.rename(columns={'region_id':'fu2_T32_region_id'})
 if kind=='paired':
  if section!='Cervilla':raise ValueError('Paired old cell result exists only for Cervilla')
  return pd.read_parquet(root/'tables/Cervilla_paired_analysis.parquet')
 if kind=='neighbors':return pd.read_csv(root/'tables/ROI_neighbor_sources.tsv.gz',sep='\t').query('section == @section')
 raise ValueError(kind)
def attach_region_evidence(table,section,id_column='region_id',root=ROOT):
 """Attach to regional table only. No original columns or labels are overwritten."""
 d=load_followup(section,root=root).set_index('region_id');ids=table[id_column]
 if ids.duplicated().any()or not ids.isin(d.index).all():raise ValueError('Unique aligned T32 IDs required')
 fields=['old_evidence','followup_evidence','calibration_status','label_change','score_zero_reason','formal_communication_ready'];out=table.copy()
 for f in fields:
  col='fu2_'+f
  if col in out:raise ValueError('Refusing overwrite: '+col)
  out[col]=d.loc[ids,f].to_numpy()
 return out
if __name__=='__main__':
 import sys
 s,kind,out=sys.argv[1:4];load_followup(s,kind).to_csv(out,sep='\t',index=False,na_rep='NA')
