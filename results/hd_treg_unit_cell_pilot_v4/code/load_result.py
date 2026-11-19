from pathlib import Path
import pandas as pd
from scipy import sparse
def load_unit(section,branch,root=None,with_counts=False):
 root=Path(root) if root else Path(__file__).resolve().parents[1]
 base=root/section
 meta=pd.read_parquet(base/(branch+'_evidence.parquet')).set_index('id',verify_integrity=True)
 members=pd.read_parquet(base/(branch+'_members.parquet'))
 out=dict(metadata=meta,members=members)
 if with_counts:
  out['counts']=sparse.load_npz(base/(branch+'.npz'));out['genes']=pd.read_csv(base/'genes.tsv',sep='\t')
  assert out['counts'].shape==(len(meta),len(out['genes']))
 return out
