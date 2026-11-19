from pathlib import Path
import pandas as pd
from scipy.sparse import load_npz

def load_result(section, counts=False, root=None):
    if section not in ('Cervilla','D1_Johnson'):
        raise ValueError(section)
    folder=(Path(root) if root is not None else Path(__file__).resolve().parents[1])/section
    result={k:pd.read_parquet(folder/(k+'.parquet')) for k in ('cells_scored','candidate_ledger','neighborhoods')}
    if counts:
        result['counts']=load_npz(folder/'counts.npz')
        result['genes']=pd.read_csv(folder/'genes.tsv',sep='\t')
        assert result['counts'].shape==(len(result['cells_scored']),len(result['genes']))
    return result
