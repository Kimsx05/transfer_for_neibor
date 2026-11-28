from pathlib import Path
import os,json,hashlib,itertools,sys
import numpy as np,pandas as pd
from scipy import stats,sparse
from scipy.spatial.distance import cdist
import pyarrow as pa
pa.set_cpu_count(1);pa.set_io_thread_count(1)
O=Path(__file__).resolve().parents[1];V=O.parent/'mp10_visium_niche_boundary_v2'
OUT=Path(os.environ.get('V3_MODEL_OUT',O)).resolve()
for f in ['classification','tables','figures/groups','review','logs','audit']:(OUT/f).mkdir(parents=True,exist_ok=True)
saved_parameters=json.loads((OUT/'classification_lock.json').read_text())['parameters'] if (OUT/'classification_lock.json').exists() else {}
GRID=os.environ.get('V3_GRID',saved_parameters.get('grid','500'));AXIS=os.environ.get('V3_AXIS',saved_parameters.get('axis','T'))
if GRID not in ['500','1000','x250','y250','xy250'] or AXIS not in ['T','I']:raise ValueError('Unknown grid or immune axis')
def read(p):return pd.read_csv(p,sep='\t')
def save(d,p):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);pd.DataFrame(d).to_csv(p,sep='\t',index=False,na_rep='NA',float_format='%.12g')
def sha(p):
    with open(p,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
M=read(O/'cache/immune_metadata.tsv').set_index('section_id')
M['tumor']=~M.source_context.eq('adjacent_normal_source')
M['unit']=M.patient_uid.fillna(M.analysis_unit_uid)
B=['all_E','E_without_I_enriched'];G=['C1','C2','C3'];seed=int(os.environ.get('V3_SEED','20261127'));area=.008660254037844386
def spot(s,cols=None):return pd.read_parquet(V/'spots'/(s+'.parquet'),columns=cols,use_threads=False)
def bh(p,n=None):
    a=np.asarray(p,float);k=np.isfinite(a);o=np.full(len(a),np.nan)
    if k.any():o[k]=np.minimum(1,stats.false_discovery_control(a[k])*(n or len(a))/k.sum())
    return o
def graph(d,valid):
    loc={(r,c,b):i for i,(r,c,b) in enumerate(zip(d.array_row,d.array_col,d.block)) if valid[i]}
    edges=[(i,j) for (r,c,b),i in loc.items() for dr,dc in [(0,2),(1,1),(1,-1),(0,-2),(-1,-1),(-1,1)] if (j:=loc.get((r+dr,c+dc,b))) is not None]
    i,j=np.array(edges,dtype=int).reshape(-1,2).T
    return sparse.csr_matrix((np.ones(len(i)),(i,j)),shape=(len(d),len(d)))
def divide(x,y):return x/y if y>0 else np.nan
def theme():
    import matplotlib as mpl
    from matplotlib import font_manager as fm
    import subprocess
    path=subprocess.check_output(['fc-match','-f','%{file}','Arial'],text=True)
    if Path(path).is_file():fm.fontManager.addfont(path)
    mpl.rcParams.update({'font.family':fm.FontProperties(fname=path).get_name(),'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.facecolor':'white'})

