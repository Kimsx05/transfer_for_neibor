from pathlib import Path
import sys,os,json,re
P=Path('/home/data/t070721/codex_workspace/Bladder Metabolism');V=P/'mp10_subtype_screen_round1_v2';O=P/'mp10_subtype_screen_round3';T=O/'tables'
sys.path=[p for p in sys.path if 'dist-packages' not in p]
sys.path[:0]=[str(P/'spatial_ecology_top20_v2/runtime'),str(P/'cell2location_setup_v1/cell2location_env/lib/python3.12/site-packages')]
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
os.environ['MPLBACKEND']='Agg';os.environ['MPLCONFIGDIR']=str(O/'logs/mplconfig')
import pandas as pd,numpy as np
EPI=['Epi_MP10_low','Epi_MP10_high']
def read(name):return pd.read_csv(T/name,sep='\t',low_memory=False)
def save(x,name):pd.DataFrame(x).to_csv(T/name,sep='\t',index=False,na_rep='NA')
def safe(x):return re.sub(r'[^a-zA-Z0-9_-]+','_',str(x))
