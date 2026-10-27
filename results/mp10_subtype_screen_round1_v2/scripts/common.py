from pathlib import Path
import os,sys,json
P=Path('/home/data/t070721/codex_workspace/Bladder Metabolism');O=P/'mp10_subtype_screen_round1_v2';V=P/'mp10_subtype_screen_round1';T=O/'tables'
sys.path=[p for p in sys.path if 'dist-packages' not in p]
sys.path[:0]=[str(P/'spatial_ecology_top20_v2/runtime'),str(P/'cell2location_setup_v1/cell2location_env/lib/python3.12/site-packages')]
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
os.environ['MPLBACKEND']='Agg';os.environ['MPLCONFIGDIR']=str(O/'logs/mplconfig')
import numpy as np,pandas as pd
TOP='HMGCS1 MSMO1 INSIG1 FDFT1 SQLE FDPS IDI1 ACSS2 ELOVL5 LDLR MVD SCD HMGCR ACAT2 SLC20A1 ACSL3 ACSL1 ADIPOR2 MXD1 NPC1'.split()
def save(x,name):pd.DataFrame(x).to_csv(T/name,sep='\t',index=False,na_rep='NA')
def read(name):return pd.read_csv(T/name,sep='\t',low_memory=False)
