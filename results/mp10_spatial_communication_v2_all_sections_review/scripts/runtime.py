from pathlib import Path
import sys,os
P=Path('/home/data/t070721/codex_workspace/Bladder Metabolism'); V=P/'mp10_spatial_communication_v1'; O=P/'mp10_spatial_communication_v2_all_sections'
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']: os.environ[k]='1'
os.environ['MPLBACKEND']='Agg';os.environ['MPLCONFIGDIR']=str(O/'logs/mplconfig')
sys.path=[p for p in sys.path if 'dist-packages' not in p]
sys.path[:0]=[str(V/'runtime'),str(P/'spatial_ecology_top20_v2/runtime'),str(P/'cell2location_setup_v1/cell2location_env/lib/python3.12/site-packages')]
