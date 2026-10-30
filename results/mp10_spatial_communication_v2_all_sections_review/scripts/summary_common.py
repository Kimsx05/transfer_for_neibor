from pathlib import Path
import os,sys,json,re
P=Path('/home/data/t070721/codex_workspace/Bladder Metabolism');V=P/'mp10_spatial_communication_v1';O=P/'mp10_spatial_communication_v2_all_sections'
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
os.environ['MPLBACKEND']='Agg';os.environ['MPLCONFIGDIR']=str(O/'logs/mplconfig')
sys.path=[p for p in sys.path if 'dist-packages' not in p]
sys.path[:0]=[str(V/'runtime'),str(P/'spatial_ecology_top20_v2/runtime'),str(P/'cell2location_setup_v1/cell2location_env/lib/python3.12/site-packages')]
import numpy as np,pandas as pd
EPI=['Epi_MP10_high','Epi_MP10_low','Epi_MP10_other'];LINEAGES=['B','Plasma','T/NK','Myeloid','Endothelial','Fibroblast','Mural','Mast']
def read(path):return pd.read_csv(path,sep='\t',low_memory=False)
def save(df,path):Path(path).parent.mkdir(parents=True,exist_ok=True);pd.DataFrame(df).to_csv(path,sep='\t',index=False,na_rep='NA')
def safe(x):return re.sub(r'[^a-zA-Z0-9_-]+','_',str(x))
def truth(x):return str(x).strip().lower() in ['true','1','yes'] if not isinstance(x,(bool,np.bool_)) else bool(x)
def unknown(x):return 'unknown' if pd.isna(x) or str(x).strip() in ['', 'None','nan','NA'] else str(x)
def atom_json(obj,path):
 path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=str)+'\n');tmp.replace(path)
def job_dirs(method,section):
 folder='commot' if method=='COMMOT' else 'spatialcellchat';variant='primary2dNN' if method=='COMMOT' else 'main';out=O/folder/section/variant;source=out
 for ref in [out/'CACHE_REFERENCE.json',out.parent/'CACHE_REFERENCE.json']:
  if ref.exists():
   j=json.loads(ref.read_text());p=next((j[k] for k in ['source_dir','cache_path','source_path','original_path','path'] if k in j),None)
   if p is None:raise ValueError(f'Unrecognized cache reference schema {ref}: {list(j)}')
   source=Path(p);break
 return out,source,variant

def section_info(section):
 inp=O/'input_review/sections'/section
 info=json.loads((inp/'input_summary.json').read_text());meta=read(inp/'meta.tsv')
 # Parent input agent may attach updated section metadata beside read-only reused inputs.
 for fn in ['section_manifest.tsv','section_status.tsv','section_index.tsv','pilot_sections.tsv']:
  f=O/'input_review'/fn
  if f.exists():
   tab=read(f)
   if 'section_id' in tab:
    sub=tab[tab.section_id==section]
    if len(sub):info.update({k:v for k,v in sub.iloc[0].items() if pd.notna(v)})
 for k in ['dataset','patient_uid','patient_id_status','analysis_unit_uid','analysis_unit_type','tissue_condition','treatment_status','score_coverage_status','hl_status']:
  if k not in info and k in meta:info[k]=meta[k].iloc[0]
 for k in ['tissue_condition','treatment_status']:info[k]=unknown(info.get(k))
 info['score_coverage_status']=info.get('score_coverage_status','complete20' if info.get('n_mp10_genes_measured',0)==20 else 'score_coverage_incomplete')
 def table(name):
  local=inp/name;globalp=O/'input_review'/name
  if local.exists():return read(local)
  tab=read(globalp);return tab[tab.section_id==section].copy() if 'section_id' in tab else tab
 support=table('type_support.tsv');crosswalk=table('spatial_label_crosswalk.tsv')
 return inp,info,meta,support,crosswalk
