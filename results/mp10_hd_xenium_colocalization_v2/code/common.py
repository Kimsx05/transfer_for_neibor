import sys,json,gzip,os,itertools
from pathlib import Path
sys.dont_write_bytecode=True
BASE=Path(__file__).resolve().parents[1]
OUT=Path(os.environ.get('MP_V2_OUTPUT',BASE))
V=BASE.parent/'mp10_hd_xenium_colocalization_v1'
sys.path.insert(1,str(V/'code'))
from base import HD,XE,read_export,tenx,put,np,pd,sparse
pd.options.mode.copy_on_write=True
ROOT=BASE.parents[1]
CFG=json.loads((OUT/'config/frozen.json').read_text())
TOP=CFG['top20']; STATES=['MP10_covered','Chol_biosynthesis_8','H05_transfer']
PAN='CD3D CD3E CD3G TRAC TRBC1 TRBC2'.split()
CORE='FOXP3 IL2RA CTLA4 IKZF2'.split()
MAJ=json.loads((V/'shared/major_markers.json').read_text())
TAX=pd.read_csv(V/'shared/taxonomy.tsv',sep='\t').drop_duplicates('label').set_index('label')
def folders(s):
 k='hd' if s in HD else 'xenium';d=OUT/k/s;d.mkdir(parents=True,exist_ok=True);return V/k/s,d
def expression(s):
 v,d=folders(s)
 if s in HD:return read_export(v/'counts')
 f=pd.read_csv(V/'shared/xenium_inputs.tsv',sep='\t').query('sample==@s and role=="expression"').path.item()
 x,ids,g=tenx(f);b=g.feature_type.eq('Gene Expression').values
 return x[:,b],g.name[b].values,ids
def labels(s):return pd.read_parquet(folders(s)[0]/'labels.parquet')
def raw(x,g,ids,keys):
 if not pd.Index(g).is_unique:
  return pd.DataFrame({k:np.asarray(x[:,g==k].sum(1)).ravel() if np.any(g==k) else np.nan for k in keys},index=ids)
 cols=pd.Index(g).get_indexer(keys);ok=cols>=0
 a=pd.DataFrame(x[:,cols[ok]].toarray(),columns=np.array(keys)[ok],index=ids)
 return a.reindex(columns=keys)
def paths(s):return OUT/'review'/s
def imagedir(s):
 root=ROOT/'highres_VisiumHD/data'
 return root/{'Cervilla':'Cervilla/visiumhd_Bladder/binned_outputs/square_002um/spatial','D1_Johnson':'GSE300692/D1_Johnson/outs/spatial','UC_SCC':'GSE299573/GSM9040999/square_008um/spatial'}[s]
def sf(s):
 p=imagedir(s)/'scalefactors_json.json'
 if not p.exists():p=imagedir(s).parents[1]/'outs/binned_outputs/square_008um/spatial/scalefactors_json.json'
 return json.loads(p.read_text())
