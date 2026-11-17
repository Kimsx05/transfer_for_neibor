import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:4]))
from pathlib import Path
import numpy as np,pandas as pd,json,yaml,h5py
import warnings
warnings.filterwarnings("ignore",category=pd.errors.PerformanceWarning)
from scipy import sparse,ndimage
from scipy.spatial import cKDTree
from scipy.stats import rankdata,spearmanr
from numba import njit
from datetime import datetime
from zoneinfo import ZoneInfo
import time
R=Path(__file__).resolve().parents[1];H=R.parents[1];P=H.parent;A=H/'analysis';OLD=A/'hd_mp10_treg_endo_precommunication_v1'
SECTIONS=['Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2']
TOP=(P/'neighborhood_mp10_calibration_v5/01_score_fidelity/top20.txt').read_text().split()
PANELS={'T':['CD3D','CD3E','CD3G','TRAC','CD2'],'Treg_candidates':['FOXP3','IL2RA','CTLA4','IKZF2','RTKN2','LRRC32','CCR8','TNFRSF18','TNFRSF4','ICOS','LAYN','MAGEH1'],'Epi':['EPCAM','KRT7','KRT19','UPK2','KRT5','KRT14','KRT6A','KRT6B','KRT17','COL17A1'],'Fib':['COL1A1','COL1A2','DCN','LUM'],'Myeloid':['LYZ','LST1','TYROBP'],'Endo':['PECAM1','VWF','CDH5'],'Tfh':['CXCL13','PDCD1','BCL6','MAF','IL21'],'Activated_T':['CD69','IL2','FOS','JUN','IFNG','TNF'],'CD8':['CD8A','CD8B','GZMK','CCL5'],'NK':['NKG7','GNLY','KLRD1']}
COMM=['TIGIT','NECTIN2','CD96','NECTIN1','PVR','CD226']
def gate():
 if os.environ.get('HD_ALLOW_DAYTIME')=='1':return
 while 7<=datetime.now(ZoneInfo('Asia/Shanghai')).hour<23:time.sleep(30)
def put(d,name):d.to_csv(R/'tables'/name,sep='\t',index=False,na_rep='NA')
def load8(s):
 if s in SECTIONS[:2]:p=H/json.loads((A/'hd_bin8_comparison_v1'/s/'input_summary.json').read_text())['matrix_path']
 elif s=='UC_SCC':p=H/'data/GSE299573/GSM9040999/square_008um/filtered_feature_bc_matrix.h5'
 else:p=A/'hd_bin8_banksy_extension_v1'/s/'counts8.h5'
 with h5py.File(p) as f:
  g=f['matrix'];x=sparse.csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=tuple(g['shape'][:])).T.tocsr();ids=g['barcodes'][:].astype(str);ge=g['features/id'][:].astype(str);sy=g['features/name'][:].astype(str)
 d=pd.read_parquet(A/'hd_epi_niche_dual_methods_v1'/s/'nodes.parquet').set_index('bin_id').loc[ids].reset_index()
 assert not d.bin_id.duplicated().any() and x.shape==(len(d),18085) and len(set(ge))==18085
 assert np.array_equal(np.asarray(x.sum(1)).ravel(),d.nCount_RNA)
 return x,d,ge,sy,p
@njit
def selected_ranks(indptr,indices,data,n_genes,select):
 out=np.empty((len(indptr)-1,len(select)),np.float32);lookup=np.full(n_genes,-1,np.int32)
 for k in range(len(select)):lookup[select[k]]=k
 for j in range(len(indptr)-1):
  lo=indptr[j];hi=indptr[j+1];vals=data[lo:hi];ix=indices[lo:hi];nz=len(vals);out[j,:]=(nz+1+n_genes)/2.
  order=np.argsort(-vals.astype(np.float64));a=0
  while a<nz:
   b=a+1
   while b<nz and vals[order[b]]==vals[order[a]]:b+=1
   rank=(a+1+b)/2.
   for z in range(a,b):
    k=lookup[ix[order[z]]]
    if k>=0:out[j,k]=rank
   a=b
 return out

def ucell(ranks,genes,maxrank):
 z=np.minimum(ranks[genes].to_numpy(dtype=float),maxrank);n=len(genes)
 return np.maximum(0,1-(z.sum(1)-n*(n+1)/2)/(n*maxrank-n*(n+1)/2))
def mp10(ranks):return np.maximum(905-ranks[TOP].to_numpy(dtype=float),0).sum(1)/sum(905-np.arange(1,21))
@njit
def rootof(parent,i):
 while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
 return i
@njit
def union_grid(grid,rows,cols,groups):
 p=np.arange(len(rows))
 for i in range(len(rows)):
  for dr,dc in [(0,-1),(-1,-1),(-1,0),(-1,1)]:
   rr=rows[i]+dr;cc=cols[i]+dc
   if rr>=0 and cc>=0 and rr<grid.shape[0] and cc<grid.shape[1]:
    j=grid[rr,cc]
    if j>=0 and groups[j]==groups[i]:
     a=rootof(p,i);b=rootof(p,j)
     if a!=b:p[max(a,b)]=min(a,b)
 for i in range(len(p)):p[i]=rootof(p,i)
 return p

def group_connected(d,keep,side):
 ii=np.flatnonzero(keep);v=d.iloc[ii];rr=v.array_row.to_numpy(int);cc=v.array_col.to_numpy(int)
 group=pd.factorize(pd.MultiIndex.from_arrays([v.component,rr//(side//8),cc//(side//8)]))[0]
 grid=np.full((d.array_row.max()+1,d.array_col.max()+1),-1,np.int32);grid[rr,cc]=np.arange(len(ii));roots=union_grid(grid,rr,cc,group);labels=pd.factorize(roots)[0]
 return ii,labels

def summarize(x,d,sy,kind,side,keep,s):
 ii,lab=group_connected(d,keep,side);ng=lab.max()+1
 agg=sparse.csr_matrix((np.ones(len(ii),np.int64),(lab,ii)),shape=(ng,len(d)));y=(agg@x).tocsr();y.eliminate_zeros()
 gs=np.asarray(y.sum(0)).ravel();expected=np.asarray(x[ii].sum(0)).ravel();assert np.array_equal(gs,expected)
 first=np.unique(lab,return_index=True)[1];base=d.iloc[ii[first]]
 ids=np.array([f'{s}_{kind}{side}_{j:06d}' for j in range(ng)])
 n=np.bincount(lab);v=d.iloc[ii]
 out=pd.DataFrame({'region_id':ids,'section':s,'kind':kind,'grid_um':side,'component':base.component.to_numpy(),'array_row_grid':base.array_row.to_numpy()//(side//8),'array_col_grid':base.array_col.to_numpy()//(side//8),'n_bins':n,'support_area_um2':n*64,'nominal_grid_fraction':n/(side//8)**2,'x_um':np.bincount(lab,weights=v.x_um)/n,'y_um':np.bincount(lab,weights=v.y_um)/n,'nCount':np.asarray(y.sum(1)).ravel(),'nFeature':np.diff(y.indptr),'tissue_status':'original_mask' if s in SECTIONS[:3] else 'unknown'})
 for c in ['reference_conflict','analysis_qc','epi_mixed_flag']:
  out[c+'_fraction']=np.bincount(lab,weights=v[c].fillna(False).to_numpy(float))/n
 for c in ['type_blca','type_hpca']:
  bg=v[c].fillna('Unknown').astype(str)
  for typ in ['Epithelial','Fibroblast','Mural','T/NK','Myeloid','Endothelial','Unknown']:
   out[c+'_'+typ.replace('/','_')+'_fraction']=np.bincount(lab,weights=(bg==typ).to_numpy(float))/n
 genes=sorted(set(sum(PANELS.values(),[])+TOP+COMM+literature_genes()))
 avail=[g for g in genes if g in sy];sel=np.array([np.where(sy==g)[0][0] for g in avail]);assert all(np.sum(sy==g)==1 for g in avail)
 raw=y[:,sel].toarray();ranks=selected_ranks(y.indptr,y.indices,y.data,y.shape[1],sel)
 for j,g in enumerate(avail):out['raw_'+g]=raw[:,j];out['rank_'+g]=ranks[:,j]
 for p,gg in PANELS.items():
  z=raw[:,[avail.index(g) for g in gg if g in avail]];out[p+'_detected']=(z>0).sum(1);out[p+'_UMI']=z.sum(1)
 out['MP10_detected']=(raw[:,[avail.index(g) for g in TOP]]>0).sum(1)
 out.to_parquet(R/s/f'{kind}{side}_measure.parquet',index=False)
 membership=pd.DataFrame({'bin_id':v.bin_id.to_numpy(),'region_id':ids[lab]});membership.to_parquet(R/s/f'{kind}{side}_members.parquet',index=False)
 # Exact counts-only thinning, fixed balanced nFeature subset, never selected by association.
 rng=np.random.default_rng(123);idx=np.unique(np.concatenate([rng.choice(z,min(150,len(z)),False) for _,z in pd.Series(np.arange(ng)).groupby(pd.cut(out.nFeature,[-1,100,300,600,1000,2000,np.inf]),observed=True)]))
 z=y[idx].copy();z.data=rng.binomial(z.data.astype(int),.5);z.eliminate_zeros();rankhalf=selected_ranks(z.indptr,z.indices,z.data,z.shape[1],sel)
 vd=pd.DataFrame(rankhalf,columns=['rank_'+g for g in avail]);vd['region_id']=ids[idx];vd['nFeature_half']=np.diff(z.indptr);vd.to_parquet(R/s/f'{kind}{side}_halfdepth.parquet',index=False)
 # A small full-background input for a native UCell numerical check; all-gene matrices not persisted.
 if s=='Cervilla' and kind=='T' and side==32:
  from scipy.io import mmwrite
  mmwrite(R/'temporary/native_UCell_check.mtx',y[idx[:80]].T);pd.Series(sy).to_csv(R/'temporary/native_UCell_genes.txt',index=False,header=False);pd.Series(ids[idx[:80]]).to_csv(R/'temporary/native_UCell_ids.txt',index=False,header=False)
 return dict(section=s,branch=f'{kind}{side}',regions=int(ng),input_UMI=int(x.sum()),included_UMI=int(y.sum()),excluded_UMI=int(x.sum()-y.sum()),per_gene_conservation=True,raw_integer=True)

def literature_genes():
 p=R/'sources/literature_signature.tsv'
 return pd.read_csv(p,sep='\t').gene.dropna().tolist() if p.exists() else []
