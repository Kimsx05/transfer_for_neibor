import os,sys
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
R=Path(__file__).resolve().parents[3];sys.path.insert(0,str(R/'spatial_ecology_top20_v2/runtime'))
import numpy as np,pandas as pd
from scipy import stats
O=R/'results/mp10_dataset_distance_heterogeneity_v1';S=R/'results/mp10_section_specificity_pearson_v1';OLD=R/'results/mp10_p95_front_programs_v1'
SEED=20261106
M=pd.read_csv(S/'full_sections.tsv',sep='\t');T=M[M.source_context.eq('tumor_study_specimen_not_region_diagnosis')].copy();DS=sorted(T.dataset_id.unique())
META=['section_id','dataset_id','analysis_unit_uid','patient_uid','source_context','n_Epi','n_P95','n_niche_labelled_Epi','fraction_front_exact_Epi']
CORE=['BUDDEN_TG_CATABOLISM_24','BUDDEN_CHOLESTEROL_SYNTHESIS_27','BUDDEN_S1P_RESPONSE_18','BUDDEN_TNF_EMT_RESPONSE_15','HALLMARK_MTORC1_SIGNALING','HALLMARK_GLYCOLYSIS','HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT','HALLMARK_UNFOLDED_PROTEIN_RESPONSE','HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY','PURAM_PEMT_TOP30','FLORES_CRISPR_PRO_INVASION_10']
def save(d,p):
 pd.DataFrame(d).to_csv(p if isinstance(p,Path) else O/'tables'/p,sep='\t',index=False,na_rep='NA',float_format='%.7g')
def read(n):return pd.read_csv(O/'tables'/n,sep='\t',low_memory=False)
def bh(p):
 p=np.asarray(p,float);o=np.flatnonzero(np.isfinite(p));o=o[np.argsort(p[o])];q=np.full(len(p),np.nan);q[o]=np.minimum.accumulate((p[o]*len(o)/np.arange(1,len(o)+1))[::-1])[::-1].clip(0,1);return q
 def_dummy=0
def transform(v,kind):return np.arctanh(np.clip(v,-.999999,.999999)) if kind=='r' else np.asarray(v,float)
def inverse(v,kind):return np.tanh(v) if kind=='r' else np.asarray(v,float)
def combos(ds):
 from itertools import combinations
 return [tuple(c) for k in range(1,len(ds)+1) for c in combinations(sorted(ds),k)]
def combo_id(ds):return 'D'+str(sum(1<<DS.index(d) for d in ds)).zfill(3)
def category(c):
 if c in ['Matrix_CAF','iCAF','rCAF','vCAF']:return 'CAF'
 if c in ['C1QC TAM','SPP1 TAM','Classical Mono','Neutrophil','cDC1','cDC2','pDC','LAMP3 DC','Mast']:return 'Myeloid'
 if 'EC' in c:return 'Endothelial'
 if c in ['SMC','Pericyte']:return 'Mural'
 if c in ['GC B','Naive/Memory B','Plasma']:return 'B_Plasma'
 return 'T_NK'
