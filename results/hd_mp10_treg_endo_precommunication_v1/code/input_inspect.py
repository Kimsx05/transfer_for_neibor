import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import anndata as ad,pandas as pd,h5py
P=Path.cwd();H=P/'highres_VisiumHD';V=H/'analysis/hd_epi_niche_dual_methods_v1'
p=P/'subtype/v2/reference_training_v2/11_final_spatial_reference/reference_v2_final_30type_regression_input.h5ad'
a=ad.read_h5ad(p,backed='r');print('REF',a.shape,a.obs.columns.tolist(),a.var.columns.tolist(),list(a.layers));print(a.obs.head().to_string());print(a.var.head().to_string());print(a.obs['spatial_reference_celltype_v2_final'].value_counts());a.file.close()
for s in ['Cervilla','D1_Johnson']:
 d=pd.read_parquet(V/s/'nodes.parquet');print(s,'nodes',d.shape,d.columns.tolist());print(d.head(2).to_string())
 a=ad.read_h5ad(H/'analysis/segmentation_assignment_v2'/s/'cells.h5ad',backed='r');print('cells',a.shape,a.obs.columns.tolist(),a.var.head().to_string(),list(a.obsm));a.file.close()
