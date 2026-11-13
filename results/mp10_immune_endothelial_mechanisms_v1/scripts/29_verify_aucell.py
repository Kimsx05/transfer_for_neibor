import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import pathlib,sys,pandas as pd,json,numpy as np
B=pathlib.Path('results/mp10_immune_endothelial_mechanisms_v1');lin=sys.argv[1];members=pd.read_csv(B/f'tables/scenic_membership_{lin}.tsv.gz',sep='\t').barcode.tolist();checks=[]
for mask in (['Full','M20-masked'] if lin=='Epi' else ['Full']):
 stem=lin if mask=='Full' else lin+'_'+mask;parts=pd.read_csv(B/f'tables/scenic_{stem}_aucell_parts.tsv',sep='\t');ids=[];columns=None
 for f in parts.file:
  p=pathlib.Path(f);z=pd.read_csv(p,sep='\t');assert z.columns[0]=='barcode';assert z.barcode.is_unique;assert np.isfinite(z.iloc[:,1:].values).all();assert (z.iloc[:,1:].values>=0).all();assert (z.iloc[:,1:].values<=1+1e-10).all();ids+=z.barcode.tolist()
  if columns is None:columns=z.columns.tolist()
  else:assert columns==z.columns.tolist()
 assert ids==members,'AUCell cell order differs from frozen full membership'
 checks.append(dict(lineage=lin,mask=mask,n_cells=len(ids),n_regulons=len(columns)-1,exact_membership_order=True,finite_unit_interval=True))
q=pd.read_csv(B/f'tables/source_normalization_audit_{lin}.tsv',sep='\t');assert (q.total_counts==q.nCount_metadata).all() and (q.normalization_max_abs_error<1e-7).all()
(B/f'qc/aucell_{lin}.json').write_text(json.dumps(checks,indent=2));print('FULL AUCELL VERIFIED',lin,checks)
