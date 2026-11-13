import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
import pathlib,sys,pandas as pd
from pyscenic.rss import regulon_specificity_scores
B=pathlib.Path('results/mp10_immune_endothelial_mechanisms_v1');lin=sys.argv[1];parts=pd.read_csv(B/f'tables/scenic_{lin}_aucell_parts.tsv',sep='\t');a=pd.concat([pd.read_csv(x,sep='\t',index_col=0) for x in parts.file]);m=pd.read_csv(B/f'tables/scenic_membership_{lin}.tsv.gz',sep='\t',index_col='barcode');m=m.loc[a.index]
if lin=='Epi':
 old=pd.read_csv('mp10_subtype_screen_round1/tables/epithelial_existing_metadata.tsv.gz',sep='\t',usecols=['cell_id','epi_allgene_cluster']).set_index('cell_id');labels=old.reindex(a.index).epi_allgene_cluster.fillna('Unavailable existing cluster').astype(str)
else:labels=m.analysis_celltype_c2l_final
rss=regulon_specificity_scores(a,labels);rss.index.name='subgroup';rss.reset_index().melt('subgroup',var_name='regulon',value_name='RSS').assign(lineage=lin).to_csv(B/f'tables/scenic_{lin}_RSS.tsv.gz',sep='\t',index=False);print('RSS COMPLETE',lin)
