import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import pathlib,csv,pyarrow as pa
pa.set_cpu_count(1);pa.set_io_thread_count(1)
B=pathlib.Path('results/mp10_immune_endothelial_mechanisms_v1');P=pathlib.Path('resources/mp10_regulatory_priors/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather')
with pa.memory_map(str(P),'r') as f:names=set(pa.ipc.open_file(f).schema.names)
rows=[]
for lin in ['Epi','T_lineage','Endothelial']:
 genes=[r['gene'] for r in csv.DictReader(open(B/f'tables/scenic_genes_{lin}.tsv'),delimiter='\t')]
 rows += [dict(lineage=lin,gene=g,in_motif_ranking=g in names,known_previous_symbol='PVRL2' if g=='NECTIN2' else '',previous_symbol_in_ranking=('PVRL2' in names) if g=='NECTIN2' else '',handling='exact frozen symbols; unmapped genes not silently added or scored') for g in genes]
with open(B/'tables/scenic_feature_motif_coverage.tsv','w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys(),delimiter='\t');w.writeheader();w.writerows(rows)
print('MOTIF FEATURE COVERAGE COMPLETE')
