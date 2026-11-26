from pathlib import Path
import os, sys
O = Path(__file__).resolve().parents[1]
R = O.parents[1]
sys.path.insert(0, str(R/'spatial_ecology_top20_v2/runtime'))
import numpy as np
import pandas as pd
S = {'GSE171351_B1':959, 'GSE238145_B17':1494, 'GSE319536_BC01':1567, 'GSE285715_T1':4167}
TOP = 'HMGCS1 MSMO1 INSIG1 FDFT1 SQLE FDPS IDI1 ACSS2 ELOVL5 LDLR MVD SCD HMGCR ACAT2 SLC20A1 ACSL3 ACSL1 ADIPOR2 MXD1 NPC1'.split()
A = R/'results/mp10_immune_endothelial_mechanisms_v1/tables'
F = R/'results/mp10_p95_front_programs_v1'
K = R/'mp10_k32_reference_identifiability_v5/scope_frozen_resume_v1'
REF = R/'subtype/v2/reference_training_v2/04_gene_harmonization/reference_v2_raw_counts_minimally_harmonized.h5ad'
M = pd.read_csv(F/'tables/input_manifest.tsv',sep='\t').set_index('section_id').loc[list(S)]
def read(p): return pd.read_csv(p,sep='\t')
def save(d,p): pd.DataFrame(d).to_csv(O/p,sep='\t',index=False,na_rep='NA',float_format='%.12g')

