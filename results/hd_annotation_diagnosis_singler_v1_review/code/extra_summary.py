from pathlib import Path
import pandas as pd,numpy as np
root=Path(__file__).resolve().parents[1]
for s in ['Cervilla','D1_Johnson']:
 o=root/s;d=pd.read_csv(o/'diagnosis.tsv.gz',sep='\t',low_memory=False);u=d.final_major.eq('Unknown');rows=[]
 for key,g in d[u].groupby('annotation_reason'):
  rows.append(dict(reason=key,n=len(g),UMI=g.nCount_RNA.sum(),fraction_unknown=len(g)/u.sum(),UMI_fraction_unknown=g.nCount_RNA.sum()/d.loc[u,'nCount_RNA'].sum(),median_nFeature=g.nFeature_RNA.median(),median_nCount=g.nCount_RNA.median()))
 pd.DataFrame(rows).to_csv(o/'original_unknown_reasons.tsv',sep='\t',index=False)
 rows=[]
 for scope,mask in [('original_known_covered',d.final_major.isin(['Epithelial','T/NK','Myeloid','Fibroblast','Endothelial','B'])),('original_known_uncovered',d.final_major.isin(['Mast','Plasma','Mural'])),('original_Unknown_BLCA_candidate_covered',u&d.singler_blca_pruned_labels.isin(['Epithelial','T/NK','Myeloid','Fibroblast','Endothelial','B'])),('original_Unknown_BLCA_candidate_uncovered',u&d.singler_blca_pruned_labels.isin(['Mast','Plasma','Mural']))]:
  for label,g in d[mask].groupby('diag_candidate_pattern'):
   rows.append(dict(scope=scope,pattern=label,n=len(g),UMI=g.nCount_RNA.sum(),fraction_scope=len(g)/mask.sum()))
 pd.DataFrame(rows).to_csv(o/'coverage_stratified_comparison.tsv',sep='\t',index=False)
 rows=[]
 for route in ['blca','hpca']:
  for (status,band),g in d.groupby(['annotation_status','diag_information_band']):
   for metric in ['delta_next','delta_median','selected_initial_score']:
    v=g[f'singler_{route}_{metric}'];rows.append(dict(route=route,original_status=status,nFeature_band=band,metric=metric,evaluated=v.notna().sum(),q10=v.quantile(.1),median=v.median(),q90=v.quantile(.9)))
 pd.DataFrame(rows).to_csv(o/'singler_score_diagnostics.tsv',sep='\t',index=False)
 # Per-candidate marker evidence remains descriptive; no new acceptance or consensus gate.
 rows=[]
 for route in ['blca','hpca']:
  for lab,g in d[u].groupby(f'singler_{route}_pruned_labels',dropna=False):
   rows.append(dict(route=route,candidate=lab,n=len(g),UMI=g.nCount_RNA.sum(),fraction_original_unknown=len(g)/u.sum(),UMI_fraction_original_unknown=g.nCount_RNA.sum()/d.loc[u,'nCount_RNA'].sum(),two_or_more_markers=(g[f'diag_{route}_candidate_marker_count']>=2).sum(),no_candidate_marker=(g[f'diag_{route}_candidate_marker_count']==0).sum(),other_group_at_least_as_strong=g[f'diag_{route}_other_group_supported'].sum(),median_nFeature=g.nFeature_RNA.median()))
 pd.DataFrame(rows).to_csv(o/'unknown_candidate_summary.tsv',sep='\t',index=False)
 print('EXTRA_SUMMARY_COMPLETE',s)
