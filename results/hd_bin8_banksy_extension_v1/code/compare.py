from pathlib import Path
import pandas as pd,numpy as np,json
r=Path(__file__).resolve().parents[1];old=r.parent/'hd_bin8_banksy_v1';base=r.parent/'hd_bin8_comparison_v1';tabs=r/'tables';tabs.mkdir(exist_ok=True)
sections=['Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2'];inputs=[];info=[];dis=[];conf=[];markers=[];states=[]
for s in sections:
 o=r/s;new=s not in ['Cervilla','D1_Johnson'];j=o/'input_summary.json'
 if new:
  inp=json.loads(j.read_text()) if j.exists() else {};states.append(dict(section=s,input_complete=bool(inp),normalization_complete=(o/'normalized.rds').exists(),annotation_complete=(o/'annotation_summary.json').exists(),lambda0_complete=(o/'l00_complete.json').exists(),lambda02_complete=(o/'l02_complete.json').exists(),lambda08_complete=(o/'l08_complete.json').exists(),validation_complete=(o/'validation.json').exists(),visual_review_complete=(o/'visual_review.md').exists(),tissue_mask='original' if s=='UC_SCC' else 'unknown'))
  if inp:inputs.append(dict(section=s,input_bins=inp['original_bins'],input_resolution_um=inp['input_resolution_um'],bins8=inp['bins8'],input_UMI=inp['input_UMI'],analysis_bins=inp['QC_bins'],analysis_UMI=inp['QC_UMI'],excluded_UMI=inp['excluded_UMI'],tissue_mask=inp['in_tissue_known'],denominator='supplied in-tissue nCount>5' if s=='UC_SCC' else 'supplied nCount>5; tissue unknown'))
 else:
  v=json.loads((old/s/'validation.json').read_text());inputs.append(dict(section=s,input_bins=np.nan,input_resolution_um=8,bins8=np.nan,input_UMI=np.nan,analysis_bins=v['n'],analysis_UMI=v['UMI'],excluded_UMI=np.nan,tissue_mask=True,denominator='supplied in-tissue nCount>5'));states.append(dict(section=s,input_complete=True,normalization_complete=True,annotation_complete=True,lambda0_complete=True,lambda02_complete=True,lambda08_complete=True,validation_complete=(o/'validation.json').exists(),visual_review_complete=(o/'visual_review.md').exists(),tissue_mask='original'))
 p=o/'comparison_metadata.tsv.gz'
 if not p.exists():continue
 d=pd.read_csv(p,sep='\t');assert d.bin_id.is_unique
 b=d.singler_blca_pruned_labels;h=d.singler_hpca_pruned_labels;both=b.notna()&h.notna();fine=b.fillna('Unknown')!=h.fillna('Unknown');coarse=d.coarse_blca!=d.coarse_hpca
 d['conflict_group']=np.select([(b=='Fibroblast')&(h=='Epithelial'),(b=='Epithelial')&(h=='Fibroblast'),(b=='Epithelial')&(h=='Epithelial'),(d.coarse_blca=='Mesenchymal')&(d.coarse_hpca=='Mesenchymal'),~both],['BLCA_Fib_HPCA_Epi','BLCA_Epi_HPCA_Fib','Agree_Epithelial','Agree_Mesenchymal','Not_both_evaluable'],default='Other')
 for band in ['All','<20','20-99','100-499','>=500']:
  q=np.ones(len(d),bool) if band=='All' else (d.information_band==band);z=d[q]
  if len(z)==0:continue
  info.append(dict(section=s,information_band=band,n=len(z),UMI=int(z.nCount_RNA.sum()),fraction=len(z)/len(d),UMI_fraction=z.nCount_RNA.sum()/d.nCount_RNA.sum(),nFeature_median=z.nFeature_RNA.median(),nCount_median=z.nCount_RNA.median(),nFeature_q10=z.nFeature_RNA.quantile(.1),nFeature_q90=z.nFeature_RNA.quantile(.9),old_epi_fib_ge2=z.epi_fib_ge2.mean(),expanded_keratin_collagen_ge2=z.ext_keratin_collagen_ge2.mean(),controls_collagen_ge2=z.ext_controls_collagen_ge2.mean(),keratin_ge2=(z.ext_Keratin_extended_detected>=2).mean(),collagen_ge2=(z.ext_Collagen_stroma_detected>=2).mean()))
  for den,mask in [('all_analysis_bins',q),('both_candidates',q&both)]:
   if not mask.any():continue
   for level,hit in [('original_major',fine),('four_coarse',coarse)]:
    yes=mask&hit;dis.append(dict(section=s,information_band=band,denominator=den,level=level,n=int(mask.sum()),UMI=int(d.loc[mask,'nCount_RNA'].sum()),disagree_n=int(yes.sum()),disagree_fraction=yes.sum()/mask.sum(),disagree_UMI=int(d.loc[yes,'nCount_RNA'].sum()),disagree_UMI_fraction=d.loc[yes,'nCount_RNA'].sum()/d.loc[mask,'nCount_RNA'].sum()))
  for group,zz in z.groupby('conflict_group'):
   conf.append(dict(section=s,information_band=band,group=group,n=len(zz),UMI=int(zz.nCount_RNA.sum()),within_band_fraction=len(zz)/len(z),within_band_UMI_fraction=zz.nCount_RNA.sum()/z.nCount_RNA.sum(),median_nFeature=zz.nFeature_RNA.median(),keratin_ge2=(zz.ext_Keratin_extended_detected>=2).mean(),controls_ge2=(zz.ext_Urothelial_controls_detected>=2).mean(),collagen_ge2=(zz.ext_Collagen_stroma_detected>=2).mean(),expanded_keratin_collagen_ge2=zz.ext_keratin_collagen_ge2.mean(),controls_collagen_ge2=zz.ext_controls_collagen_ge2.mean(),old_epi_fib_ge2=zz.epi_fib_ge2.mean()))
 p=o/'supplement_marker_statistics.tsv'
 if p.exists():t=pd.read_csv(p,sep='\t');t['section']=s;markers.append(t)
for name,rows in [('input_QC_summary',inputs),('completion_status',states),('information_strata',info),('reference_disagreement_strata',dis),('reference_conflict_strata',conf)]:pd.DataFrame(rows).to_csv(tabs/f'{name}.tsv',sep='\t',index=False)
if markers:pd.concat(markers).to_csv(tabs/'expanded_marker_evidence.tsv.gz',sep='\t',index=False)
# Descriptive continuity comes from already-computed original or new models, no refitting.
v=[]
for s in sections:
 p=(old if s in ['Cervilla','D1_Johnson'] else r)/s/'spatial_structure_metrics.tsv'
 if p.exists():t=pd.read_csv(p,sep='\t');t['section']=s;v.append(t)
if v:pd.concat(v).to_csv(tabs/'spatial_structure_comparison.tsv',sep='\t',index=False)
print('COMPARISON_UPDATED',len(info),'information rows; incomplete sections retained in status table')
# One bounded descriptive sensitivity analysis, triggered by observed very different information distributions.
# No bin labels/counts are changed and no model is refit; common information strata use 20-gene width.
main=['Cervilla','D1_Johnson','UC_SCC'];frames={}
if all((r/s/'comparison_metadata.tsv.gz').exists() for s in main):
 for s in main:
  d=pd.read_csv(r/s/'comparison_metadata.tsv.gz',sep='\t');d['feature_stratum20']=np.minimum((d.nFeature_RNA//20).astype(int),50);b=d.singler_blca_pruned_labels;h=d.singler_hpca_pruned_labels;d['both']=b.notna()&h.notna();d['fine_disagree']=b.fillna('Unknown')!=h.fillna('Unknown');d['coarse_disagree']=d.coarse_blca!=d.coarse_hpca;frames[s]=d
 counts=pd.concat({s:d.groupby('feature_stratum20').size() for s,d in frames.items()},axis=1).fillna(0);common=counts.min(axis=1);common=common[common>0];weights=common/common.sum();out=[]
 for s,d in frames.items():
  metrics=['fine_disagree','coarse_disagree','ext_keratin_collagen_ge2','ext_controls_collagen_ge2','epi_fib_ge2']
  for metric in metrics:
   rate=d.groupby('feature_stratum20')[metric].mean().reindex(weights.index);out.append(dict(section=s,metric=metric,standardized_fraction=float((rate*weights).sum()),strata=len(weights),weight_definition='minimum all-analysis stratum size across Cervilla,D1,UC; normalized; nFeature width20, >=1000 pooled',model_refit=False))
 pd.DataFrame(out).to_csv(tabs/'information_standardized_sensitivity.tsv',sep='\t',index=False)
 weights.rename('weight').to_csv(tabs/'information_standardized_weights.tsv',sep='\t')
