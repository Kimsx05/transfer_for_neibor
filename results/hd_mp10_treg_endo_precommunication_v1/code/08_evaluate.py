from common import *
from scipy.stats import spearmanr
from datetime import datetime
DIC=pd.read_csv(P/'mp10_sc_c2l_final_annotation_v1/tables/final_label_dictionary.tsv',sep='\t');MAJOR=dict(zip(DIC.celltype,DIC.parent))
d=pd.read_parquet(R/'tables/Cervilla_cell_measured.parquet');p=pd.read_csv(R/'tables/Cervilla_cell_predictions.tsv.gz',sep='\t');d=d.merge(p,on='cell_id',validate='one_to_one');d['fine_major']=d.fine_pruned.map(MAJOR);d['type_supported']='Unknown';d['annotation_reason']='reference_or_marker_insufficient';d['epi_mixed_flag']=(d.Epi_detected>=2)&((d.Fib_detected>=2)|(d.T_detected>=2)|(d.Endo_detected>=2));d['T_crosslineage_mixed']=(d.T_detected>=1)&((d.Epi_detected>=2)|(d.Fib_detected>=2))
lookup={'Epithelial':'Epi','Endothelial':'Endo','Fibroblast':'Fib','Myeloid':'Myeloid','B':'B','Plasma':'Plasma','Mast':'Mast','Mural':'Mural'}
for major,k in lookup.items():
 valid=d.fine_major.eq(major)&(d[k+'_detected']>= (2 if major in ['Epithelial','Endothelial'] else 1));d.loc[valid,'type_supported']=major;d.loc[valid,'annotation_reason']='reference plus unsmoothed lineage-panel support'
tline=d.fine_major.eq('T/NK')&(d.T_detected>=1)&~d.T_crosslineage_mixed;d.loc[tline,'type_supported']='T_subtype_unresolved';d.loc[tline,'annotation_reason']='T lineage supported; fine subtype not yet established';nk=d.fine_pruned.eq('NK')&(d.NK_detected>=2)&(d.T_detected==0);d.loc[nk,'type_supported']='NK'
# Freeze conservative target eligibility from sample-heldout depth validation, not spatial associations.
v=pd.read_csv(R/'tables/annotation_validation.tsv',sep='\t');acceptable=set(v[(v.depth==300)&(v.precision>=.7)&(v.recall>=.4)].type)
for typ in ['Treg','Tfh/CXCL13-like','CD4 Helper/Conventional','CD4 Naive/Memory-like','CD8 Exhausted','CD8 GZMK/Effector-memory-like']:
 if typ not in acceptable:continue
 ok=tline&d.fine_pruned.eq(typ)&(d.nCount_RNA>=300)
 if typ=='Treg':ok&=(d.Treg_detected>=2)&((d.raw_FOXP3>0)|(d.raw_IL2RA>0))
 elif typ.startswith('Tfh'):ok&=(d.Tfh_detected>=2)
 elif typ.startswith('CD8'):ok&=(d.CD8_detected>=1)
 else:ok&=d.CD4_detected>=2
 d.loc[ok,'type_supported']=typ;d.loc[ok,'annotation_reason']='sample-heldout target performance at 300 UMI, no TIGIT; multiple lineage evidence; no detected cross-lineage mixture'
d['annotation_confidence']=np.where(d.type_supported.eq('Unknown'),'unresolved',np.where(d.type_supported.eq('T_subtype_unresolved'),'lineage_only','supported_candidate'))
d['epi_anchor']=d.type_supported.eq('Epithelial');d['mp10_epithelial_eligible']=d.epi_anchor&~d.epi_mixed_flag;d['mp10_single_score_eligible']=d.mp10_epithelial_eligible&(d.MP10_detected>=5)&(d.nFeature_RNA>=905);d['analysis_unit']='cell';d.to_parquet(R/'tables/Cervilla_cell_annotation.parquet',index=False)
put(d.groupby(['type_supported','annotation_confidence'],dropna=False).agg(n=('cell_id','size'),UMI=('nCount_RNA','sum'),median_UMI=('nCount_RNA','median'),median_nFeature=('nFeature_RNA','median')).reset_index(),'cell_annotation_summary.tsv')
rows=[];roi=d[d.roi.ne('')]
rows.append(dict(unit='cell',scope='fixed_ROI_union',n=len(roi),fitted=int(roi.fine_pruned.notna().sum()),Treg_candidates=int(roi.fine_prediction.eq('Treg').sum()),Treg_multigene=int(roi.type_supported.eq('Treg').sum()),Endo_multigene=int(roi.type_supported.eq('Endothelial').sum()),Epi_eligible=int(roi.mp10_epithelial_eligible.sum()),Epi_score_eligible=int(roi.mp10_single_score_eligible.sum()),mixed_Epi_fraction=float(roi.loc[roi.epi_anchor,'epi_mixed_flag'].mean()),median_nFeature=float(roi.nFeature_RNA.median())))
for unit in ['bin8','bin16']:
 q=pd.read_csv(R/'tables'/f'{unit}_ROI.tsv.gz',sep='\t');a=pd.read_csv(R/'tables'/f'{unit}_rctd_results.tsv.gz',sep='\t');w=pd.read_csv(R/'tables'/f'{unit}_rctd_weights.tsv.gz',sep='\t');q=q.merge(a,on='bin_id',how='left',validate='one_to_one');types=sorted(DIC.celltype);W=np.zeros((len(q),len(types)));wi=w.set_index('bin_id').reindex(q.bin_id);wc=wi.columns.tolist()
 # spacexr 2.2.1 weights_doublet columns are first_type and second_type; never mistake them for fixed identities.
 if set(wc)=={'first_type','second_type'}:
  for column in wc:
   labels=q[column].map({t.replace('/','_'):t for t in types})
   for t in types:W[:,types.index(t)]+=np.where(labels.eq(t),wi[column].fillna(0),0)
 else:
  for column in wc:
   label={t.replace('/','_'):t for t in types}.get(column,column)
   if label in types:W[:,types.index(label)]=wi[column].fillna(0)
 reliable=q.spot_class.isin(['singlet','doublet_certain']).to_numpy();tot=W.sum(1);W=np.divide(W,tot[:,None],out=np.zeros_like(W),where=tot[:,None]>0)
 # Installed spacexr get_beta_doublet semantics: singlet is first_type=1; uncertain/reject are not complete proportions.
 pairW=W.copy();W[~reliable]=np.nan
 sing=q.spot_class.eq('singlet').to_numpy();W[sing]=0
 for t in types:W[sing & q.first_type.eq(t.replace('/','_')),types.index(t)]=1
 native=pd.read_csv(R/'tables'/f'{unit}_native_selected_weights.tsv.gz',sep='\t').set_index('bin_id').reindex(q.bin_id)
 for t in types:
  c=t.replace('/','_')
  if c in native:assert np.allclose(W[reliable,types.index(t)],native.loc[reliable,c],atol=1e-12)
 con=pd.DataFrame(W,columns=['contribution_'+x for x in types]);con.insert(0,'bin_id',q.bin_id);con['reliable_fit']=reliable;con.to_parquet(R/'tables'/f'{unit}_ROI_contributions.parquet',index=False)
 tr=(W[:,types.index('Treg')]>=.5)&reliable;en=W[:,[i for i,t in enumerate(types)if MAJOR[t]=='Endothelial']].sum(1);ep=W[:,types.index('Epithelial')];mixed=(q.Epi_detected>=2)&(q.Fib_detected>=2)
 rows.append(dict(unit=unit,scope='fixed_ROI_union',n=len(q),fitted=int(q.spot_class.notna().sum()),reliable_fit=int(reliable.sum()),Treg_candidates=int(tr.sum()),Treg_multigene=int((tr&(q.T_detected>=1)&(q.Treg_detected>=2)&((q.raw_FOXP3>0)|(q.raw_IL2RA>0))).sum()),Endo_multigene=int(((en>=.5)&reliable&(q.Endo_detected>=2)).sum()),Epi_eligible=int(((ep>=.8)&reliable&~mixed).sum()),Epi_score_eligible=int(((ep>=.8)&reliable&~mixed&(q.MP10_detected>=5)&(q.nFeature_RNA>=905)).sum()),mixed_Epi_fraction=float(mixed[ep>=.8].mean()),median_nFeature=float(q.nFeature_RNA.median())))
put(pd.DataFrame(rows),'unit_comparison.tsv')
score=[]
for unit in ['cell','bin8','bin16']:
 a=pd.read_csv(R/'tables'/f'{unit}_score_seed_validation.tsv.gz',sep='\t');b=pd.read_csv(R/'tables'/f'{unit}_score_depth_validation.tsv',sep='\t')
 for scope,sel in [('all',np.ones(len(b),bool)),('Top20_detected_ge5',b.MP10_detected>=5)]:
  z=b[sel];score.append(dict(unit=unit,scope=scope,n=len(z),half_depth_spearman=spearmanr(z.full,z.half_depth).statistic if len(z)>3 else np.nan,full_nCount_dependence_note='see per-unit direct score/QC plot; association not inspected'))
 for mode in ['auc_original','auc_observed']:
  z=a.pivot(index='id',columns='seed',values=mode);score.append(dict(unit=unit,scope=mode,n=len(z),seed123_124_spearman=spearmanr(z[123],z[124]).statistic,median_seed_range=float((z.max(1)-z.min(1)).median())))
put(pd.DataFrame(score),'mp10_score_validation.tsv')
# Assess target availability independently from an association outcome.
eval=[]
for s in ['Cervilla','D1_Johnson']:
 z=pd.read_parquet(R/'tables'/f'{s}_bin8_measurement.parquet');z=z[z.analysis_qc]
 for target,mask,note in [('T_lineage',z.T_detected>=2,'two raw T-lineage panel genes; not a subtype identity'),('Treg_multigene_raw',(z.T_detected>=1)&(z.Treg_detected>=2)&((z.raw_FOXP3>0)|(z.raw_IL2RA>0)),'raw marker co-detection only, not a final label'),('Endothelial',z.Endo_detected>=2,'multi-gene raw endothelial support'),('Epi_MP10',(z.Epi_detected>=2)&(z.MP10_detected>=5),'raw detectability; purity and depth must still be assessed')]:
  eval.append(dict(section=s,target=target,analysis_bins=len(z),raw_supported_bins=int(mask.sum()),status='partially_evaluable',reason=note+('; D1 light check only, no new full annotation or associations'if s!='Cervilla'else'')))
put(pd.DataFrame(eval),'section_target_evaluability.tsv');print(pd.DataFrame(rows).to_string(index=False));print('CELL_TYPES',d.type_supported.value_counts().to_dict());print('SCORE',pd.DataFrame(score).to_string(index=False))
