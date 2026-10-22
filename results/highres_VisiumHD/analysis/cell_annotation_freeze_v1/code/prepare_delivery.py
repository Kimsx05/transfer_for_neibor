# FROZEN_VERSION_GUARD: reproduce in a new version directory; never overwrite v1.
from pathlib import Path as _FreezePath
if (_FreezePath(__file__).resolve().parents[1] / 'metadata/FROZEN.json').exists():
    raise RuntimeError('Frozen v1 is immutable; create v2 for changes.')
"""Assemble one final raw-count object and aligned labels; freeze marker written only after validation."""
from pathlib import Path
import anndata as ad,pandas as pd,numpy as np,json,os
from scipy import sparse
from marker_rules import MAJOR,SUBTYPE
root=Path(__file__).resolve().parents[1];assert not (root/'metadata/FROZEN.json').exists()
a=ad.read_h5ad(root/'cell_annotation.h5ad');o=pd.read_parquet(root/'metadata/computed_cell_labels.parquet');assert a.obs_names.equals(o.index) and o.index.is_unique
reviews=json.loads((root/'metadata/morphology_review_objects.json').read_text())
for col in ['qc_flags','reconstruction_flag']:o[col]=o[col].astype(object)
for r in reviews:
 cid=r['cell_id'];o.loc[cid,'reconstruction_flag']='suspected_merge_or_atypical_large_nucleus';o.loc[cid,'annotation_status']='mixed_or_ambiguous';o.loc[cid,'subtype_label']='mixed_or_ambiguous';o.loc[cid,'unresolved_reason']='morphology_suspected_merge_or_atypical_large_nucleus; native image detail insufficient for definitive rejection';o.loc[cid,'label_evidence']='object_specific_HE_review';o.loc[cid,'qc_flags']=str(o.loc[cid,'qc_flags'])+';morphology_review_uncertain'
 assert not o.loc[cid,'model_eligible'], 'Reviewed large-object flag should already exclude it from structure core'
z=np.load(root/'metadata/expression_embeddings.npz');assert np.array_equal(z['cell_ids'],o.index);cl=pd.read_parquet(root/'metadata/expression_clusters.parquet').set_index('cell_id');assert cl.index.equals(o.index)
o['expression_cluster']=cl.expression_cluster.to_numpy();o['scvi_training_core']=z['core_mask'];o['model_role']=np.where(o.model_eligible,'structure_core','inference_only');o['eligible_to_define_subtypes']=o.model_eligible&~o.low_information&~o.annotation_status.isin(['Unknown','mixed_or_ambiguous'])
# Keep direct observed marker evidence for each model candidate, including unresolved cells.
gm=pd.read_csv(root/'metadata/gene_mapping_Cervilla.tsv',sep='\t');canon=gm.reference_symbol.fillna('').to_numpy();a.var['reference_symbol']=canon
for field,pred,rules in [('observed_major_markers','major_prediction',MAJOR),('observed_subtype_markers','hierarchical_prediction',SUBTYPE)]:
 strings=np.repeat('',a.n_obs).astype(object)
 for label,markers in rules.items():
  select=np.flatnonzero(o[pred].eq(label));found=[];cols=[]
  for gene in markers:
   ids=np.flatnonzero(canon==gene)
   if len(ids):cols.append(np.asarray(a.X[select][:,ids].sum(1)).ravel()>0);found.append(gene)
  if cols:
   hits=np.stack(cols,axis=1)
   strings[select]=[';'.join(np.array(found)[row]) for row in hits]
 o[field]=strings
 if field=='observed_subtype_markers':
  o['subtype_marker_genes_detected']=[0 if not x else len(x.split(';')) for x in strings]
# Annotation evidence from one or very few cells is explicitly marked, never deleted for rarity.
support=o[o.annotation_status.eq('subtype_supported')].subtype_label.value_counts();o['subtype_support_context']=o.subtype_label.map(lambda x:'isolated_candidate_fewer_than_5_supported_cells' if x in support and support[x]<5 else 'multiple_supported_cells' if x in support else 'unresolved_or_broad')
a.obs=o;a.obsm['X_scVI']=z['X_scVI'];a.obsm['X_umap']=z['X_umap'];a.uns['freeze_status']='frozen_subtypes_v1_computational_pending_user_review';a.uns['annotation_version']='frozen_subtypes_v1';a.uns['annotation_rules_path']=str(root/'metadata/label_rules.json');a.uns['scvi_config']=json.loads((root/'metadata/scvi_config.json').read_text());a.uns['annotation_limitations']='Partial image coverage; native 1.643um/pixel; no pathology ground truth; scRNA source-bundle calibration is not HD cross-platform accuracy; Unknown retained';a.uns['unknown_retention']='all valid cells passing raw count/gene floor kept regardless of label confidence';a.uns['model_role_note']='Only 1 section entered scVI; core defined before marker-conflict review; unconfirmed marker conflicts remain visible and never define supported subtypes'
# H5AD cannot serialize None in some nested metadata versions: encode explicit absent batch.
a.uns['scvi_config']['batch_key']='none_single_section'
tmp=root/'cell_annotation.tmp.h5ad';a.write_h5ad(tmp,compression='gzip');tmp.replace(root/'cell_annotation.h5ad');o.to_parquet(root/'frozen_cell_labels.parquet',index=True,compression='zstd')
counts=o.annotation_status.value_counts();summary=json.loads((root/'metadata/annotation_summary.json').read_text());summary.update(annotation_status=counts.to_dict(),major_labels=o.major_label.value_counts().to_dict(),subtype_supported=support.to_dict(),morphology_review_uncertain_retained=len(reviews),low_information_retained=int(o.low_information.sum()),full_sections_frozen=0,partial_sections_frozen=1,blocked_sections=2);(root/'metadata/annotation_summary.json').write_text(json.dumps(summary,indent=2))
dictionary=pd.read_csv(root/'metadata/reference_label_dictionary.tsv',sep='\t');rows=[]
for _,r in dictionary.iterrows():
 lab=r.reference_celltype_v2;cand=o.hierarchical_prediction.eq(lab);hit=o.annotation_status.eq('subtype_supported')&o.subtype_label.eq(lab);n=int(hit.sum());rows.append(dict(subtype=lab,major=r.major_label,reference_n=int(r.n_cells),individual_model_candidates=int(o.individual_prediction.eq(lab).sum()),hierarchical_candidates=int(cand.sum()),subtype_supported=n,low_information_supported=int((hit&o.low_information).sum()),candidate_unresolved_n=int((cand&~hit).sum()),candidate_low_information_fraction=float(o.loc[cand,'low_information'].mean()) if cand.any() else None,interpretation='reference_broad_only' if lab in ['Epithelial','Plasma','Mast'] else 'supported_but_fewer_than_5_cells_review_individually' if 0<n<5 else 'multiple_supported_cells' if n else 'candidates_unresolved' if cand.any() else 'no_supported_detection_not_proof_of_absence'))
pd.DataFrame(rows).to_csv(root/'metadata/subtype_support_summary.tsv',sep='\t',index=False)
pd.crosstab(o.major_label,o.annotation_status).to_csv(root/'metadata/major_by_annotation_status.tsv',sep='\t');pd.crosstab(o.low_information,o.annotation_status).to_csv(root/'metadata/low_information_by_status.tsv',sep='\t')
qc=json.loads((root/'metadata/reconstruction_qc_summary.json').read_text());sections=[dict(dataset_id='Cervilla',section_id='visiumhd_Bladder',candidate=True,bin_usable=True,cell_input_status='qualified_with_native_resolution_and_partial_image_limits',completion='partial_ROI',donor_id='unknown',reconstructed=qc['reconstructed'],retained=a.n_obs,low_information=int(o.low_information.sum()),subtype_supported=int(counts.get('subtype_supported',0)),major_only=int(counts.get('major_only',0)),mixed_or_ambiguous=int(counts.get('mixed_or_ambiguous',0)),Unknown=int(counts.get('Unknown',0)),blocking_gap='original HE not recovered; 1274 in-tissue bins outside available image'),dict(dataset_id='D1_Johnson',section_id='D1_Johnson',candidate=True,bin_usable=True,cell_input_status='blocked',completion='blocked',donor_id='unknown',blocking_gap='Original JohnsonA1_D1_Scan1.qptiff and reliable mapping or provenance-verified raw cell counts+segmentation not publicly located; 5.59um/pixel preview inadequate'),dict(dataset_id='UC_SCC',section_id='UC_SCC',candidate=True,bin_usable=True,cell_input_status='blocked',completion='blocked',donor_id='unknown',blocking_gap='Only 8um counts/coordinates publicly listed; no 2um or provenance-verified author raw cell counts+mapping')];pd.DataFrame(sections).to_csv(root/'metadata/section_status.tsv',sep='\t',index=False,na_rep='NA')
print(json.dumps(summary,indent=2),flush=True)
