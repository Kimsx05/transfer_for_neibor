"""Complete a compact evidence schema without changing any historical estimate or q.
The exact initial M0 export is preserved. No new model is fitted here.
"""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import pathlib,pandas as pd,numpy as np,json,hashlib
P=pathlib.Path.cwd();B=P/'results/mp10_hypothesis_followup_v1';O=P/'results/mp10_immune_endothelial_mechanisms_v1';f=B/'tables/evidence_matrix.tsv';frozen=B/'tables/evidence_matrix_M0_original.tsv'
if not frozen.exists():frozen.write_bytes(f.read_bytes())
d=pd.read_csv(frozen,sep='\t',low_memory=False);candidate=pd.read_csv(O/'tables/candidate_de_effects.tsv',sep='\t');candidate['evidence_type']='historical_candidate_DE';candidate['source_path']=str(O/'tables/candidate_de_effects.tsv');d=pd.concat([d,candidate],ignore_index=True)
plan=pd.read_csv(O/'tables/model_plan.tsv',sep='\t').set_index('model_id')
for col in ['kind','target','min_cells','patient','exposure','subset','omit_dataset','only_dataset','primary']:
 mapped=d.model_id.map(plan[col])
 if col in d:d[col]=d[col].fillna(mapped)
 else:d[col]=mapped
for col in ['variable','effect','P','q','FDR_family','status','reason','receiver','source','TF','gene','regulon','sample_model_p','sample_model_q','activity_effect','AUCell_effect','estimate','p','logFC','FDR','PValue','n_specimens','n_recorded_patient_units','n_datasets']:
 if col not in d:d[col]=np.nan
def coalesce(cols):return d[cols].bfill(axis=1).iloc[:,0]
d['variable']=coalesce(['variable','gene','source','TF','regulon']);d['effect']=coalesce(['effect','logFC','activity_effect','AUCell_effect','estimate']);d['P']=coalesce(['P','sample_model_p','PValue','p']);d['q']=coalesce(['q','sample_model_q','FDR']);d['sender']='Epi';d['receiver']=d.receiver.fillna(pd.Series(np.where(d.kind.eq('Epi'),'Epi',d.target),index=d.index));d['FDR_family']=d.FDR_family.fillna(pd.Series(np.where(d.evidence_type.eq('historical_candidate_DE'),'historical full-gene BH per model','historical per-model regulon or abundance design family; see source'),index=d.index))
d['statistical_unit']=np.where(d.kind.eq('Epi'),'within-specimen paired RNA profile','dataset-adjusted recorded patient unit; identities provisional');d['sensitivity']=d[['min_cells','patient','exposure','subset','omit_dataset','only_dataset']].fillna('NA').astype(str).agg(';'.join,axis=1);d['evidence_nature']='same scRNA data: association or rank enrichment, not independent validation';d['shared_data']='v1 RNA counts and nested/overlapping cohorts';d['shared_genes']=d['mask'].fillna('not_applicable') if 'mask' in d else 'see resource targets';d['gap']='protein/kinase/flux/function/direction not measured; P and q retain original family';d['status']=d.status.fillna(pd.Series(np.where(d.P.notna(),'computed_historical','not_evaluable'),index=d.index))
# Every row has an unambiguous identifier, including abundance rows without a DE model ID.
d['model_id']=d.model_id.fillna(pd.Series(['historical_abundance_row_'+str(i) for i in range(len(d))]));d['evidence_id']=['E%06d'%i for i in range(1,len(d)+1)];d['CI_low']=np.nan;d['CI_high']=np.nan
cols=['evidence_id','sender','receiver','variable','source_path','model_id','statistical_unit','n_specimens','n_recorded_patient_units','n_datasets','effect','CI_low','CI_high','P','q','joint_main_q','FDR_family','sensitivity','evidence_type','evidence_nature','shared_data','shared_genes','status','gap']
d[cols].to_csv(f,sep='\t',index=False)
axis=pd.read_csv(O/'tables/key_axis_evidence.tsv',sep='\t');axis['source_path']=str(O/'tables/key_axis_evidence.tsv');axis.to_csv(B/'tables/reused_contact_endothelial_axis.tsv',sep='\t',index=False)
lit=pd.read_csv(B/'tables/literature_evidence.tsv',sep='\t');records=json.loads((B/'logs/literature_source_records.json').read_text());maps={x['id']:x for x in records};lit['species_metadata']=lit.PMID.astype(str).map(lambda x:';'.join(y.get('descriptorName','') for y in maps[x].get('meshHeadingList',{}).get('meshHeading',[]) if y.get('descriptorName','') in ['Humans','Mice','Animals','Mice, Inbred C57BL','Mice, Knockout']) or 'not fully specified in retrieved abstract/metadata')
lit['intervened_cell']= ['cell-intrinsic mTORC1/S6K/SREBP model','BLCA cancer cells: FGFR3/SCD','Treg-specific Rptor deletion','gastric cancer-cell RHOA/PI3K environment','tumour Treg SCAP/FASN','tumour Treg CD36','tumour Treg metabolic programs','Treg mitochondrial complex III']
lit['cancer_scope']=['not cancer-type-specific','bladder cancer','immune homeostasis, not an Epi-extrinsic cancer experiment','gastric cancer','tumour models; full model list not inferred from abstract','tumour models; full model list not inferred from abstract','murine tumours plus human tumour observations','immune homeostasis / inflammatory phenotype']
lit.to_csv(B/'tables/literature_evidence.tsv',sep='\t',index=False)
(B/'config/evidence_schema_completion.md').write_text('The original M0 evidence_matrix_M0_original.tsv is retained unchanged. This delivery schema completion adds canonical sender/receiver, units, historical candidate DE and source fields from already-frozen v1 tables. No estimate, P, q or hypothesis candidate was selected or refitted. Full original historical tables remain authoritative for family-specific details. CI not present in the original export is NA, not manufactured from a different inferential method. The primary-paper retrieval level is explicitly retained; unavailable full-text details are not inferred.\n')
print('Evidence schema rows',len(d),'candidate rows',len(candidate))
