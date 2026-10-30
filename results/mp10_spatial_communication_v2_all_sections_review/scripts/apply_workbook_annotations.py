#!/usr/bin/python3
"""Apply user-authorized clinical annotation only; never touch model inputs or v1."""
import csv,json,sys,os,hashlib,argparse
from pathlib import Path
P=Path('/home/data/t070721/codex_workspace/Bladder Metabolism');O=P/'mp10_spatial_communication_v2_all_sections';I=O/'input_review'
sys.path.insert(0,str(P/'cell2location_setup_v1/cell2location_env/lib/python3.12/site-packages'))
import openpyxl
SRC=Path('/home/data/t070721/.codex/attachments/7106e7cf-0ddd-4288-8dac-f82f9c4f1f86/空转+ScRNA样本汇总.xlsx')
REV='user_workbook_annotation_20261027_v1';p=argparse.ArgumentParser();p.add_argument('--inputs-only',action='store_true');args=p.parse_args()
def read(p):
 with p.open(newline='') as f:r=csv.DictReader(f,delimiter='\t');return list(r.fieldnames),list(r)
def write(p,fields,rows):
 assert p.is_relative_to(O) and not p.is_symlink(),str(p)
 tmp=p.with_name(p.name+f'.tmp.{os.getpid()}')
 with tmp.open('w',newline='') as f:w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rows)
 tmp.replace(p)
def js(p,x):
 assert p.is_relative_to(O) and not p.is_symlink(),str(p)
 tmp=p.with_name(p.name+f'.tmp.{os.getpid()}');tmp.write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n');tmp.replace(p)
w=openpyxl.load_workbook(SRC,data_only=True,read_only=True);sheet=w['UC空转样本'];group='';records=[]
for ix,row in enumerate(sheet.iter_rows(values_only=True),1):
 if ix==1:continue
 if row[0] is not None:group=str(row[0])
 if row[1] is not None:records.append(dict(row=ix,group=group,sample=str(row[1]),sample_type=row[3],stage=row[9]))
_,cohort=read(I/'section_status.tsv');assert len(cohort)==63
cd=[r for r in records if r['group'].startswith('CancerDiscovery')];nh=[r for r in records if r['group']=='科室自测数据']
assert len(cd)==13 and sum(x['sample_type']=='MIBC' for x in cd)==12 and sum(x['sample_type']=='NMIBC' for x in cd)==1
assert len(nh)==7 and all(x['sample_type']=='MIBC' for x in nh)
updates={};audit=[]
for r in cohort:
 sec=r['section_id']
 if not sec.startswith(('NPJ_','CanDis_','Inhousedata_')):continue
 z=dict(annotation_revision=REV,tissue_evidence_file=str(SRC),clinical_annotation_source_sheet='UC空转样本',patient_identity_change='none; previous verification status retained',specimen_site='unknown',tumor_involvement_verified='unknown',sample_pathologic_stage='unknown')
 if sec.startswith('NPJ_'):
  h=[x for x in records if x['group'].startswith('NPJ ') and x['sample']==sec.removeprefix('NPJ_')];assert len(h)==1 and h[0]['sample_type']=='NMIBC' and h[0]['stage']=='pT1';h=h[0]
  z.update(tissue_condition='tumor',sample_disease_subtype='NMIBC',clinical_disease_background='NMIBC',sample_pathologic_stage='pT1',tissue_evidence_level='user_workbook_direct_sample_alias',clinical_annotation_source_rows=str(h['row']),tissue_evidence='User-supplied UC spatial sample sheet: same dataset and sample alias explicitly NMIBC/pT1; patient S-alias remains unverified',annotation_limitations='Sample type is annotated; published patient identity and detailed anatomical site remain unverified')
 elif sec.startswith('CanDis_'):
  z.update(tissue_condition='tumor',sample_disease_subtype='unknown',clinical_disease_background='bladder_cancer',tissue_evidence_level='user_workbook_cohort_uniform_tumor_classification',clinical_annotation_source_rows=','.join(str(x['row']) for x in cd),tissue_evidence='Documented CanDis source cohort has 13 workbook samples, all MIBC or NMIBC; broad tumor classification is invariant to unresolved local-to-public mapping',annotation_limitations='Cohort evidence only: 12 MIBC and 1 NMIBC; do not map local Sample numbers by row order or assign individual subtype, stage, patient or treatment')
 else:
  z.update(tissue_condition='unknown',sample_disease_subtype='unknown',clinical_disease_background='MIBC',tissue_evidence_level='user_workbook_cohort_diagnosis_only_specimen_unresolved',clinical_annotation_source_rows=','.join(str(x['row']) for x in nh),tissue_evidence='Workbook RMYY-P1..P6 and RMYY-LN1 all carry MIBC diagnosis; local T-to-RMYY mapping missing; diagnosis alone does not establish specimen tumor involvement',annotation_limitations='Retain unknown tissue condition; cannot identify LN section or its tumor involvement; no patient identity promotion')
 updates[sec]=z
 audit.append(dict(section_id=sec,previous_tissue_condition=r['tissue_condition'],**z))
assert len(updates)==27
snapshot=I/'workbook_annotation_before_20261027.tsv'
if not snapshot.exists():write(snapshot,list(audit[0]),audit)
write(I/'workbook_annotations.tsv',['section_id']+list(next(iter(updates.values()))),[dict(section_id=k,**v) for k,v in updates.items()])
changed=[]
def patch_table(p,append_fields=True):
 if not p.exists():return
 fields,rows=read(p)
 if not rows:return
 idcol='section_id' if 'section_id' in fields else 'section' if 'section' in fields else None
 if idcol is None:return
 oldfields=fields[:];newkeys=list(next(iter(updates.values()))) if append_fields else [k for k in next(iter(updates.values())) if k in fields]
 fields=fields+[k for k in newkeys if k not in fields]
 changed_rows=0
 for row in rows:
  z=updates.get(row[idcol])
  if not z:continue
  before=row.copy()
  for k in newkeys:row[k]=z[k]
  assert all(row[k]==before[k] for k in oldfields if k not in newkeys)
  changed_rows+=1
 if changed_rows:write(p,fields,rows);changed.append(dict(file=str(p.relative_to(O)),updated_rows=changed_rows,non_annotation_fields_unchanged=True))
for fn in ['section_status.tsv','all_sections.tsv','tissue_treatment_evidence.tsv','spot_quality_diagnostics.tsv','geometry_denominator_diagnostics.tsv']:
 patch_table(I/fn)
for sec,z in updates.items():
 d=I/'sections'/sec;j=json.loads((d/'input_summary.json').read_text());before=j.copy();j.update(z);assert all(j[k]==v for k,v in before.items() if k not in z);js(d/'input_summary.json',j)
 # These are independent v2 spot metadata; counts, coordinates, scores, Q and support are untouched.
 path=d/'meta.tsv';fields,rows=read(path);oldfields=fields[:]
 for k in z:
  if k not in fields:fields.append(k)
 for row in rows:
  before=row.copy();row.update(z);assert all(row[k]==before[k] for k in oldfields if k not in z)
 write(path,fields,rows);changed.append(dict(file=str(path.relative_to(O)),updated_rows=len(rows),non_annotation_fields_unchanged=True))
if not args.inputs_only:
 # Completed model-derived tables contain clinical context but the numeric inference stays frozen.
 for folder in ['commot','spatialcellchat']:
  for p in (O/folder).glob('*/*/*.tsv'):
   if (p.parent/'POSTPROCESS_COMPLETED.json').exists() and p.name in ['epi_state_type_contrasts.tsv','lineage_shares.tsv','figure_index.tsv','batch_geometry_diagnostics.tsv']:
    patch_table(p)
if (I/'workbook_unknown27_annotation_audit.tsv').exists():
 fields,rows=read(I/'workbook_unknown27_annotation_audit.tsv')
 for row in rows:row['applied_to_primary_metadata']='True';row['current_tissue_condition']=updates[row['section_id']]['tissue_condition']
 write(I/'workbook_unknown27_annotation_audit.tsv',fields,rows)
_,after=read(I/'section_status.tsv');counts={k:sum(r['tissue_condition']==k for r in after) for k in ['tumor','normal','unknown']};assert counts==dict(tumor=53,normal=3,unknown=7)
for row in after:
 sec=row['section_id'];source=next(x for x in cohort if x['section_id']==sec)
 for k in ['patient_uid','patient_id_status','analysis_unit_type','treatment_status','n_spots','hl_status','score_coverage_status']:assert row[k]==source[k],(sec,k)
js(I/'WORKBOOK_ANNOTATIONS_APPLIED.json',dict(revision=REV,source=str(SRC),source_sha256=hashlib.sha256(SRC.read_bytes()).hexdigest(),n_annotated_sections=27,n_tissue_resolved_from_original_unknown=20,tissue_counts=counts,n_patient_identities_changed=0,n_treatment_statuses_changed=0,inference_recomputed=False,scientific_input_values_changed=False,files=changed))
readme=I/'README.md';text=readme.read_text();mark='## 用户样本表补充标注（2026-10-27）'
if mark not in text:
 text+='\n'+mark+'\n\n用户授权采用《空转+ScRNA样本汇总.xlsx》：NPJ 7片按直接sample别名标tumor、NMIBC、pT1；CanDis 13片按全队列均为MIBC/NMIBC的证据标tumor，逐片亚型不指定；Inhouse 7片仅补MIBC临床背景，组织条件unknown（T与RMYY对应及LN受累未知）。最终53 tumor、3 normal、7 unknown。证据层级、Excel行号和限制见workbook_annotations.tsv；患者核验、治疗、Q、评分及网络均未改变。\n'
 readme.write_text(text)
print(json.dumps(dict(annotation_revision=REV,tissue_counts=counts,metadata_files_updated=len(changed)),ensure_ascii=False))
