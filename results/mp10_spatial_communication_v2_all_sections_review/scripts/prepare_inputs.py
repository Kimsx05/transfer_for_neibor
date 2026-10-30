#!/usr/bin/env python3
"""Read-only reuse of accepted spatial assets; no model fitting or old-file writes."""
import os,sys,json,time,gzip,traceback,shutil
from pathlib import Path
P=Path('/home/data/t070721/codex_workspace/Bladder Metabolism')
sys.path=[p for p in sys.path if 'dist-packages' not in p]
sys.path[:0]=[str(P/'spatial_ecology_top20_v2/runtime'),str(P/'cell2location_setup_v1/cell2location_env/lib/python3.12/site-packages')]
for k in ('ARROW_NUM_THREADS','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import numpy as np,pandas as pd,h5py
from scipy import sparse,io
from scipy.io import _fast_matrix_market
_fast_matrix_market.PARALLELISM=1
from threadpoolctl import threadpool_limits
_thread_limits=threadpool_limits(limits=1)
from scipy.spatial import cKDTree
V=P/'mp10_spatial_communication_v1';R=P/'mp10_spatial_communication_v2_all_sections';O=R/'input_review';O.mkdir(parents=True,exist_ok=True)
A=P/'spatial_ecology_top20_v2/00_inputs';reuse=['GSE171351_A1','Inhousedata_T4469479'];start=time.time()
def tsv(df,p):
 p=Path(p);tmp=p.with_suffix(p.suffix+'.tmp');pd.DataFrame(df).to_csv(tmp,sep='\t',index=False,na_rep='NA');tmp.replace(p)
def js(obj,p):
 p=Path(p);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2,ensure_ascii=False,default=lambda x:x.item() if hasattr(x,'item') else str(x))+'\n');tmp.replace(p)
def strings(g):
 if isinstance(g,h5py.Group):
  if 'values' in g:return strings(g['values'])
  cats=strings(g['categories']);return [cats[i] for i in g['codes'][()]]
 return [v.decode() if isinstance(v,bytes) else str(v) for v in g[()]]
manifest=pd.read_csv(A/'section_manifest.tsv',sep='\t').set_index('section_id');cm=pd.read_csv(A/'c2l_manifest.tsv',sep='\t').set_index('section_id')
spots=pd.read_parquet(A/'frozen_spot_inputs.parquet');mapping=pd.read_csv(A/'major_map.tsv',sep='\t').set_index('subtype');mapping['major']=mapping.major.replace({'T_NK':'T/NK'})
usage=pd.read_csv(P/'cell2location_final_validation_gate_v1/09_final_status/DOWNSTREAM_USAGE_TABLE.csv').set_index('factor')
ids_path=P/'spatial_ecology_top20_v2/patient_identity_audit_20261002/revised_by_user_v2/section_patient_and_analysis_unit_crosswalk.tsv'
ident=pd.read_csv(ids_path,sep='\t').set_index('section_id')
genesig=pd.read_csv(A/'frozen_signature_manifest.tsv',sep='\t').set_index('MP').loc['EpiMP10','gene_list'].split(';')
mother=pd.read_csv(V/'input_review/candidate_sections.tsv',sep='\t')
assert len(mother)==63 and mother.section_id.is_unique
secs=['CanDis_2026_Sample11','GSE238145_B10']+reuse+[x for x in mother.sort_values('n_spots').section_id if x not in reuse+['CanDis_2026_Sample11','GSE238145_B10']]
manifest=manifest.loc[secs]
assert set(manifest.index)==set(mother.section_id)
tsv(mother.assign(v1_selected=mother.selected,selected=True,selection_basis='all 63 mother-list ordinary Visium sections; no pilot-selection filter'),O/'candidate_sections.tsv')
diffs=[]
for sec in secs:
 old=mother.set_index('section_id').loc[sec]
 changed=[]
 if int(old.n_spots)!=int(manifest.loc[sec,'n_spots']):changed.append('n_spots')
 if str(old.dataset)!=str(manifest.loc[sec,'dataset']):changed.append('dataset')
 diffs.append(dict(section_id=sec,change=';'.join(changed) if changed else 'none',v1_n_spots=int(old.n_spots),final_n_spots=int(manifest.loc[sec,'n_spots']),final_input_index=str(A/'section_manifest.tsv'),identity_index=str(ids_path)))
tsv(diffs,O/'input_index_revision_audit.tsv')
lattice_source=P/'neighborhood_v1/00_inputs_and_audit/lattice_recovery'
lattice_coverage=pd.read_csv(lattice_source/'04_coordinate_coverage_by_section.csv').set_index('section')
lattice_topology=pd.read_csv(lattice_source/'06_topology_qc_by_section.csv').set_index('section')
assert lattice_coverage.loc[secs,'status'].eq('PASS').all() and lattice_topology.loc[secs,'status'].eq('PASS').all()
source_dir=P/'spatial_ecology_top20_v2/patient_identity_audit_20261002/sources'
conditions=[]
for sec in secs:
 ds=str(manifest.loc[sec,'dataset']);tissue='unknown';evidence='No verified section-level tissue metadata available';src=''
 if ds in ['GSE171351','GSE238145','GSE319536']:
  tissue='tumor';src=str(source_dir/(ds+'_family.soft'));evidence='All source GEO spatial samples explicitly annotated bladder cancer/tumor/carcinoma; inherited source series and audited sample mapping'
 elif ds=='GSE285715':
  title=sec.removeprefix('GSE285715_');records=(source_dir/'GSE285715_family.soft').read_text().split('^SAMPLE = ')[1:]
  matched=[r for r in records if ('!Sample_title = '+title+'\n') in r]
  assert len(matched)==1
  source=[l for l in matched[0].splitlines() if l.startswith('!Sample_source_name_ch1 = ')][0]
  tissue='normal' if source.endswith('Normal tissue adjacent to bladder cancer') else ('tumor' if source.endswith('= bladder cancer') else 'unknown')
  src=str(source_dir/'GSE285715_family.soft');evidence=source+'; exact GEO title match; no patient pairing inferred'
 conditions.append(dict(section_id=sec,tissue_condition=tissue,treatment_status='unknown',tissue_evidence=evidence,tissue_evidence_file=src,treatment_evidence='No verified section-level treatment metadata supplied in frozen input/identity tables; unknown preserved'))
cond=pd.DataFrame(conditions).set_index('section_id');tsv(cond.reset_index(),O/'tissue_treatment_evidence.tsv')
tsv(ident.loc[secs].reset_index(),O/'section_identity_crosswalk.tsv')
lowid=usage.index[usage.status.eq('LOW_IDENTIFIABILITY')].tolist()
index=[];pilot=[];sup=[];gencov=[];params=[];checks=[];cross=[];candidates=[]
failures=[]
for order,sec in enumerate(secs,1):
 try:
  t0=time.time();d=O/'sections'/sec;d.mkdir(parents=True,exist_ok=True);r=manifest.loc[sec];c=cm.loc[sec];identrow=ident.loc[sec];condition=cond.loc[sec];ssup=[];scross=[];sgencov=[]
  if (d/'READY.json').exists():
   if (d/'FAILED.json').exists():(d/'FAILED.json').unlink()
   pilot.append(json.loads((d/'input_summary.json').read_text()));sup.extend(pd.read_csv(d/'type_support.tsv',sep='\t').to_dict('records'));cross.extend(pd.read_csv(d/'spatial_label_crosswalk.tsv',sep='\t').to_dict('records'));gencov.extend(pd.read_csv(d/'gene_coverage.tsv',sep='\t').to_dict('records'));print(sec,'CACHED_READY',flush=True);continue
  with h5py.File(r.raw_input_path,'r') as h:
   barcodes=strings(h['obs/_index']);genes=strings(h['var/_index']);Xg=h['X'];shape=tuple(Xg.attrs['shape']);X=sparse.csr_matrix((Xg['data'][()],Xg['indices'][()],Xg['indptr'][()]),shape=shape)
   assert str(Xg.attrs['encoding-type'])=='csr_matrix';xy=h['obsm/spatial'][()]
  assert X.shape==(len(barcodes),len(genes)) and len(set(genes))==len(genes) and len(set(barcodes))==len(barcodes)
  assert np.isfinite(X.data).all() and (X.data>=0).all() and (X.data==np.round(X.data)).all()
  s=spots[spots.section_id.eq(sec)].set_index('barcode').loc[barcodes].copy();assert np.max(np.abs(xy-s[['x','y']].values))<1e-10
  mean=pd.read_csv(c.mean_path,index_col=0,float_precision='round_trip');assert list(mean.index)==barcodes
  rawQ=mean.div(mean.sum(axis=1),axis=0);assert np.isfinite(rawQ.values).all();assert set(rawQ.columns)==set(mapping.index)
  assert np.max(abs(rawQ.Epithelial.values-s.epi_fraction_gate.values))<1e-12
  score=s.EpiMP10_Top20_AUCell;anchors=rawQ.Epithelial.ge(.25);assert np.array_equal(anchors.values,s.analysis_eligible_legacy.values.astype(bool));threshold=float(score[anchors&np.isfinite(score)].median()) if (anchors&np.isfinite(score)).any() else float('nan');high=anchors&score.gt(threshold);low=anchors&score.lt(threshold);state=np.where(high,'high',np.where(low,'low','other'))
  uid=[sec+'::'+b for b in barcodes]
  meta=pd.DataFrame(dict(spot_id=uid,barcode=barcodes,section_id=sec,dataset=r.dataset,record_sample_id=sec,patient_uid=identrow.patient_uid,patient_id_status=identrow.patient_id_status,analysis_unit_uid=identrow.analysis_unit_uid,analysis_unit_type=identrow.analysis_unit_type,MP10_frozen=score.values,MP10_threshold_new=threshold,epi_anchor=anchors.values,epi_state=state,epi_fraction_frozen=rawQ.Epithelial.values,library_counts=np.asarray(X.sum(axis=1)).ravel(),n_genes_detected=X.getnnz(axis=1)))
  q=rawQ.drop(columns=['Epithelial']+lowid).copy();q['T_NK_unresolved']=rawQ[lowid].sum(axis=1)
  for st,mask in [('high',high.values),('low',low.values),('other',~(high.values|low.values))]:q['Epi_MP10_'+st]=rawQ.Epithelial.values*mask
  # Fixed display/inference schema includes every modeled mass; no dominant type or .1 cutoff.
  q=q[['Epi_MP10_high','Epi_MP10_low','Epi_MP10_other']+[x for x in q.columns if not x.startswith('Epi_MP10_')]]
  qmajor=pd.DataFrame(index=q.index)
  for st in ['high','low','other']:qmajor['Epi_MP10_'+st]=q['Epi_MP10_'+st]
  majororder=['B','Plasma','T/NK','Myeloid','Fibroblast','Mural','Endothelial','Mast']
  for major in majororder:qmajor[major]=rawQ[mapping.index[mapping.major.eq(major)]].sum(axis=1)
  assert np.max(abs(q.sum(axis=1)-1))<1e-12 and np.max(abs(qmajor.sum(axis=1)-1))<1e-12
  # Recovered actual hexagonal array indices; nominal standard Visium pitch, not misleading image scalefactor.
  coord=pd.DataFrame(dict(spot_id=uid,x_um=s.array_col.values*50.,y_um=s.array_row.values*np.sqrt(3)*50.,x_image_px=xy[:,0],y_image_px=xy[:,1],array_row=s.array_row.values,array_col=s.array_col.values))
  nn=cKDTree(coord[['x_um','y_um']]).query(coord[['x_um','y_um']],k=2)[0][:,1];pixnn=cKDTree(xy).query(xy,k=2)[0][:,1]
  lattice=np.column_stack([s.array_col.values*.5,s.array_row.values*np.sqrt(3)/2,np.ones(len(s))]);affine=np.linalg.lstsq(lattice,xy,rcond=None)[0];rms=float(np.sqrt(np.mean((lattice@affine-xy)**2)))
  assert np.isfinite(coord[['x_um','y_um']].values).all() and (nn>0).all()
  # Missing hard-QC neighbors never expand nominal platform pitch or range.
  assert np.isfinite(rms), 'nonfinite affine residual'
  arr=s[['array_row','array_col']].values
  assert np.max(abs(arr-np.round(arr)))==0 and len(np.unique((arr[:,0]+arr[:,1])%2))==1 and len(np.unique(arr,axis=0))==len(arr), 'invalid hexagonal array indices'
  # Image coordinates have source-dependent pixel rounding; this diagnostic never rescales true lattice positions.
  affine_relative_rms=rms/float(np.median(pixnn))
  assert lattice_coverage.loc[sec,'n_common_spots']==len(s) and lattice_topology.loc[sec,'n_spots']==len(s)
  # True recovered VisiumV1 array coordinates are authoritative; image affine fit is diagnostic only.
  for field in ['tissue_condition','treatment_status']:meta[field]=condition[field]
  tsv(meta,d/'meta.tsv');tsv(coord,d/'coords.tsv');tsv([dict(gene=g) for g in genes],d/'genes.tsv');tsv([dict(spot_id=u,barcode=b) for u,b in zip(uid,barcodes)],d/'spots.tsv')
  for name,df in [('Q_original_30',rawQ),('Q_fine',q),('Q_major',qmajor),('abundance_mean_original_30',mean)]:
   z=df.copy();z.insert(0,'spot_id',uid);tsv(z.reset_index(drop=True),d/(name+'.tsv'))
  if sec in reuse:
   old=V/'input_review/sections'/sec
   for name in ['counts.mtx.gz','counts_spot_gene.npz','log1p_spot_gene.npz','log1p.mtx.gz']:
    target=d/name
    if not target.exists():target.symlink_to(old/name)
   # Q and coordinates were reconstructed only as an equivalence check; use v1 files after exact numeric agreement.
   for name in ['Q_original_30.tsv','Q_fine.tsv','Q_major.tsv','abundance_mean_original_30.tsv','coords.tsv','genes.tsv','spots.tsv']:
    olddf=pd.read_csv(old/name,sep='\t');newdf=pd.read_csv(d/name,sep='\t')
    pd.testing.assert_frame_equal(olddf,newdf,check_exact=False,atol=1e-12,rtol=1e-12)
    (d/name).unlink();(d/name).symlink_to(old/name)
  else:
   with gzip.open(d/'counts.mtx.gz','wb',compresslevel=1) as handle:io.mmwrite(handle,X.T,field='integer',comment='gene x spot; complete original hard-QC positions, all measured genes')
   sparse.save_npz(d/'counts_spot_gene.npz',X)
  normerr=0.0
  for level,df in [('fine',q),('major',qmajor)]:
   for name in df:
    if name.startswith('Epi_MP10_'):
     ab=mean.Epithelial.values*(state==name.removeprefix('Epi_MP10_'));status='STATE_ASSOCIATED_SPOT_EPITHELIAL_COMPONENT'
    elif level=='major':ab=mean[mapping.index[mapping.major.eq(name)]].sum(axis=1).values;status='MAJOR_FROM_ACCEPTED_MEANS'
    elif name=='T_NK_unresolved':ab=mean[lowid].sum(axis=1).values;status='LOW_IDENTIFIABILITY_MASS_PRESERVED_UNRESOLVED'
    else:ab=mean[name].values;status=usage.loc[name,'status']
    eligible=bool((ab>=.01).sum()>=10 and ab.sum()>1)
    ssup.append(dict(section_id=sec,dataset=r.dataset,level=level,node=name,n_positive_Q=int((df[name]>0).sum()),n_abundance_ge_0_01=int((ab>=.01).sum()),total_expected_abundance=float(ab.sum()),effective_spot_mass=float(df[name].sum()),max_Q=float(df[name].max()),n_dominant_spots=int((df.idxmax(axis=1)==name).sum()),n_dominant_original_30=int((rawQ.idxmax(axis=1)==name).sum()) if name in rawQ else np.nan,type_evaluable=eligible,support_rule='accepted validation: >=10 spots mean abundance>=0.01 and total abundance>1; no spot weights masked',identifiability_status=status,composition_weight='posterior_mean / sum_all_original_30_means'))
  for gene in genesig:sgencov.append(dict(section_id=sec,dataset=r.dataset,gene_set='MP10_frozen_Top20',gene=gene,measured=gene in genes,detected_spots=int(X[:,genes.index(gene)].getnnz()) if gene in genes else np.nan,not_measured_is_zero=False))
  for original in rawQ:
   out='Epi_MP10_high|Epi_MP10_low|Epi_MP10_other' if original=='Epithelial' else ('T_NK_unresolved' if original in lowid else original)
   scross.append(dict(section_id=sec,original_factor=original,major=mapping.loc[original,'major'],analysis_fine_label=out,identifiability_status=usage.loc[original,'status'],downstream_usage=usage.loc[original,'downstream_usage'],reason=usage.loc[original,'reason'],mapping_source=str(A/'major_map.tsv'),reliability_source=str(P/'cell2location_final_validation_gate_v1/09_final_status/DOWNSTREAM_USAGE_TABLE.csv'),mass_action='split within epi only' if original=='Epithelial' else ('merge only three unidentifiable labels; no redistribution' if original in lowid else 'retained')))
  # Retain absent fine labels as NE in crosswalk, never invented columns.
  for absent,major in [('Tip EC','Endothelial'),('Cycling EC','Endothelial')]:scross.append(dict(section_id=sec,original_factor=absent,major=major,analysis_fine_label=np.nan,identifiability_status='NE_NOT_IN_SPATIAL_REFERENCE',downstream_usage='NE',reason='No independently estimated composition column in accepted 30-factor reference',mass_action='not fabricated'))
  nsum=dict(section_id=sec,dataset=r.dataset,preparation_order=order,n_spots=len(s),n_measured_genes=len(genes),n_mp10_genes_measured=sum(g in genes for g in genesig),n_epi_anchors=int(anchors.sum()),n_epi_high=int(high.sum()),n_epi_low=int(low.sum()),n_epi_equal_median=int((anchors&score.eq(threshold)).sum()),n_epi_other_spots=int((state=='other').sum()),threshold=threshold,n_original_factors=30,n_fine_columns=len(q.columns),n_major_columns=len(qmajor.columns),observational_unit='mixed Visium spot',coordinate_units='um; nominal platform array geometry',median_nn_um=float(np.median(nn)),median_nn_image_px=float(np.median(pixnn)),image_lattice_affine_rms_px=rms,patient_uid=identrow.patient_uid,patient_id_status=identrow.patient_id_status,analysis_unit_uid=identrow.analysis_unit_uid,selection_basis='all 63 ordinary Visium mother-list sections; complete original hard-QC background',prepared_dir=str(d),raw_input_path=r.raw_input_path,mean_path=c.mean_path)
  epi_support={x['node']:x['type_evaluable'] for x in ssup if x['level']=='major'}
  hl_ok=epi_support['Epi_MP10_high'] and epi_support['Epi_MP10_low']
  nsum.update(network_valid=True,input_status='ready',hl_status='evaluable' if hl_ok else 'network_only',score_coverage_status='complete_20_of_20' if nsum['n_mp10_genes_measured']==20 else 'incomplete_'+str(nsum['n_mp10_genes_measured'])+'_of_20',tissue_condition=condition.tissue_condition,treatment_status=condition.treatment_status,tissue_evidence=condition.tissue_evidence,tissue_evidence_file=condition.tissue_evidence_file,treatment_evidence=condition.treatment_evidence,analysis_unit_type=identrow.analysis_unit_type,hl_reason='both epi states pass original abundance support' if hl_ok else 'insufficient Epi high/low original abundance support; do not lower threshold',native_pitch_um=100,secreted_range_um=200,ecm_range_um=200,contact_range_um=105,reused_v1=sec in reuse,image_lattice_relative_rms=affine_relative_rms,coordinate_audit='authoritative array coverage/topology PASS; image residual diagnostic only',image_affine_warning=bool(affine_relative_rms>=.05))
  nsum={k:None if isinstance(v,float) and np.isnan(v) else v for k,v in nsum.items()}
  pilot.append(nsum);js(nsum,d/'input_summary.json')
  sup.extend(ssup);cross.extend(scross);gencov.extend(sgencov)
  tsv(ssup,d/'type_support.tsv');tsv(scross,d/'spatial_label_crosswalk.tsv');tsv(sgencov,d/'gene_coverage.tsv')
  for field in ['hl_status','score_coverage_status']:meta[field]=nsum[field]
  tsv(meta,d/'meta.tsv')
  for key,value,unit,source in [('nominal_spot_diameter',55,'um','standard Visium manufacturer platform specification; not inferred from invalid spot scale_factor'),('nominal_center_pitch',100,'um','actual recovered ordinary Visium lattice + manufacturer platform specification'),('secreted_range',200,'um','fixed 2 x nominal platform pitch; not recalibrated after QC'),('contact_range',105,'um','fixed nominal platform pitch + 5%; no claim of membrane contact'),('ecm_range',200,'um','predeclared 2 x d_NN; SpatialCellChat native ECM shares interaction.range with secreted; fixed before first inference results'),('first_section_secreted_sensitivity',300,'um','3 x d_NN, first selected section only; SpatialCellChat changes ECM and secreted together by native API; COMMOT changes secreted only'),('expression_normalization_target',10000,'counts per spot','full measured gene library normalization once; log1p'),('epi_anchor_threshold',.25,'fraction','frozen accepted c2l mean epithelial / total mean'),('MP10_HL_threshold',threshold,'frozen AUCell','new per-section anchors-only strict median; no overwritten original score'),('spot_weight_threshold',0,'Q','none: all positive posterior means preserved; no 0.1 cutoff')]:params.append(dict(section_id=sec,parameter=key,value=value,unit=unit,source=source))
  for role,path in [('raw_counts_expression',r.raw_input_path),('posterior_mean_abundance',c.mean_path),('frozen_score',r.score_path),('frozen_spots_anchors_and_lattice',A/'frozen_spot_inputs.parquet'),('patient_identity',ids_path),('accepted_major_map',A/'major_map.tsv')]:index.append(dict(section_id=sec,dataset=r.dataset,role=role,path=str(path),exists=Path(path).is_file(),read_only=True,observation_unit='spot',platform='ordinary_Visium',measured_expression_layer='X raw counts',coordinate_reference='true array row/col; original lowres image pixels also preserved',patient_id_status=identrow.patient_id_status,id_alignment='exact barcode bijection; expression,coordinates,score,abundance same order'))
  checks.extend([dict(section_id=sec,check='barcode_bijection_and_order',value=len(s),pass_check=True),dict(section_id=sec,check='fine_Q_mass_max_error',value=float(abs(q.sum(axis=1)-rawQ.sum(axis=1)).max()),pass_check=True),dict(section_id=sec,check='major_Q_mass_max_error',value=float(abs(qmajor.sum(axis=1)-rawQ.sum(axis=1)).max()),pass_check=True),dict(section_id=sec,check='TME_shared_Q_unchanged',value=float(abs(q['Treg']-rawQ['Treg']).max()),pass_check=True),dict(section_id=sec,check='measured_counts_nonnegative_integer',value=X.nnz,pass_check=True),dict(section_id=sec,check='nominal_nn_um',value=float(np.median(nn)),pass_check=True),dict(section_id=sec,check='lattice_image_affine_rms_px',value=rms,pass_check=True)])
  if (d/'FAILED.json').exists():(d/'FAILED.json').unlink()
  js(dict(status='READY',section_id=sec,seconds=time.time()-t0,network_valid=True,hl_status=nsum['hl_status'],score_coverage_status=nsum['score_coverage_status'],reused_v1=sec in reuse),d/'READY.json')
  tsv(pilot,O/'section_status.tsv');tsv(sup,O/'type_support.tsv');tsv(cross,O/'spatial_label_crosswalk.tsv')
  print(sec,'READY',nsum,'seconds',time.time()-t0,flush=True)
 except Exception as exc:
  detail=dict(section_id=sec,status='failed',reason=str(exc),traceback=traceback.format_exc(),prepared_dir=str(O/'sections'/sec),recovery='rerun scripts/prepare_inputs.py; completed READY inputs reused')
  with (R/'logs/input_attempt_history.jsonl').open('a') as h:h.write(json.dumps(detail,ensure_ascii=False)+'\n')
  failures.append(detail);js(detail,O/'sections'/sec/'FAILED.json');tsv(failures,O/'input_failures.tsv');print(sec,'FAILED',traceback.format_exc(),flush=True)
# Rebuild lightweight provenance/check/parameter tables from every READY summary on resume.
index=[];params=[];checks=[]
for row in pilot:
 sec=row['section_id'];r=manifest.loc[sec];c=cm.loc[sec]
 for role,path in [('raw_counts_expression',r.raw_input_path),('posterior_mean_abundance',c.mean_path),('frozen_score',r.score_path),('frozen_spots_anchors_and_lattice',A/'frozen_spot_inputs.parquet'),('patient_identity',ids_path),('accepted_major_map',A/'major_map.tsv')]:
  index.append(dict(section_id=sec,dataset=r.dataset,role=role,path=str(path),exists=Path(path).is_file(),read_only=True,observation_unit='spot',platform='ordinary_Visium',patient_id_status=row['patient_id_status'],id_alignment='exact barcode bijection; counts/coords/Q/score same order'))
 for key,value,unit in [('nominal_center_pitch',100,'um'),('secreted_range',200,'um'),('ecm_range',200,'um'),('contact_range',105,'um'),('expression_normalization_target',10000,'all measured counts per spot before log1p'),('epi_anchor_threshold',.25,'fraction'),('MP10_HL_threshold',row['threshold'],'frozen AUCell'),('spot_weight_threshold',0,'no mask'),('support_min_spots',10,'spots with posterior mean abundance >=0.01'),('support_min_mean_abundance',.01,'original posterior mean abundance'),('support_total_abundance_strict_min',1,'strict >; no Q clipping'),('distance_numeric_epsilon',1e-6,'um floating point comparison only')]:params.append(dict(section_id=sec,parameter=key,value=value,unit=unit,source='v1 frozen definition; nominal platform pitch not changed after hard QC'))
 for name,value in [('barcode_bijection_and_order',row['n_spots']),('fine_Q_mass_max_error',0.0),('major_Q_mass_max_error',0.0),('TME_shared_Q_unchanged',0.0),('lattice_image_affine_rms_px',row['image_lattice_affine_rms_px']),('nominal_nn_um',row['median_nn_um'])]:checks.append(dict(section_id=sec,check=name,value=value,pass_check=True,note='READY requires explicit assertions; Q numerical residual below1e-12'))
for name,rows in [('input_index',index),('all_sections',pilot),('spatial_label_crosswalk',cross),('type_support',sup),('gene_coverage',gencov),('parameters',params),('input_checks',checks)]:tsv(rows,O/(name+'.tsv'))
tsv(pilot,O/'section_status.tsv')
tsv(pd.DataFrame(failures,columns=['section_id','status','reason','traceback','prepared_dir','recovery']),O/'input_failures.tsv')
js(dict(status='READY' if not failures else 'COMPLETE_WITH_FAILURES',n_ready=len(pilot),failures=failures,sections=secs,seconds=time.time()-start,blas_threads=1,new_annotations=False,new_c2l_training=False),O/'PREPARED.json')

(O/'README.md').write_text("# 全部普通 Visium 输入审计\n\n母清单63片全部准入，v1 selected仅保存为历史标记。冻结RNA counts、AUCell、accepted posterior mean、身份索引均只读。两个试点counts/Q/coordinates经逐元素等价验证后软链接v1；meta/support为本轮独立文件。每片READY.json为队列入口；失败记录FAILED.json并继续其他切片。\n\nQ以全部30参考类型归一化；Epi严格中位数状态；支持为原始mean abundance≥0.01的spot至少10且总abundance>1，不截断Q。固定100μm平台节距，200/105μm通讯范围；QC缺失邻居不改变节距。所有全片通过QC位置保留。没有足够H/L仍network_only，19/20评分覆盖单独汇总。\n\n组织条件来自已存GEO SOFT明确注释并记录tissue_treatment_evidence.tsv；未核验组织/治疗保留unknown，不由文件名前缀猜测。身份沿用正式revised_by_user_v2交叉表，患者未核验不视独立患者。\n\n恢复：/usr/bin/python3 scripts/prepare_inputs.py，跳过已READY切片。无大文件hash，无母RNA对象复制。\n")

# Lightweight independent verification and quality/geometry audit from persisted inputs.
audit=[];quality=[];geometry=[];summaries=[]
for sec in secs:
 d=O/'sections'/sec
 if not (d/'READY.json').exists():continue
 info=json.loads((d/'input_summary.json').read_text());meta=pd.read_csv(d/'meta.tsv',sep='\t');co=pd.read_csv(d/'coords.tsv',sep='\t')
 qs={lev:pd.read_csv(d/('Q_'+lev+'.tsv'),sep='\t').set_index('spot_id') for lev in ['original_30','major','fine']}
 orig=qs['original_30'];ab=pd.read_csv(d/'abundance_mean_original_30.tsv',sep='\t').set_index('spot_id')
 assert list(orig.index)==list(meta.spot_id)==list(co.spot_id)
 assert np.max(abs(orig.values-ab.div(ab.sum(axis=1),axis=0).values))<1e-12
 a=co[['array_row','array_col']].values
 arr_integer=float(np.max(abs(a-np.round(a))));parities=np.unique((a[:,0]+a[:,1])%2)
 assert arr_integer==0 and len(parities)==1 and len(np.unique(a,axis=0))==len(a)
 epi_cols=['Epi_MP10_high','Epi_MP10_low','Epi_MP10_other']
 for lev in ['major','fine']:
  q=qs[lev]
  err=float(np.max(abs(q.sum(axis=1)-1)));epierr=float(np.max(abs(q[epi_cols].sum(axis=1)-orig.Epithelial)))
  assert err<1e-12 and epierr<1e-12
  audit.extend([dict(section_id=sec,check=lev+'_Q_mass_max_error',value=err,pass_check=True),dict(section_id=sec,check=lev+'_Epi_split_mass_max_error',value=epierr,pass_check=True)])
 assert (meta.library_counts>0).all() and len(meta)==int(manifest.loc[sec,'n_spots'])
 audit.extend([dict(section_id=sec,check='array_integer_max_error',value=arr_integer,pass_check=True),dict(section_id=sec,check='array_parity_classes',value=len(parities),pass_check=True),dict(section_id=sec,check='full_hard_QC_positions_preserved',value=len(meta),pass_check=True),dict(section_id=sec,check='raw_mean_Q_max_error',value=float(np.max(abs(orig.values-ab.div(ab.sum(axis=1),axis=0).values))),pass_check=True)])
 rel=float(info['image_lattice_affine_rms_px']/info['median_nn_image_px'])
 info.update(image_lattice_relative_rms=rel,coordinate_audit='authoritative VisiumV1 array coverage and topology PASS; image affine residual diagnostic only',image_affine_warning=rel>=.05,lattice_source=str(lattice_source/'03_recovered_coordinates_by_section.csv.gz'))
 if info['hl_status']=='evaluable':code='BOTH_STATES_SUPPORTED';reason='both Epi states satisfy original abundance support'
 elif info['n_epi_anchors']==0:code='NO_EPI_ANCHORS';reason='no spots meet Epi Q>=0.25; H/L NE'
 elif info['threshold']==0 and info['n_epi_low']==0:code='STRICT_MEDIAN_ZERO_NO_LOW';reason='anchor median AUCell=0; no strict-low anchors; median ties remain other; H/L NE'
 else:code='STATE_SUPPORT_INSUFFICIENT';reason='one or both Epi states fail original abundance support; H/L NE'
 info.update(hl_reason_code=code,hl_reason=reason,hl_evaluable=info['hl_status']=='evaluable')
 js(info,d/'input_summary.json');summaries.append(info)
 for state in ['high','low','other']:
  mask=meta.epi_state.eq(state).values;w=qs['major']['Epi_MP10_'+state].values;den=w.sum()
  quality.append(dict(section_id=sec,dataset=info['dataset'],tissue_condition=info['tissue_condition'],treatment_status=info['treatment_status'],epi_state=state,n_spots=int(mask.sum()),epi_Q_mass=float(den),median_UMI=float(meta.loc[mask,'library_counts'].median()) if mask.any() else None,median_nFeature=float(meta.loc[mask,'n_genes_detected'].median()) if mask.any() else None,mean_epi_Q=float(orig.loc[mask,'Epithelial'].mean()) if mask.any() else None,epi_Q_weighted_UMI=float(np.dot(w,meta.library_counts)/den) if den>0 else None,epi_Q_weighted_nFeature=float(np.dot(w,meta.n_genes_detected)/den) if den>0 else None,n_measured_genes=info['n_measured_genes'],n_mp10_genes_measured=info['n_mp10_genes_measured'],hl_status=info['hl_status'],score_coverage_status=info['score_coverage_status'],median_nn_um=info['median_nn_um'],image_lattice_relative_rms=rel))
 tree=cKDTree(co[['x_um','y_um']].values)
 for radius in [105,200]:
  pairs=tree.query_pairs(radius+1e-6,output_type='ndarray');rows=np.concatenate([pairs[:,0],pairs[:,1]]);cols=np.concatenate([pairs[:,1],pairs[:,0]])
  reach=sparse.csr_matrix((np.ones(len(rows)),(rows,cols)),shape=(len(meta),len(meta)))
  for lev in ['major','fine']:
   q=qs[lev];v=q.values;D=np.outer(v.sum(axis=0),v.sum(axis=0))-v.T@v;R=v.T@(reach@v)
   for state in ['high','low']:
    epi='Epi_MP10_'+state;i=q.columns.get_loc(epi)
    for j,target in enumerate(q.columns):
     if target.startswith('Epi_MP10_'):continue
     geometry.append(dict(section_id=sec,dataset=info['dataset'],level=lev,epi_state=state,target=target,radius_um=radius,n_directed_reachable_spot_pairs=len(rows),n_all_directed_cross_spot_pairs=len(meta)*(len(meta)-1),all_cross_spot_Q_pair_mass=float(D[i,j]),reachable_cross_spot_Q_pair_mass=float(R[i,j]),reachable_fraction_of_full_Q_denominator=float(R[i,j]/D[i,j]) if D[i,j]>0 else None,interpretation='geometry-only potential weighted pair mass; no expression or inferred strength; distance graph symmetric'))
 tsv([r for r in quality if r['section_id']==sec],d/'quality_diagnostics.tsv')
 tsv([r for r in geometry if r['section_id']==sec],d/'geometry_denominator_diagnostics.tsv')
tsv(audit,O/'input_checks.tsv');tsv(quality,O/'spot_quality_diagnostics.tsv');tsv(geometry,O/'geometry_denominator_diagnostics.tsv');tsv(summaries,O/'section_status.tsv');tsv(summaries,O/'all_sections.tsv')
for name in ['04_coordinate_coverage_by_section.csv','06_topology_qc_by_section.csv','08_upstream_source_trace.md']:
 shutil.copyfile(lattice_source/name,O/name)
(O/'README.md').write_text((O/'README.md').read_text()+"\n图像仿射残差仅作诊断。权威阵列索引来自 surviving VisiumV1@coordinates$row/$col，经严格section+barcode连接；63片coverage/topology原审计均PASS。保留每片相对残差及warning，不以任意图像像素阈值删片。geometry_denominator_diagnostics.tsv为固定半径图下潜在Q配对质量，不是推断通讯。spot_quality_diagnostics.tsv提供Epi状态内UMI/nFeature及Q权重质量诊断。\n")
print('FINAL_INPUT_AUDIT',len(summaries),'READY',len(audit),'assertion checks',flush=True)

final_failures=[json.loads((O/'sections'/sec/'FAILED.json').read_text()) for sec in secs if (O/'sections'/sec/'FAILED.json').exists() and not (O/'sections'/sec/'READY.json').exists()]
missing=[sec for sec in secs if not (O/'sections'/sec/'READY.json').exists()]
js(dict(status='completed' if not missing else 'completed_with_blocked_inputs',n_total=len(secs),n_ready=len(summaries),n_failed=len(final_failures),failed_sections=final_failures,missing_ready=missing,n_network_only=sum(r['hl_status']=='network_only' for r in summaries),n_incomplete_mp10=sum(r['n_mp10_genes_measured']<20 for r in summaries),audit_checks_passed=len(audit),recovery='rerun scripts/prepare_inputs.py; existing READY sections reused',seconds=time.time()-start),O/'INPUT_PREPARATION_COMPLETED.json')

# Reapply user-authorized clinical metadata after any input-cache recovery.
if (R/'scripts/apply_workbook_annotations.py').exists():
 import subprocess
 subprocess.run([sys.executable,str(R/'scripts/apply_workbook_annotations.py'),'--inputs-only'],check=True)
