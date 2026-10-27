#!/usr/bin/env python3
"""Read-only reuse of accepted spatial assets; no model fitting or old-file writes."""
import os,sys,json,time,gzip
from pathlib import Path
P=Path('/home/data/t070721/codex_workspace/Bladder Metabolism')
sys.path=[p for p in sys.path if 'dist-packages' not in p]
sys.path[:0]=[str(P/'spatial_ecology_top20_v2/runtime'),str(P/'cell2location_setup_v1/cell2location_env/lib/python3.12/site-packages')]
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import numpy as np,pandas as pd,h5py
from scipy import sparse,io
from scipy.spatial import cKDTree
R=P/'mp10_spatial_communication_v1';O=R/'input_review';O.mkdir(parents=True,exist_ok=True)
A=P/'spatial_ecology_top20_v2/00_inputs';secs=['GSE171351_A1','Inhousedata_T4469479'];start=time.time()
def tsv(df,p):pd.DataFrame(df).to_csv(p,sep='\t',index=False,na_rep='NA')
def js(obj,p):Path(p).write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n')
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
lowid=usage.index[usage.status.eq('LOW_IDENTIFIABILITY')].tolist()
index=[];pilot=[];sup=[];gencov=[];params=[];checks=[];cross=[];candidates=[]
# Candidate ranking frozen before any communication inference: valid-input complete20 and sufficient anchors; minimum size per required dataset.
gc=pd.read_csv(A/'gene_coverage.tsv',sep='\t');ga=pd.read_csv(A/'gate_audit.tsv',sep='\t').set_index('section_id')
for sec,r in manifest.iterrows():
 cov=gc[gc.section_id.eq(sec)&gc.MP.eq('EpiMP10')].iloc[0]
 candidates.append(dict(section_id=sec,dataset=r.dataset,n_spots=r.n_spots,n_epi_anchors=ga.loc[sec,'n_eligible'],mp10_genes_available=cov.n_available,raw_expression_exists=Path(r.raw_input_path).is_file(),mean_abundance_exists=Path(cm.loc[sec,'mean_path']).is_file(),selected=sec in secs,selection_basis='minimum spot count among complete20 input-complete public sections and separately Inhouse; preserve full QC-passed section'))
tsv(candidates,O/'candidate_sections.tsv')
for order,sec in enumerate(secs,1):
 t0=time.time();d=O/'sections'/sec;d.mkdir(parents=True,exist_ok=True);r=manifest.loc[sec];c=cm.loc[sec];identrow=ident.loc[sec]
 with h5py.File(r.raw_input_path,'r') as h:
  barcodes=strings(h['obs/_index']);genes=strings(h['var/_index']);Xg=h['X'];shape=tuple(Xg.attrs['shape']);X=sparse.csr_matrix((Xg['data'][()],Xg['indices'][()],Xg['indptr'][()]),shape=shape)
  assert str(Xg.attrs['encoding-type'])=='csr_matrix';xy=h['obsm/spatial'][()]
 assert X.shape==(len(barcodes),len(genes)) and len(set(genes))==len(genes) and len(set(barcodes))==len(barcodes)
 assert np.isfinite(X.data).all() and (X.data>=0).all() and (X.data==np.round(X.data)).all()
 s=spots[spots.section_id.eq(sec)].set_index('barcode').loc[barcodes].copy();assert np.max(np.abs(xy-s[['x','y']].values))<1e-10
 mean=pd.read_csv(c.mean_path,index_col=0,float_precision='round_trip');assert list(mean.index)==barcodes
 rawQ=mean.div(mean.sum(axis=1),axis=0);assert np.isfinite(rawQ.values).all();assert set(rawQ.columns)==set(mapping.index)
 assert np.max(abs(rawQ.Epithelial.values-s.epi_fraction_gate.values))<1e-12
 score=s.EpiMP10_Top20_AUCell;anchors=s.analysis_eligible_legacy.astype(bool);threshold=float(score[anchors].median());high=anchors&score.gt(threshold);low=anchors&score.lt(threshold);state=np.where(high,'high',np.where(low,'low','other'))
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
 assert abs(np.median(nn)-100)<1e-8 and rms<.05
 tsv(meta,d/'meta.tsv');tsv(coord,d/'coords.tsv');tsv([dict(gene=g) for g in genes],d/'genes.tsv');tsv([dict(spot_id=u,barcode=b) for u,b in zip(uid,barcodes)],d/'spots.tsv')
 for name,df in [('Q_original_30',rawQ),('Q_fine',q),('Q_major',qmajor),('abundance_mean_original_30',mean)]:
  z=df.copy();z.insert(0,'spot_id',uid);tsv(z.reset_index(drop=True),d/(name+'.tsv'))
 with gzip.open(d/'counts.mtx.gz','wb',compresslevel=1) as handle:io.mmwrite(handle,X.T,field='integer',comment='gene x spot; measured counts; original complete hard-QC section; ordering genes.tsv and meta.tsv')
 sparse.save_npz(d/'counts_spot_gene.npz',X)
 # Once-only measured expression normalization. Both methods may directly reuse this layer.
 total=np.asarray(X.sum(axis=1)).ravel();norm=X.astype(np.float64).multiply(10000/total[:,None]).tocsr();norm.data=np.log1p(norm.data);sparse.save_npz(d/'log1p_spot_gene.npz',norm)
 with gzip.open(d/'log1p.mtx.gz','wb',compresslevel=1) as handle:io.mmwrite(handle,norm.T,field='real',precision=9,comment='gene x spot; natural log1p(10000*count/total counts); full measured gene background')
 normerr=0.0
 for level,df in [('fine',q),('major',qmajor)]:
  for name in df:
   if name.startswith('Epi_MP10_'):
    ab=mean.Epithelial.values*(state==name.removeprefix('Epi_MP10_'));status='STATE_ASSOCIATED_SPOT_EPITHELIAL_COMPONENT'
   elif level=='major':ab=mean[mapping.index[mapping.major.eq(name)]].sum(axis=1).values;status='MAJOR_FROM_ACCEPTED_MEANS'
   elif name=='T_NK_unresolved':ab=mean[lowid].sum(axis=1).values;status='LOW_IDENTIFIABILITY_MASS_PRESERVED_UNRESOLVED'
   else:ab=mean[name].values;status=usage.loc[name,'status']
   eligible=bool((ab>=.01).sum()>=10 and ab.sum()>1)
   sup.append(dict(section_id=sec,dataset=r.dataset,level=level,node=name,n_positive_Q=int((df[name]>0).sum()),n_abundance_ge_0_01=int((ab>=.01).sum()),total_expected_abundance=float(ab.sum()),effective_spot_mass=float(df[name].sum()),max_Q=float(df[name].max()),n_dominant_spots=int((df.idxmax(axis=1)==name).sum()),n_dominant_original_30=int((rawQ.idxmax(axis=1)==name).sum()) if name in rawQ else np.nan,type_evaluable=eligible,support_rule='accepted validation: >=10 spots mean abundance>=0.01 and total abundance>1; no spot weights masked',identifiability_status=status,composition_weight='posterior_mean / sum_all_original_30_means'))
 for gene in genesig:gencov.append(dict(section_id=sec,dataset=r.dataset,gene_set='MP10_frozen_Top20',gene=gene,measured=gene in genes,detected_spots=int(X[:,genes.index(gene)].getnnz()) if gene in genes else np.nan,not_measured_is_zero=False))
 for original in rawQ:
  out='Epi_MP10_high|Epi_MP10_low|Epi_MP10_other' if original=='Epithelial' else ('T_NK_unresolved' if original in lowid else original)
  cross.append(dict(section_id=sec,original_factor=original,major=mapping.loc[original,'major'],analysis_fine_label=out,identifiability_status=usage.loc[original,'status'],downstream_usage=usage.loc[original,'downstream_usage'],reason=usage.loc[original,'reason'],mapping_source=str(A/'major_map.tsv'),reliability_source=str(P/'cell2location_final_validation_gate_v1/09_final_status/DOWNSTREAM_USAGE_TABLE.csv'),mass_action='split within epi only' if original=='Epithelial' else ('merge only three unidentifiable labels; no redistribution' if original in lowid else 'retained')))
 # Retain absent fine labels as NE in crosswalk, never invented columns.
 for absent,major in [('Tip EC','Endothelial'),('Cycling EC','Endothelial')]:cross.append(dict(section_id=sec,original_factor=absent,major=major,analysis_fine_label=np.nan,identifiability_status='NE_NOT_IN_SPATIAL_REFERENCE',downstream_usage='NE',reason='No independently estimated composition column in accepted 30-factor reference',mass_action='not fabricated'))
 nsum=dict(section_id=sec,dataset=r.dataset,pilot_order=order,n_spots=len(s),n_measured_genes=len(genes),n_mp10_genes_measured=sum(g in genes for g in genesig),n_epi_anchors=int(anchors.sum()),n_epi_high=int(high.sum()),n_epi_low=int(low.sum()),n_epi_equal_median=int((anchors&score.eq(threshold)).sum()),n_epi_other_spots=int((state=='other').sum()),threshold=threshold,n_original_factors=30,n_fine_columns=len(q.columns),n_major_columns=len(qmajor.columns),observational_unit='mixed Visium spot',coordinate_units='um; nominal platform array geometry',median_nn_um=float(np.median(nn)),median_nn_image_px=float(np.median(pixnn)),image_lattice_affine_rms_px=rms,patient_uid=identrow.patient_uid,patient_id_status=identrow.patient_id_status,analysis_unit_uid=identrow.analysis_unit_uid,selection_basis='input complete; full20 MP10; minimum computational size per selected dataset; no communication outcome used',prepared_dir=str(d),raw_input_path=r.raw_input_path,mean_path=c.mean_path)
 pilot.append(nsum);js({k:None if isinstance(v,float) and np.isnan(v) else v for k,v in nsum.items()},d/'input_summary.json')
 for key,value,unit,source in [('nominal_spot_diameter',55,'um','standard Visium manufacturer platform specification; not inferred from invalid spot scale_factor'),('nominal_center_pitch',100,'um','actual recovered ordinary Visium lattice + manufacturer platform specification'),('secreted_range',200,'um','predeclared 2 x measured d_NN'),('contact_range',105,'um','predeclared nearest-neighbor d_NN + 5% tolerance; does not prove membrane contact'),('ecm_range',200,'um','predeclared 2 x d_NN; SpatialCellChat native ECM shares interaction.range with secreted; fixed before first inference results'),('first_section_secreted_sensitivity',300,'um','3 x d_NN, first selected section only; SpatialCellChat changes ECM and secreted together by native API; COMMOT changes secreted only'),('expression_normalization_target',10000,'counts per spot','full measured gene library normalization once; log1p'),('epi_anchor_threshold',.25,'fraction','frozen accepted c2l mean epithelial / total mean'),('MP10_HL_threshold',threshold,'frozen AUCell','new per-section anchors-only strict median; no overwritten original score'),('spot_weight_threshold',0,'Q','none: all positive posterior means preserved; no 0.1 cutoff')]:params.append(dict(section_id=sec,parameter=key,value=value,unit=unit,source=source))
 for role,path in [('raw_counts_expression',r.raw_input_path),('posterior_mean_abundance',c.mean_path),('frozen_score',r.score_path),('frozen_spots_anchors_and_lattice',A/'frozen_spot_inputs.parquet'),('patient_identity',ids_path),('accepted_major_map',A/'major_map.tsv')]:index.append(dict(section_id=sec,dataset=r.dataset,role=role,path=str(path),exists=Path(path).is_file(),read_only=True,observation_unit='spot',platform='ordinary_Visium',measured_expression_layer='X raw counts',coordinate_reference='true array row/col; original lowres image pixels also preserved',patient_id_status=identrow.patient_id_status,id_alignment='exact barcode bijection; expression,coordinates,score,abundance same order'))
 checks.extend([dict(section_id=sec,check='barcode_bijection_and_order',value=len(s),pass_check=True),dict(section_id=sec,check='fine_Q_mass_max_error',value=float(abs(q.sum(axis=1)-rawQ.sum(axis=1)).max()),pass_check=True),dict(section_id=sec,check='major_Q_mass_max_error',value=float(abs(qmajor.sum(axis=1)-rawQ.sum(axis=1)).max()),pass_check=True),dict(section_id=sec,check='TME_shared_Q_unchanged',value=float(abs(q['Treg']-rawQ['Treg']).max()),pass_check=True),dict(section_id=sec,check='measured_counts_nonnegative_integer',value=X.nnz,pass_check=True),dict(section_id=sec,check='nominal_nn_um',value=float(np.median(nn)),pass_check=True),dict(section_id=sec,check='lattice_image_affine_rms_px',value=rms,pass_check=True)])
 print(sec,'READY',nsum,'seconds',time.time()-t0,flush=True)
for name,rows in [('input_index',index),('pilot_sections',pilot),('spatial_label_crosswalk',cross),('type_support',sup),('gene_coverage',gencov),('parameters',params),('input_checks',checks)]:tsv(rows,O/(name+'.tsv'))
js(dict(status='READY',sections=secs,seconds=time.time()-start,blas_threads=1,new_annotations=False,new_c2l_training=False),O/'PREPARED.json')
