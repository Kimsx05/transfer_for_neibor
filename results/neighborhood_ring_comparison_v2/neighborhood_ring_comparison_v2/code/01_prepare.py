from common import *
import importlib.util,platform,shutil
import matplotlib.pyplot as plt
from scipy.sparse.csgraph import connected_components
A=ROOT/'neighborhood_v1/00_inputs_and_audit';C=ROOT/'spatial_ecology_top20_v2/03_copykat_geography'
spec=importlib.util.spec_from_file_location('hex',ROOT/'neighborhood_v1/code/visium_hex_graph.py');hex=importlib.util.module_from_spec(spec);spec.loader.exec_module(hex)

def rings(coords):
 r1,r12=hex.build_hex_graph(coords);r2=r12-r1;r2.eliminate_zeros();assert not r1.multiply(r2).nnz;assert not r2.diagonal().any();return r1,r2

def main():
 if (OUT/'03_features/COMPLETE.json').exists():log('M00-M03 checkpoint exists');return
 contract=json.loads((ROOT/'spatial_ecology_top20_v2/05_07_shared_K12/v1/analysis_contract.json').read_text());types=contract['factor_order'];mapping=pd.read_csv(A/'metadata/subtype_to_major_mapping.csv');tab(OUT/'00_preflight/type_mapping.tsv',mapping)
 f=pd.read_csv(A/'frozen_input_manifest.csv');f=f[f.layer.eq('subtype_q05')].sort_values('section');assert len(f)==63
 identity=pd.read_csv(ROOT/'spatial_ecology_top20_v2/05_07_shared_K12/v1/section_identity_overlay.tsv',sep='\t').set_index('section_id');assert identity.index.is_unique
 tab(OUT/'00_preflight/section_identity.tsv',identity.reset_index())
 lattice=pd.read_csv(OUT/'00_preflight/full_image_tissue_coordinates.csv.gz');lattice=lattice[lattice.section.isin(f.section)&lattice.tissue.eq(1)].copy();legacy=pd.read_csv(A/'metadata/spot_manifest.csv').set_index(['section','barcode'])
 cnv=pd.read_parquet(C/'spot_geography.parquet',columns=['section_id','barcode','x','y','native_copykat_label','profile_available','spot_MSQ','baseline_confidence','mixture_nonconvergence','histology_region','region_consensus','nCount','nFeature']).set_index(['section_id','barcode']);assert cnv.index.is_unique
 assert set(f.section)==set(identity.index)==set(lattice.section.unique())
 coords=pd.DataFrame([(r,c) for r in range(-4,5) for c in range(-8,9) if (r+c)%2==0],columns=['array_row','array_col']);g1,g2=rings(coords);ci=coords.index[(coords.array_row==0)&(coords.array_col==0)][0];assert g1[ci].sum()==6 and g2[ci].sum()==12
 # Missing observation leaves graph intact; absent tissue node removes only that node's paths.
 missing=np.ones(len(coords));missing[g1[ci].indices[0]]=0;assert (g1@missing)[ci]==5
 hole=coords.drop(g1[ci].indices[0]).reset_index(drop=True);h1,h2=rings(hole);assert h1.max()<=1 and h2.max()<=1
 js(OUT/'02_graphs/implementation_checks.json',dict(hex_6_12=True,disjoint=True,center_excluded=True,missing_observation_retains_topology=True,no_knn_gap_bridging=True,per_section_graphs=True))
 audits=[];allmeta=[];qs1=[];qs2=[];rois=[];errors=[];assets=[]
 for row in f.itertuples():
  s=row.section;ds=str(identity.loc[s,'dataset']);log('prepare '+s)
  dest=OUT/'01_regions'/ds;dest.mkdir(exist_ok=True)
  co=lattice[lattice.section.eq(s)].set_index('barcode').sort_index();assert co.index.is_unique
  qpath=Path(row.path);actual=sha(qpath);assert actual==row.sha256
  raw=pd.read_csv(qpath,index_col=0);raw=raw.drop(columns=['section','dataset'],errors='ignore');assert set(raw.columns)==set(types);raw=raw.reindex(index=co.index,columns=types)
  v=raw.to_numpy(float);valid=np.isfinite(v).all(1)&(v>=0).all(1)&(v.sum(1)>0)
  for i in np.flatnonzero(~valid):errors.append(dict(section_id=s,barcode=co.index[i],reason='missing_nonfinite_negative_or_zero_total'))
  p=np.zeros_like(v);p[valid]=v[valid]/v[valid].sum(1,keepdims=True)
  mp=ROOT/'cell2location_full_run/v2_30type/02_sections'/s/'abundance_mean.csv.gz';mean=pd.read_csv(mp,index_col=0);mean.columns=[c.replace('meanscell_abundance_w_sf_','').replace('means_cell_abundance_w_sf_','') for c in mean.columns];assert set(types).issubset(mean.columns),mean.columns
  mean=mean.reindex(index=co.index,columns=types);mv=mean.to_numpy(float);mvok=np.isfinite(mv).all(1)&(mv>=0).all(1)&(mv.sum(1)>0);epi=np.full(len(co),np.nan);epi[mvok]=mv[mvok,types.index('Epithelial')]/mv[mvok].sum(1)
  r1,r2=rings(co);n1=np.asarray(r1.sum(1)).ravel();n2=np.asarray(r2.sum(1)).ravel();nv1=r1@valid.astype(int);nv2=r2@valid.astype(int)
  q1=np.full_like(p,np.nan);q2=np.full_like(p,np.nan);np.divide(r1@p,nv1[:,None],out=q1,where=nv1[:,None]>0);np.divide(r2@p,nv2[:,None],out=q2,where=nv2[:,None]>0)
  common=valid&mvok&(epi>=.25)&(nv1>=4)&(nv2>=8);full=(n1==6)&(n2==12)&(nv1==6)&(nv2==12)
  ev=cnv.loc[s].reindex(co.index);prediction=ev.native_copykat_label.fillna('missing').astype(str)
  region=np.full(len(co),'uncertain_or_conflicting',dtype=object);region[(epi<.25)&~prediction.eq('aneuploid')]='non_epithelial_region';region[(epi>=.25)&prediction.eq('aneuploid')]='tumor_supported_epithelium'
  uid=ds+'||'+s+'||'+co.index.astype(str)
  meta=pd.DataFrame(dict(spot_uid=uid,dataset_id=ds,section_id=s,barcode=co.index,patient_uid=identity.loc[s,'patient_uid'],patient_identity_status=identity.loc[s,'patient_id_status'],source_patient_identity_known=identity.loc[s,'source_patient_identity_known'],array_row=co.array_row.to_numpy(),array_col=co.array_col.to_numpy(),x=ev.x.to_numpy(),y=ev.y.to_numpy(),epithelial_fraction=epi,epithelial_fraction_source='v2_30type posterior mean; Epithelial / all30',composition_valid=valid,center_qc=mvok,ring1_geometric_n=n1,ring2_geometric_n=n2,ring1_valid_n=nv1,ring2_valid_n=nv2,ring1_coverage=nv1/6,ring2_coverage=nv2/12,complete_rings=full,common_anchor=common,region_auto=region,region_final='',region_review_status='pending_review',region_version='cnv_candidate_v1',boundary_status='unknown_pending_region_review',boundary_type='unknown',depth_valid=False,depth_hops=np.nan,truncation_reason=np.where((n1<6)|(n2<12),'capture_edge_or_tissue_gap','none_detected'),copykat_prediction_original=prediction.to_numpy(),cnv_score_original=ev.spot_MSQ.to_numpy(),cnv_score_definition_status='mean(native_final_CNA_profile_squared); not log-ratio; scale parity unverified',cnv_profile_available=ev.profile_available.to_numpy(),pathology_label_existing=ev.histology_region.fillna('unknown').to_numpy(),pathology_evidence_source='no accepted region annotation in Module03',nCount=ev.nCount.to_numpy(),nFeature=ev.nFeature.to_numpy(),baseline_confidence=ev.baseline_confidence.to_numpy(),mixture_nonconvergence=ev.mixture_nonconvergence.to_numpy()))
  meta['physical_position_source']='Seurat VisiumV1 image coordinates; tissue==1, independently of composition QC';meta['patient_stratum']=meta.patient_uid.fillna('patient_unresolved');meta['main_inclusion_reason']=np.where(common,'included',np.where(epi<.25,'below_epi_gate',np.where(~valid|~mvok,'invalid_composition','insufficient_ring_coverage')))
  meta['old_RCTD_anchor']=legacy.loc[(s,),:].reindex(co.index).is_epi025_anchor.fillna(False).astype(bool).to_numpy()
  meta['cnv_candidate_tumor_evidence']=prediction.eq('aneuploid').to_numpy() # includes low-Epi spots; no Epi tumor contour
  meta['candidate_inner_edge_to_unresolved']=meta.cnv_candidate_tumor_evidence.to_numpy()&((r1@(~meta.cnv_candidate_tumor_evidence.to_numpy()).astype(int))>0)
  # ROI = connected evidence-consistent components, preserving small islands.
  meta['roi_id']=''
  for label in sorted(set(region)):
   ix=np.flatnonzero(region==label);nc,cc=connected_components(r1[ix][:,ix],directed=False)
   for k in range(nc):
    ids=ix[cc==k];rid=f'{s}__{label}__{k+1:04d}';meta.loc[ids,'roi_id']=rid
    question=('上皮与非整倍体证据相符；形态是否支持肿瘤区域，边缘是否受缺口影响？' if label=='tumor_supported_epithelium' else '是否正常上皮、肿瘤或其他组织？diploid/缺失CNV不能单独裁定正常。' if label=='uncertain_or_conflicting' else '低上皮区域的组织身份是什么；是否仍含肿瘤证据？')
    rois.append(dict(roi_id=rid,dataset_id=ds,section_id=s,n_spots=len(ids),n_common_anchors=int(common[ids].sum()),candidate_label=label,array_row_min=int(co.array_row.iloc[ids].min()),array_row_max=int(co.array_row.iloc[ids].max()),array_col_min=int(co.array_col.iloc[ids].min()),array_col_max=int(co.array_col.iloc[ids].max()),mean_epi=float(np.nanmean(epi[ids]) if np.isfinite(epi[ids]).any() else np.nan),aneuploid_n=int(prediction.iloc[ids].eq('aneuploid').sum()),diploid_n=int(prediction.iloc[ids].eq('diploid').sum()),edge_or_gap_n=int(((n1<6)|(n2<12))[ids].sum()),priority='high' if label=='uncertain_or_conflicting' else 'small_island' if len(ids)<=3 else 'audit',evidence_figure=str(dest/(s+'__evidence.png')),spot_evidence=str(dest/(s+'__regions.parquet')),question=question,impact='改变区域分层/候选界面；不改变共同anchor或六方案标签',reviewer='',review_date='',reviewed_label='',review_comment=''))
  meta.to_parquet(dest/(s+'__regions.parquet'),index=False)
  gd=OUT/'02_graphs'/ds;gd.mkdir(exist_ok=True);sparse.save_npz(gd/(s+'__R1.npz'),r1);sparse.save_npz(gd/(s+'__R2.npz'),r2);tab(gd/(s+'__nodes.tsv'),meta[['spot_uid','barcode','array_row','array_col']])
  fig,axes=plt.subplots(1,4,figsize=(20,5));xx=co.array_col.to_numpy();yy=co.array_row.to_numpy()*np.sqrt(3)
  panels=[(epi,'Epi fraction (mean)', 'viridis'),(prediction.map({'aneuploid':2,'diploid':1}).fillna(0),'CopyKAT: 0 unknown / 1 diploid / 2 aneuploid','viridis'),(ev.spot_MSQ,'CNA mean square; descriptive','magma'),(pd.Series(region).map({'uncertain_or_conflicting':0,'non_epithelial_region':1,'tumor_supported_epithelium':2}),'Candidate: 0 uncertain / 1 non-Epi / 2 tumor-supported','viridis')]
  for ax,(values,title,cmap) in zip(axes,panels):
   sc=ax.scatter(xx,yy,c=values,s=3,cmap=cmap);ax.set_aspect('equal');ax.invert_yaxis();ax.set_title(title,fontsize=9);fig.colorbar(sc,ax=ax,shrink=.6);ax.set_xlabel('array_col');ax.set_ylabel('sqrt(3) array_row')
  fig.suptitle(f'{ds} / {s} | common anchors={common.sum()} | cnv_candidate_v1 | NOT pathology-confirmed');fig.tight_layout();fig.savefig(dest/(s+'__evidence.png'),dpi=130);plt.close(fig)
  allmeta.append(meta);qs1.append(q1[common].astype('float32'));qs2.append(q2[common].astype('float32'))
  sp=C/'by_dataset'/ds/s;profiles=list(sp.glob('final_CNA_profile_original_barcodes.npz'));assets.append(dict(section_id=s,dataset_id=ds,prediction_path=str(sp/'native_prediction.tsv'),prediction_exists=(sp/'native_prediction.tsv').exists(),score_path=str(C/'spot_geography.parquet'),profile_container=str(profiles[0]) if profiles else '',profile_container_exists=bool(profiles),profile_not_loaded_this_round=True,spots_profile_available=int(ev.profile_available.fillna(False).sum())))
  audits.append(dict(section_id=s,dataset_id=ds,n_tissue_nodes=len(co),n_q05=int(valid.sum()),n_missing_composition=int((~valid).sum()),n_valid=int(valid.sum()),n_epi=int((epi>=.25).sum()),n_common=int(common.sum()),n_complete=int((full&common).sum()),old_RCTD_n=int(meta.old_RCTD_anchor.sum()),gate_overlap=int((meta.old_RCTD_anchor&(epi>=.25)).sum()),q05_path=str(qpath),q05_sha256=actual,mean_path=str(mp),mean_sha256=sha(mp),q05_manifest_match=True,status='PASS'))
 master=pd.concat(allmeta,ignore_index=True);assert master.spot_uid.is_unique;anchors=master[master.common_anchor].reset_index(drop=True);q1=np.concatenate(qs1);q2=np.concatenate(qs2);assert len(q1)==len(anchors);assert np.allclose(q1.sum(1),1) and np.allclose(q2.sum(1),1)
 ref,alloc=sample_hierarchy(anchors,20000,123);anchors['reference']=False;anchors.loc[ref,'reference']=True;master['reference']=master.spot_uid.isin(anchors.loc[ref,'spot_uid']);master.to_parquet(OUT/'03_features/anchor_master.parquet',index=False);anchors.to_parquet(OUT/'03_features/common_anchors.parquet',index=False)
 np.savez_compressed(OUT/'03_features/features.npz',q1=q1,q2=q2,X1=np.sqrt(q1),X2=np.sqrt(q2),X12=np.concatenate([np.sqrt(q1),np.sqrt(q2)],axis=1)/np.sqrt(2),reference=ref)
 tab(OUT/'03_features/reference_anchor_ids.tsv',anchors.loc[ref,['spot_uid','dataset_id','patient_uid','section_id']]);tab(OUT/'03_features/reference_allocation.tsv',alloc);tab(OUT/'00_preflight/coverage63.tsv',audits);tab(OUT/'00_preflight/composition_anomalies.tsv',errors);tab(OUT/'01_regions/cnv_assets.tsv',assets);tab(OUT/'01_regions/ROI_review.tsv',rois)
 js(OUT/'03_features/feature_columns.json',dict(q1=types,q2=types,X1=types,X2=types,X12=['R1_'+x for x in types]+['R2_'+x for x in types]))
 config=dict(summary='frozen neighborhood_v1 subtype_q05 regenerated from accepted v2_30type, 1000 draws',gate='accepted newer c2l posterior mean Epithelial/all30 >=0.25',user_confirmation='2026-10-12 q05 environment + c2l mean gate',ring1_min=4,ring2_min=8,seed=123,reference_n=len(ref),factor_order=types,high_epi_threshold=.7,high_epi_basis='prospective composition-only sensitivity; not malignant cutoff',high_epi_n=int((anchors.epithelial_fraction>=.7).sum()),region_status='pending_review',MP10_values_loaded=False)
 js(OUT/'00_preflight/input_contract.json',config)
 js(OUT/'03_features/COMPLETE.json',dict(n_sections=len(audits),n_datasets=anchors.dataset_id.nunique(),n_source_patient_uids=anchors.loc[anchors.source_patient_identity_known,'patient_uid'].nunique(),n_common=len(anchors),n_reference=len(ref),n_full=int(anchors.complete_rings.sum()),n_roi=len(rois),sha_features=sha(OUT/'03_features/features.npz'),sha_common=sha(OUT/'03_features/common_anchors.parquet')))
 log('M00-M03 complete')
if __name__=='__main__':main()
