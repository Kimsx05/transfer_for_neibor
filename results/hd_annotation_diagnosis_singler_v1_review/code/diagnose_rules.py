from pathlib import Path
import pandas as pd,numpy as np,json
root=Path(__file__).resolve().parents[1];old=root.parent/'hd_major_annotation_v1'
classes=np.array(['Epithelial','T/NK','Myeloid','Fibroblast','Endothelial','Mast','B','Plasma','Mural'])
for section in ['Cervilla','D1_Johnson']:
 out=root/section;out.mkdir(exist_ok=True);t=old/section/'tables'
 if (out/'diagnosis.tsv.gz').exists():print('USING_EXISTING_DIAGNOSIS',section,flush=True);continue
 m=pd.read_csv(t/'analysis_metadata.tsv.gz',sep='\t',low_memory=False).set_index('cell_id');e=pd.read_csv(t/'marker_evidence.tsv.gz',sep='\t').set_index('cell_id');a=pd.read_csv(t/'annotations.tsv.gz',sep='\t').set_index('cell_id');assert m.index.equals(e.index) and m.index.equals(a.index)
 D=e[['detected_'+k for k in classes]].values.copy();M=e[['strength_'+k for k in classes]].values.copy();D[e.plasma_core_detected.values==0,7]=0;M[e.plasma_core_detected.values==0,7]=0
 top=pd.Index(classes).get_indexer(a.raw_prediction);rr=np.arange(len(m));best=M.argmax(1);strong=(D>=2)&(M>=M.max(1)[:,None]*.7);mesench=strong[:,3]&strong[:,8]&(strong.sum(1)==2);conflict=(strong.sum(1)>=2)&~mesench
 bands=pd.read_csv(t/'confidence_support_bands.tsv',sep='\t').set_index('major').loc[classes];refgood=(a.reference_score_max.values>=bands.score_q10.values[top])&(a.reference_score_margin.values>=bands.margin_q10.values[top]);coherent=a.cluster_reference_fraction.values;E=e[['cluster_enriched_detected_'+k for k in classes]].values
 def evaluate(remove=None):
  match=(best==top) if remove!='marker_agreement' else np.ones(len(m),bool);multi=(D[rr,top]>=2) if remove!='main_marker_count' else np.ones(len(m),bool);nc=~conflict if remove!='conflict' else np.ones(len(m),bool);rg=refgood if remove!='reference_band' else np.ones(len(m),bool)
  accept=match&multi&nc;single=(D[rr,top]==1)&match&nc&rg&(coherent>=.8)
  f=np.where(accept|single,classes[top],'Unknown');status=np.where(accept&rg,'supported',np.where(accept|single,'low_confidence','unknown'))
  override=(f=='Unknown')&nc&(D[rr,best]>=2)&(E[rr,best]>=2)&(best!=top)&(D[rr,top]==0);f[override]=classes[best[override]];status[override]='low_confidence';return f,status,accept,single,override
 f,status,accept,single,override=evaluate();assert np.array_equal(f,a.final_major) and np.array_equal(status,a.annotation_status);assert np.array_equal(conflict,a.marker_conflict)
 keep=['section_id','object_id','nCount_RNA','nFeature_RNA','n_bins','assigned_bin_area_um2','qc_keep_analysis','geometry_review','cluster_projected','sketch_projection_score','nucleus_x_um','nucleus_y_um','nucleus_x_fullres','nucleus_y_fullres','UMAP1','UMAP2','final_major','annotation_status','annotation_reason','raw_prediction','reference_score_max','reference_score_margin']
 z=m[keep].copy();z=z.rename(columns={'section_id':'section'});unknown=f=='Unknown'
 fail={'marker_agreement':best!=top,'main_marker_count':D[rr,top]<2,'conflict':conflict,'reference_band':~refgood}
 # Vectorized conjunction has no recorded short-circuit log. This is the explicit written order in agree=match & multi & ~conflict.
 first=np.select([best!=top,D[rr,top]<2,conflict],['marker_agreement','main_marker_count','conflict'],default='none')
 z['diag_first_rejection_reconstructed']=np.where(unknown,first,'not_finally_rejected');z['diag_reason_provenance']='original annotation_reason retained; first-failure and rule flags reconstructed from frozen code/caches';z['diag_original_reason']=a.annotation_reason
 z['diag_primary_accept']=accept;z['diag_single_marker_rescue']=single;z['diag_cluster_marker_revision']=override;z['diag_detected_predicted_markers']=D[rr,top];z['diag_detected_best_markers']=D[rr,best];z['diag_marker_top']=classes[best];z['diag_plasma_core_missing']=e.plasma_core_detected.values==0
 for k,v in fail.items():z['diag_fail_'+k]=v
 z['diag_n_main_rule_failures']=np.stack([fail[k] for k in ['marker_agreement','main_marker_count','conflict']],1).sum(1)
 z['diag_rule_overlap']=np.array(['|'.join(k for k in ['marker_agreement','main_marker_count','conflict'] if fail[k][i]) or 'none' for i in rr])
 z['diag_information_band']=pd.cut(z.nFeature_RNA,[-np.inf,19,49,99,199,499,np.inf],labels=['<20','20-49','50-99','100-199','200-499','>=500'])
 # Intermediate evidence table is compact and will be replaced by the combined delivery table after SingleR.
 z.to_csv(out/'diagnosis.tsv.gz',sep='\t',index_label='cell_id')
 rows=[]
 for scope,sel in [('all',np.ones(len(m),bool)),('original_Unknown',unknown)]:
  for k,v in fail.items():
   hit=sel&v;rows.append(dict(section=section,scope=scope,rule=k,n=int(hit.sum()),UMI=int(m.nCount_RNA.values[hit].sum()),fraction_scope=float(hit.sum()/sel.sum()),is_direct_main_rejection_rule=k!='reference_band'))
 pd.DataFrame(rows).to_csv(out/'rule_failure_counts.tsv',sep='\t',index=False)
 z[unknown].groupby('diag_first_rejection_reconstructed').agg(n=('object_id','size'),UMI=('nCount_RNA','sum')).to_csv(out/'first_rejection_reconstructed.tsv',sep='\t')
 z.groupby(['diag_rule_overlap','annotation_status']).agg(n=('object_id','size'),UMI=('nCount_RNA','sum')).to_csv(out/'rule_overlaps.tsv',sep='\t')
 rows=[]
 for k in fail:
  alt,st,_,_,_=evaluate(k);changed=alt!=f;recovered=unknown&(alt!='Unknown');rows.append(dict(section=section,removed_rule=k,labels_changed=int(changed.sum()),status_changed=int((st!=status).sum()),original_unknown_get_label=int(recovered.sum()),original_unknown_UMI_get_label=int(m.nCount_RNA.values[recovered].sum()),unknown_remaining=int((alt=='Unknown').sum()),thresholds='original cached support bands held fixed; diagnostic only'))
 pd.DataFrame(rows).to_csv(out/'single_rule_ablation.tsv',sep='\t',index=False)
 rows=[]
 for stage,keep in [('raw_prediction',np.ones(len(m),bool)),('main_marker_accept',accept),('plus_single_marker_rescue',accept|single),('plus_cluster_marker_revision',~unknown)]:rows.append(dict(stage=stage,labelled=int(keep.sum()),unknown=int((~keep).sum()),labelled_UMI=int(m.nCount_RNA.values[keep].sum())))
 pd.DataFrame(rows).to_csv(out/'rule_transitions.tsv',sep='\t',index=False)
 z.groupby(['diag_information_band','annotation_status'],observed=True).agg(n=('object_id','size'),UMI=('nCount_RNA','sum')).to_csv(out/'information_strata.tsv',sep='\t')
 rows=[]
 for state,g in z.groupby('annotation_status'):
  for field in ['nFeature_RNA','nCount_RNA','n_bins']:
   d=g[field].describe(percentiles=[.1,.25,.5,.75,.9]).to_dict();rows.append(dict(status=state,metric=field,**d))
 pd.DataFrame(rows).to_csv(out/'QC_distributions.tsv',sep='\t',index=False)
 cr=z.groupby('cluster_projected').agg(n=('object_id','size'),UMI=('nCount_RNA','sum'),unknown=('final_major',lambda x:(x=='Unknown').sum()),unknown_UMI=('nCount_RNA',lambda x:x[z.loc[x.index,'final_major']=='Unknown'].sum()),median_nFeature=('nFeature_RNA','median'),major_raw=('raw_prediction',lambda x:x.value_counts().index[0]),main_original_rejection=('annotation_reason',lambda x:x[z.loc[x.index,'final_major']=='Unknown'].value_counts().index[0] if (z.loc[x.index,'final_major']=='Unknown').any() else 'none'));cr['unknown_fraction']=cr.unknown/cr.n;cr.to_csv(out/'cluster_diagnosis.tsv',sep='\t')
 manifest=dict(section=section,objects=len(m),UMI=int(m.nCount_RNA.sum()),unknown=int(unknown.sum()),unknown_UMI=int(m.nCount_RNA.values[unknown].sum()),original_labels_and_status_reconstructed_exactly=True,GetTransferPredictions_used=False,score_filter=None,original_transfer_level='nine major classes directly, not subtype transfer',original_marker_assay='RNA data and counts, full analysis cells, not HVG/scale.data',first_failure_semantics='reconstructed written conjunction order: marker agreement, marker count, conflict; only final Unknown; not historical execution-time logging')
 (out/'diagnosis_manifest.json').write_text(json.dumps(manifest,indent=2));print(manifest,flush=True)
