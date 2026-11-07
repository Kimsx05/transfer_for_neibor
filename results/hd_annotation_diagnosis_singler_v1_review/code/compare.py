from pathlib import Path
import pandas as pd,numpy as np,json,itertools
root=Path(__file__).resolve().parents[1];old=root.parent/'hd_major_annotation_v1';classes=['Epithelial','T/NK','Myeloid','Fibroblast','Endothelial','Mast','B','Plasma','Mural'];hpclasses=['Epithelial','T/NK','Myeloid','Fibroblast','Endothelial','B']
summ=[]
for s in ['Cervilla','D1_Johnson']:
 out=root/s;z=pd.read_csv(out/'diagnosis.tsv.gz',sep='\t',low_memory=False).set_index('cell_id');e=pd.read_csv(old/s/'tables/marker_evidence.tsv.gz',sep='\t').set_index('cell_id');assert z.index.equals(e.index);n=len(z);umi=z.nCount_RNA.sum();unk=z.final_major.eq('Unknown');B=z.singler_blca_pruned_labels;H=z.singler_hpca_pruned_labels;bp=B.notna();hp=H.notna()
 z['diag_candidate_pattern']=np.select([bp&hp&B.eq(H),bp&hp&~B.eq(H),bp&~hp,~bp&hp],['both_nonpruned_agree','both_nonpruned_disagree','BLCA_only_nonpruned','HPCA_only_nonpruned'],default='neither_nonpruned')
 z['diag_hpca_original_coverage']=np.where(z.final_major.eq('Unknown'),'unknown_original_identity',np.where(z.final_major.isin(hpclasses),'covered','not_covered'))
 z['diag_hpca_blca_candidate_coverage']=np.where(~bp,'no_nonpruned_BLCA_candidate',np.where(B.isin(hpclasses),'covered','not_covered'))
 D=e[['detected_'+k for k in classes]].values;M=e[['strength_'+k for k in classes]].values;rr=np.arange(n)
 effectiveD=D.copy();effectiveM=M.copy();effectiveD[e.plasma_core_detected.values==0,7]=0;effectiveM[e.plasma_core_detected.values==0,7]=0
 z['diag_any_panel_marker_detected']=effectiveD.max(1)>0;z['diag_any_group_two_markers']=effectiveD.max(1)>=2;z['diag_all_marker_strength_zero']=effectiveM.max(1)==0;z['diag_marker_top_tied']=(effectiveM==effectiveM.max(1)[:,None]).sum(1)>1
 for route,pred in [('blca',B),('hpca',H)]:
  ix=pd.Index(classes).get_indexer(pred);valid=ix>=0;det=np.full(n,np.nan);strength=np.full(n,np.nan);det[valid]=D[rr[valid],ix[valid]];strength[valid]=M[rr[valid],ix[valid]]
  other=M.copy();other[rr[valid],ix[valid]]=-np.inf;otherbest=other.argmax(1);otherdet=D[rr,otherbest];otherstrong=(otherdet>=2)&(other.max(1)>=strength)&valid
  z['diag_'+route+'_candidate_marker_count']=det;z['diag_'+route+'_other_strong_marker_group']=np.where(valid,np.array(classes)[otherbest],None);z['diag_'+route+'_other_group_supported']=otherstrong
  z['diag_'+route+'_marker_evidence']=np.select([~valid,otherstrong,det>=2,det==1,det==0],['no_nonpruned_candidate','other_lineage_group_at_least_as_strong','two_or_more_candidate_markers','one_candidate_marker','no_candidate_marker'],default='missing')
  z['diag_'+route+'_conflicts_with_original_known']=valid&~unk&pred.ne(z.final_major)
  # These are evidence descriptions only, never acceptance criteria or automatic revisions.
 z.to_csv(out/'diagnosis.tsv.gz',sep='\t',index_label='cell_id')
 methods={'Seurat_raw':z.raw_prediction,'Original_final':z.final_major,'BLCA_raw':z.singler_blca_labels.fillna('Not_assessable'),'BLCA_pruned':B.fillna('Unknown'),'HPCA_raw':z.singler_hpca_labels.fillna('Not_assessable'),'HPCA_pruned':H.fillna('Unknown')}
 rows=[]
 for method,lab in methods.items():
  tmp=pd.DataFrame({'label':lab,'UMI':z.nCount_RNA.values});g=tmp.groupby('label').agg(n=('UMI','size'),UMI=('UMI','sum')).reset_index();g['method']=method;g['fraction_all']=g.n/n;g['UMI_fraction_all']=g.UMI/umi;rows.append(g)
 pd.concat(rows).to_csv(out/'method_composition.tsv',sep='\t',index=False)
 rows=[]
 for (ka,va),(kb,vb) in itertools.combinations(methods.items(),2):
  tmp=pd.DataFrame({'from_label':va,'to_label':vb,'UMI':z.nCount_RNA});g=tmp.groupby(['from_label','to_label']).agg(n=('UMI','size'),UMI=('UMI','sum')).reset_index();g['from_method']=ka;g['to_method']=kb;g['fraction_all']=g.n/n;g['UMI_fraction_all']=g.UMI/umi;g['fraction_source']=g.n/g.groupby('from_label').n.transform('sum');rows.append(g)
 pd.concat(rows).to_csv(out/'label_crosstabs.tsv',sep='\t',index=False)
 rows=[]
 for scope,mask in [('all',np.ones(n,bool)),('original_Unknown',unk.values),('original_known',~unk.values)]:
  for k,g in z[mask].groupby(['diag_candidate_pattern','diag_hpca_blca_candidate_coverage'],dropna=False):rows.append(dict(scope=scope,pattern=k[0],hpca_blca_coverage=k[1],n=len(g),UMI=int(g.nCount_RNA.sum()),fraction_scope=len(g)/sum(mask),UMI_fraction_scope=g.nCount_RNA.sum()/z.nCount_RNA.values[mask].sum()))
 pd.DataFrame(rows).to_csv(out/'candidate_patterns.tsv',sep='\t',index=False)
 rows=[]
 for route in ['blca','hpca']:
  for (label,ev),g in z[unk].groupby([f'singler_{route}_pruned_labels',f'diag_{route}_marker_evidence'],dropna=False):rows.append(dict(route=route,label=label,evidence=ev,n=len(g),UMI=int(g.nCount_RNA.sum()),median_nFeature=g.nFeature_RNA.median(),median_nCount=g.nCount_RNA.median(),median_bins=g.n_bins.median()))
 pd.DataFrame(rows).to_csv(out/'unknown_candidate_marker_evidence.tsv',sep='\t',index=False)
 z[unk].groupby(['cluster_projected','singler_blca_pruned_labels','singler_hpca_pruned_labels'],dropna=False).agg(n=('object_id','size'),UMI=('nCount_RNA','sum'),median_nFeature=('nFeature_RNA','median')).to_csv(out/'unknown_candidates_by_cluster.tsv',sep='\t')
 rows=[]
 for band,g in z.groupby('diag_information_band',observed=True):
  u=g.final_major.eq('Unknown');rows.append(dict(information_band=band,n=len(g),original_unknown=int(u.sum()),BLCA_nonpruned=int(g.singler_blca_pruned_labels.notna().sum()),HPCA_nonpruned=int(g.singler_hpca_pruned_labels.notna().sum()),unknown_BLCA_candidate=int((u&g.singler_blca_pruned_labels.notna()).sum()),unknown_HPCA_candidate=int((u&g.singler_hpca_pruned_labels.notna()).sum())))
 pd.DataFrame(rows).to_csv(out/'information_candidate_counts.tsv',sep='\t',index=False)
 item=dict(section=s,n=n,UMI=int(umi),original_unknown=int(unk.sum()),unknown_UMI=int(z.loc[unk,'nCount_RNA'].sum()),unknown_UMI_fraction=float(z.loc[unk,'nCount_RNA'].sum()/umi),unknown_BLCA_nonpruned=int((unk&bp).sum()),unknown_HPCA_nonpruned=int((unk&hp).sum()),unknown_both_agree=int((unk&bp&hp&B.eq(H)).sum()),unknown_both_disagree=int((unk&bp&hp&~B.eq(H)).sum()),unknown_neither_nonpruned=int((unk&~bp&~hp).sum()),unknown_BLCA_not_HPCA_covered=int((unk&bp&~B.isin(hpclasses)).sum()))
 for route,pred in [('blca',B),('hpca',H)]:
  hit=unk&pred.notna();item['unknown_'+route+'_candidate_UMI']=int(z.loc[hit,'nCount_RNA'].sum());item['unknown_'+route+'_two_marker_support']=int((hit&(z['diag_'+route+'_candidate_marker_count']>=2)).sum());item[route+'_known_conflict']=int(z['diag_'+route+'_conflicts_with_original_known'].sum())
 (out/'comparison_summary.json').write_text(json.dumps(item,indent=2));summ.append(item);print(item,flush=True)
pd.DataFrame(summ).to_csv(root/'summary.tsv',sep='\t',index=False)
