from pathlib import Path
import pandas as pd,numpy as np,json,sys
r=Path(__file__).resolve().parents[1];s=sys.argv[1];o=r/s;b=r.parent/'hd_bin8_comparison_v1'/s
su=pd.read_csv(o/'cluster_summary.tsv',sep='\t',dtype={'cluster':str});mk=pd.read_csv(o/'cluster_markers.tsv.gz',sep='\t',dtype={'cluster':str});sr=pd.read_csv(o/'cluster_SingleR.tsv',sep='\t',dtype={'cluster':str});comp=pd.read_csv(o/'cluster_reference_composition.tsv',sep='\t',dtype={'cluster':str});over=json.loads((o/'reviewed_region_overrides.json').read_text()) if (o/'reviewed_region_overrides.json').exists() else {}
rows=[]
for _,z in su.iterrows():
 get=lambda k:float(z.get('Plasma_core_supported_ge2_fraction',0) if k=='Plasma' else z.get('detected_'+k+'_ge2_fraction',0))
 strength={'Epithelial':get('Epithelial'),'Mesenchymal':max(get('Fibroblast'),get('Mural')),'Immune':max(get(k) for k in ['T/NK','Myeloid','B','Plasma','Mast']),'Endothelial':get('Endothelial')};best=max(strength,key=strength.get)
 sub=comp[(comp.branch==z.branch)&(comp.cluster==z.cluster)&(comp.level=='coarse')]
 tops={ref:sub[sub.reference==ref].sort_values('n',ascending=False).iloc[0].label for ref in ['blca','hpca']}
 if z.nFeature_median<20:
  label='Low-information / unresolved';reason='Median nFeature <20 (descriptive information state; no bins removed)'
 elif get('Epithelial')>=.35 and get('Fibroblast')>=.35 and z.epi_fib_ge2_fraction>=.2:
  label='Epithelial-mesenchymal mixed signal';reason='Both original panels supported across bins and co-detected within >=20% of bins; not automatically a boundary or hybrid cell type'
 elif get('Epithelial')>=.35 and get('Fibroblast')>=.35:
  label='Epithelial-stromal composite region';reason='Both original panels frequent across the cluster but fewer bins co-detect both; juxtaposed components, not a uniform cell type'
 elif strength[best]>=.2:
  label=best+'-enriched';reason='Largest original panel detection support, supported in >=20% of bins; cluster-level enrichment, not uniform bin identity'
 elif tops['blca']==tops['hpca'] and tops['blca']!='Unknown' and strength[tops['blca']]>=.05:
  label=tops['blca']+'-associated (sparse support)';reason='Concordant dominant coarse reference composition with sparse original marker support; tentative regional interpretation'
 else:
  label='Unresolved expression region';reason='Sparse or competing original marker evidence; references do not justify a uniform region type'
 key=f'{z.branch}:{z.cluster}';review='evidence_summary_only'
 if key in over:label=over[key]['interpretation'];reason=over[key]['reason'];review='reviewed_markers_and_spatial'
 top=mk[(mk.branch==z.branch)&(mk.cluster==z.cluster)].head(12).symbol.astype(str).str.cat(sep=', ')
 zz=sr[(sr.branch==z.branch)&(sr.cluster==z.cluster)];pred={q.route:q.label for _,q in zz.iterrows()}
 rows.append(dict(branch=z.branch,cluster=z.cluster,interpretation=label,interpretation_reason=reason,review_status=review,n=int(z.n),UMI=int(z.UMI),nFeature_median=z.nFeature_median,low_info_fraction=z.low_info_fraction,edge16_fraction=z.edge16_fraction,long_neighbor80_fraction=z.long_neighbor80_fraction,epi_fib_ge2_fraction=z.epi_fib_ge2_fraction,dominant_BLCA_coarse=tops['blca'],dominant_HPCA_coarse=tops['hpca'],cluster_SingleR_BLCA=pred.get('blca','not_tested'),cluster_SingleR_HPCA=pred.get('hpca','not_tested'),top_original_markers=top,**{f'{k}_panel_fraction':v for k,v in strength.items()}))
regions=pd.DataFrame(rows);regions.to_csv(o/'cluster_domain_interpretation.tsv',sep='\t',index=False)
d=pd.read_csv(o/'metadata.tsv.gz',sep='\t',dtype={f'banksy_{t}_cluster':str for t in ['l00','l02','l08']});old=pd.read_csv(b/'analysis_metadata.tsv.gz',sep='\t');assert d.bin_id.equals(old.bin_id)
for tag in ['l00','l02','l08']:
 mp=regions[regions.branch==tag].set_index('cluster').interpretation;d[f'banksy_{tag}_region']=d[f'banksy_{tag}_cluster'].map(mp);assert d[f'banksy_{tag}_region'].notna().all()
d.to_csv(o/'metadata.tsv.gz',sep='\t',index=False)
tr=json.loads((b/'transform.json').read_text());sc=tr['local_to_model_scale_xy'];cr=tr['fullres_crop_xyxy'];mp=tr['native_mpp'];mx=(d.x_um/mp-cr[0]+.5)*sc[0]-.5;my=(d.y_um/mp-cr[1]+.5)*sc[1]-.5
rr=[]
for _,roi in pd.read_csv(o/'fixed_rois.tsv',sep='\t').iterrows():
 mask=(mx>=roi.x)&(mx<roi.X)&(my>=roi.y)&(my<roi.Y)
 for tag in ['l00','l02','l08']:
  for cl,z in d[mask].groupby(f'banksy_{tag}_cluster'):
   idx=z.index;rr.append(dict(roi=roi.roi,branch=tag,cluster=cl,interpretation=z[f'banksy_{tag}_region'].iloc[0],n=len(z),ROI_n=int(mask.sum()),UMI=int(old.loc[idx,'nCount_RNA'].sum()),nFeature_median=old.loc[idx,'nFeature_RNA'].median(),forward_conflict_n=(z.reference_conflict_group=='BLCA_Fibroblast_HPCA_Epithelial').sum(),edge16_fraction=(z.native_edge_distance_um<=16).mean(),long_neighbor80_fraction=(z.neighbor_max_um>80).mean()))
pd.DataFrame(rr).to_csv(o/'ROI_cluster_summary.tsv',sep='\t',index=False)
print(s,'INTERPRETATIONS',regions[regions.branch=='l08'][['cluster','interpretation','n','nFeature_median','top_original_markers']].to_string(index=False))
