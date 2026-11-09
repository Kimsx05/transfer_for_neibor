from pathlib import Path
import pandas as pd,numpy as np,json,sys
root=Path(__file__).resolve().parents[1];hd=root.parents[1];sections=sys.argv[1:] or ['Cervilla','D1_Johnson'];bands=['<20','20-99','100-499','>=500'];covered=['Epithelial','T/NK','Myeloid','Fibroblast','Endothelial','B']
def save(d,p):d.to_csv(p,sep='\t',index=False)
def band(z):return pd.cut(z,[-np.inf,19,99,499,np.inf],labels=bands)
def metrics(d):
 n=len(d);u=d.nCount_RNA.sum();b=d.singler_blca_pruned_labels;h=d.singler_hpca_pruned_labels;both=b.notna()&h.notna();cov=both&b.isin(covered);dis=both&(b!=h);ef=both&(((b=='Epithelial')&(h=='Fibroblast'))|((b=='Fibroblast')&(h=='Epithelial')))
 z=dict(n=n,UMI=int(u),median_nFeature=d.nFeature_RNA.median(),median_UMI=d.nCount_RNA.median(),both_evaluable=int(both.sum()),reference_disagree=int(dis.sum()),reference_disagree_UMI=int(d.loc[dis,'nCount_RNA'].sum()),covered_both=int(cov.sum()),covered_disagree=int((dis&cov).sum()),epi_fib_ge2=int(d.epi_fib_ge2.sum()),epi_fib_ge2_UMI=int(d.loc[d.epi_fib_ge2,'nCount_RNA'].sum()),marker_conflict=int(d.marker_conflict.sum()))
 for r in ['blca','hpca']:
  valid=d[f'singler_{r}_labels'].notna();p=d[f'singler_{r}_pruned_labels'].notna();z[r+'_evaluated']=int(valid.sum());z[r+'_pruned']=int((valid&~p).sum());z[r+'_candidate_ge2']=int((d[f'singler_{r}_marker_count']>=2).sum())
 for k in ['reference_disagree','epi_fib_ge2','marker_conflict']:z[k+'_fraction']=z[k]/n if n else np.nan
 z['epi_fib_reference_disagree']=int(ef.sum());z['epi_fib_reference_disagree_UMI']=int(d.loc[ef,'nCount_RNA'].sum());z['covered_disagree_fraction']=z['covered_disagree']/z['covered_both'] if z['covered_both'] else np.nan
 z['reference_disagree_UMI_fraction']=z['reference_disagree_UMI']/u if u else np.nan;z['epi_fib_ge2_UMI_fraction']=z['epi_fib_ge2_UMI']/u if u else np.nan
 return z
for s in sections:
 o=root/s;d=pd.read_csv(o/'analysis_metadata.tsv.gz',sep='\t',low_memory=False);p=pd.read_csv(o/'predictions.tsv.gz',sep='\t');d=d.drop(columns=[k for k in p.columns if k!='bin_id' and k in d],errors='ignore').merge(p,on='bin_id',validate='one_to_one');d['cell_id']=d.bin_id;d['unit']='bin8';d['information_band']=band(d.nFeature_RNA)
 c=pd.read_csv(hd/'analysis/hd_annotation_diagnosis_singler_v1'/s/'diagnosis.tsv.gz',sep='\t',low_memory=False);ev=pd.read_csv(hd/'analysis/hd_major_annotation_v1'/s/'tables/marker_evidence.tsv.gz',sep='\t',usecols=['cell_id','detected_Epithelial','detected_Fibroblast']);c=c.merge(ev,on='cell_id',validate='one_to_one');c['epi_fib_ge2']=(c.detected_Epithelial>=2)&(c.detected_Fibroblast>=2);c['marker_conflict']=c.diag_fail_conflict;c['unit']='cell8';c['information_band']=band(c.nFeature_RNA)
 for r in ['blca','hpca']:c[f'singler_{r}_marker_count']=c[f'diag_{r}_candidate_marker_count']
 uni=[];st=[];comp=[];cross=[]
 for un,z in [('cell8',c),('bin8',d)]:
  uni.append(dict(section=s,unit=un,**metrics(z)))
  for b in bands:st.append(dict(section=s,unit=un,information_band=b,**metrics(z[z.information_band==b])))
  for r in ['blca','hpca']:
   for typ in ['labels','pruned_labels']:
    k=f'singler_{r}_{typ}';zz=z.assign(label=z[k].fillna('Not evaluable / pruned'));t=zz.groupby('label',observed=True).agg(n=('label','size'),UMI=('nCount_RNA','sum')).reset_index();t['fraction']=t.n/len(z);t['UMI_fraction']=t.UMI/z.nCount_RNA.sum();t['unit']=un;t['reference']=r;t['stage']=typ;comp.append(t)
  zz=z.assign(BLCA=z.singler_blca_pruned_labels.fillna('Not evaluable / pruned'),HPCA=z.singler_hpca_pruned_labels.fillna('Not evaluable / pruned'));t=zz.groupby(['BLCA','HPCA']).agg(n=('cell_id','size'),UMI=('nCount_RNA','sum')).reset_index();t['unit']=un;t['fraction']=t.n/len(z);t['within_BLCA']=t.n/t.groupby('BLCA').n.transform('sum');t['HPCA_covers_BLCA']=t.BLCA.isin(covered);cross.append(t)
 save(pd.DataFrame(uni),o/'unit_comparison.tsv');save(pd.DataFrame(st),o/'information_strata.tsv');save(pd.concat(comp),o/'method_composition.tsv');save(pd.concat(cross),o/'reference_crosstabs.tsv')
 # Fixed ROI tables use physical bounds and original cell IDs, not nearest-centroid pairing.
 rx=pd.read_csv(o/'roi_expression.tsv.gz',sep='\t');rstat=[]
 for un,z in [('cell8',c),('bin8',d)]:
  rr=rx[rx.unit==un][['object_id','roi']];zz=z.drop(columns=['roi'],errors='ignore').merge(rr,left_on='cell_id',right_on='object_id',validate='one_to_one')
  for ro,g in zz.groupby('roi'):rstat.append(dict(unit=un,roi=ro,**metrics(g)))
 save(pd.DataFrame(rstat),o/'ROI_reference_comparison.tsv')
 rr=[]
 for (un,ro,an),z in rx.groupby(['unit','roi','analysis']):
  for b in ['all']+bands:
   zz=z if b=='all' else z[band(z.nFeature)==b]
   if not len(zz):continue
   rr.append(dict(unit=un,roi=ro,analysis=an,information_band=b,n=len(zz),UMI=zz.nCount.sum(),epi_UMI=zz.raw_Epithelial.sum(),collagen_UMI=zz.raw_Fibroblast.sum(),immune_UMI=zz.raw_Immune.sum(),endothelial_UMI=zz.raw_Endothelial.sum(),both_any=int((zz.signal_category=='Both panels').sum()),both_ge2=int(zz.epi_fib_ge2.sum()),both_ge2_fraction=float(zz.epi_fib_ge2.mean())))
 save(pd.DataFrame(rr),o/'ROI_signal_strata.tsv')
 z4=rx[rx.unit=='cell4'];z8=rx[rx.unit=='cell8'];pp=z4.merge(z8,on=['object_id','roi'],suffixes=('_4','_8'),validate='one_to_one');rows=[]
 for ro,z in pp.groupby('roi'):
  rows.append(dict(roi=ro,matched_objects=len(z),UMI_4=z.nCount_4.sum(),UMI_8=z.nCount_8.sum(),epi_fib_ge2_4=int(z.epi_fib_ge2_4.sum()),epi_fib_ge2_8=int(z.epi_fib_ge2_8.sum()),gained_ge2=int((~z.epi_fib_ge2_4&z.epi_fib_ge2_8).sum()),lost_ge2=int((z.epi_fib_ge2_4&~z.epi_fib_ge2_8).sum()),counts_decreased=int((z.nCount_8<z.nCount_4).sum())))
 save(pd.DataFrame(rows),o/'assignment4_vs8_ROI.tsv')
 coverage0=pd.read_csv(o/'bin_metadata.tsv.gz',sep='\t');coverage0=coverage0[coverage0.qc_keep&coverage0.roi.notna()].groupby('roi').agg(bins=('bin_id','size'),native8_UMI=('nCount_RNA','sum'),supplied2_UMI=('common_2um_UMI','sum'),unassigned2_UMI=('cell_unassigned_UMI','sum'));coverage0['unassigned_fraction']=coverage0.unassigned2_UMI/coverage0.supplied2_UMI;save(coverage0.reset_index(),o/'ROI_coverage_accounting.tsv')

 if (o/'rctd_ROI_results.tsv.gz').exists():
  rt=pd.read_csv(o/'rctd_ROI_results.tsv.gz',sep='\t');rt=rt.merge(d[['bin_id','singler_blca_pruned_labels','singler_hpca_pruned_labels','epi_fib_ge2']],on='bin_id',how='left');save(rt.groupby(['roi','rctd_status'],dropna=False).agg(n=('bin_id','size'),UMI=('nCount_RNA','sum')).reset_index(),o/'RCTD_ROI_status_summary.tsv');valid=rt.first_type.notna();rt['epi_fib_pair']=valid&((rt.first_type=='Epithelial')&(rt.second_type=='Fibroblast')|(rt.first_type=='Fibroblast')&(rt.second_type=='Epithelial'));rt['epi_fib_reference_disagreement']=((rt.singler_blca_pruned_labels=='Epithelial')&(rt.singler_hpca_pruned_labels=='Fibroblast'))|((rt.singler_blca_pruned_labels=='Fibroblast')&(rt.singler_hpca_pruned_labels=='Epithelial'))
  zz=[]
  for ro,z in rt.groupby('roi'):
   for sub in ['all','epi_fib_reference_disagreement','epi_fib_ge2']:
    zz0=z if sub=='all' else z[z[sub]==True];zz.append(dict(roi=ro,subset=sub,n=len(zz0),fitted=int(zz0.first_type.notna().sum()),epi_fib_pair=int(zz0.epi_fib_pair.sum()),epi_fib_doublet_certain=int((zz0.epi_fib_pair&(zz0.rctd_status=='doublet_certain')).sum()),epi_fib_doublet_uncertain=int((zz0.epi_fib_pair&(zz0.rctd_status=='doublet_uncertain')).sum())))
  save(pd.DataFrame(zz),o/'RCTD_mixture_evidence.tsv')
 if (o/'author_sensitivity_predictions.tsv.gz').exists():
  coverage=pd.read_csv(o/'marker_coverage.tsv',sep='\t');missing=set(pd.read_csv(o/'author_SCT_missing_genes.tsv',sep='\t').symbol);coverage['SCT_retained']=coverage.present&~coverage.symbol.isin(missing);save(coverage,o/'author_SCT_marker_coverage.tsv');a=pd.read_csv(o/'author_sensitivity_predictions.tsv.gz',sep='\t');subset=d[d.bin_id.isin(a.bin_id)];composition=subset[[k for k in subset if k.startswith('detected_')]].ge(2).sum().reset_index();composition.columns=['panel','bins_ge2_markers'];save(composition,o/'author_subset_panel_composition.tsv');a['information_band']=band(a.nFeature);t=a.groupby(['author_LogNormalize_major','author_SCT_data_major'],dropna=False).agg(n=('bin_id','size'),UMI=('nCount','sum')).reset_index();save(t,o/'author_major_crosstab.tsv');rows=[]
  for b,z in a.groupby('information_band',observed=True):
   v=z.author_LogNormalize_major.notna()&z.author_SCT_data_major.notna();rows.append(dict(band=b,n=len(z),both_evaluable=int(v.sum()),changed=int((v&(z.author_LogNormalize_major!=z.author_SCT_data_major)).sum()),fine_changed=int((v&(z.author_LogNormalize_labels!=z.author_SCT_data_labels)).sum()),Log_pruned=int((z.author_LogNormalize_labels.notna()&z.author_LogNormalize_pruned.isna()).sum()),SCT_pruned=int((z.author_SCT_data_labels.notna()&z.author_SCT_data_pruned.isna()).sum())))
  save(pd.DataFrame(rows),o/'author_information_strata.tsv')
 print('COMPARE_COMPLETE',s)
