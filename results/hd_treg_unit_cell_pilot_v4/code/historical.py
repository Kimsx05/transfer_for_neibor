from base import *
rows=[];paired=[]
for s in sections:
 root=R/s;roi=rd(R/'roi_manifest.tsv').query('section==@s')
 b=pd.read_parquet(V/s/'bin_evidence.parquet',columns=['bin_id','array_row','array_col','x_um','y_um','eligible_pool','analysis_qc','nCount_RNA','nFeature_RNA'])
 b['roi']=roi_for(b.x_um,b.y_um,roi);b['id']=b.array_row.astype(str)+'_'+b.array_col.astype(str)
 t=pd.read_parquet(root/'T8_from2.parquet')[['id','roi','nCount','nFeature','n_supplied_usable']]
 pair=t.merge(b[['id','bin_id','roi','nCount_RNA','nFeature_RNA','analysis_qc']],how='left',on='id',suffixes=('_from2','_native'));pair['section']=s;paired.append(pair)
 bb=b[b.roi.ne('')];rows.append(dict(section=s,unit='historical_native8',units=len(bb),nCount=bb.nCount_RNA.sum(),eligible=bb.analysis_qc.sum(),FOXP3_supported=np.nan,strict=np.nan,eligibility_definition='native8 analysis_qc'))
 for k in ['T16','T32']:
  d=pd.read_parquet(V/s/(k+'_regions.parquet'),columns=['region_id','x_um','y_um','measurement_eligible','nCount','evidence','raw_FOXP3']);d=d[roi_for(d.x_um,d.y_um,roi)!='']
  rows.append(dict(section=s,unit='historical_'+k,units=len(d),nCount=d.nCount.sum(),eligible=d.measurement_eligible.sum(),FOXP3_supported=(d.measurement_eligible&d.raw_FOXP3.gt(0)).sum(),strict=d.evidence.eq('Treg_like_multigene').sum(),eligibility_definition='historical measured T lineage'))
put(pd.DataFrame(rows),'historical_unit_baselines.tsv');pd.concat(paired).to_parquet(R/'native8_supplied2_comparison.parquet',index=False)
