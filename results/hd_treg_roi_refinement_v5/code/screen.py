from common5 import *
def run(s):
 out=R/s;out.mkdir(exist_ok=True);d,audit=inputdata(s);gid,sy=genedata(audit);panel=T+reg;ix=pd.Index(sy);parts=[]
 for a,b,x in slabs(audit):
  if not np.array_equal(np.asarray(x.sum(axis=1)).ravel(),d.raw_counts.to_numpy()[a:b]):raise ValueError('Barcode/count alignment')
  genes=np.array([np.asarray(x[:,sy==g].sum(axis=1)).ravel()for g in panel]).T
  j=np.flatnonzero(genes.sum(1));z=d.iloc[a+j][['source_bin_index','array_row','array_col','x','y','nucleus_label','cell_label','usable']].copy();z[panel]=genes[j];parts.append(z)
 raw=pd.concat(parts,ignore_index=True);raw.to_parquet(out/'raw_candidates.parquet',index=False)
 cells=qc(s);old=pd.read_csv(A/'hd_annotation_diagnosis_singler_v1'/s/'diagnosis.tsv.gz',sep='\t',usecols=['cell_id','final_major','singler_hpca_labels','singler_blca_labels','raw_prediction'])
 hints=old.drop(columns='cell_id').isin(['T/NK','Treg','T']).any(axis=1);o=cells.merge(old[hints],on='cell_id');o['source']='historical_major'
 if s=='Cervilla':
  f=pd.read_parquet(O/'tables/Cervilla_cell_annotation.parquet',columns=['cell_id','fine_prediction']).reset_index(drop=True);f=f[ist(f.fine_prediction)];f=cells.merge(f,on='cell_id');f['source']='historical_fine';o=pd.concat([o,f])
 for k in ['cell_nucleus','cell_conservative','existing_cell8','T8_from2','T16_from2']:
  v=pd.read_parquet(U/s/(k+'_evidence.parquet'));v=v[ist(v.predicted)];v=v.rename(columns={'x_um':'x','y_um':'y'});v['source']='v4_'+k;v['historical_unit_id']=v.id.astype(str);v['id']=v.nid if 'nid'in v else pd.NA;o=pd.concat([o,v[['source','id','x','y','historical_unit_id','predicted']]])
 for k in ['T16','T32']:
  v=pd.read_parquet(V/s/(k+'_regions.parquet'),columns=['region_id','Treg_main','measurement_eligible','raw_FOXP3','x_um','y_um']);v=v[v.measurement_eligible|v.Treg_main.gt(0)|v.raw_FOXP3.gt(0)].rename(columns={'x_um':'x','y_um':'y','region_id':'historical_unit_id'});v['source']='historical_'+k;o=pd.concat([o,v])
 o.reset_index(drop=True).to_parquet(out/'historical_candidates.parquet',index=False)
 print(s,len(d),len(raw),len(o),flush=True)
if __name__=='__main__':
 for s in sections:run(s)
