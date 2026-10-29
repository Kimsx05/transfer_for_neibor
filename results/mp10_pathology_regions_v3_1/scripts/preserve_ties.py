from revise import *
from concurrent.futures import ProcessPoolExecutor

def fix(sec):
 r=M.loc[sec];p=OUT/'by_dataset'/r.dataset_id/sec;d=pd.read_parquet(p/'spot_regions_v3_1.parquet');un=pd.read_csv(p/'unknown_region_boundary.tsv.gz',sep='\t');adj,_=graph(d);tc,tcs=components(adj,np.ones(len(d),bool));changed=[]
 for k,ix in enumerate(tcs):
  ix=ix[d.signed_distance_um.iloc[ix].notna().to_numpy()];u=un[un.component==k]
  if not len(ix) or not len(u):continue
  q,jj,_=nearest(d[['x_um','y_um']].to_numpy()[ix],segs(u));bad=q<=abs(d.signed_distance_um.iloc[ix].to_numpy())+1e-6
  for i,j in zip(ix[bad],jj[bad]):d.loc[i,'signed_distance_um']=np.nan;d.loc[i,'exact_distance_status']=u.unknown_type.iloc[j];d.loc[i,'distance_status']=u.unknown_type.iloc[j];changed.append(i)
 if changed:
  for w in [100,200,300]:
   for i in changed:d.loc[i,f'band_status_{w}um']='upper_bound_confirms_margin' if d.abs_distance_upper_bound_um.iloc[i]<=w else 'lower_bound_confirms_non_margin'
  d['band_status']=d.band_status_200um
  old=d.v3_signed_distance_um.notna();new=d.signed_distance_um.notna();d['exact_distance_change_reason']=np.where(~old&new,'technical_gap_geometry_support_and_no_closer_unknown_front',np.where(old&~new,'more_conservative_all_unknown_region_model','unchanged_estimability'))
  d.to_parquet(p/'spot_regions_v3_1.parquet',index=False);core=pd.read_csv(p/'spot_labels.tsv.gz',sep='\t');d[list(core.columns)].to_csv(p/'spot_labels.tsv.gz',sep='\t',index=False)
  co=json.loads((p/'coordinate_scale.json').read_text());aff=np.array(co['array_to_image_affine']);keys=np.c_[3*d.array_col,3*d.array_row].astype(int)[:,None,:]+VERT;poly=np.concatenate([keys[:,:,:1]/6,keys[:,:,1:]/(2*np.sqrt(3)),np.ones((len(d),6,1))],axis=2)@aff
  ed=pd.read_csv(p/'lesion_boundary_edges.tsv.gz',sep='\t');holes=pd.read_csv(p/'internal_gap_assessment.tsv',sep='\t').to_dict('records');plot(sec,d,poly,np.asarray(Image.open(r.image_path).convert('RGB')),ed,un,holes,aff,p)
 info=json.loads((p/'complete.json').read_text());old=d.v3_signed_distance_um.notna();new=d.signed_distance_um.notna();info.update(exact_v3_1=int(new.sum()),exact_restored=int((~old&new).sum()),exact_newly_censored=int((old&~new).sum()),exact_unestimated=int((~new).sum()),band_only_restored=int((~new&d.distance_band_200um.isin(list(BP)[:4])&~d.v3_distance_band_200um.isin(['Tumor interior','Inner margin','Outer margin','Outside'])).sum()));(p/'complete.json').write_text(json.dumps(info,indent=2));print(sec,len(changed),flush=True)
if __name__=='__main__':
 with ProcessPoolExecutor(max_workers=2) as pool:list(pool.map(fix,list(M.index)))
