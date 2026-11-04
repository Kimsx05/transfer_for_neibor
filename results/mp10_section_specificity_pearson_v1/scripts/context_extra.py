from context import *
def extra(sec):
 if (O/'tables'/f'{sec}_context_extra.tsv').exists():return sec
 d=pd.read_parquet(O/'spots'/f'{sec}.parquet');c=pd.read_parquet(O/'centers'/f'{sec}_new5.parquet').set_index('barcode').reindex(d.barcode);q=d['query'].to_numpy(bool);tum=d.tumor_epi.to_numpy(bool);df=d.d_front.to_numpy(float);ex=np.column_stack([api.z(df),api.z(df)**2]);fr=tum&(df>=-200)&(df<=0);out=[];curves=[];rings=[]
 for ring in RINGS:
  y=c[[ring+'__'+p for p in NEW]].to_numpy(float)
  for adj in [False,True]:out+=api.fit(d,d.MP10_P95.to_numpy(float),y,NEW,tum,'B4_position_'+ring,'B_neighborhood',adjust=adj,extra=ex,binary=True,blocksize=2000)
 for adj in [False,True]:out+=api.fit(d,d.pEMT_high_front.to_numpy(float),d[NEW].to_numpy(float),NEW,fr,'A_active_front_support','A_front',adjust=adj,binary=True,extra=ex)
 for ref in ['P95','control_base','P95_matched_base','control_front','P95_matched_front']:
  dist=d['d_MP10' if ref=='P95' else 'd_'+ref].to_numpy(float);bins=np.where(dist==0,0,np.ceil(dist/100));bins[dist>1000]=11;bins[~np.isfinite(dist)]=np.nan;hop=d['hops_'+ref].to_numpy(float)
  for pop,mask in [('tissue',d.valid_tissue.to_numpy()),('Epi',q)]:
   for b in np.unique(bins[np.isfinite(bins)&mask]):
    ix=mask&(bins==b)
    for p in NEW:curves.append(dict(section_id=sec,axis='d_MP10',reference=ref,population=pop,bin=int(b),distance=np.nanmedian(dist[ix]),program=p,n=ix.sum(),mean=d.loc[ix,p].mean() if ix.sum()>=5 else np.nan))
   for ring,(lo,hi) in zip(RINGS,api.RB):
    ix=mask&(hop>=lo)&(hop<=hi)
    for p in NEW:rings.append(dict(section_id=sec,reference=ref,population=pop,ring=ring,program=p,n=ix.sum(),mean=d.loc[ix,p].mean() if ix.sum()>=5 else np.nan,distance=np.nanmedian(dist[ix]) if ix.any() else np.nan))
 a=pd.read_csv(O/'tables'/f'{sec}_new5_context.tsv.gz',sep='\t');b=pd.DataFrame(out);b['section_id']=sec;save(pd.concat([a,b],ignore_index=True),O/'tables'/f'{sec}_new5_context.tsv.gz');save(curves,O/'tables'/f'{sec}_new5_reference_curves.tsv.gz');save(rings,O/'tables'/f'{sec}_new5_reference_rings.tsv.gz');oldcenter=pd.read_parquet(OLD/'centers'/f'{sec}.parquet');newcenter=pd.read_parquet(O/'centers'/f'{sec}_new5.parquet');extra_cols=[x for x in oldcenter if '__' not in x and x not in newcenter];newcenter=newcenter.merge(oldcenter[['section_id','barcode']+extra_cols],on=['section_id','barcode'],validate='1:1');newcenter.to_parquet(O/'centers'/f'{sec}_new5.parquet',index=False);save([dict(section_id=sec,status='complete')],O/'tables'/f'{sec}_context_extra.tsv');return sec
RINGS=['R0','R1','R2','R3_R5','R6_R10']
if __name__=='__main__':
 with ProcessPoolExecutor(max_workers=2) as ex:
  for s in ex.map(extra,M.section_id):print(s,flush=True)
