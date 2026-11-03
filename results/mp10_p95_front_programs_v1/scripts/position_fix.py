from analyze import *
def work(sec):
 d=pd.read_parquet(O/'spots'/f'{sec}.parquet');c=pd.read_parquet(O/'centers'/f'{sec}.parquet').set_index('barcode').reindex(d.barcode);df=d.d_front.to_numpy(float);ex=np.column_stack([z(df),z(df)**2]);out=[]
 for ring in RINGS:
  y=c[[ring+'__'+p for p in P]].to_numpy(float)
  for adj in [False,True]:out+=fit(d,d.MP10_P95.to_numpy(float),y,P,d.tumor_epi.to_numpy(bool),'B4_position_'+ring,'B_neighborhood',adjust=adj,extra=ex,binary=True,blocksize=2000)
 ny=np.column_stack([d.accepted_label.eq(k).fillna(False).to_numpy(float) for k in range(9)]);lab=d['query'].to_numpy(bool)&d.accepted_label.notna().to_numpy()
 for adj in [False,True]:out+=fit(d,d[RAW].to_numpy(float),ny,['R12_'+str(k) for k in range(9)],lab,'D_continuous_all_Epi_position','D_MP10_niche',adjust=adj,raw_y=True,extra=np.column_stack([df,df**2]))
 p=O/'tables'/f'{sec}_effects.tsv.gz';a=pd.read_csv(p,sep='\t');a=a[~a.analysis.str.startswith('B4_position_')&a.analysis.ne('D_continuous_all_Epi_position')];b=pd.DataFrame(out);b['section_id']=sec;save(pd.concat([a,b],ignore_index=True),p);return sec
if __name__=='__main__':
 with ProcessPoolExecutor(max_workers=4) as ex:
  for s in ex.map(work,M.section_id):print(s,flush=True)
