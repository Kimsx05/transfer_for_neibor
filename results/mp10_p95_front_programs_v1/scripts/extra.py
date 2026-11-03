from analyze import *
def extra(sec):
 d=pd.read_parquet(O/'spots'/f'{sec}.parquet');q=d['query'].to_numpy(bool);df=d.d_front.to_numpy(float);fr=d.tumor_epi.to_numpy(bool)&(df>=-200)&(df<=0);mp=d[RAW].to_numpy(float);Y=d[P].to_numpy(float);out=[]
 for adj in [False,True]:
  out+=fit(d,df,mp[:,None],[RAW],q,'C1_all_Epi_continuous','C_MP10_front',adjust=adj,quadratic=True)
  names=[p for p in P if not p.startswith('PURAM_PEMT')]
  out+=fit(d,d.pEMT_high_front.to_numpy(float),d[names].to_numpy(float),names,fr,'A_active_front_support','A_front',adjust=adj,binary=True,extra=np.column_stack([df,df**2]))
  lab=q&d.accepted_label.notna().to_numpy();ny=np.column_stack([d.accepted_label.eq(k).fillna(False).to_numpy(float) for k in range(9)])
  out+=fit(d,mp,ny,['R12_'+str(k) for k in range(9)],lab,'D_continuous_all_Epi_position','D_MP10_niche',adjust=adj,raw_y=True,extra=np.column_stack([df,df**2]))
 p=O/'tables'/f'{sec}_effects.tsv.gz';a=pd.read_csv(p,sep='\t');a=a[~a.analysis.isin(['C1_all_Epi_continuous','A_active_front_support','D_continuous_all_Epi_position'])];out=pd.DataFrame(out);out['section_id']=sec;save(pd.concat([a,out],ignore_index=True),p)
 return sec
if __name__=='__main__':
 with ProcessPoolExecutor(max_workers=4) as ex:
  for sec in ex.map(extra,M.section_id):print(sec,flush=True)
