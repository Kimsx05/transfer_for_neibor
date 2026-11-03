from analyze import *
def run_bands(sec):
 d=pd.read_parquet(O/'spots'/f'{sec}.parquet');df=d.d_front.to_numpy(float);q=d.tumor_epi.to_numpy(bool);fr=q&(df>=-200)&(df<=0);effects=[];desc=[]
 for contrast,mask in [('front_vs_near',q&(df>=-500)&(df<=0)),('front_vs_deep',q&((df< -500)|fr))]:
  for adjust in [False,True]:
   effects+=fit(d,fr.astype(float),d[P].to_numpy(float),P,mask,'A_'+contrast,'A_front',adjust=adjust,binary=True)
   effects+=fit(d,fr.astype(float),d[[RAW]].to_numpy(float),[RAW],mask,'C1_'+contrast,'C_MP10_front',adjust=adjust,binary=True)
   effects+=fit(d,fr.astype(float),d[['MP10_P95']].to_numpy(float),['P95_probability'],mask,'C2_'+contrast,'C_MP10_front',adjust=adjust,binary=True,raw_y=True)
 for pop,mask in [('tumor_Epi',q),('P95_Epi',d.MP10_P95.to_numpy(bool)),('outside_tissue',d.lesion_location.eq('outside').to_numpy()&d.Epithelial_fraction.notna().to_numpy())]:
  for band,g in d[mask].groupby('front_band'):
   for p in [RAW]+P:
    v=g[p].dropna();desc.append(dict(section_id=sec,population=pop,band=band,program=p,n=len(v),mean=v.mean() if len(v)>=5 else np.nan,median=v.median() if len(v)>=5 else np.nan,q25=v.quantile(.25) if len(v)>=5 else np.nan,q75=v.quantile(.75) if len(v)>=5 else np.nan))
 a=pd.read_csv(O/'tables'/f'{sec}_effects.tsv.gz',sep='\t');a=a[~a.analysis.str.endswith(('front_vs_near','front_vs_deep'))];b=pd.DataFrame(effects);b['section_id']=sec;save(pd.concat([a,b],ignore_index=True),O/'tables'/f'{sec}_effects.tsv.gz');save(desc,O/'tables'/f'{sec}_front_program_bands.tsv.gz');return sec
if __name__=='__main__':
 with ProcessPoolExecutor(max_workers=4) as ex:
  for sec in ex.map(run_bands,M.section_id):print(sec,flush=True)
 a=pd.concat([pd.read_csv(O/'tables'/f'{s}_front_program_bands.tsv.gz',sep='\t') for s in M.section_id],ignore_index=True);save(a,O/'tables/front_program_bands.tsv.gz')
