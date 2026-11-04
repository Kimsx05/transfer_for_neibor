from core import *
def main():
 out=[]
 for sec in M.section_id:
  d=pd.read_parquet(O/'spots'/f'{sec}.parquet');normal=M.set_index('section_id').loc[sec,'source_context']=='adjacent_normal_source';groups={'normal_Epi':d['query'].to_numpy()} if normal else {'P95':d.MP10_P95.to_numpy(),'nonP95':(d.tumor_epi&~d.MP10_P95).to_numpy()};m=pd.read_csv(OLD/'tables'/f'{sec}_matching.tsv.gz',sep='\t')
  if 'matching' in m:
   mm=m[m.matching.eq('base')];groups['P95_base_matched']=d.barcode.isin(mm.high_barcode).to_numpy();groups['control_base_matched']=d.barcode.isin(mm.control_barcode).to_numpy()
  for t in CAND:
   for pop,mask in groups.items():
    x=d.loc[mask,'d_candidate__'+t].dropna().to_numpy();v=np.quantile(x,[0,.05,.25,.5,.75,.95,1]) if len(x) else np.repeat(np.nan,7);r=dict(section_id=sec,candidate_type=t,population=pop,n=len(x),n_no_same_component_reference=int(np.sum(mask)-len(x)),**dict(zip(['min_um','q05_um','q25_um','median_um','q75_um','q95_um','max_um'],v)));r.update({name:int(k.sum()) for name,k in [('n_zero',x==0),('n_0_200',(x>0)&(x<=200)),('n_200_500',(x>200)&(x<=500)),('n_500_1000',(x>500)&(x<=1000)),('n_above1000',x>1000)]});out.append(r)
 a=pd.DataFrame(out).merge(M[['section_id','dataset_id','analysis_unit_uid','patient_uid','source_context']],on='section_id',validate='m:1');save(a,O/'tables/candidate_distance_distributions.tsv.gz');print('Candidate distance distributions complete',flush=True)
if __name__=='__main__':main()
