from plots import *
def main():
 e=pd.read_csv(O/'tables/effects.tsv.gz',sep='\t').merge(M[['section_id','dataset_id','source_context','analysis_unit_uid']],on='section_id')
 e=e[e.source_context.str.startswith('tumor')&e.analysis_unit_uid.notna()]
 rows=[]
 for axis,g in e.groupby('axis'):
  cols={}
  for label,kind,pr,col in [('raw','distance',RAW,'r'),('adjusted','distance',RAW,'adj'),('model_spot','matched','spot_MP10_model_universe','r'),('model_epi','matched','epi_MP10_model_universe','r')]:
   z=g[g.kind.eq(kind)&g.program.eq(pr)].dropna(subset=[col]).copy()
   z['z']=np.arctanh(z[col].clip(-.999999,.999999))
   u=z.groupby(['dataset_id','analysis_unit_uid']).z.mean().groupby('dataset_id').mean()
   cols[label]=np.tanh(u.mean()) if len(u) else np.nan
  rows.append(dict(axis=axis,**cols))
 a=pd.DataFrame(rows);write(a,O/'tables/raw_adjusted_model_comparison.tsv')
 heat(a.set_index('axis'),'P95/P90 and interface / raw, adjusted, matched model layers','sensitivity_comparison',.3)
if __name__=='__main__':main()
