from run import *
def extra(sec):
 d=pd.read_parquet(O/'spots'/f'{sec}.parquet');rows=[];region=[]
 axes={k:d[k+'_P95_distance'].to_numpy() for k in list(TARGETS)[:7]};axes['local_interface']=d.local_interface_um.to_numpy()
 for col in ['accepted_label','local_niche_id']:
  for label,g in d[d['query']].groupby(col):
   if len(g)<100:continue
   for axis,x in axes.items():
    fit=stats(g,x[g.index],g[PROG].to_numpy(),True)
    if fit:
     for j,p in enumerate(PROG):rows.append(dict(section_id=sec,niche_version=col,niche=label,axis=axis,program=p,**{k:(v[j] if isinstance(v,np.ndarray) else v) for k,v in fit.items()}))
 for axis,x in axes.items():
  bins=pd.cut(x,[-np.inf,-500,-200,0,100,200,500,1000,np.inf])
  for label,g in d.groupby(bins,observed=True):
   for pop,sub in [('tissue',g[g.center_qc]),('epi_query',g[g['query']])]:
    vals=sub[PROG];q=vals.quantile([.25,.5,.75]).T;q.columns=['q25','median','q75'];q['n']=vals.count();q['program']=q.index
    q=q.reset_index(drop=True).assign(section_id=sec,axis=axis,bin=str(label),population=pop);region.append(q)
 for name,items in [('niche_interaction',pd.DataFrame(rows)),('program_band',pd.concat(region))]:write(items,O/'tables'/f'{sec}_{name}.tsv.gz')
 return sec
if __name__=='__main__':
 with ProcessPoolExecutor(max_workers=6) as p:
  for s in p.map(extra,M.section_id):print(s,flush=True)
