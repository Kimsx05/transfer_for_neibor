from base_data import *
for pl,samples in [('hd',HD),('xenium',XE)]:
 pieces={};lineages={};records=[]
 for s in samples:
  a=dataset(s);z=pd.read_parquet(paths(s)/'analysis_index.parquet');assert z.index.equals(a.obs_names)
  x=norm(a.X)
  for cluster,idx in z.groupby('lineage_cluster',observed=True).indices.items():
   if cluster=='unclustered':continue
   par=cluster.rsplit('.',1)[0];v=(len(idx),x[idx].sum(0).A1.astype('f8'),x[idx].getnnz(0).astype('f8'))
   def add(store,k):
    if k not in store:store[k]=list(v[:1])+[b.copy() for b in v[1:]]
    else:
     for j in range(3):store[k][j]+=v[j]
   add(pieces,(s,cluster));add(lineages,(s,par))
  if s==samples[0]:gg=a.var_names.copy()
  else:assert gg.equals(a.var_names)
 for (_,cluster) in sorted(set((0,c) for s,c in pieces)):
  par=cluster.rsplit('.',1)[0];delta=[];frac=[];weights=[];positive=[];names=[]
  for s in samples:
   if (s,cluster) not in pieces:continue
   n,xx,h=pieces[s,cluster];N,total,hit=lineages[s,par]
   if n<5 or N-n<5:continue
   v=xx/n-(total-xx)/(N-n);delta.append(v);frac.append(h/n);weights.append(n);positive.append(v>0);names.append(s)
  if not delta:continue
  av=np.mean(delta,axis=0);det=np.mean(frac,axis=0);rev=np.sum(positive,axis=0);order=np.flatnonzero((det>=.05)&(av>0));order=order[np.argsort(-av[order])[:12]]
  for rank,j in enumerate(order,1):records.append((cluster,par,rank,gg[j],av[j],det[j],rev[j],len(names),int(sum(weights))))
 out=pd.DataFrame(records,columns=['cluster','lineage','rank','gene','mean_within_section_log_expression_delta','mean_cluster_detection_fraction','sections_positive','sections_evaluable','objects'])
 write(out,pl+'_lineage_marker_ranks.tsv')

