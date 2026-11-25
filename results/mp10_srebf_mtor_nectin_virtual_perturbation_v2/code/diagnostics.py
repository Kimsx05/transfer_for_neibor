from common import *
T=Path('tables');P0=rd(T/'program_members.tsv');cc=rd(T/'control_matching.tsv');focus=E+K
po=rd(T/'panel_members.tsv').set_index('gene');spec=rd(T/'control_specificity.tsv')
# Compare like with like, excluding a control's own KO row.
a=rd('cache/jobs/clean/heat.tsv.gz');measured=a[~a.is_KO]
for k in K:
 ids=spec.ko==k;ctrl=measured[measured.ko.isin(cc[cc.ko==k].control)].pivot(index='gene',columns='ko',values='distance')
 idx=spec[ids].set_index('gene');vv=ctrl.reindex(idx.index)
 spec.loc[ids,'n_controls_lower']=vv.lt(idx.distance,axis=0).sum(1).values
save(spec,T/'control_specificity.tsv')
pr=[];aud=[];cases=[];mm=[]
for d in sorted(Path('cache/jobs').iterdir()):
 state=rd(d/'status.tsv').iloc[0]
 if state.status!='complete':continue
 co=rd(d/'coverage.tsv.gz').set_index('gene');sources=rd(d/'KO_status.tsv').set_index('ko')
 for branch in ['heat']+(['manifold'] if d.name=='clean' else []):
  zz=rd(d/(branch+'.tsv.gz'));zz=zz[zz.ko.isin(K)]
  for k in K:
   q=zz[(zz.ko==k)&~zz.is_KO].set_index('gene')
   for g in focus:
    if g==k:continue
    why='';eff=g in q.index;reach=bool(q.at[g,'reachable']) if eff else False
    if not eff:
     if sources.at[k,'status']=='NE':why=sources.at[k,'reason']
     elif g not in co.index:why='not in local panel'
     else:why='zero variance' if co.at[g,'variance']==0 else 'low detection cells/patients'
    elif not reach:why='no directed WT path'
    aud.append(dict(analysis=d.name,branch=branch,ko=k,gene=g,native_statistic_available=eff,directed_interpretation='NE' if why else 'EVALUABLE',reason=why))
   for (n,mask),v in P0.groupby(['program','mask']):
    z=q.reindex(v.gene);e=z[z.reachable.eq(True)];good=len(e)>=3 and 2*len(e)>=len(z)
    pr.append(dict(analysis=d.name,branch=branch,ko=k,program=n,mask=mask,n_frozen=len(z),n_native=z.distance.count(),n_directed_reachable=len(e),path_gated_status='EVALUABLE' if good else 'NE',reachable_median_rank=e.rank_fraction.median() if good else np.nan,reachable_negative_fraction=(e.directionScore<0).mean() if good else np.nan))
 if d.name=='clean':
  ma=rd(d/'manifold.tsv.gz');hc=rd(d/'heat.tsv.gz')
  for k,r in cc[cc.match_order==1][['ko','control']].values:
   for name,obj in [('heat',hc),('manifold',ma)]:
    l=obj[(obj.ko==k)&~obj.is_KO];c=obj[(obj.ko==r)&~obj.is_KO];z=l.merge(c,on='gene',suffixes=('_main','_control'));mm.append(z[['gene','rank_main','rank_control','rank_fraction_main','rank_fraction_control','p.adj_main','p.adj_control']].assign(ko=k,control=r,branch=name))
save(pd.DataFrame(aud),T/'endpoint_path_coverage.tsv');save(pd.DataFrame(pr),T/'program_path_sensitivity.tsv');save(pd.concat(mm),T/'manifold_matched_comparison.tsv.gz')
# Non-restricted top intersections alongside comparisons restricted to common support.
hist=rd(V/'qc/jobs/main_effects.tsv.gz');mp=rd(T/'gene_mapping.tsv').set_index('feature_id').gene
hist=hist.assign(ko=hist.perturbed_gene.map(mp),gene=hist.gene.map(mp))
v={'v1':hist,'full':rd('cache/jobs/full/heat.tsv.gz'),'clean':a,'manifold':rd('cache/jobs/clean/manifold.tsv.gz')}
for left,right in [('v1','full'),('v1','clean'),('full','clean'),('clean','manifold')]:
 for k in K:
  g=[z[(z.ko==k)&(z.gene!=k)] for z in [v[left],v[right]]]
  if any(t.empty for t in g):continue
  for n in (20,50,100):
   u,w=[set(t.nsmallest(n,'rank').gene.dropna()) for t in g];cases.append(dict(left=left,right=right,ko=k,top_n=n,n_left=len(u),n_right=len(w),intersection=len(u&w),union=len(u|w)))
save(pd.DataFrame(cases),T/'literal_top_overlap.tsv')

