from base6 import *
from score_stats import score,check
from scipy.stats import spearmanr
from scipy.sparse.csgraph import connected_components
metrics=['H05_fraction','H10_fraction','mean_MP10']
def calc(s):
 o=R/s;d=pd.read_parquet(o/'cells_pre_MP10.parquet');sy=rd(V5/s/'genes.tsv').symbol.to_numpy();X=sparse.load_npz(o/'counts.npz');assert np.array_equal(np.asarray(X.sum(1)).ravel(),d.nCount)
 put(check(s),o/'historical_score_check.tsv')
 d['score_ok']=d.epi&d.nFeature.ge(50);ok=d.score_ok;v=score(X[ok.to_numpy()],sy);d['MP10_cell_score']=np.nan;d.loc[ok,'MP10_cell_score']=v
 th=[]
 for key,q in [('H05',.95),('H10',.90)]:
  cut=np.quantile(v,q,method='linear');d[key]=pd.array([pd.NA]*len(d),dtype='boolean');d.loc[ok,key]=v>cut
  th.append(dict(section=s,label=key,epi=d.epi.sum(),scored=len(v),n_high=(v>cut).sum(),cutoff=cut,ties=(v==cut).sum(),fraction=(v>cut).mean(),zero=(v==0).sum(),n_scores=len(np.unique(v)),depth_rho=spearmanr(v,d.loc[ok,'nCount']).statistic,median_UMI=d.loc[ok,'nCount'].median(),median_nFeature=d.loc[ok,'nFeature'].median()))
 put(pd.DataFrame(th),o/'thresholds.tsv');d.to_parquet(o/'cells_scored.parquet',index=False)
 native=pd.read_parquet(V/s/'bin_evidence.parquet',columns=['x_um','y_um','component','eligible_pool','in_tissue_original'])
 c=d[d.extended|d.control];rows=[];link=[]
 for comp,cc in c.groupby('component'):
  e=d[d.epi&d.component.eq(comp)].reset_index(drop=True);p=native[native.component.eq(comp)&native.eligible_pool];allp=native[native.component.eq(comp)&native.in_tissue_original.eq(1)]
  xy=cc[['x','y']].to_numpy();tree=cKDTree(e[['x','y']]);bc=cKDTree(p[['x_um','y_um']]);bt=cKDTree(allp[['x_um','y_um']])
  for r in [50,200]:
   for z,ix,cv,tv in zip(cc.itertuples(),tree.query_ball_point(xy,r),bc.query_ball_point(xy,r,return_length=True),bt.query_ball_point(xy,r,return_length=True)):
    neigh=e.iloc[ix];ev=neigh[neigh.score_ok];n=len(ev);cov=min(1,cv*64/(np.pi*r*r));row={k:getattr(z,k) for k in ['id','roi','x','y','A','B','main','strict','extended','control','nonTreg_type','nCount','component']}
    row.update(radius=r,coverage=cov,observed_in_tissue_coverage=min(1,tv*64/(np.pi*r*r)),epi=len(neigh),scored=n,scoreable_fraction=n/len(neigh) if len(neigh) else np.nan,eligible=(comp>=0 and cov>=.5 and n>0),mean_MP10=ev.MP10_cell_score.mean())
    for key in ['H05','H10']:
     row[key+'_n']=ev[key].sum();row[key+'_fraction']=ev[key].mean() if n else np.nan
    rows.append(row);link.append((z.id,r,neigh.id.tolist(),ev.id.tolist()))
 nb=pd.DataFrame(rows);nb.to_parquet(o/'neighborhoods.parquet',index=False)
 pd.DataFrame(link,columns=['id','radius','all_epi_ids','scored_epi_ids']).to_parquet(o/'members.parquet',index=False)
 out=[];match=[];rng=np.random.default_rng(123)
 for r,alls in nb.groupby('radius'):
  a=alls[alls.eligible].copy();features=np.column_stack([np.log1p(a.scored),np.log1p(a.nCount),a.coverage]);scale=features.std(axis=0);scale[scale==0]=1;a[['f1','f2','f3']]=features/scale
  for ct in ['other_T','CD4']:
   controls=a[a.control & (a.nonTreg_type.str.startswith('CD4') if ct=='CD4' else True)].copy();f=['f1','f2','f3']
   m={}
   for z in a[a.extended].itertuples():
    dx=np.hypot(controls.x-z.x,controls.y-z.y);candidate=controls[controls.component.eq(z.component)&dx.le(1000)].copy();candidate['distance_um']=dx[candidate.index];candidate['cov_distance']=np.linalg.norm(candidate[f].to_numpy()-np.array([z.f1,z.f2,z.f3]),axis=1)
    nn=candidate.sort_values(['cov_distance','distance_um','id']).head(3)
    if len(nn):
     m[z.id]=nn;match.extend(dict(case=z.id,control=k.id,radius=r,control_pool=ct,distance_um=k.distance_um,cov_distance=k.cov_distance,epi_opportunity=k.scored,case_epi=z.scored,case_coverage=z.coverage,control_coverage=k.coverage,T_UMI=k.nCount,case_UMI=z.nCount) for k in nn.itertuples())
   for grp in ['A','B','main','strict','extended']:
    cases=a[a[grp]];paired=cases[cases.id.isin(m)].reset_index(drop=True);n=len(paired);refs=set.union(*(set(m[i].id) for i in paired.id)) if n else set()
    # Dependence components join cases in the same 400um block, overlapping circles or sharing controls (including controls' circles).
    deps=np.eye(n,dtype=bool);supports=[]
    for z in paired.itertuples():
     nodes=np.r_[[[z.x,z.y]],m[z.id][['x','y']].to_numpy()];supports.append((nodes,set(map(tuple,np.floor(nodes/400).astype(int))),set(m[z.id].id)))
    from scipy.spatial.distance import cdist
    for j,(p,bl,ids) in enumerate(supports):
     for k,(p2,bl2,ids2) in enumerate(supports[:j]):
      deps[j,k]=deps[k,j]=bool(bl&bl2 or ids&ids2 or cdist(p,p2).min()<=2*r)
    ng,clusters=connected_components(sparse.csr_matrix(deps),directed=False)
    for metric in metrics:
     dif=np.array([getattr(z,metric)-m[z.id][metric].mean() for z in paired.itertuples()]);mean=dif.mean() if n else np.nan;lo=hi=np.nan
     if ng>=5:
      sums=np.bincount(clusters,weights=dif);nums=np.bincount(clusters);boot=rng.integers(ng,size=(2000,ng));lo,hi=np.quantile(sums[boot].sum(1)/nums[boot].sum(1),[.025,.975])
     out.append(dict(section=s,radius=r,group=grp,controls=ct,metric=metric,total_cases=alls[grp].sum(),eligible_cases=len(cases),matched=n,control_used=len(refs),control_eligible=len(controls),dependency_blocks=ng,matched_case_mean=paired[metric].mean(),matched_control_mean=np.mean([m[i][metric].mean() for i in paired.id]) if n else np.nan,difference=mean,lower=lo,upper=hi,reason='no case' if not alls[grp].any() else 'no eligible pair' if not n else 'fewer than five independent block groups' if ng<5 else 'block cluster bootstrap'))
 pd.DataFrame(match,columns=['case','control','radius','control_pool','distance_um','cov_distance','epi_opportunity','case_epi','case_coverage','control_coverage','T_UMI','case_UMI']).to_parquet(o/'matches.parquet',index=False);put(pd.DataFrame(out),o/'effects.tsv')
if __name__=='__main__':
 for s in sections:calc(s)
