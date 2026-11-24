from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
O=Path(__file__).resolve().parents[1]; C=Path.cwd()/'cache/mp10_public_resources_v1'
def output(d,name): d.to_csv(O/name,index=False,sep='\t',na_rep='NA')
def read(name):return pd.read_csv(O/name,sep='\t')
sets={r.module:r.genes.split() for r in read('scripts/modules.tsv').itertuples()}; top=sets['MP10']
def adjust(d,group):
 d['q']=np.nan
 for _,g in d.dropna(subset=['p']).groupby(group):
  d.loc[g.index,'q']=stats.false_discovery_control(g.p) if hasattr(stats,'false_discovery_control') else np.minimum.accumulate((g.p.sort_values().values*len(g)/np.arange(1,len(g)+1))[::-1])[::-1][np.argsort(np.argsort(g.p.values))].clip(0,1)
 return d
def associations(data,platform):
 # One observation per patient/lineage/drug, mean across derived models with both assays.
 groups=['endpoint','lineage','unit','eligible']
 d=data.dropna(subset=['score','response','lineage']).groupby(groups,as_index=False).agg(score=('score','mean'),response=('response','mean'),models=('name',lambda s:'|'.join(map(str,s))))
 output(d,'plot_data/'+platform+'_units.tsv')
 rng=np.random.default_rng(627); rows=[]; omissions=[];dirs=[]
 for endpoint,x in d.groupby('endpoint'):
  for lin,y in x.groupby('lineage'):
   if len(y)>=3 and min(y.score.nunique(),y.response.nunique())>1:
    r,p=stats.spearmanr(y.score,y.response);dirs.append(dict(endpoint=endpoint,lineage=lin,n=len(y),rho=r,p=p))
  for cohort,y in [('BLCA',x[x.eligible]),('pan_lineage',x.groupby('lineage').filter(lambda a:len(a)>=3 and min(a.score.nunique(),a.response.nunique())>1))]:
   if len(y)<3:continue
   a,b=y.score.values,y.response.values
   if cohort=='BLCA':
    r,p=stats.spearmanr(a,b); idx=rng.integers(0,len(a),(2000,len(a)))
    ar=stats.rankdata(a[idx],axis=1);br=stats.rankdata(b[idx],axis=1);ar-=ar.mean(1,keepdims=True);br-=br.mean(1,keepdims=True); den=np.sqrt((ar**2).sum(1)*(br**2).sum(1)); boot=np.divide((ar*br).sum(1),den,out=np.full(len(idx),np.nan),where=den>0)
    for j,row in enumerate(y.itertuples()):
     omissions.append(dict(endpoint=endpoint,omit=row.models,rho=stats.spearmanr(np.delete(a,j),np.delete(b,j))[0]))
   else:
    # Percentile ranks and subtraction of lineage means. Stratified model bootstrap re-ranks each lineage.
    ar=[];br=[]; chunks=[]
    for l,g in y.groupby('lineage'):
     aa,bb=g.score.values,g.response.values
     def rr(z):return stats.rankdata(z,axis=-1)/z.shape[-1]-(z.shape[-1]+1)/(2*z.shape[-1])
     ax,bx=rr(aa),rr(bb);ar.extend(ax);br.extend(bx)
     idx=rng.integers(0,len(g),(1000,len(g)));rx,ry=rr(aa[idx]),rr(bb[idx]);chunks.append(np.stack([(rx*ry).sum(1),(rx**2).sum(1),(ry**2).sum(1)]))
    num,s1,s2=np.sum(chunks,axis=0);boot=num/np.sqrt(s1*s2)
    ar,br=np.array(ar),np.array(br);r=np.corrcoef(ar,br)[0,1];df=len(y)-y.lineage.nunique()-1;p=2*stats.t.sf(abs(r)*np.sqrt(df/(1-r*r)),df)
   lo,hi=np.nanquantile(boot,[.025,.975]);rows.append(dict(endpoint=endpoint,cohort=cohort,n=len(y),lineages=y.lineage.nunique(),effect=r,lo=lo,hi=hi,p=p,score_min=min(a),score_max=max(a),response_min=min(b),response_max=max(b)))
 output(adjust(pd.DataFrame(rows),'cohort'),'tables/'+platform+'_associations.tsv');output(pd.DataFrame(omissions),'tables/'+platform+'_influence.tsv');output(adjust(pd.DataFrame(dirs),'endpoint'),'tables/'+platform+'_within_lineage.tsv')

