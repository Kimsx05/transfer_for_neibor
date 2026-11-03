import os,sys
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'spatial_ecology_top20_v2/runtime'))
import numpy as np,pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr,norm
from concurrent.futures import ProcessPoolExecutor
O=ROOT/'mp10_paper_distance_subtypes_v2';F=ROOT/'mp10_paper_distance_full_v1'
M=pd.read_csv(F/'tables/input_manifest.tsv',sep='\t')
SETS=pd.read_csv(F/'tables/gene_sets.tsv.gz',sep='\t');PROG=sorted(SETS.gs_name.unique());RAW='EpiMP10_Top20_AUCell'
TARGETS={k:[k] for k in ['Matrix_CAF','iCAF','rCAF','vCAF','C1QC TAM','SPP1 TAM','Treg']}
TARGETS.update(CAF=['Matrix_CAF','iCAF','rCAF','vCAF'],TAM=['C1QC TAM','SPP1 TAM'])
def write(d,p):pd.DataFrame(d).to_csv(p,sep='\t',index=False,na_rep='NA')
def nearest(d,mask):
 out=np.full(len(d),np.nan);xy=d[['x_um','y_um']].to_numpy()
 for _,idx in d.groupby('tissue_component_id').groups.items():
  t=idx[mask[idx]]
  if len(t):out[idx]=cKDTree(xy[t]).query(xy[idx])[0]
 return out
def interface(d,r):
 old=pd.read_parquet(Path(r.output_dir)/'spot_region_labels.parquet').set_index('barcode').reindex(d.barcode)
 author=d.author_panel_class_at_center
 ep=(d.Epithelial_fraction>=.5)&(old.marker_epithelial.to_numpy()>old.marker_fibroblast.to_numpy())
 st=(d.CAF_fraction>=.2)&(old.marker_fibroblast.to_numpy()>old.marker_epithelial.to_numpy())
 ep &= author.isin(['Tumor Cell','Urothelium'])|old.tissue_region.isin(['Tumor','Normal urothelium']).to_numpy()|old.detected_epithelial.ge(2).to_numpy()
 st &= author.isin(['Tumor Stroma','Normal Stroma'])|old.tissue_region.eq('Non-tumor stromal region').to_numpy()|old.detected_fibroblast.ge(2).to_numpy()
 ep &= ~author.isin(['Tumor Stroma','Normal Stroma']);st &= ~author.isin(['Tumor Cell','Urothelium'])
 side=np.where(ep,1,np.where(st,-1,0));side[~d.center_qc.to_numpy()]=0
 xy=d[['x_um','y_um']].to_numpy();lookup=dict(zip(zip(d.array_row,d.array_col),d.index));edges=[];adj=[[] for _ in range(len(d))]
 for i,row in d.iterrows():
  for dr,dc in [(0,2),(1,1),(1,-1),(-1,-1),(-1,1),(0,-2)]:
   j=lookup.get((row.array_row+dr,row.array_col+dc))
   if j is None or d.tissue_component_id[i]!=d.tissue_component_id[j]:continue
   adj[i].append(j)
   if i<j and side[i]*side[j]==-1:
    v=xy[j]-xy[i];mid=(xy[i]+xy[j])/2;perp=np.array([-v[1],v[0]])/np.linalg.norm(v)*100/np.sqrt(3)/2
    edges.append([i,j,*(mid-perp),*(mid+perp),d.tissue_component_id[i],d.lesion_location[i]=='inside' and d.lesion_location[j]=='inside'])
 dist=np.full(len(d),np.nan);hops=np.full(len(d),np.nan)
 for comp,idx in d.groupby('tissue_component_id').groups.items():
  e=[a for a in edges if a[6]==comp]
  if not e:continue
  e=np.array([a[:6] for a in e],float);a=e[:,2:4];v=e[:,4:6]-a
  for chunk in np.array_split(idx,max(1,len(idx)//100)):
   delta=xy[chunk,None,:]-a;u=np.clip((delta*v).sum(2)/(v*v).sum(1),0,1)
   dist[chunk]=np.linalg.norm(delta-u[:,:,None]*v,axis=2).min(1)
  front=set(e[:,:2].astype(int).ravel());hops[list(front)]=1
  second={j for i in front for j in adj[i]}-front;hops[list(second)]=2
 d['local_side']=side;d['local_interface_um']=dist*np.where(side==0,np.nan,side);d['local_unsigned_um']=dist;d['local_ring']=hops
 d['local_reason']=np.where(~np.isfinite(dist),'no_supported_interface_in_component',np.where(side==0,'side_uncertain','candidate_geometry_distance'))
 write(edges,O/'tables'/f'{r.section_id}_edges.tsv')
 return d,dict(section_id=r.section_id,n_edges=len(edges),n_internal_edges=sum(e[7] for e in edges),n_components=len(set(e[6] for e in edges)),n_spots=np.isfinite(d.local_interface_um).sum(),n_queries=(np.isfinite(d.local_interface_um)&d['query']).sum())
def stats(d,x,Y,interaction=False):
 ok=d['query'].to_numpy()&np.isfinite(x)&np.isfinite(Y).all(1)&d.nCount_RNA.notna().to_numpy()
 n=ok.sum();p=Y.shape[1]
 if n<20:return None
 y=Y[ok];xx=np.asarray(x)[ok];mp=d[RAW].to_numpy()[ok];epi=d.Epithelial_fraction.to_numpy()[ok];depth=np.log1p(d.nCount_RNA.to_numpy()[ok])
 if np.std(xx)<1e-10:return None
 z=lambda v:(v-v.mean(0))/np.where(v.std(0)>1e-12,v.std(0),1)
 y=z(y);xx=z(xx);mp=z(mp)
 C=np.column_stack([np.ones(n),z(epi),z(depth),pd.get_dummies(d.loc[ok,'lesion_location'],drop_first=True,dtype=float)])
 resid=lambda v:v-C@np.linalg.lstsq(C,v,rcond=None)[0]
 yr=resid(y);xr=resid(xx)
 cor=(xx@y)/n;partial=(xr@yr)/np.maximum(np.linalg.norm(xr)*np.linalg.norm(yr,axis=0),1e-15)
 X=np.column_stack([C,mp,xx,xx*mp]) if interaction else np.column_stack([C,xx])
 if np.linalg.matrix_rank(X)<X.shape[1]: # remove redundant covariates only
  keep=[]
  for j in range(X.shape[1]- (3 if interaction else 1)):
   if np.linalg.matrix_rank(X[:,keep+[j]])>len(keep):keep.append(j)
  X=X[:,keep+list(range(X.shape[1]-(3 if interaction else 1),X.shape[1]))]
 if np.linalg.matrix_rank(X)<X.shape[1]:return None
 inv=np.linalg.pinv(X.T@X);b=inv@X.T@y;err=y-X@b
 block=(d.loc[ok,'tissue_component_id'].astype(str)+'|'+(d.loc[ok,'x_um']//500).astype(str)+'|'+(d.loc[ok,'y_um']//500).astype(str))
 ids,labels=pd.factorize(block);B=len(labels)
 if B<8:return dict(n=n,blocks=B,r=cor,adj=partial,beta=b[-1],low=np.full(p,np.nan),high=np.full(p,np.nan),p=np.full(p,np.nan))
 influence=(X@inv[-1])[:,None]*err;g=np.zeros((B,p));np.add.at(g,ids,influence)
 rng=np.random.default_rng(20261103);boots=rng.choice([-1,1],size=(399,B))@g
 se=g.std(0,ddof=1)*np.sqrt(B)
 pv=(1+(np.abs(boots)>=np.abs(b[-1])).sum(0))/400
 return dict(n=n,blocks=B,r=cor,adj=partial,beta=b[-1],low=b[-1]-np.quantile(boots,.975,axis=0),high=b[-1]-np.quantile(boots,.025,axis=0),p=pv)
def run(sec):
 r=M.set_index('section_id').loc[sec];r.section_id=sec
 d=pd.read_parquet(F/'spots'/f'{sec}.parquet');q05=pd.read_csv(r.c2l_q05_path,index_col=0).reindex(d.barcode)
 d,ic=interface(d,r);targets=[];results=[];curves=[];strata=[];niches=[]
 means=d.filter(regex='^mean__');den=means.sum(1);valid=d.center_qc & d.lesion_location.ne('quality_excluded')&means.notna().all(1)&den.gt(0)
 axes={}
 for target,factors in TARGETS.items():
  ab=d[['mean__'+f for f in factors]].sum(1);v=ab/den
  lineage=d.CAF_abundance if target in list(TARGETS)[:4] else d.TAM_abundance if 'TAM' in target else d.T_NK_fraction*den
  if target not in ['CAF','TAM']:
   d[target+'_within_lineage']=ab/lineage.where(lineage>0)
   d[target+'_lineage_denominator_lt1']=lineage<1
  for q in [95,90]:
   th=v[valid].quantile(q/100);mask=(valid & v.ge(th)&v.gt(0)).to_numpy() if v[valid].nunique()>1 else np.zeros(len(d),bool)
   name=target+'_P'+str(q);dist=nearest(d,mask);axes[name]=dist;d[name+'_distance']=dist;d[name+'_target']=mask
   rec=dict(section_id=sec,target=target,quantile=q,n_eligible=valid.sum(),threshold=th,n_target=mask.sum(),fraction=mask.sum()/valid.sum(),abundance_median=ab[mask].median(),abundance_q25=ab[mask].quantile(.25),abundance_q75=ab[mask].quantile(.75),composition_median=v[mask].median(),n_target_components=d.loc[mask,'tissue_component_id'].nunique(),n_query=d['query'].sum(),n_pairs=(d['query']&np.isfinite(dist)).sum(),zero_fraction=(dist[d['query']]==0).sum()/max(1,np.isfinite(dist[d['query']]).sum()),no_target_query=(d['query']&~np.isfinite(dist)).sum(),low_abundance=ab[mask].median()<1)
   for f in factors:rec['marginal_q05_'+f]=q05.loc[mask,f].median()
   targets.append(rec)
 axes.update(local_interface=d.local_interface_um.to_numpy(),observed_outer=d.distance_to_observed_outer_candidate_um.to_numpy(),strict_outer=d.signed_distance_um.to_numpy())
 Y=d[[RAW]+PROG].to_numpy()
 for axis,x in axes.items():
  for kind in ['distance','interaction']:
   fit=stats(d,x,Y,kind=='interaction')
   for j,name in enumerate([RAW]+PROG):
    if kind=='interaction' and name==RAW:continue
    z=dict(section_id=sec,axis=axis,program=name,kind=kind)
    if fit:
     z.update({k:(v[j] if isinstance(v,np.ndarray) else v) for k,v in fit.items()})
     ok=d['query']&np.isfinite(x);z['rho']=spearmanr(np.asarray(x)[ok],Y[ok,j]).statistic
    results.append(z)
  # Actual distance bins, all tissue composition and Epi-only programs
  bins=pd.cut(x,[-np.inf,-1000,-500,-200,0,100,200,500,1000,np.inf])
  for b,idx in d.groupby(bins,observed=True).groups.items():
   z=d.loc[idx];query=z[z['query']]
   curves.append(dict(section_id=sec,axis=axis,bin=str(b),n_tissue=len(z),n_query=len(query),**{k:query[k].median() if len(query)>=5 else np.nan for k in [RAW]+PROG},**{k:z['fraction__'+k].mean() for k in list(TARGETS)[:7]}))
   for version,col,eligible in [('R12','accepted_label',z.query if False else z['epi_anchor']),('local','local_niche_id',valid.loc[idx])]:
    lab=z[eligible&z[col].notna()];n=len(lab)
    for label,g in lab.groupby(col):
     gq=g[g['query']]
     niches.append(dict(section_id=sec,axis=axis,bin=str(b),version=version,niche=label,eligible=eligible.sum(),labelled=n,n=len(g),fraction=len(g)/n,n_MP10=len(gq),MP10=gq[RAW].median() if len(gq)>=5 else np.nan))
 fit=stats(d,d[RAW].to_numpy(),d[PROG].to_numpy())
 for j,name in enumerate(PROG):
  z=dict(section_id=sec,axis='MP10',program=name,kind='program')
  if fit:z.update({k:(v[j] if isinstance(v,np.ndarray) else v) for k,v in fit.items()})
  results.append(z)
 for col in ['lesion_location','distance_band_200um','accepted_label','local_niche_id']:
  for label,g in d[d['query']].groupby(col):
   if len(g)<20:continue
   for name in PROG:
    strata.append(dict(section_id=sec,stratum=col,label=label,program=name,n=len(g),r=g[[RAW,name]].corr().iloc[0,1]))
 # Matched model layers
 for axis,x in axes.items():
  ok=d['query']&d.epi_score_valid.eq(True)&d.spot_MP10_model_universe.notna()&d.epi_MP10_model_universe.notna()&np.isfinite(x)
  for name in ['spot_MP10_model_universe','epi_MP10_model_universe']:
   results.append(dict(section_id=sec,axis=axis,program=name,kind='matched',n=ok.sum(),r=np.corrcoef(np.asarray(x)[ok],d.loc[ok,name])[0,1] if ok.sum()>=20 else np.nan))
 d.to_parquet(O/'spots'/f'{sec}.parquet',index=False)
 for data,name in [(targets,'targets'),(results,'effects'),(curves,'curves'),(strata,'strata'),(niches,'niche'),([ic],'interface')]:write(data,O/'tables'/f'{sec}_{name}.tsv.gz')
 return sec
if __name__=='__main__':
 with ProcessPoolExecutor(max_workers=8) as pool:
  for s in pool.map(run,M.section_id):print(s,flush=True)
