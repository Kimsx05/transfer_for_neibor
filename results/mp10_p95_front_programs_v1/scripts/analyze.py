import os,sys
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
R=Path(__file__).resolve().parents[3];O=R/'results/mp10_p95_front_programs_v1';F=R/'mp10_paper_distance_full_v1'
sys.path.insert(0,str(R/'spatial_ecology_top20_v2/runtime'))
import numpy as np,pandas as pd
import warnings
warnings.filterwarnings("ignore",category=pd.errors.PerformanceWarning)
from scipy import sparse,stats
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree
from concurrent.futures import ProcessPoolExecutor
M=pd.read_csv(O/'tables/input_manifest.tsv',sep='\t');G=pd.read_csv(O/'tables/gene_sets.tsv',sep='\t');P=sorted(G.gs_name.unique());RAW='EpiMP10_Top20_AUCell';PE='PURAM_PEMT_TOP30'
RINGS=['R0','R1','R2','R3_R5','R6_R10'];RB=[(0,0),(1,1),(2,2),(3,5),(6,10)]
def save(a,path):pd.DataFrame(a).to_csv(path,sep='\t',index=False,na_rep='NA')
def z(a):
 a=np.asarray(a,float);sd=np.nanstd(a,axis=0);return (a-np.nanmean(a,axis=0))/np.where(sd>1e-12,sd,1)
def frontband(x):
 return np.select([~np.isfinite(x),x< -500,x< -200,x<=0,x<=200,x<=500],['unavailable','deep','near_inside','front','outer_200','outer_500'],default='outer_far')
def graph(d,valid):
 look={(r.array_row,r.array_col):i for i,r in d.iterrows() if valid[i]};a=[];b=[]
 for (rr,cc),i in look.items():
  for dr,dc in [(0,2),(0,-2),(1,1),(1,-1),(-1,1),(-1,-1)]:
   j=look.get((rr+dr,cc+dc))
   if j is not None and pd.notna(d.tissue_component_id[i]) and d.tissue_component_id[i]==d.tissue_component_id[j]:a.append(i);b.append(j)
 return sparse.csr_matrix((np.ones(len(a)),(a,b)),shape=(len(d),len(d)))
def near(d,mask,valid):
 out=np.full(len(d),np.nan)
 for _,ids in d[valid].groupby('tissue_component_id').groups.items():
  t=ids[mask[ids]]
  if len(t):out[ids]=cKDTree(d.loc[t,['x_um','y_um']]).query(d.loc[ids,['x_um','y_um']])[0]
 return out

def fit(d,x,Y,names,mask,tag,family,adjust=True,binary=False,raw_y=False,extra=None,blocksize=500,quadratic=False):
 x=np.asarray(x,float);Y=np.asarray(Y,float);C=np.ones((len(d),1))
 if extra is not None:
  ev=np.asarray(extra,float);mask=np.asarray(mask,bool)&(np.isfinite(ev).all(1) if ev.ndim==2 else np.isfinite(ev))
 if adjust:
  C=np.column_stack([C,d.Epithelial_fraction.to_numpy(float),np.log1p(d.nCount_RNA.to_numpy(float))])
 if extra is not None and adjust:C=np.column_stack([C,extra])
 out=[]
 # Per-program completeness, retaining empty or constant programs explicitly.
 patterns={}
 for j in range(len(names)):
  ok=np.asarray(mask,bool)&np.isfinite(x)&np.isfinite(Y[:,j])&np.isfinite(C).all(1)
  patterns.setdefault(ok.tobytes(),[ok,[]])[1].append(j)
 for ok,js in patterns.values():
  n=ok.sum();base=dict(analysis=tag,family=family,adjustment='adjusted' if adjust else 'raw',n=n)
  if n<20 or np.unique(x[ok]).size<2:
   out.extend([dict(base,program=names[j],reason='insufficient_pairs_or_constant_predictor') for j in js]);continue
  xx=x[ok] if binary else z(x[ok]);y=Y[ok][:,js];sd=y.std(0);yc=y if raw_y else z(y)
  c=C[ok];c[:,1:]=z(c[:,1:]);keep=[]
  for j in range(c.shape[1]):
   if np.linalg.matrix_rank(c[:,keep+[j]])>len(keep):keep.append(j)
  c=c[:,keep];X=np.column_stack([c,xx,xx**2]) if quadratic else np.column_stack([c,xx]);k=X.shape[1]
  if np.linalg.matrix_rank(X)<k:
   out.extend([dict(base,program=names[j],reason='rank_deficient') for j in js]);continue
  inv=np.linalg.inv(X.T@X);coef=inv@X.T@yc;err=yc-X@coef
  groups=d.loc[ok,'tissue_component_id'].astype(str)+'|'+(d.loc[ok,'x_um']//blocksize).astype(str)+'|'+(d.loc[ok,'y_um']//blocksize).astype(str)
  bid,labs=pd.factorize(groups);B=len(labs);m=2 if quadratic else 1
  influence=(X@inv[:,-m:])[:,:,None]*err[:,None,:];sums=np.zeros((B,m,len(js)));np.add.at(sums,bid,influence)
  cov=np.einsum('bip,bjp->ijp',sums,sums)*(B/max(B-1,1))*((n-1)/max(n-k,1));se=np.sqrt(np.maximum(np.diagonal(cov,axis1=0,axis2=1).T,0));b=coef[-m:];ci=stats.t.ppf(.975,B-1) if B>=8 else np.nan
  xr=xx-c@np.linalg.lstsq(c,xx,rcond=None)[0];yr=yc-c@np.linalg.lstsq(c,yc,rcond=None)[0]
  rr=(z(xx)@z(y))/n;ar=xr@yr/np.maximum(np.linalg.norm(xr)*np.linalg.norm(yr,axis=0),1e-15)
  for jj,j in enumerate(js):
   if sd[jj]<1e-12:out.append(dict(base,program=names[j],reason='constant_response',blocks=B));continue
   pp=2*stats.t.sf(abs(b[0,jj]/max(se[0,jj],1e-15)),B-1) if B>=8 else np.nan
   overall=stats.f.sf(b[:,jj]@np.linalg.pinv(cov[:,:,jj])@b[:,jj]/m,m,B-1) if B>=8 else np.nan
   out.append(dict(base,program=names[j],r=rr[jj],partial_r=ar[jj],effect=b[0,jj],low=b[0,jj]-ci*se[0,jj],high=b[0,jj]+ci*se[0,jj],p=pp,overall_p=overall,quadratic_beta=b[1,jj] if quadratic else np.nan,blocks=B,block_um=blocksize,reason='ok' if B>=8 else 'effect_only_fewer_than_8_spatial_blocks'))
 return out

def match(d,hi,pool,position=False):
 rng=np.random.default_rng(20261103);h=np.flatnonzero(hi);rng.shuffle(h);available=set(np.flatnonzero(pool));pairs=[]
 cov=z(np.column_stack([np.log1p(d['mean__Epithelial']),np.log1p(d.nCount_RNA)]))
 for i in h:
  cand=np.array([j for j in available if d.tissue_component_id[j]==d.tissue_component_id[i]],int)
  if position:
   cand=cand[(d.front_band.iloc[cand].to_numpy()==d.front_band[i])&np.isfinite(d.d_front.iloc[cand])&np.isfinite(d.d_front[i])&(np.abs(d.d_front.iloc[cand].to_numpy()-d.d_front[i])<=200)]
  if not len(cand) or not np.isfinite(cov[i]).all():continue
  ds=((cov[cand]-cov[i])**2).sum(1);j=cand[np.argmin(ds)]
  if not np.isfinite(ds.min()):continue
  available.remove(j);pairs.append(dict(high_index=i,control_index=j,high_barcode=d.barcode[i],control_barcode=d.barcode[j],match_distance=np.sqrt(ds.min()),epi_log_difference=np.log1p(d['mean__Epithelial'][i])-np.log1p(d['mean__Epithelial'][j]),depth_log_difference=np.log1p(d.nCount_RNA[i])-np.log1p(d.nCount_RNA[j]),front_difference=d.d_front[i]-d.d_front[j],matching='front' if position else 'base'))
 return pairs

def run(sec):
 done=O/'tables'/f'{sec}_coverage.tsv'
 if done.exists():return sec+' cached'
 d=pd.read_parquet(F/'spots'/f'{sec}.parquet').reset_index(drop=True);r=M.set_index('section_id').loc[sec]
 new=pd.read_csv(O/'scores'/f'{sec}_programs.tsv.gz',sep='\t');d=d.merge(new,on=['section_id','barcode'],how='left',validate='1:1')
 assert len(d)==r.n_spots
 means=d.filter(regex='^mean__');valid=(means.ge(0)&np.isfinite(means)).all(1)&means.sum(1).gt(0)&d.center_qc&d.lesion_location.ne('quality_excluded')
 q=d['query'].to_numpy(bool);tumor=q&d.lesion_location.eq('inside').to_numpy()&(r.source_context!='adjacent_normal_source')
 df=d.signed_distance_um.to_numpy(float);d['d_front']=df;d['front_band']=frontband(df)
 th=d.loc[tumor,RAW].quantile(.95);hi=tumor&d[RAW].ge(th).to_numpy() if np.isfinite(th) else np.zeros(len(d),bool)
 p75=d.loc[tumor,PE].quantile(.75);fr=tumor&(df>=-200)&(df<=0);active=fr&d[PE].ge(p75).to_numpy()
 d['tumor_epi']=tumor;d['MP10_P95']=hi;d['pEMT_high_front']=active;d['d_MP10']=near(d,hi,valid);d['d_MP10_reason']=np.where(~valid,'ineligible_tissue',np.where(np.isfinite(d.d_MP10),'ok','no_P95_in_component'))
 adj=graph(d,valid);refs={'P95':hi};pairs=[]
 for pos in [False,True]:
  pair=match(d,hi,tumor&~hi,pos);pairs.extend(pair);m=np.zeros(len(d),bool)
  if pair:m[[p['control_index'] for p in pair]]=True
  name='control_front' if pos else 'control_base';refs[name]=m;d['reference_'+name]=m;d['d_'+name]=near(d,m,valid)
  hm=np.zeros(len(d),bool)
  if pair:hm[[p['high_index'] for p in pair]]=True
  refs['P95_matched_front' if pos else 'P95_matched_base']=hm
 for name,m in refs.items():
  if m.any():hop=dijkstra(adj,directed=False,indices=np.flatnonzero(m),unweighted=True,min_only=True);hop[~np.isfinite(hop)]=np.nan;hop[~valid]=np.nan
  else:hop=np.full(len(d),np.nan)
  d['hops_'+name]=hop
  if name.startswith('P95_matched'):d['d_'+name]=near(d,m,valid)
 Y=d[P].to_numpy(float);mp=d[RAW].to_numpy(float);effects=[];curves=[];counts=[];niches=[];ringrows=[]
 def add(x,y,names,mask,tag,family,**kw):
  for a in [False,True]:effects.extend(fit(d,x,y,names,mask,tag,family,adjust=a,**kw))
 # All Epi centers: true shortest-path rings. No complete-radius requirement.
 qi=np.flatnonzero(q);hh=dijkstra(adj,directed=False,indices=qi,unweighted=True,limit=10)
 phys=np.linalg.norm(d.loc[qi,['x_um','y_um']].to_numpy()[:,None,:]-d[['x_um','y_um']].to_numpy()[None,:,:],axis=2)
 center=d.loc[qi,['section_id','barcode',RAW,'d_front','accepted_label','tumor_epi','MP10_P95']].copy();ringmats={}
 for ring,(lo,up) in zip(RINGS,RB):
  ii,jj=np.where((hh>=lo)&(hh<=up)&valid.to_numpy()[None,:]);w=sparse.csr_matrix((np.ones(len(ii)),(ii,jj)),shape=(len(qi),len(d)))
  denom=np.asarray(w@np.isfinite(Y).astype(float));val=(w@np.nan_to_num(Y))/np.where(denom>0,denom,np.nan)
  full=np.full((len(d),len(P)),np.nan);full[qi]=val;ringmats[ring]=full
  nc=np.asarray(w.sum(1)).ravel();avg_epi=np.asarray(w@d.Epithelial_fraction.fillna(0).to_numpy()).ravel()/np.where(nc>0,nc,np.nan);extra=np.full(len(d),np.nan);extra[qi]=avg_epi
  center[ring+'_n_neighbors']=nc;center[ring+'_mean_distance_um']=np.bincount(ii,weights=phys[ii,jj],minlength=len(qi))/np.where(nc>0,nc,np.nan)
  for j,p in enumerate(P):center[ring+'__'+p]=val[:,j]
  add(mp,full,P,q,'B1_'+ring,'B_center' if ring=='R0' else 'B_neighborhood',extra=None if ring=='R0' else extra,blocksize=500 if ring=='R0' else 2000)
  # P95 vs non-P95 at comparable front position, all available tumors.
  if ring in RINGS:
   ex=np.column_stack([z(df),z(df)**2])
   add(hi.astype(float),full,P,tumor,'B4_position_'+ring,'B_neighborhood',extra=ex,binary=True,blocksize=2000)
  for typ in ['base','front']:
   pp=[a for a in pairs if a['matching']==typ]
   if not pp:continue
   h=np.array([a['high_index'] for a in pp]);c=np.array([a['control_index'] for a in pp]);diff=full[h]-full[c]
   for j,p in enumerate(P):
    v=diff[:,j];ok=np.isfinite(v);sd=np.nanstd(full[tumor,j]);n=ok.sum()
    effects.append(dict(analysis='B4_matched_'+typ+'_'+ring,family='B_neighborhood',adjustment='matched',program=p,n=n,effect=np.nanmean(v)/sd if n>=5 and sd>1e-12 else np.nan,raw_difference=np.nanmean(v) if n else np.nan,reason='paired_section_effect_hierarchical_inference' if n>=5 else 'insufficient_pairs'))
 center.to_parquet(O/'centers'/f'{sec}.parquet',index=False);del hh,phys
 # A and C: strict formal front; signed distance, quadratic curves, bands.
 allnames=[RAW]+P;allY=d[allnames].to_numpy(float)
 add(df,Y,P,tumor,'A_front_distance','A_front',quadratic=True)
 add(fr.astype(float),Y,P,tumor&np.isfinite(df),'A_front_vs_internal','A_front',binary=True)
 add(df,mp[:,None],[RAW],tumor,'C1_continuous','C_MP10_front',quadratic=True)
 add(fr.astype(float),mp[:,None],[RAW],tumor&np.isfinite(df),'C1_front_vs_internal','C_MP10_front',binary=True)
 add(fr.astype(float),hi.astype(float)[:,None],['P95_probability'],tumor&np.isfinite(df),'C2_P95_enrichment','C_MP10_front',binary=True,raw_y=True)
 add(df,Y,P,hi,'C3_P95_state','C_MP10_front',quadratic=True)
 add(active.astype(float),mp[:,None],[RAW],fr,'C_active_front_MP10','C_MP10_front',binary=True,extra=np.column_stack([df,df**2]))
 # niche probabilities, all nine labels. No missing to zero assignment.
 labelled=q&d.accepted_label.notna().to_numpy();ny=np.column_stack([d.accepted_label.eq(k).fillna(False).to_numpy(float) for k in range(9)]);nn=['R12_'+str(k) for k in range(9)]
 add(fr.astype(float),ny,nn,tumor&labelled&np.isfinite(df),'D_front_enrichment','D_niche_front',binary=True,raw_y=True)
 add(active.astype(float),ny,nn,fr&labelled,'D_active_front','D_niche_front',binary=True,raw_y=True,extra=np.column_stack([df,df**2]))
 add(mp,ny,nn,tumor&labelled&np.isfinite(df),'D_continuous_MP10_position','D_MP10_niche',raw_y=True,extra=np.column_stack([df,df**2]))
 add(hi.astype(float),ny,nn,tumor&labelled&np.isfinite(df),'D_P95_position','D_MP10_niche',binary=True,raw_y=True,extra=np.column_stack([df,df**2]))
 for typ in ['base','front']:
  pp=[a for a in pairs if a['matching']==typ and labelled[a['high_index']] and labelled[a['control_index']]]
  for k in range(9):
   dv=np.array([ny[a['high_index'],k]-ny[a['control_index'],k] for a in pp])
   effects.append(dict(analysis='D_matched_'+typ,family='D_MP10_niche',adjustment='matched',program=nn[k],n=len(dv),effect=dv.mean() if len(dv)>=5 else np.nan,reason='paired_proportion_difference_hierarchical_inference'))
 # P95 unique target distance fields, per-reference curves and formal effects.
 for name,m in refs.items():
  dist=d['d_MP10' if name=='P95' else 'd_'+name].to_numpy(float)
  for target,mask in [('tissue',valid.to_numpy()),('Epi',q)]:
   if name=='P95':add(dist,Y,P,mask&np.isfinite(dist)&(dist<=1000),'B3_'+target+'_distance','B_neighborhood',quadratic=True,blocksize=2000)
   bins=np.where(dist==0,0,np.ceil(dist/100));bins[dist>1000]=11;bins[~np.isfinite(dist)]=np.nan
   for b in np.unique(bins[np.isfinite(bins)&mask]):
    sel=mask&(bins==b);n=sel.sum()
    for j,p in enumerate(P):curves.append(dict(axis='d_MP10',reference=name,population=target,bin=int(b),distance=float(np.nanmedian(dist[sel])),program=p,n=n,mean=np.nanmean(Y[sel,j]) if n>=5 else np.nan,median=np.nanmedian(Y[sel,j]) if n>=5 else np.nan))
   hop=d['hops_'+name].to_numpy(float)
   for ring,(lo,up) in zip(RINGS,RB):
    sel=mask&(hop>=lo)&(hop<=up);n=sel.sum()
    for j,p in enumerate(P):ringrows.append(dict(reference=name,population=target,ring=ring,program=p,n=n,mean=np.nanmean(Y[sel,j]) if n>=5 else np.nan,distance=np.nanmedian(dist[sel]) if n else np.nan))
 # Front band and 100um curves, separate outside tissue and query/P95.
 for pop,mask in [('tumor_Epi',tumor),('P95_Epi',hi),('outside_tissue',valid.to_numpy()&d.lesion_location.eq('outside').to_numpy())]:
  bins=np.floor(df/100)*100
  for b in np.unique(bins[np.isfinite(bins)&mask]):
   sel=mask&(bins==b);n=sel.sum()
   for j,p in enumerate(allnames):curves.append(dict(axis='d_front',reference='v3.4',population=pop,bin=b,distance=np.nanmedian(df[sel]),program=p,n=n,mean=np.nanmean(allY[sel,j]) if n>=5 else np.nan,median=np.nanmedian(allY[sel,j]) if n>=5 else np.nan))
 for band in ['deep','near_inside','front','outer_200','outer_500','outer_far','unavailable','pEMT_high_front','other_front']:
  sel=q&(d.front_band.eq(band).to_numpy()) if band not in ['pEMT_high_front','other_front'] else (active if band=='pEMT_high_front' else fr&~active)
  lab=sel&labelled;nt=(sel&tumor).sum();nh=(sel&hi).sum()
  counts.append(dict(band=band,n_eligible=sel.sum(),n_tumor_epi=nt,n_P95=nh,P95_fraction=nh/nt if nt else np.nan,n_labelled=lab.sum(),n_unlabelled=(sel&~labelled).sum()))
  for k in range(9):
   sub=lab&d.accepted_label.eq(k).fillna(False).to_numpy(bool)
   niches.append(dict(band=band,niche=k,n_eligible=sel.sum(),n_labelled=lab.sum(),n=sub.sum(),fraction=sub.sum()/lab.sum() if lab.sum() else np.nan,n_MP10=sub.sum(),MP10_median=d.loc[sub,RAW].median() if sub.sum()>=5 else np.nan))
 # Minimal useful spot table plus scores. Existing full data are referenced, not recopied.
 keep=['section_id','barcode','spot_uid','x','y','x_um','y_um','array_row','array_col','tissue_component_id','lesion_location','distance_band_200um','exact_distance_status','distance_to_observed_outer_candidate_um','observed_outer_contact','query','epi_anchor','common_anchor','E0','score_eligible',RAW,'spot_MP10_model_universe','epi_MP10_model_universe','epi_score_valid','accepted_label','niche_version','nCount_RNA','Epithelial_fraction','CAF_fraction','TAM_fraction','Treg_fraction','mean__Epithelial','tumor_epi','MP10_P95','pEMT_high_front','d_front','front_band','d_MP10','d_MP10_reason']+[c for c in d if c.startswith(('reference_control','d_control','d_P95_matched','hops_'))]+P
 d[keep].to_parquet(O/'spots'/f'{sec}.parquet',index=False)
 for a,name in [(effects,'effects'),(curves,'curves'),(ringrows,'rings'),(counts,'bands'),(niches,'niche'),(pairs,'matching')]:
  a=pd.DataFrame(a);a.insert(0,'section_id',sec);save(a,O/'tables'/f'{sec}_{name}.tsv.gz')
 cov=dict(section_id=sec,n_spots=len(d),n_query=q.sum(),n_tumor_epi=tumor.sum(),n_P95=hi.sum(),P95_threshold=th,P95_fraction=hi.sum()/max(1,tumor.sum()),pEMT_P75=p75,n_pEMT_high_front=active.sum(),n_front=fr.sum(),n_query_front_distance=(q&np.isfinite(df)).sum(),n_query_MP10_distance=(q&np.isfinite(d.d_MP10)).sum(),n_R12_query=labelled.sum(),n_match_base=sum(p['matching']=='base' for p in pairs),n_match_front=sum(p['matching']=='front' for p in pairs),status='complete')
 save([cov],done);return sec+' complete'
if __name__=='__main__':
 secs=sys.argv[1:] or M.section_id.tolist()
 with ProcessPoolExecutor(max_workers=4) as ex:
  for s in ex.map(run,secs):print(s,flush=True)
