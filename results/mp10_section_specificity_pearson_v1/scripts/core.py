import os,sys
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
R=Path(__file__).resolve().parents[3];sys.path.insert(0,str(R/'spatial_ecology_top20_v2/runtime'))
import numpy as np,pandas as pd
from scipy import sparse,stats
from scipy.spatial import cKDTree
from concurrent.futures import ProcessPoolExecutor
O=R/'results/mp10_section_specificity_pearson_v1';OLD=R/'results/mp10_p95_front_programs_v1';F=R/'mp10_paper_distance_full_v1';SUB=R/'mp10_paper_distance_subtypes_v2'
M=pd.read_csv(O/'tables/input_manifest.tsv',sep='\t');G=pd.read_csv(O/'tables/gene_sets.tsv',sep='\t');P=sorted(G.gs_name.unique());PAN=G.groupby('gs_name').panel.first().to_dict();CAND=pd.read_csv(O/'tables/candidate_mapping.tsv',sep='\t').candidate_type.tolist();RAW='EpiMP10_Top20_AUCell';B=1999
NEW=sorted(G[G.gs_name.str.startswith('BUDDEN_')].gs_name.unique())
def save(a,p):pd.DataFrame(a).to_csv(p,sep='\t',index=False,na_rep='NA')
def bh(p):
 p=np.asarray(p,float);o=np.flatnonzero(np.isfinite(p));o=o[np.argsort(p[o])];q=np.full(len(p),np.nan);q[o]=np.minimum.accumulate((p[o]*len(o)/np.arange(1,len(o)+1))[::-1])[::-1].clip(0,1);return q

def prepare(sec):
 d=pd.read_parquet(OLD/'spots'/f'{sec}.parquet');f=pd.read_parquet(F/'spots'/f'{sec}.parquet');means=f.filter(regex='^mean__');den=means.sum(1);valid=(means.ge(0)&np.isfinite(means)).all(1)&den.gt(0)&f.center_qc&f.lesion_location.ne('quality_excluded');d['valid_tissue']=valid.to_numpy()
 assert d.barcode.equals(f.barcode)
 sub=pd.read_parquet(SUB/'spots'/f'{sec}.parquet');assert sub.barcode.equals(d.barcode)
 targetrows=[]
 for t in CAND:
  v=means['mean__'+t]/den;th=v[valid].quantile(.95);mask=(valid&v.ge(th)&v.gt(0)).to_numpy();key=t+'_P95';oldcol=key+'_distance'
  if oldcol in sub and np.array_equal(mask,sub[key+'_target']):dist=sub[oldcol].to_numpy();origin='reused_subtypes_v2'
  else:
   dist=np.full(len(d),np.nan)
   for _,ix in d[valid].groupby('tissue_component_id').groups.items():
    ti=ix[mask[ix]]
    if len(ti):dist[ix]=cKDTree(d.loc[ti,['x_um','y_um']]).query(d.loc[ix,['x_um','y_um']],workers=1)[0]
   origin='new_same_component_exact_Euclidean'
  d['d_candidate__'+t]=dist;d['candidate_P95__'+t]=mask;d['mean__'+t]=means['mean__'+t].to_numpy();d['fraction__'+t]=v.where(valid).to_numpy()
  targetrows.append(dict(section_id=sec,candidate_type=t,n_valid_tissue=valid.sum(),relative_P95_threshold=th,n_reference=mask.sum(),actual_reference_fraction=mask.sum()/max(1,valid.sum()),absolute_abundance_median=means.loc[mask,'mean__'+t].median(),absolute_abundance_q25=means.loc[mask,'mean__'+t].quantile(.25),absolute_abundance_q75=means.loc[mask,'mean__'+t].quantile(.75),relative_abundance_median=v[mask].median(),n_reference_components=d.loc[mask,'tissue_component_id'].nunique(),n_query_distance=(d['query']&np.isfinite(dist)).sum(),n_query_no_reference=(d['query']&~np.isfinite(dist)).sum(),origin=origin,reason='ok' if mask.any() else 'no_positive_reference'))
 # Score joins never change query, P95 or region labels.
 ams=O/'scores'/f'{sec}_SCT_AddModuleScore.tsv.gz';uc=O/'scores'/f'{sec}_new5_UCell.tsv.gz'
 if ams.exists():d=d.merge(pd.read_csv(ams,sep='\t'),on=['section_id','barcode'],how='left',validate='1:1')
 else:
  for p in P:d['AMS_SCT__'+p]=np.nan
 if uc.exists():d=d.merge(pd.read_csv(uc,sep='\t'),on=['section_id','barcode'],how='left',validate='1:1')
 else:
  for p in NEW:d[p]=np.nan
 d.to_parquet(O/'spots'/f'{sec}.parquet',index=False);save(targetrows,O/'tables'/f'{sec}_targets.tsv.gz');return d

def block_axis(d,seed):
 lab=d.tissue_component_id.astype(str)+'|'+(d.x_um//500).astype(str)+'|'+(d.y_um//500).astype(str);ids,levels=pd.factorize(lab);K=len(levels);H=sparse.csr_matrix((np.ones(len(d)),(ids,np.arange(len(d)))),shape=(K,len(d)));rng=np.random.default_rng(seed);W=rng.multinomial(K,np.repeat(1/K,K),size=B).astype(float);return H,W,ids

def calc(d,x,Y,mask,H,W,boot=True):
 x=np.asarray(x,float);Y=np.asarray(Y,float);mask=np.asarray(mask,bool)
 if mask.ndim==1:mask=mask[:,None]
 ok=mask&np.isfinite(x)[:,None]&np.isfinite(Y);n=ok.sum(0);xx=np.where(ok,x[:,None],0);yy=np.where(ok,Y,0)
 # Center globally before forming sufficient statistics to reduce cancellation.
 mx=xx.sum(0)/np.maximum(n,1);my=yy.sum(0)/np.maximum(n,1);xx=np.where(ok,x[:,None]-mx,0);yy=np.where(ok,Y-my,0)
 moments=[np.asarray(H@v) for v in [ok.astype(float),xx,yy,xx*xx,yy*yy,xx*yy]];ss=[v.sum(0) for v in moments];r=from_sums(ss);r[n<2]=np.nan
 pp=np.full(len(n),np.nan);fine=(n>2)&np.isfinite(r);pp[fine]=2*stats.t.sf(np.abs(r[fine])*np.sqrt((n[fine]-2)/np.maximum(1-r[fine]**2,1e-16)),n[fine]-2);pp[(n==2)&np.isfinite(r)]=1
 nb=(moments[0]>0).sum(0)
 if boot:
  # A single shared multinomial weight matrix preserves repeated sampled blocks.
  v=W@np.concatenate(moments,axis=1);boots=from_sums(np.split(v,6,axis=1));valid=np.isfinite(boots).sum(0);lo=np.full(len(n),np.nan);hi=lo.copy();zse=lo.copy()
  stable=(nb>=8)&(valid>=1600)&np.isfinite(r)
  for j in np.flatnonzero(stable):lo[j],hi[j]=np.nanquantile(boots[:,j],[.025,.975]);zse[j]=np.nanstd(np.arctanh(np.clip(boots[:,j],-.999999,.999999)),ddof=1)
 else:boots=None;valid=np.zeros(len(n),int);lo=np.full(len(n),np.nan);hi=lo.copy();zse=lo.copy();stable=np.zeros(len(n),bool)
 sx=np.sqrt(np.maximum(ss[3]/np.maximum(n,1),0));sy=np.sqrt(np.maximum(ss[4]/np.maximum(n,1),0));zx=xx/np.where(sx>0,sx,1);zy=yy/np.where(sy>0,sy,1)
 infl=np.where(ok,(zx*zy-r[None,:]*(zx*zx+zy*zy)/2)/np.maximum(n,1),0);cl=np.asarray(H@np.nan_to_num(infl))/np.maximum(1-r*r,1e-10)
 return dict(n=n,r=r,p=pp,lo=lo,hi=hi,blocks=nb,valid_replicates=valid,zse=zse,boots=boots,cluster_z_influence=cl,moments_n=moments[0],stable=stable,ok=ok)
def from_sums(s):
 n,x,y,xx,yy,xy=s;vxx=np.maximum(xx-x*x/np.maximum(n,1),0);vyy=np.maximum(yy-y*y/np.maximum(n,1),0);den=np.sqrt(vxx*vyy);return np.divide(xy-x*y/np.maximum(n,1),den,out=np.full_like(den,np.nan),where=(den>1e-14)&(n>=2)).clip(-1,1)
def partial(d,x,Y,mask,front=False):
 C=np.column_stack([np.ones(len(d)),d.Epithelial_fraction,np.log1p(d.nCount_RNA)])
 if front:C=np.column_stack([C,d.d_front])
 mask=np.asarray(mask,bool)
 if mask.ndim==1:mask=np.repeat(mask[:,None],Y.shape[1],axis=1)
 out=np.full(Y.shape[1],np.nan);nn=np.zeros(Y.shape[1],int);patterns={}
 for j in range(Y.shape[1]):
  ok=mask[:,j]&np.isfinite(x)&np.isfinite(Y[:,j])&np.isfinite(C).all(1);patterns.setdefault(ok.tobytes(),[ok,[]])[1].append(j)
 for ok,js in patterns.values():
  nn[js]=ok.sum()
  if ok.sum()<5:continue
  c=C[ok];xr=x[ok]-c@np.linalg.lstsq(c,x[ok],rcond=None)[0];y=Y[ok][:,js];yr=y-c@np.linalg.lstsq(c,y,rcond=None)[0];den=np.linalg.norm(xr)*np.linalg.norm(yr,axis=0);out[js]=np.divide(xr@yr,den,out=np.full(len(js),np.nan),where=den>1e-12)
 return out,nn

def run(sec):
 if (O/'tables'/f'{sec}_status.tsv').exists():return sec+' cached'
 d=prepare(sec);meta=M.set_index('section_id').loc[sec];seed=20261104+M.section_id.tolist().index(sec);H,W,ids=block_axis(d,seed);q=d['query'].to_numpy(bool);tum=d.tumor_epi.to_numpy(bool);high=d.MP10_P95.to_numpy(bool);low=tum&~high;normal=meta.source_context=='adjacent_normal_source'
 A=d[['AMS_SCT__'+p for p in P]].to_numpy(float);U=d[P].to_numpy(float);pair=np.isfinite(A)&np.isfinite(U);A=np.where(pair,A,np.nan);U=np.where(pair,U,np.nan);Y=np.column_stack([A,U]);names=P+P;methods=['AddModuleScore_SCT_data']*len(P)+['UCell_RNA']*len(P)
 rows=[];deltas=[];close=[];mpres=[]
 try:matches=pd.read_csv(OLD/'tables'/f'{sec}_matching.tsv.gz',sep='\t')
 except pd.errors.EmptyDataError:matches=pd.DataFrame()
 for t in CAND:
  x=d['d_candidate__'+t].to_numpy(float);nref=d['candidate_P95__'+t].sum();querypop=q if normal else tum
  shared=np.repeat(querypop[:,None],Y.shape[1],axis=1)&np.isfinite(Y)&np.isfinite(x)[:,None]
  if normal:ranges={'full_range':shared,'fixed_0_1000':shared&(x[:,None]<=1000)}
  else:
   hmin=np.min(np.where(shared&high[:,None],x[:,None],np.inf),axis=0);hmax=np.max(np.where(shared&high[:,None],x[:,None],-np.inf),axis=0);lmin=np.min(np.where(shared&low[:,None],x[:,None],np.inf),axis=0);lmax=np.max(np.where(shared&low[:,None],x[:,None],-np.inf),axis=0);lower=np.maximum(hmin,lmin);upper=np.minimum(hmax,lmax);common=shared&(x[:,None]>=lower)&(x[:,None]<=upper)
   ranges={'full_range':shared,'common_support':common,'common_support_0_1000':common&(x[:,None]<=1000)}
  for rang,base in ranges.items():
   groups=[('normal_Epi',q)] if normal else [('P95',high),('nonP95',low),('all_internal_Epi',tum)]
   fits={}
   for pop,mask in groups:
    mask=base&mask[:,None];fit=calc(d,x,Y,mask,H,W);fits[pop]=fit
    part,pn=partial(d,x,Y,mask);fm=mask&np.isfinite(d.d_front.to_numpy(float))[:,None];p0,fn=partial(d,x,Y,fm);p1,_=partial(d,x,Y,fm,True)
    for j,p in enumerate(names):
     ok=fit['ok'][:,j];xx=x[ok];uni=np.unique(xx);gap=np.diff(uni).max() if len(uni)>1 else np.nan;nn=int(fit['n'][j]);rr=fit['r'][j]
     reason='no_candidate_reference' if not nref else 'no_common_support_or_pairs' if nn<2 else 'constant_distance_or_score' if not np.isfinite(rr) else 'very_sparse_n_lt10' if nn<10 else 'sparse_n_lt30' if nn<30 else 'ok'
     rows.append(dict(section_id=sec,candidate_type=t,program=p,program_family=PAN[p],score_method=methods[j],population=pop,range_type=rang,distance_min_um=np.min(xx) if nn else np.nan,distance_max_um=np.max(xx) if nn else np.nan,distance_sd_um=np.std(xx) if nn else np.nan,largest_distance_gap_um=gap,gap_over_200um=bool(gap>200),n_spots=nn,n_candidate_reference=nref,fraction_distance_zero=np.mean(xx==0) if nn else np.nan,pearson_r=rr,p_pearson=fit['p'][j],spatial_ci_low=fit['lo'][j],spatial_ci_high=fit['hi'][j],spatial_fisher_z_se=fit['zse'][j],n_spatial_blocks=int(fit['blocks'][j]),n_bootstrap_valid=int(fit['valid_replicates'][j]),spatial_ci_status='stable_by_prespecified_coverage' if fit['stable'][j] else 'insufficient_blocks_or_degenerate',partial_r_epi_depth=part[j],n_partial=pn[j],partial_r_front_subset_before=p0[j],partial_r_front_subset_after=p1[j],n_front_subset=fn[j],status='computed' if np.isfinite(rr) else 'not_estimable',reason=reason))
   if not normal:
    h,l=fits['P95'],fits['nonP95'];dr=h['r']-l['r'];db=h['boots']-l['boots'];nrep=np.isfinite(db).sum(0);ciok=h['stable']&l['stable']&(nrep>=1600);inf=h['cluster_z_influence']-l['cluster_z_influence'];nb=((h['moments_n']+l['moments_n'])>0).sum(0);se=np.sqrt((inf*inf).sum(0)*nb/np.maximum(nb-1,1));zd=np.arctanh(h['r'].clip(-.999999,.999999))-np.arctanh(l['r'].clip(-.999999,.999999));pv=2*stats.t.sf(np.abs(zd)/np.maximum(se,1e-15),np.maximum(nb-1,1));pv[(~ciok)|(np.abs(h['r'])>.999999)|(np.abs(l['r'])>.999999)]=np.nan
    if rang=='full_range':pv[:]=np.nan
    for j,p in enumerate(names):
     lo,hi=np.nanquantile(db[:,j],[.025,.975]) if ciok[j] else [np.nan,np.nan]
     deltas.append(dict(section_id=sec,candidate_type=t,program=p,program_family=PAN[p],score_method=methods[j],range_type=rang,common_min_um=lower[j] if rang!='full_range' and np.isfinite(lower[j]) else np.nan,common_max_um=min(upper[j],1000) if rang=='common_support_0_1000' else upper[j] if rang=='common_support' else np.nan,n_P95=h['n'][j],n_nonP95=l['n'][j],r_P95=h['r'][j],r_nonP95=l['r'][j],delta_r=dr[j],spatial_ci_low=lo,spatial_ci_high=hi,p_spatial_fisher_z=pv[j],n_blocks_P95=h['blocks'][j],n_blocks_nonP95=l['blocks'][j],n_blocks_union=nb[j],n_bootstrap_valid=nrep[j],spatial_method='joint_500um_block_bootstrap_1999; clustered_Fisher_z_difference_t',status='computed' if np.isfinite(dr[j]) else 'not_estimable',reason='descriptive_full_ranges' if rang=='full_range' else 'ok' if ciok[j] else 'insufficient_common_support_or_spatial_blocks'))
  # Comparison four: all frozen Epi, separately labelled raw Pearson and partial correlation.
  xx=d[RAW].to_numpy(float)[:,None];f=calc(d,x,xx,q,H,W);pa,pn=partial(d,x,xx,q);mpres.append(dict(section_id=sec,candidate_type=t,n=f['n'][0],pearson_r=f['r'][0],p_pearson=f['p'][0],spatial_ci_low=f['lo'][0],spatial_ci_high=f['hi'][0],n_blocks=f['blocks'][0],partial_r_epi_depth=pa[0],n_partial=pn[0],distance_min_um=np.nanmin(x[q]) if np.isfinite(x[q]).any() else np.nan,distance_max_um=np.nanmax(x[q]) if np.isfinite(x[q]).any() else np.nan,origin='same_definition_recomputed_raw_Pearson_and_new_spatial_CI'))
  # Comparison one: distance summaries, and spatial inference for mean difference.
  if not normal:
   hv=x[high&np.isfinite(x)];lv=x[low&np.isfinite(x)];rec=dict(section_id=sec,candidate_type=t,n_P95=len(hv),n_nonP95=len(lv),median_P95_um=np.median(hv) if len(hv) else np.nan,median_nonP95_um=np.median(lv) if len(lv) else np.nan,median_difference_um=np.median(hv)-np.median(lv) if len(hv) and len(lv) else np.nan,zero_fraction_P95=np.mean(hv==0) if len(hv) else np.nan,zero_fraction_nonP95=np.mean(lv==0) if len(lv) else np.nan)
   bootmeans=[];ifs=[];nbs=[]
   for mask in [high,low]:
    ok=mask&np.isfinite(x);nn=ok.sum();sums=np.asarray(H@np.where(ok,x,0)).ravel();cnt=np.asarray(H@ok.astype(float)).ravel();bc=W@cnt;bootmeans.append(np.divide(W@sums,bc,out=np.full(B,np.nan),where=bc>0));mean=x[ok].mean() if nn else np.nan;ifs.append(np.asarray(H@np.where(ok,(x-mean)/max(1,nn),0)).ravel());nbs.append((cnt>0).sum())
   if len(hv) and len(lv):
    diff=hv.mean()-lv.mean();db=bootmeans[0]-bootmeans[1];ci=np.nanquantile(db,[.025,.975]) if min(nbs)>=8 else [np.nan,np.nan];inf=ifs[0]-ifs[1];K=(inf!=0).sum();se=np.linalg.norm(inf)*np.sqrt(K/max(1,K-1));pv=2*stats.t.sf(abs(diff)/se,K-1) if min(nbs)>=8 and se>0 else np.nan;rec.update(mean_difference_um=diff,spatial_ci_low=ci[0],spatial_ci_high=ci[1],p_spatial=pv)
   if 'matching' in matches:
    pairs=matches[matches.matching.eq('base')];look=d.set_index('barcode')['d_candidate__'+t];a=look.reindex(pairs.high_barcode).to_numpy();b=look.reindex(pairs.control_barcode).to_numpy();ok=np.isfinite(a)&np.isfinite(b);rec.update(n_matched=ok.sum(),matched_mean_difference_um=np.mean(a[ok]-b[ok]) if ok.any() else np.nan,matched_median_pair_difference_um=np.median(a[ok]-b[ok]) if ok.any() else np.nan,matched_median_P95_um=np.median(a[ok]) if ok.any() else np.nan,matched_median_control_um=np.median(b[ok]) if ok.any() else np.nan,matched_inference='descriptive_paired_distances')
   close.append(rec)
 for data,name in [(rows,'pearson'),(deltas,'delta'),(close,'closeness'),(mpres,'MP10_candidate')]:
  a=pd.DataFrame(data)
  if len(a):a=a.merge(M[['section_id','dataset_id','analysis_unit_uid','patient_uid','source_context','eligible_for_patient_balanced_analysis']],on='section_id',validate='m:1')
  a.to_parquet(O/'tables'/f'{sec}_{name}.parquet',index=False)
 save([dict(section_id=sec,n_spots=len(d),n_Epi=q.sum(),n_internal=tum.sum(),n_P95=high.sum(),n_AMS_valid=d.filter(regex='^AMS_SCT__').notna().all(1).sum(),n_candidates=len(CAND),n_P95_estimable=sum(np.isfinite(r['pearson_r']) for r in rows if r['population']=='P95' and r['range_type']=='full_range' and r['score_method']=='AddModuleScore_SCT_data'),status='complete')],O/'tables'/f'{sec}_status.tsv');return sec+' complete'
if __name__=='__main__':
 secs=sys.argv[1:] or M.section_id.tolist()
 with ProcessPoolExecutor(max_workers=4) as ex:
  for sec in ex.map(run,secs):print(sec,flush=True)
