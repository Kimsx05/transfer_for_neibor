from core import *
import importlib.util
from scipy.sparse.csgraph import dijkstra
spec=importlib.util.spec_from_file_location('frozen_front_api',OLD/'scripts/analyze.py');api=importlib.util.module_from_spec(spec);spec.loader.exec_module(api)

def mean_boot(H,W,values,mask):
 mask=np.asarray(mask,bool);n=mask.sum();v=np.asarray(values,float);cnt=np.asarray(H@mask.astype(float)).ravel();s=np.asarray(H@np.where(mask[:,None],v,0));bn=W@cnt;b=np.divide(W@s,bn[:,None],out=np.full((B,v.shape[1]),np.nan),where=bn[:,None]>0);mean=v[mask].mean(0) if n else np.full(v.shape[1],np.nan);inf=np.asarray(H@np.where(mask[:,None],(v-mean)/max(n,1),0));return mean,b,inf,(cnt>0).sum()

def run_context(sec):
 if (O/'tables'/f'{sec}_context_done.tsv').exists():return sec+' cached'
 d=pd.read_parquet(O/'spots'/f'{sec}.parquet');q=d['query'].to_numpy(bool);tum=d.tumor_epi.to_numpy(bool);hi=d.MP10_P95.to_numpy(bool);front=tum&(d.d_front>=-200)&(d.d_front<=0);Y=d[NEW].to_numpy(float);effects=[];bands=[];niches=[];curves=[]
 H,W,ids=block_axis(d,20261104+M.section_id.tolist().index(sec));lab=q&d.accepted_label.notna().to_numpy();ny=np.column_stack([d.accepted_label.eq(k).fillna(False).to_numpy(float) for k in range(9)])
 for scope,base in [('all_Epi',q),('internal_Epi',tum),('front_only',np.asarray(front)),('near_inside',tum&d.front_band.eq('near_inside')),('deep',tum&d.front_band.eq('deep'))]:
  for pop,mask in [('all',base),('P95',base&hi),('nonP95',base&~hi)]:
   marked=mask&lab;den=marked.sum()
   for k in range(9):niches.append(dict(section_id=sec,scope=scope,population=pop,niche=k,n_eligible=np.sum(mask),n_labelled=den,n_unlabelled=np.sum(mask&~lab),fraction_unlabelled=np.sum(mask&~lab)/max(1,np.sum(mask)),n_niche=np.sum(marked&d.accepted_label.eq(k).fillna(False).to_numpy()),fraction=ny[marked,k].mean() if den else np.nan))
 # Truly front-only P95 versus nonP95, not the old across-band matched comparison.
 fm=np.asarray(front)&lab;h=mean_boot(H,W,ny,fm&hi);l=mean_boot(H,W,ny,fm&~hi);delta=h[0]-l[0];boots=h[1]-l[1];nrep=np.isfinite(boots).sum(0);inf=h[2]-l[2];K=((inf!=0).any(1)).sum();se=np.sqrt((inf*inf).sum(0)*K/max(1,K-1));nfront=[]
 for k in range(9):
  stable=min(h[3],l[3])>=8 and nrep[k]>=1600;ci=np.nanquantile(boots[:,k],[.025,.975]) if stable else [np.nan,np.nan];pv=2*stats.t.sf(abs(delta[k])/se[k],K-1) if stable and se[k]>0 else np.nan
  nfront.append(dict(section_id=sec,niche=k,n_P95=np.sum(fm&hi),n_nonP95=np.sum(fm&~hi),fraction_P95=h[0][k],fraction_nonP95=l[0][k],difference_percentage_points=delta[k]*100,ci_low_percentage_points=ci[0]*100,ci_high_percentage_points=ci[1]*100,p_spatial=pv,blocks_P95=h[3],blocks_nonP95=l[3],n_bootstrap_valid=nrep[k],status='ok' if stable else 'descriptive_insufficient_blocks_or_groups'))
 # Old overlap-neighborhood implementation; only new five programs are scored here.
 adj=api.graph(d,d.valid_tissue);qi=np.flatnonzero(q);hops=dijkstra(adj,directed=False,indices=qi,unweighted=True,limit=10);center=d.loc[qi,['section_id','barcode']].copy();ringdata={}
 for ring,(lo,up) in zip(api.RINGS,api.RB):
  ii,jj=np.where((hops>=lo)&(hops<=up)&d.valid_tissue.to_numpy()[None,:]);w=sparse.csr_matrix((np.ones(len(ii)),(ii,jj)),shape=(len(qi),len(d)));n=np.asarray(w@np.isfinite(Y));val=(w@np.nan_to_num(Y))/np.where(n>0,n,np.nan);full=np.full((len(d),len(NEW)),np.nan);full[qi]=val;ringdata[ring]=full
  cnt=np.asarray(w.sum(1)).ravel();ep=np.asarray(w@d.Epithelial_fraction.fillna(0).to_numpy()).ravel()/np.where(cnt>0,cnt,np.nan);extra=np.full(len(d),np.nan);extra[qi]=ep
  for j,p in enumerate(NEW):center[ring+'__'+p]=val[:,j]
  for adjust in [False,True]:effects+=api.fit(d,d[RAW].to_numpy(),full,NEW,q,'B1_'+ring,'B_center' if ring=='R0' else 'B_neighborhood',adjust=adjust,extra=None if ring=='R0' else extra,blocksize=500 if ring=='R0' else 2000)
 try:pairs=pd.read_csv(OLD/'tables'/f'{sec}_matching.tsv.gz',sep='\t')
 except pd.errors.EmptyDataError:pairs=pd.DataFrame()
 if 'matching' in pairs:
  for kind,g in pairs.groupby('matching'):
   ix=d.reset_index().set_index('barcode')['index'];h=ix.reindex(g.high_barcode).to_numpy(int);c=ix.reindex(g.control_barcode).to_numpy(int)
   for ring,values in ringdata.items():
    for j,p in enumerate(NEW):
     a=values[h,j]-values[c,j];n=np.isfinite(a).sum();sd=np.nanstd(values[tum,j]);effects.append(dict(analysis='B4_matched_'+kind+'_'+ring,family='B_neighborhood',adjustment='matched',program=p,n=n,effect=np.nanmean(a)/sd if n>=5 and sd>1e-12 else np.nan,raw_difference=np.nanmean(a) if n else np.nan,reason='paired_section_effect_hierarchical_inference'))
 center.to_parquet(O/'centers'/f'{sec}_new5.parquet',index=False);del hops
 for adjust in [False,True]:
  for pop,mask in [('tissue',d.valid_tissue.to_numpy()),('Epi',q)]:effects+=api.fit(d,d.d_MP10.to_numpy(),Y,NEW,mask&(d.d_MP10<=1000),'B3_'+pop+'_distance','B_neighborhood',adjust=adjust,quadratic=True,blocksize=2000)
  effects+=api.fit(d,d.d_front.to_numpy(),Y,NEW,tum,'A_front_distance','A_front',adjust=adjust,quadratic=True)
  for tag,mask in [('internal',tum&d.d_front.notna()),('near',tum&d.d_front.between(-500,0)),('deep',tum&((d.d_front< -500)|front))]:effects+=api.fit(d,np.asarray(front,float),Y,NEW,mask,'A_front_vs_'+tag,'A_front',adjust=adjust,binary=True)
  effects+=api.fit(d,d.d_front.to_numpy(),Y,NEW,hi,'C3_P95_state','C_MP10_front',adjust=adjust,quadratic=True)
 for pop,mask in [('tumor_Epi',tum),('P95_Epi',hi),('all_tissue',d.valid_tissue.to_numpy())]:
  for band,g in d[mask].groupby('front_band'):
   for p in NEW:
    v=g[p].dropna();bands.append(dict(section_id=sec,population=pop,band=band,program=p,n=len(v),mean=v.mean() if len(v)>=5 else np.nan,median=v.median() if len(v)>=5 else np.nan,q25=v.quantile(.25) if len(v)>=5 else np.nan,q75=v.quantile(.75) if len(v)>=5 else np.nan))
  for niche,g in d[mask&d.accepted_label.notna().to_numpy()].groupby('accepted_label'):
   for p in NEW:bands.append(dict(section_id=sec,population=pop,band='niche_'+str(niche),program=p,n=len(g),mean=g[p].mean() if len(g)>=5 else np.nan,median=g[p].median() if len(g)>=5 else np.nan))
 for axis in ['d_MP10','d_front']:
  x=d[axis].to_numpy();bb=np.ceil((x-1e-8)/100)*100 if axis=='d_MP10' else np.floor(x/100)*100
  for pop,mask in [('tissue',d.valid_tissue.to_numpy()),('Epi',q),('P95',hi)]:
   for v in np.unique(bb[mask&np.isfinite(bb)]):
    ok=mask&(bb==v)
    for p in NEW:curves.append(dict(section_id=sec,axis=axis,population=pop,bin=v,distance_um=np.median(x[ok]),program=p,n=ok.sum(),mean=d.loc[ok,p].mean() if ok.sum()>=5 else np.nan))
 effects=pd.DataFrame(effects);effects['section_id']=sec
 for a,n in [(effects,'new5_context'),(bands,'new5_bands'),(niches,'niche_counts'),(nfront,'front_only_niche'),(curves,'new5_curves')]:save(a,O/'tables'/f'{sec}_{n}.tsv.gz')
 save([dict(section_id=sec,status='complete')],O/'tables'/f'{sec}_context_done.tsv');return sec
if __name__=='__main__':
 with ProcessPoolExecutor(max_workers=3) as ex:
  for s in ex.map(run_context,M.section_id):print(s,flush=True)
