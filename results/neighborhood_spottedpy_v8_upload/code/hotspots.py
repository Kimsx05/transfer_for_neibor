from common import *
import warnings,ast,types,inspect,resource,importlib,functools
warnings.filterwarnings('ignore',category=FutureWarning)
from esda import G_Local
from libpysal.weights import WSP

_CRAND=importlib.import_module('esda.crand')
_ORIGINAL_VEC=_CRAND.vec_permutations
_CACHED_VEC=functools.lru_cache(maxsize=8)(_ORIGINAL_VEC)
def glocal(y,w,keep=False,star=False):
 _CRAND.vec_permutations=_CACHED_VEC
 try:
  with warnings.catch_warnings():
   warnings.simplefilter('ignore');return G_Local(y,WSP(w.copy()).to_W(silence_warnings=True),transform='R',permutations=9999,star=star,keep_simulations=keep,n_jobs=1,seed=123)
 finally:_CRAND.vec_permutations=_ORIGINAL_VEC


def run_variable(sec,var,k):
 key,name,domain=var;dest=OUT/'02_hotspots'/sec;dest.mkdir(exist_ok=True);file=dest/f'{key}_k{k}.parquet';state=dest/f'{key}_k{k}.json'
 if state.exists() and file.exists() and file.with_name(file.stem+'_components.parquet').exists():
  s=json.loads(state.read_text());assert s['fingerprint']==fingerprint() and file.exists();return s
 t,a,g=load_section(sec);geom=np.load(OUT/'01_geometry'/sec/f'{domain}_k{k}_axis.npz');ids=geom['node_index'];v=t[name].to_numpy(float) if domain=='SOURCE' else a['q05'][:,TYPES.index(name)];y=v[ids];status='CALCULABLE' if len(y)>=30 and len(np.unique(y))>=3 else 'INSUFFICIENT_POSITIONS' if len(y)<30 else 'INSUFFICIENT_DISTINCT_VALUES';start=time.time();z=p=gs=np.full(len(y),np.nan)
 if status=='CALCULABLE':
  w=sparse.load_npz(OUT/'01_geometry'/sec/f'{domain}_k{k}_W.npz');res=glocal(y,w);gs=res.Gs;z=res.Zs;p=res.p_sim.astype(float)
 p2=np.minimum(1,2*p);q=bh(p2);df=pd.DataFrame(dict(node_index=ids,spot_uid=t.spot_uid.to_numpy()[ids],value=y,G=gs,Z=z,p_sim_native=p,p_two_sided=p2,q_local=q,geometry_local=geom['geometry_local'],kth_distance=geom['kth_distance'],n_within_R3=geom['n_within_R3']));components=[];coverage=[]
 for track,cut in [('SOURCE_P05',p<.05),('LOCAL_BH05',q<.05)]:
  for direction,sgn in [('hot',z>0),('cold',z<0)]:
   raw=ids[cut&sgn];lab,formal,rec=induced_components(g,raw);df[f'{track}_{direction}_raw']=cut&sgn;df[f'{track}_{direction}']=formal[ids];df[f'{track}_{direction}_component']=lab[ids]
   for r in rec:
    nodes=np.flatnonzero(lab==r['component_id']);vals=v[nodes];r.update(track=track,direction=direction,value_mean=vals.mean(),value_median=np.median(vals),value_q25=np.quantile(vals,.25),value_q75=np.quantile(vals,.75),mean_percentile=np.mean(np.searchsorted(np.sort(y),vals,side='right')/len(y)),physical_component=int(t.physical_component.iloc[nodes[0]]));components.append(r)
   nh=int(formal.sum());coverage.append(dict(track=track,direction=direction,n_raw=len(raw),n_formal=nh,n_components=sum(r['formal'] for r in rec),fraction=nh/len(ids) if len(ids) else np.nan,single_component=sum(r['formal'] for r in rec)==1,large_fraction_warning=nh/len(ids)>.8 if len(ids) else False,mean_absolute_relative_composition=float(v[formal].mean()) if domain=='CELL' and nh else np.nan,LOW_ABSOLUTE_RELATIVE_COMPOSITION=bool(domain=='CELL' and nh and v[formal].mean()<.001)))
 df['pvalue_convention']=pd.Series('NATIVE_DIRECTED_SMALLER_TAIL; TWO_SIDED=MIN(1,2*NATIVE)',index=df.index,dtype='string');df.to_parquet(file,index=False,compression='zstd');pd.DataFrame(components).to_parquet(dest/f'{key}_k{k}_components.parquet',index=False);s=dict(section_id=sec,variable=key,variable_name=name,domain=domain,k=k,status=status,n_supported=len(ids),n_distinct=len(np.unique(y)),n_local_calculable=int(np.isfinite(p).sum()),seconds=time.time()-start,fingerprint=fingerprint(),coverage=coverage,rss_max_MB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024);js(str(state.relative_to(OUT)),s);return s

def source_parity(y,w,xy,uid):
 import anndata
 from typing import Union,Tuple,List
 tree=ast.parse((OUT/'00_audit/source/hotspot_helper.py').read_text());fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='calculate_hotspots_with_hotspots_numbered');captured=[]
 def capture(y,ww,**kw):
  r=G_Local(y,ww,n_jobs=1,keep_simulations=False,**kw);captured.append(r);return r
 fake=types.SimpleNamespace(cg=types.SimpleNamespace(KDTree=lambda x:x),weights=types.SimpleNamespace(KNN=lambda x,k:WSP(w.copy()).to_W(silence_warnings=True)))
 ns=dict(np=np,pd=pd,AnnData=anndata.AnnData,Union=Union,Tuple=Tuple,List=List,libpysal=fake,G_Local=capture,find_connected_components=lambda h,a:(h,0))
 exec(compile(ast.Module(body=[fn],type_ignores=[]),'locked_hotspot_function','exec'),ns);obs=pd.DataFrame(dict(scores=y,array_row=xy[:,0],array_col=xy[:,1],batch='pilot'),index=uid);ad=anndata.AnnData(obs=obs)
 with warnings.catch_warnings():
  warnings.simplefilter('ignore');out=ns[fn.name](ad,score_column='scores',neighbours_param=10,permutation=9999,seed_number=123)
 return captured[0],out

def pilot():
 if done('M02_PILOT_PASSED'):return
 rows=[];timings=[]
 # Native tail and simulation-retention parity on bounded nondegenerate synthetic sample.
 from libpysal.weights import KNN
 rng=np.random.default_rng(73);xy=rng.normal(size=(64,2));y=rng.uniform(.01,1,64);w=KNN.from_array(xy,k=10).sparse
 r=glocal(y,w,True);s=glocal(y,w,False);assert np.array_equal(r.p_sim,s.p_sim) and np.array_equal(r.Gs,s.Gs) and np.array_equal(r.Zs,s.Zs)
 counts=(r.sim>=r.Gs).sum(axis=0);native=(np.minimum(counts,9999-counts)+1)/10000
 assert np.allclose(native,r.p_sim,atol=3e-8,rtol=0);rows.append(dict(check='native smaller tail; kept versus discarded simulations',passed=True,max_error=float(abs(native-r.p_sim).max())))
 # Distance fixtures retain zero, NA, and component isolation.
 xx=np.array([[0.,0.],[1.,0.],[2.,0.],[0.,.1],[0.,.2]]);cc=np.array([0,0,0,1,1]);d=nearest(xx,cc,np.arange(5),np.array([0]));assert np.array_equal(d[:3],[0,1,2]) and np.isnan(d[3:]).all();assert np.isnan(nearest(xx,cc,np.arange(5),np.array([],int))).all();assert np.allclose(d,nearest(xx+100,cc,np.arange(5),np.array([0])),equal_nan=True);rows.append(dict(check='shared zero, empty target NA, disconnected target NA, translation',passed=True))
 # Locked source helper on small one-component case, including zero-valued target membership.
 from scipy.spatial import distance_matrix
 import anndata
 tree=ast.parse((OUT/'00_audit/source/hotspot_helper.py').read_text());fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='calculateDistancesHelper');ns=dict(np=np,pd=pd,distance_matrix=distance_matrix);exec(compile(ast.Module(body=[fn],type_ignores=[]),'locked_distance_helper','exec'),ns)
 obs=pd.DataFrame(dict(array_row=[0.,1.,2.],array_col=[0.,0.,0.],source=[0.,1.,2.],target=[0.,np.nan,np.nan]),index=['a','b','c']);ad=anndata.AnnData(obs=obs);ad.obsm['spatial']=obs[['array_row','array_col']].to_numpy();orig=ns[fn.name](ad,['source'],['target'],'pilot',False,False);assert np.allclose(orig.min_distance,[0,1,2]);rows.append(dict(check='locked source distance versus KDTree including zero-valued membership',passed=True))
 for pr in pd.read_csv(OUT/'00_audit/pilot_sections.tsv',sep='\t').itertuples():
  t,a,g=load_section(pr.section_id);axis=np.load(OUT/'01_geometry'/pr.section_id/'SOURCE_k10_axis.npz');ids=axis['node_index'];ww=sparse.load_npz(OUT/'01_geometry'/pr.section_id/'SOURCE_k10_W.npz')
  for var in VARIABLES[:2]:
   start=time.time();y=t[var[1]].to_numpy(float)[ids];direct=glocal(y,ww);original,out=source_parity(y,ww,a['xy_pitch'][ids],t.spot_uid.to_numpy()[ids]);assert np.array_equal(direct.Gs,original.Gs) and np.array_equal(direct.Zs,original.Zs) and np.array_equal(direct.p_sim,original.p_sim);raw=(direct.Zs>0)&(direct.p_sim<.05);assert np.array_equal(out.obs.scores_hot.notna().to_numpy(),raw);star=glocal(y,ww,star=True);hs=(star.Zs>0)&(star.p_sim<.05);rows.append(dict(check='locked source adapter parity and deterministic repeat',dataset_id=pr.dataset_id,section_id=pr.section_id,variable=var[0],passed=True,Gstar_max_stat_difference=np.max(abs(direct.Gs-star.Gs)),G_hot_fraction=raw.mean(),Gstar_hot_fraction=hs.mean(),G_vs_Gstar_hot_Jaccard=(raw&hs).sum()/(raw|hs).sum() if (raw|hs).any() else np.nan));timings.append(dict(section_id=pr.section_id,variable=var[0],n=len(ids),three_calls_seconds=time.time()-start));run_variable(pr.section_id,var,10);log('pilot '+pr.section_id+' '+var[0])
 tab('00_audit/pilot_checks.tsv',rows);tab('00_audit/pilot_runtime.tsv',timings);pernode=sum(x['three_calls_seconds'] for x in timings)/3/sum(x['n'] for x in timings);cov=pd.read_csv(OUT/'01_geometry/coordinate_and_graph_audit.tsv',sep='\t');nodes=3*(2*cov.n_E0.sum()+29*cov.n_valid_composition.sum());js('00_audit/runtime_estimate.json',dict(sequential_seconds=pernode*nodes,two_worker_hours=pernode*nodes/2/3600,peak_MB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,interpretation='pilot extrapolation; actual variable runtime logged',all_checks_pass=True));checkpoint('M02_PILOT_PASSED')

def section_worker(task):
 sec,stage=task
 work=[(v,10) for v in VARIABLES[:2]] if stage=='SOURCE_K10' else [(v,10) for v in VARIABLES[2:]] if stage=='CELL_K10' else [(v,k) for k in [6,18] for v in VARIABLES]
 for var,k in work:run_variable(sec,var,k)
 log(stage+' complete '+sec);return sec

def full():
 assert done('M02_PILOT_PASSED')
 from concurrent.futures import ProcessPoolExecutor
 with ProcessPoolExecutor(max_workers=2) as pool:
  for stage in ['SOURCE_K10','CELL_K10','SENSITIVITY_K6_K18']:
   for sec in pool.map(section_worker,[(s,stage) for s in COV.section_id]):pass
   checkpoint('ALL_'+stage+'_COMPLETE')
 records=[];cover=[]
 for sec in COV.section_id:
  for k in KS:
   for var in VARIABLES:
    s=json.loads((OUT/'02_hotspots'/sec/f'{var[0]}_k{k}.json').read_text());rr={x:v for x,v in s.items() if x!='coverage'};records.append(rr)
    for r in s['coverage']:cover.append(dict(**rr,**r))
 tab('02_hotspots/variable_execution_audit.tsv',records);tab('02_hotspots/hotspot_coverage.tsv',cover);assert len(records)==63*31*3;checkpoint('M03_ALL_HOTSPOTS_COMPLETE')
if __name__=='__main__':pilot();full()

def materialize_pilot_point_audit():
 path=OUT/'00_audit/pilot_G_vs_Gstar_points.parquet'
 if path.exists():return
 rows=[]
 for pr in pd.read_csv(OUT/'00_audit/pilot_sections.tsv',sep='\t').itertuples():
  w=sparse.load_npz(OUT/'01_geometry'/pr.section_id/'SOURCE_k10_W.npz')
  for var in VARIABLES[:2]:
   d=pd.read_parquet(OUT/'02_hotspots'/pr.section_id/f'{var[0]}_k10.parquet');r=glocal(d.value.to_numpy(),w,star=True);z=d[['spot_uid','value','G','Z','p_sim_native']].copy();z['dataset_id']=pr.dataset_id;z['section_id']=pr.section_id;z['variable']=var[0];z['Gstar']=r.Gs;z['Z_Gstar']=r.Zs;z['p_sim_native_Gstar']=r.p_sim;z['G_raw_hot']=(d.Z>0)&(d.p_sim_native<.05);z['Gstar_raw_hot']=(r.Zs>0)&(r.p_sim<.05);rows.append(z)
 pd.concat(rows,ignore_index=True).to_parquet(path,index=False,compression='zstd');log('Seven pilot G/Gstar point audit stored; not used in primary inference')
