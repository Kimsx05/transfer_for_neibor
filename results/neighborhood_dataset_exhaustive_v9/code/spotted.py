from core import *
from stats import *
import warnings,importlib,functools
from esda import G_Local
from libpysal.weights import WSP
_CRAND=importlib.import_module('esda.crand');_ORIGINAL_VEC=_CRAND.vec_permutations;_CACHED_VEC=functools.lru_cache(maxsize=2)(_ORIGINAL_VEC)
extract(V8/'code/hotspots.py',['glocal'],globals())

def section(sec,labels,master,scorebase):
 p=V8/'01_geometry'/sec;t=pd.read_parquet(p/'nodes.parquet');a=np.load(p/'values.npz');g=sparse.load_npz(p/'R1.npz');assert np.array_equal(a['spot_uid'],t.spot_uid)
 take=pd.Index(master.spot_uid).get_indexer(t.spot_uid);valid=take>=0;acc=np.full(len(t),-1,np.int16);acc[valid]=labels[take[valid]];e0=np.zeros(len(t),bool);e0[valid]=scorebase[take[valid]]&(acc[valid]>=0)
 # Existing geometry contains the frozen scores for all common anchors, not only old E0.
 assert np.isfinite(t.loc[e0,LAYERS]).all().all()
 t['accepted_label']=acc;t['E0']=e0;return t,a,g

def geometry(t,a,g,k):
 xy=a['xy_pitch'];comp=t.physical_component.to_numpy();mask=t.E0.to_numpy();counts=np.bincount(comp[mask],minlength=int(comp.max())+1);ids=np.flatnonzero(mask&(counts[comp]>=k+1));index=np.full(len(t),-1,int);index[ids]=np.arange(len(ids));nn=np.empty((len(ids),k),int);dist=np.full(len(ids),np.nan);within=np.zeros(len(ids),int)
 g2=(g.astype(np.int32)@g.astype(np.int32)).tocsr();g2.data[:]=1;g2=g2-g2.multiply(g);g2.setdiag(0);g2.eliminate_zeros();g3=g2@g;g3.data[:]=1;g3=g3-g3.multiply(g)-g3.multiply(g2);g3.setdiag(0);g3.eliminate_zeros();reach=g+g2+g3
 # Match v8 arithmetic: physical pixel subtraction THEN pitch normalization in tie sorting.
 rawxy=t[['imagecol','imagerow']].to_numpy(float);pitch=float(t.pitch_px.iloc[0])
 for cc in np.unique(comp[ids]):
  local=ids[comp[ids]==cc];tree=cKDTree(xy[local]);dd,_=tree.query(xy[local],k=k+1,workers=1)
  for j,node in enumerate(local):
   cand=local[tree.query_ball_point(xy[node],dd[j,-1]+1e-10)];cand=cand[cand!=node];ds=np.linalg.norm((rawxy[cand]-rawxy[node])/pitch,axis=1);order=np.lexsort((cand,np.round(ds,12)));choose=cand[order[:k]];row=index[node];nn[row]=index[choose];dist[row]=np.linalg.norm((rawxy[choose[-1]]-rawxy[node])/pitch);within[row]=len(set(reach[node].indices).intersection(choose))
 w=sparse.csr_matrix((np.ones(len(ids)*k),(np.repeat(np.arange(len(ids)),k),nn.ravel())),shape=(len(ids),len(ids)))
 return ids,w,(within==k)&(dist<=3+1e-12)

def source(sec,t,a,g,k,folder):
 dest=folder/'domains';dest.mkdir(exist_ok=True);file=dest/f'{sec}_k{k}.npz';meta=dest/f'{sec}_k{k}.json'
 key=digest(dict(input_hash=sha(OUT/'audit/input_manifest.tsv'),uid=arrhash(t.spot_uid.to_numpy(dtype=str)),E0=arrhash(t.E0.to_numpy()),scores=arrhash(t[LAYERS].to_numpy(float)),k=k,config=config_hash(),policy='current_accepted&base_score;full_component_W;native_P05;2native_BH_all;R1min5',versions=json.loads((V8/'00_audit/source_versions.json').read_text())))
 if meta.exists() and file.exists():
  m=json.loads(meta.read_text());assert m['key']==key and m['sha256']==sha(file),'DOMAIN_CACHE_INVALID';return dict(np.load(file)),m
 guard(folder.name+':source:'+sec,2*1024**2);start=time.time();ids,w,local=geometry(t,a,g,k);result=dict(node_index=ids.astype(np.int32),geometry_local=local,E0=np.packbits(t.E0.to_numpy()));records=[]
 oldids=np.load(V8/'01_geometry'/sec/f'SOURCE_k{k}_axis.npz')['node_index'];oldE0=pd.read_parquet(V8/'01_geometry'/sec/'nodes.parquet',columns=['E0']).E0.to_numpy();same_domain=np.array_equal(t.E0,oldE0) and np.array_equal(ids,oldids)
 for li,layer in enumerate(LAYERS):
  y=t[layer].to_numpy(float)[ids];eligible=len(y)>=30 and len(np.unique(y))>=3
  if same_domain:
   old=pd.read_parquet(V8/'02_hotspots'/sec/f'{["MP10_SPOT","MP10_EPI"][li]}_k{k}.parquet');assert np.array_equal(old.node_index,ids) and np.array_equal(old.value,y,equal_nan=True)
   oldw=sparse.load_npz(V8/'01_geometry'/sec/f'SOURCE_k{k}_W.npz');assert (w!=oldw).nnz==0
   z=old.Z.to_numpy();p=old.p_sim_native.to_numpy();reason='EXACT_OLD_DOMAIN_VALUES_WEIGHTS_REUSE'
  elif eligible:
   rr=glocal(y,w);z=rr.Zs;p=rr.p_sim.astype(float);reason='NEW_ACCEPTED_DOMAIN_RECOMPUTED'
  else:z=p=np.full(len(ids),np.nan);reason='INSUFFICIENT_POSITIONS_OR_DISTINCT'
  result[f'E{li}']=np.isfinite(z)&np.isfinite(p);q=bh(np.minimum(1,2*p));coverage=[]
  for ti,(track,cut) in enumerate(zip(TRACKS,[p<.05,q<.05])):
   lab,formal,rec=induced_components(g,ids[cut&(z>0)]);result[f'H{li}_{ti}']=formal[ids];result[f'C{li}_{ti}']=lab[ids].astype(np.int32)
   coverage.append(dict(track=track,n_hot=int(formal.sum()),n_raw=int((cut&(z>0)).sum()),n_components=sum(r['formal'] for r in rec)))
  records.append(dict(layer=layer,reason=reason,n_supported=len(ids),n_calculable=int(result[f'E{li}'].sum()),coverage=coverage))
 npz(file,**result);info=dict(key=key,domain_hash=arrhash(t.E0.to_numpy()),sha256=sha(file),seconds=time.time()-start,records=records);js(meta,info);return result,info

def effects(sec,t,a,g,src,K,k=10,domains=['MAIN']):
 meta=json.loads((V8/'01_geometry'/sec/'metadata.json').read_text());rows=[];n=len(t);xy=a['xy_pitch'];physical=t.physical_component.to_numpy();ids=src['node_index'];cell={}
 for j in NON:
  old=pd.read_parquet(V8/'02_hotspots'/sec/f'CELL_{j:02d}_k{k}.parquet',columns=['node_index']+[f'{tr}_hot' for tr in TRACKS]);cell[j]=old
 for ti,track in enumerate(TRACKS):
  targets=[('NICHE',f'niche{i}',np.flatnonzero(t.accepted_label.to_numpy()==i)) for i in range(K)]
  targets += [('CELL',TYPES[j],cell[j].node_index.to_numpy()[cell[j][track+'_hot'].to_numpy()]) for j in NON]
  for li,layer in enumerate(LAYERS):
   E=ids[src[f'E{li}']];H=np.zeros(n,bool);H[ids[src[f'H{li}_{ti}']]]=True;comps=np.full(n,-1,int);comps[ids]=src[f'C{li}_{ti}'];local=np.zeros(n,bool);local[ids]=src['geometry_local']
   for kind,target,tid in targets:
    tm=np.zeros(n,bool);tm[tid]=True;d=nearest(xy,physical,E,tid);reason='NO_LOCAL_SOURCE_SUPPORT' if not len(E) else ('TARGET_NICHE_ABSENT' if kind=='NICHE' else 'NO_TARGET_HOTSPOT') if not len(tid) else 'NO_SOURCE_HOTSPOT' if not H.any() else 'READY'
    for domain in domains:
     dom=np.ones(len(E),bool) if domain=='MAIN' else local[E];ee=E[dom];dd=d[dom];hh=H.copy()
     if domain!='MAIN':hh[~local]=False
     reach=np.isfinite(dd);base=dict(dataset_id=meta['dataset_id'],section_id=sec,analysis_unit_uid=meta['analysis_unit_uid'],known_patient=meta['known_patient'],layer=layer,k=k,track=track,target_kind=kind,target_id=target,analysis_domain=domain,n_source_supported=len(ee),n_source_hot=int(hh[ee].sum()),source_formal_components=len(np.unique(comps[hh])),n_E_unreachable=int((~reach).sum()),source_domain_hash=arrhash(t.E0.to_numpy()))
     for metric,mask in [('O',np.ones(len(ee),bool)),('D',reach),('D_out',reach&~tm[ee])]:
      at=ee[mask];values=tm[at].astype(float) if metric=='O' else dd[mask];why='TARGET_UNREACHABLE' if metric!='O' and reason=='READY' and not reach.any() else reason;rr=metric_record(metric,at,hh,values,tm,t,comps,len(tid),why);rows.append(dict(**base,**rr))
 return pd.DataFrame(rows)

def run(cid):
 folder=OUT/'combinations'/cid
 if valid(folder,'SPOTTED_MAIN_DONE'):return
 assert valid(folder,'NICHE_FROZEN');lab=np.load(folder/'labels.npz');master=pd.read_parquet(OUT/'shared/master_anchors.parquet');labels=np.full(len(master),-1,np.int16);labels[lab['anchor_index']]=lab['accepted'];score=np.load(OUT/'shared/base_domains.npz')['base_score'];K=json.loads((folder/'selection.json').read_text())['K'];members=pd.read_csv(OUT/'combinations.tsv',sep='\t',dtype={'combo_id':str}).set_index('combo_id').loc[cid,'datasets'].split('|');cov=pd.read_csv(V2/'00_preflight/coverage63.tsv',sep='\t');parts=[];domainaudit=[]
 for sec in cov[cov.dataset_id.isin(members)].section_id:
  dest=folder/'spotted_sections'/f'{sec}.parquet';state=dest.with_suffix('.json')
  if dest.exists() and state.exists():
   st=json.loads(state.read_text());assert st['label_hash']==sha(folder/'labels.npz') and st['sha256']==sha(dest);parts.append(pd.read_parquet(dest));continue
  t,a,g=section(sec,labels,master,score);src,info=source(sec,t,a,g,10,folder);df=effects(sec,t,a,g,src,K);table(dest,df);js(state,dict(label_hash=sha(folder/'labels.npz'),sha256=sha(dest),domain_hash=info['domain_hash']));parts.append(df);domainaudit.append(dict(section_id=sec,n_E0=int(t.E0.sum()),n_old_E0=int(pd.read_parquet(V8/'01_geometry'/sec/'nodes.parquet',columns=['E0']).E0.sum()),domain_hash=info['domain_hash'],seconds=info['seconds']))
  log(cid+' spotted '+sec)
 domainaudit=[]
 for sec in cov[cov.dataset_id.isin(members)].section_id:
  meta=json.loads((folder/'domains'/f'{sec}_k10.json').read_text());t,_,_=section(sec,labels,master,score);assert arrhash(t.E0.to_numpy())==meta['domain_hash'];old=pd.read_parquet(V8/'01_geometry'/sec/'nodes.parquet',columns=['E0']);domainaudit.append(dict(section_id=sec,n_E0=int(t.E0.sum()),n_old_E0=int(old.E0.sum()),domain_hash=meta['domain_hash'],seconds=meta['seconds'],source_decisions='|'.join(x['reason'] for x in meta['records'])))
 s=pd.concat(parts,ignore_index=True);u,r=summarize_spotted(s,K,'REFIT');table(folder/'spotted_units.parquet',u);table(folder/'spotted_tests.parquet',r);table(folder/'domain_audit.tsv',domainaudit)
 mark(folder,'SPOTTED_MAIN_DONE',[folder/'spotted_units.parquet',folder/'spotted_tests.parquet']+list((folder/'spotted_sections').glob('*'))+list((folder/'domains').glob('*')))
if __name__=='__main__':run(sys.argv[1])
