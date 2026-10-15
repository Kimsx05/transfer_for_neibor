from common import *
from relations import family,KEY

def validate(final=False):
 records=[]
 def check(name,condition,detail=''):
  records.append(dict(check=name,passed=bool(condition),detail=detail));assert condition,(name,detail)
 manifest=pd.read_csv(OUT/'00_audit/input_manifest.tsv',sep='\t');import pyarrow.parquet as pq
 uidcounts=[]
 for rr in manifest.itertuples():
  if str(rr.path).endswith('.parquet'):
   pf=pq.ParquetFile(rr.path)
   if 'spot_uid' in pf.schema.names:uidcounts.append(dict(path=rr.path,n_UID_rows=pf.metadata.num_rows,purpose=rr.purpose))
 tab('00_audit/input_uid_counts.tsv',uidcounts)
 check('All frozen old input SHA256 unchanged',all(Path(r.path).exists() and sha(r.path)==r.sha256 for r in manifest.itertuples()),f'{len(manifest)} input files')
 check('Pilot passed',done('M02_PILOT_PASSED'));check('Gstar pilot point audit stored', (OUT/'00_audit/pilot_G_vs_Gstar_points.parquet').exists())
 check('All variables all scales completed',done('M03_ALL_HOTSPOTS_COMPLETE'))
 audit=pd.read_csv(OUT/'02_hotspots/variable_execution_audit.tsv',sep='\t');check('63 sections x31 variables x3 k',len(audit)==5859 and not audit.duplicated(['section_id','variable','k']).any())
 checks=0;pointrows=0;components=[];pointmanifest=[];cov=[]
 for cr in COV.itertuples():
  sec=cr.section_id;t,a,g=load_section(sec);support=t[['spot_uid','section_id','physical_component','E0','composition_valid']].copy();
  for domain,mask in [('SOURCE',t.E0.to_numpy()),('CELL',t.composition_valid.to_numpy())]:
   for kk in KS:
    ax=np.load(OUT/'01_geometry'/sec/f'{domain}_k{kk}_axis.npz');st=np.full(len(t),'NOT_IN_FROZEN_E0' if domain=='SOURCE' else 'COMPOSITION_MISSING',dtype=object);st[mask]='LOCAL_SUPPORT_INSUFFICIENT';st[ax['node_index']]='SUPPORTED';support[f'{domain}_k{kk}_status']=st
  support.to_parquet(OUT/'01_geometry'/sec/'local_support_audit.parquet',index=False,compression='zstd');check('UID axis and unique '+sec,t.spot_uid.is_unique and np.array_equal(t.spot_uid,a['spot_uid']))
  shuffled=t.sample(frac=1,random_state=71).set_index('spot_uid').loc[t.spot_uid];check('Explicit UID restore after shuffle '+sec,np.array_equal(shuffled[LAYERS],t[LAYERS],equal_nan=True))
  try:shuffled.loc[['__BAD_UID__']];bad=False
  except KeyError:bad=True
  check('Unknown UID rejected '+sec,bad)
  check('Pixel-to-pitch identity and scale invariance '+sec,np.allclose(a['xy_pitch'],t[['imagecol','imagerow']].to_numpy()/t.pitch_px.iloc[0]) and np.allclose(a['xy_pitch'],(t[['imagecol','imagerow']].to_numpy()*.3)/(t.pitch_px.iloc[0]*.3)));check('No score imputation '+sec,np.isfinite(t.loc[t.E0,LAYERS]).all().all())
  for k in KS:
   sourceuids=None
   for key,name,domain in VARIABLES:
    path=OUT/'02_hotspots'/sec/f'{key}_k{k}.parquet';d=pd.read_parquet(path);pointrows+=len(d);geom=np.load(OUT/'01_geometry'/sec/f'{domain}_k{k}_axis.npz');checkok=np.array_equal(d.spot_uid,geom['spot_uid']);assert checkok
    if domain=='SOURCE':
     if sourceuids is None:sourceuids=d.spot_uid.to_numpy()
     else:assert np.array_equal(sourceuids,d.spot_uid)
    finite=np.isfinite(d.p_sim_native);assert (d.loc[finite,'p_sim_native']>0).all() and (d.loc[finite,'p_sim_native']<=1).all();assert np.allclose(d.p_two_sided,np.minimum(1,2*d.p_sim_native),equal_nan=True);assert np.allclose(d.q_local,bh(d.p_two_sided),equal_nan=True)
    for track,cut in [('SOURCE_P05',d.p_sim_native<.05),('LOCAL_BH05',d.q_local<.05)]:
     for direction,sign in [('hot',d.Z>0),('cold',d.Z<0)]:
      assert np.array_equal(d[f'{track}_{direction}_raw'],cut&sign);raw=d.node_index.to_numpy()[cut&sign];_,formal,rec=induced_components(g,raw);assert np.array_equal(d[f'{track}_{direction}'],formal[d.node_index])
    ss=json.loads(path.with_suffix('.json').read_text());assert ss['fingerprint']==fingerprint();bg=d.value.to_numpy();ssbase={k:v for k,v in ss.items() if k!='coverage'}
    for rr in ss['coverage']:cov.append(dict(dataset_id=cr.dataset_id,**ssbase,**rr,full_tissue_background_mean=float(np.nanmean(a['q05'][:,TYPES.index(name)])) if domain=='CELL' else float(t.loc[t.E0,name].mean()),background_domain='component-supported variable positions',background_mean=float(np.mean(bg)) if len(bg) else np.nan,background_median=float(np.median(bg)) if len(bg) else np.nan,background_q25=float(np.quantile(bg,.25)) if len(bg) else np.nan,background_q75=float(np.quantile(bg,.75)) if len(bg) else np.nan,p_min_resolution=.0001,pvalue_convention='native directed smaller tail; p_two=min(1,2*p_native)',n_tissue=len(t),n_E0=int(t.E0.sum())))
    cp=pd.read_parquet(path.with_name(path.stem+'_components.parquet'))
    if len(cp):cp=cp.assign(dataset_id=cr.dataset_id,section_id=sec,variable=key,variable_name=name,domain=domain,k=k);components.append(cp)
    pointmanifest.append(dict(dataset_id=cr.dataset_id,section_id=sec,variable=key,variable_name=name,k=k,n_rows=len(d),path=str(path.relative_to(OUT)),sha256=sha(path),pvalue_convention='native directed smaller tail; p_two=min(1,2*p_native)',q_family='all finite local positions of section/variable/k'))
    if 'pvalue_convention' not in d:
     d['pvalue_convention']=pd.Series('NATIVE_DIRECTED_SMALLER_TAIL; TWO_SIDED=MIN(1,2*NATIVE)',index=d.index,dtype='string');tmp=path.with_suffix('.parquet.tmp');d.to_parquet(tmp,index=False,compression='zstd');tmp.replace(path);pointmanifest[-1]['sha256']=sha(path)
    checks+=1
 check('All 5859 local outputs/UIDs/P/BH/induced components validated',checks==5859,f'{pointrows} point statistic rows')
 tab('02_hotspots/point_hotspot_statistics_manifest.tsv',pointmanifest);tab('02_hotspots/hotspot_coverage_all.tsv',cov);pd.concat(components,ignore_index=True).to_csv(OUT/'02_hotspots/hotspot_components.tsv.gz',sep='\t',index=False,compression='gzip')
 # Full locked calculateDistances wrapper, with explicit no max-distance fill.
 import ast,anndata
 from scipy.spatial import distance_matrix
 tree=ast.parse((OUT/'00_audit/source/hotspot_helper.py').read_text());fns=[x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name in ['calculateDistances','calculateDistancesHelper']];ns=dict(np=np,pd=pd,distance_matrix=distance_matrix);exec(compile(ast.Module(body=fns,type_ignores=[]),'locked_distance_API','exec'),ns)
 ad=anndata.AnnData(obs=pd.DataFrame(dict(array_row=[0.,1.,2.],array_col=[0.,0.,0.],source=[0.,1.,2.],target=[0.,np.nan,np.nan],batch=['pilot']*3),index=['a','b','c']));original=ns['calculateDistances'](ad,['source'],['target'],False,False,False);check('Full locked distance API retains zero membership and agrees',np.allclose(original.min_distance,[0,1,2]));ad.obs['target']=np.nan;empty=ns['calculateDistances'](ad,['source'],['target'],False,False,False);check('Full locked API empty-target max-fill explicitly disabled',empty.empty)
 check('Exact/no-tie and approximate/tie signed-rank branches',signed([1,2,3,4,5,6])['actual_method'].startswith('exact') and signed([1,1,1,1,1,0])['actual_method'].startswith('approx'));check('Numerical nearzero audit',signed([1e-14,-1e-14,0,0,0])['n_rounded_nearzero']==2)
 # Deliberately disconnected induced graph; no rebuilt KNN bridge.
 g=sparse.csr_matrix(np.array([[0,1,0,0],[1,0,1,0],[0,1,0,1],[0,0,1,0]]));_,fo,rec=induced_components(g,np.array([0,3]),2);check('Induced hot graph does not bridge omitted tissue nodes',not fo.any() and len(rec)==2)
 x=np.random.default_rng(9).normal(size=(40,2));c=np.arange(40)%3;target=np.arange(0,40,4);dd=nearest(x,c,np.arange(40),target);brute=np.array([np.linalg.norm(x[target[c[target]==c[i]]]-x[i],axis=1).min() for i in range(40)]);check('KDTree vs brute within-component distances',np.allclose(dd,brute));check('Translation preserves metric',np.allclose(dd,nearest(x+123,c,np.arange(40),target)));check('Empty target returns only NA',np.isnan(nearest(x,c,np.arange(40),np.array([],int))).all());check('All-zero signed rank returns P1 only sufficient N',signed([0]*5)['p']==1 and np.isnan(signed([0]*4)['p']));check('Signed rank insufficient nonzero remains NA',np.isnan(signed([0,0,0,1,2])['p']))
 s=pd.read_parquet(OUT/'07_synthesis/all_section_effects.parquet');u=pd.read_parquet(OUT/'07_synthesis/all_unit_effects.parquet');r=pd.read_parquet(OUT/'07_synthesis/all_scope_tests.parquet');check('Complete section target metric scale domain rows',len(s)==63*2*2*38*3*4);check('Shrinkage identity audited',s.shrinkage_error.abs().max()<1e-8);check('Outside domain excludes target on both sides',s.loc[s.metric.eq('D_out'),'n_target_in_E'].eq(0).all());check('Unknown source effects remain NA',s.loc[s.status.eq('NO_LOCAL_SOURCE_SUPPORT'),'effect'].isna().all());check('Target absent remains NA',s.loc[s.status.isin(['TARGET_NICHE_ABSENT','NO_TARGET_HOTSPOT','TARGET_UNREACHABLE']),'effect'].isna().all());check('All policy scopes',set(r.scope)==set(SCOPES));check('Unit UID not duplicated within tests',not u.duplicated(KEY+['analysis_unit_uid']).any())
 for scope in SCOPES:
  ss=s if scope=='ALL_DATASETS' else s[s.dataset_id.eq(scope)];expected=ss[ss.eligible].groupby(KEY[1:]+['analysis_unit_uid']).effect.mean();actual=u[u.scope.eq(scope)].set_index(KEY[1:]+['analysis_unit_uid']).effect.reindex(expected.index);assert np.allclose(expected,actual,equal_nan=True)
 check('Equal section aggregation independently recomputed',True)
 main=r[r.k.eq(10)&r.analysis_domain.eq('MAIN')]
 for key,ii in main.groupby(['scope','track','subset','family']).groups.items():assert len(ii)==FAMILY_N[key[-1]] and np.allclose(r.loc[ii,'q_scope'],bh(r.loc[ii,'p'],FAMILY_N[key[-1]]),equal_nan=True)
 check('Fixed planned BH families per scope/track/subset',True);check('Sensitivity k never provides competing P',r.loc[r.k.ne(10),'p'].isna().all())
 matched=pd.read_csv(OUT/'05_sensitivity/two_score_layers_matched_units.tsv',sep='\t');mr=pd.read_csv(OUT/'05_sensitivity/two_score_layers_matched_unit_rows.tsv',sep='\t');check('Both layers matched table equal paired rows',not mr[['effect_spot','effect_epi']].isna().any().any())
 check('All nine niche and29 cell targets represented',s[s.target_kind.eq('NICHE')].target_id.nunique()==9 and s[s.target_kind.eq('CELL')].target_id.nunique()==29)
 if final:
  check('Figures complete checkpoint',done('M07_FIGURES_COMPLETE'))
  for sec in COV.section_id:
   ds=COV.set_index('section_id').loc[sec,'dataset_id'];base=OUT/'06_figures'/ds/'sections'/sec
   for f in ['MP10_overview','cell_targets_1','cell_targets_2','cell_targets_3']:
    for ext in ['.png','.pdf']:assert (base/(f+ext)).exists()
  check('All 63 spatial source and29 target pages present',True)
  from PIL import Image
  fp=next((OUT/'06_figures').rglob('*.png'));im=Image.open(fp);check('Transparent 300dpi PNG',im.mode=='RGBA' and min(im.getextrema()[-1])==0 and abs(im.info['dpi'][0]-300)<1)
 tab('reports/validation_results.tsv',records);js('reports/validation_summary.json',dict(passed=len(records),failed=0,point_rows=pointrows,variables=checks,final=final));return records
if __name__=='__main__':validate(final=True)
