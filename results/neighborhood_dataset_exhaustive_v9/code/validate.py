from core import *
from stats import *
from spotted import section,geometry,source,effects,glocal,_CACHED_VEC,_CRAND,_ORIGINAL_VEC

def validate():
 if (OUT/'validation_summary.json').exists():return
 checks=[]
 def check(name,condition,details=None):
  assert condition,name;checks.append(dict(check=name,passed=True,details=details));log('VALID '+name)
 c=pd.read_csv(OUT/'combinations.tsv',sep='\t',dtype={'combo_id':str});check('99_57_42',len(c)==99 and c.contains_inhouse.sum()==57 and (~c.contains_inhouse).sum()==42 and c.combo_id.is_unique)
 from scipy.stats import false_discovery_control
 p=np.array([.001,.02,.03,np.nan,.9,0.,1.]);ok=np.isfinite(p)
 for meth,by in [('fdr_bh',False),('fdr_by',True)]:check('planned_missing_'+meth,np.allclose(correction(p,by)[ok],false_discovery_control(np.where(ok,p,1),method='by' if by else 'bh')[ok],atol=1e-15,rtol=1e-14))
 check('scientific_NA_stays_NA',np.isnan(correction(p)[3]))
 check('all_zero_and_insufficient',signed([0]*5)['p']==1 and np.isnan(signed([1]*4)['p']))
 xx=np.array([[0.,0.],[1,0],[2,0],[0,.1]]);cc=np.array([0,0,0,1]);d=nearest(xx,cc,np.arange(4),np.array([0]));check('zero_distance_disconnected_empty_target',np.array_equal(d[:3],[0,1,2]) and np.isnan(d[3]) and np.isnan(nearest(xx,cc,np.arange(4),np.array([],int))).all())
 master=pd.read_parquet(OUT/'shared/master_anchors.parquet');b=np.load(OUT/'shared/base_domains.npz');m=pd.read_parquet(V7/'01_ring123/anchor_metadata.parquet');sec='Inhousedata_T4398343';t,a,g=section(sec,b['old_accepted'],master,b['base_score']);old=pd.read_parquet(V8/'01_geometry'/sec/'nodes.parquet');check('E0_reconstructed_old',np.array_equal(t.E0,old.E0));ids,w,local=geometry(t,a,g,10);oldw=sparse.load_npz(V8/'01_geometry'/sec/'SOURCE_k10_W.npz');axis=np.load(V8/'01_geometry'/sec/'SOURCE_k10_axis.npz');check('exact_v8_geometry',np.array_equal(ids,axis['node_index']) and (w!=oldw).nnz==0 and np.array_equal(local,axis['geometry_local']))
 folder=OUT/'audit/validation';folder.mkdir(exist_ok=True);src,info=source(sec,t,a,g,10,folder);src2,info2=source(sec,t,a,g,10,folder);check('source_compact_cache_parity',all(np.array_equal(src[k],src2[k]) for k in src))
 tchanged=t.copy();at=np.flatnonzero(tchanged.E0)[0];tchanged.loc[at,'E0']=False
 try:source(sec,tchanged,a,g,10,folder);raise AssertionError('cache accepted changed domain')
 except AssertionError as e:check('changed_source_mask_invalidates_cache',str(e)=='DOMAIN_CACHE_INVALID')
 y=t[LAYERS[0]].to_numpy(float)[ids];r=glocal(y,w);oldr=pd.read_parquet(V8/'02_hotspots'/sec/'MP10_SPOT_k10.parquet');check('real_section_direct_G_p_parity',np.array_equal(r.Gs,oldr.G) and np.array_equal(r.Zs,oldr.Z) and np.array_equal(r.p_sim.astype(float),oldr.p_sim_native))
 memo=_CACHED_VEC(10,len(y),9999,123);direct=_ORIGINAL_VEC(10,len(y),9999,123);check('bounded_permutation_memo_bitwise',np.array_equal(memo,direct));_CACHED_VEC.cache_clear()
 eff=effects(sec,t,a,g,src,9);old_eff=pd.read_parquet(V8/'03_niche_geometry'/f'{sec}__all_effects.parquet');old_eff=old_eff[old_eff.k.eq(10)&old_eff.analysis_domain.eq('MAIN')];keys=['layer','k','track','target_kind','target_id','metric','analysis_domain'];z=eff.merge(old_eff,on=keys,validate='one_to_one',suffixes=('_new','_old'));check('all_real_section_effects_exact',len(z)==456 and np.array_equal(z.effect_new,z.effect_old,equal_nan=True) and np.array_equal(z.eligible_new,z.eligible_old))
 # Unit key must collapse repeated section records even if dataset differs.
 tiny=eff.iloc[:1].copy();tiny['eligible']=True;tiny['effect']=1.;other=tiny.copy();other['dataset_id']='OTHER';other['section_id']='TECHNICAL_SECTION';other['effect']=3.;pair=pd.concat([tiny,other]);uid=pair.groupby('analysis_unit_uid').effect.mean();check('cross_dataset_unit_merge',len(uid)==1 and uid.iloc[0]==2.)
 # Input q1/q2 provenance and old domain reproduction, no imputation.
 q=np.load(V7/'01_ring123/ring123_arrays.npz');f=np.load(V2/'03_features/features.npz');check('ring_float32_provenance_tolerance',np.max(abs(q['q'][:,:2]-np.stack([f['q1'],f['q2']],1)))<1e-6)
 from rings import section_rings
 take=np.flatnonzero(m.section_id.eq(sec));meta={k:json.loads((V8/'01_geometry'/sec/'metadata.json').read_text())[k] for k in ['dataset_id','section_id','analysis_unit_uid','known_patient']};cr,ar=section_rings(m.iloc[take].reset_index(drop=True),q['q'][take],q['center_mean'][take],9,meta);oldar=pd.read_parquet(V7/'04_MP10_association/section_cache'/f'{sec}.parquet');oldar=oldar[oldar.analysis_set.eq('MAIN')].rename(columns={'accepted_label':'niche','ring':'endpoint','score_layer':'layer'});zz=ar.merge(oldar,on=['niche','celltype','endpoint','layer','method'],suffixes=('_new','_old'),validate='one_to_one');check('v7_partial_correlation_parity',np.allclose(zz.rho_new,zz.rho_old,rtol=0,atol=1e-12,equal_nan=True),dict(n=len(zz),max_abs=float(np.nanmax(abs(zz.rho_new-zz.rho_old)))))
 temp=temp_path('validation','small.tmp');temp.write_text('regenerable');keep=OUT/'audit/validation/keep.txt';keep.write_text('retain');clean('validation');check('owned_temp_cleanup_only',not temp.exists() and keep.exists() and (V8/'01_geometry'/sec/'nodes.parquet').exists())
 dest=OUT/'audit/validation/checkpoint';dest.mkdir(exist_ok=True);f=dest/'small.tsv';table(f,[{'x':1}]);mark(dest,'TEST_STAGE',[f]);check('atomic_checkpoint_resume',valid(dest,'TEST_STAGE'))
 js(OUT/'validation_summary.json',dict(status='PASSED',checks=checks,v8_earlier_evidence=str(V8/'reports'),time=time.time()))
if __name__=='__main__':validate()
