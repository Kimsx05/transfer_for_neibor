from core import *
import itertools,platform,importlib.metadata

def prepare():
 if (OUT/'audit/input_manifest.tsv').exists():
  m=pd.read_csv(OUT/'audit/input_manifest.tsv',sep='\t');assert all(Path(r.path).exists() and sha(r.path)==r.sha256 for r in m.itertuples()),'INPUT_CHANGED';return
 q=quota();free=shutil.disk_usage(OUT).free;available=min(free,q.get('remaining_bytes',free));budget=min(3*1024**3,int(.4*available))
 cfg=dict(run_id='neighborhood_dataset_exhaustive_v9_20261018',datasets=DATASETS,seed=123,budget_bytes=budget,cache_budget_bytes=int(.25*budget),temp_budget_bytes=int(.15*budget),figure_budget_bytes=int(.2*budget),reserve_bytes=min(2*1024**3,max(256*1024**2,int(.2*available))),startup_effective_available_bytes=available,quota=q,combination_workers=1,section_workers=1,numerical_threads=1,primary='REFIT|SOURCE_P05|k10|ALL_POLICY_UNITS|NICHE|O_D|two_layers',ring_pairs='all K x 29 nonEpi; prior v7 manual list cannot be transferred',selection='v3 audited positive silhouettes, nondegenerate, moderate K then embedding silhouette and accepted fraction',source_domain='base geometry 4/8/12 & fullTop20 & epi valid & two finite scores & CURRENT accepted',frozen_source='v8 original E0; separate control family',k=10,permutations=9999,SpottedPy_commit='9dc364f203b7617b6701aa8fb579cbbb35357b26')
 js(OUT/'config.json',cfg)
 rows=[]
 for n in range(3,8):
  for ids in itertools.combinations(range(7),n):
   bits=''.join('1' if i in ids else '0' for i in range(7));ds=[DATASETS[i] for i in ids];rows.append(dict(combo_id=bits,datasets='|'.join(ds),n_datasets=n,contains_inhouse=5 in ids,combo_hash=digest(ds),status='REGISTERED'))
 rows.sort(key=lambda r:(r['n_datasets'],r['combo_id']));rows=[next(r for r in rows if r['n_datasets']==7)]+[r for r in rows if r['n_datasets']!=7]
 for i,r in enumerate(rows):r['execution_order']=i
 c=pd.DataFrame(rows);assert len(c)==99 and c.contains_inhouse.sum()==57 and c.combo_id.is_unique
 assert c.groupby('n_datasets').size().to_dict()=={3:35,4:35,5:21,6:7,7:1};table(OUT/'combinations.tsv',c)
 cov=pd.read_csv(V2/'00_preflight/coverage63.tsv',sep='\t');assert sorted(cov.dataset_id.unique())==sorted(DATASETS) and len(cov)==63
 paths=[V2/'03_features'/f for f in ['features.npz','common_anchors.parquet','feature_columns.json','COMPLETE.json']]+[V2/'00_preflight/coverage63.tsv',V3/'01_features/metadata.parquet',V3/'03_main/SAMPLE_R12_LEIDEN/labels.parquet',V7/'01_ring123/anchor_metadata.parquet',V7/'01_ring123/ring123_arrays.npz',V8/'07_synthesis/all_section_effects.parquet',V8/'00_audit/source_versions.json',ROOT/'neighborhood_niche_celltype_screen_v6/00_audit/analysis_unit_policy.tsv',ROOT/'neighborhood_mp10_calibration_v5/00_audit_and_spec/identity_and_independence.tsv',ROOT/'spatial_ecology_top20_v2/05_07_shared_K12/v1/section_identity_overlay.tsv']
 for d in [V2,V3,V7,V8]:paths+=list((d/'code').glob('*.py'))+list((d/'code').glob('*.R'));paths+=list((d/'reports').glob('METHODS.md'))
 paths+=list((V8/'01_geometry').glob('*/*'))+list((V8/'02_hotspots').glob('*/*.parquet'))
 paths=[p for p in paths if p.is_file()];log('hashing '+str(len(paths))+' read-only input files')
 table(OUT/'audit/input_manifest.tsv',[dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(set(paths))])
 m=pd.read_parquet(V2/'03_features/common_anchors.parquet');a=pd.read_parquet(V7/'01_ring123/anchor_metadata.parquet');assert m.spot_uid.is_unique and np.array_equal(m.spot_uid,a.spot_uid)
 pol=pd.read_csv(ROOT/'neighborhood_niche_celltype_screen_v6/00_audit/analysis_unit_policy.tsv',sep='\t');over=pd.read_csv(ROOT/'spatial_ecology_top20_v2/05_07_shared_K12/v1/section_identity_overlay.tsv',sep='\t');identity=cov[['dataset_id','section_id']].merge(over.drop(columns=['dataset_id'],errors='ignore'),on='section_id',how='left',validate='one_to_one');table(OUT/'audit/full_section_identity.tsv',identity)
 assert (m.epithelial_fraction>=.25).all();assert len(m)==129632
 contract=json.loads((V2/'03_features/COMPLETE.json').read_text());assert sha(V2/'03_features/features.npz')==contract['sha_features'] and sha(V2/'03_features/common_anchors.parquet')==contract['sha_common']
 table(OUT/'shared/master_anchors.parquet',m)
 base=(a[['R1_valid_n','R2_valid_n','R3_valid_n']].to_numpy()>=[4,8,12]).all(1);sc=base&a.epi_score_valid.to_numpy()&a.n_signature_available_model.eq(20).to_numpy()&a.n_signature_available_raw.eq(20).to_numpy()&np.isfinite(a[LAYERS]).all(1).to_numpy()
 assert np.array_equal(sc&a.accepted_label.ge(0),a.score_eligible)
 npz(OUT/'shared/base_domains.npz',base_geometry=base,base_score=sc,old_accepted=a.accepted_label.to_numpy(np.int16))
 identity_notes={'total_sections':len(cov),'sections_with_common_anchors':m.section_id.nunique(),'old_scored_sections':a.loc[a.score_eligible,'section_id'].nunique(),'inhouse':identity[identity.dataset_id.eq(INHOUSE)].to_dict('records'),'B17_B18':identity[identity.section_id.str.contains('B17|B18')].to_dict('records'),'policy':'Existing traceable unit mapping retained. Unknown inhouse patient identities remain independent sample policy, not confirmed patients. Cross-dataset unit aggregation uses analysis_unit_uid alone.'};js(OUT/'audit/identity_audit.json',identity_notes)
 resources={}
 for name,cmd in [('df',['df','-B1',str(OUT),str(OUT/'tmp')]),('inode',['df','-i',str(OUT)]),('memory',['free','-b']),('processes',['ps','-u','t070721','-o','pid,etime,pcpu,rss,comm']),('mounts',['cat','/proc/self/cgroup','/proc/self/limits'])]:resources[name]=subprocess.check_output(cmd,text=True)
 resources['cgroup_files']={str(p):p.read_text() for root in [Path('/sys/fs/cgroup')] for p in root.glob('**/memory.limit_in_bytes') if str(os.getuid()) in str(p)}
 js(OUT/'audit/resources_at_start.json',resources);js(OUT/'audit/availability.json',dict(inputs_complete=True,all_seven_present=True,n_anchors=len(m),n_reference_original=20000,quota=q))
 table(OUT/'cleanup_manifest.tsv',pd.DataFrame(columns=['time','path','bytes','stage','reason']))
 log('Preflight complete; 99 combinations registered')
if __name__=='__main__':prepare()
