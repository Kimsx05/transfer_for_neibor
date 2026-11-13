from pathlib import Path
import pandas as pd,numpy as np,json,yaml,hashlib
R=Path(__file__).resolve().parents[1];V=R.parent/'hd_epi_niche_dual_methods_v1';E=R.parent/'hd_bin8_banksy_extension_v1';sections=['Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2']
if (R/'tables/T_v2.tsv').exists():print('FROZEN_ALREADY');raise SystemExit
for x in ['tables','models','figures','logs']: (R/x).mkdir(exist_ok=True)
a=pd.concat([pd.read_parquet(V/s/'anchor_manifest.parquet')[['section','bin_id','x_um','y_um','component','epi_anchor','epi_marker_detected_n','epi_mixed_flag','analysis_qc']] for s in sections],ignore_index=True);a['key']=a.section+'::'+a.bin_id;assert a.key.is_unique and len(a)==525552
a['common_evaluable']=True;over=[];manifest=[]
for ref in ['blca','hpca']:
 for r in [50,200]:
  p=V/f'{ref}_r{r}'/'labels.parquet';d=pd.read_parquet(p);assert np.array_equal(a.key,d.section+'::'+d.bin_id);a.common_evaluable &= d.evaluable
  manifest.append(dict(role='frozen_v1_labels',path=str(p),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns))
old=pd.read_csv(V/'blca_r50/training_anchor_ids.tsv',sep='\t');keys=old.section+'::'+old.bin_id;T=a.set_index('key').loc[keys].copy();T['old_key']=keys.to_numpy();T['replacement_reason']='unchanged';rng=np.random.default_rng(123)
for key,row in T[~T.common_evaluable].iterrows():
 pool=a[a.common_evaluable&(a.section==row.section)&~a.key.isin(T.index)];same=pool[(pool.x_um//400==row.x_um//400)&(pool.y_um//400==row.y_um//400)];pool=same if len(same) else pool
 chosen=pool.sort_values('key').iloc[rng.integers(len(pool))];T=T.drop(key);new=chosen.drop('key').to_dict();new.update(old_key=key,replacement_reason='not_in_all_four_S_same_block' if len(same) else 'not_in_all_four_S_same_section');T.loc[chosen.key]=new
T=T.sort_values(['section','bin_id']);assert len(T)==80000 and T.common_evaluable.all();assert (T.groupby('section').size()==20000).all()
for ref in ['blca','hpca']:
 for r in [50,200]:
  t=pd.read_csv(V/f'{ref}_r{r}'/'training_anchor_ids.tsv',sep='\t');k=set(t.section+'::'+t.bin_id);over.append(dict(reference=ref,radius=r,old_n=len(k),shared_with_T_v2=len(k&set(T.index)),removed=len(k-set(T.index)),added=len(set(T.index)-k)))
T.reset_index().rename(columns={'index':'key'}).to_csv(R/'tables/T_v2.tsv',sep='\t',index=False);a.to_parquet(R/'tables/frozen_anchors.parquet',index=False);pd.DataFrame(over).to_csv(R/'tables/training_overlap.tsv',sep='\t',index=False)
manifest += [dict(role='node_pool',path=str(V/s/'nodes.parquet'),bytes=(V/s/'nodes.parquet').stat().st_size,mtime_ns=(V/s/'nodes.parquet').stat().st_mtime_ns) for s in sections];pd.DataFrame(manifest).to_csv(R/'input_manifest.tsv',sep='\t',index=False)
config=dict(seed=123,seeds=[123,456,789],sections=sections,training_sections=sections[:4],radii_um=[50,200],radius_tolerance_um=1e-7,exclude_center=True,minimum_known=10,threads=4,night_window='23:00-07:00 Asia/Shanghai',training_source='v1 BLCA50; deterministic common-S replacement',block_um=400,K=list(range(2,11)),resolutions=[.005,.01,.02,.05,.1,.2,.4],knn_candidates=[15,30,50],composition_PCA='sqrt(P), centered, all numerical nonzero PCs, no whitening/no zscore',Harmony=dict(theta=1,lamb=5,max_iter_harmony=10,nclust=50,fit_scope='all_five_sections_transductive'),MENDER=dict(ring_width_um=25,include_self=False,normalization='native normalize_total then log1p',neighbor_graph='frequency-weighted identical-profile contraction adaptation'),CellCharter=dict(T_expr_per_reference_section=50000,HVG=3000,latent=10,max_epochs=100,early_stopping_patience=10,batch_size=1024,query_epochs=50,query_adaptation='SCVI.load_query_data MPBC2',GMM_K=list(range(2,11)),AutoK_max_runs=3),selection='simplest stable interpretable candidate; retain K2; no scalar cross-method score',buffer='test spatial blocks hash mod 7=0; remove training centers within 2R + tolerance; fixed representation conditional validation')
(R/'config.yaml').write_text(yaml.safe_dump(config,sort_keys=False));rows=[]
for method in ['hellinger_kmeans','adapted_leiden']:
 for ref in ['blca','hpca']:
  for r in [50,200]:
   for rep in ['NONE','WEAK']:rows.append(dict(method=method,reference=ref,radius=r,representation=rep,status='not_started',niche_fit='T_v2',representation_fit='T_v2' if rep=='NONE' else 'all_five_sections_transductive'))
for ref in ['blca','hpca']:
 for r in [50,200]:rows.append(dict(method='MENDER_adapted',reference=ref,radius=r,representation='ring_counts',status='not_started',niche_fit='T_v2',representation_fit='T_v2'))
for rep in ['scVI','scVI_without_Top20']:
 for r in [50,200]:rows.append(dict(method='CellCharter',reference='both_for_interpretation_only',radius=r,representation=rep,status='not_started',niche_fit='T_v2',representation_fit='reference_T_expr_plus_MPBC2_query_adaptation'))
pd.DataFrame(rows).to_csv(R/'run_matrix.tsv',sep='\t',index=False)
(R/'type_mapping.json').write_text((V/'type_mapping.json').read_text());print(pd.DataFrame(over).to_string(index=False));print('FROZEN',len(a),len(T),T.replacement_reason.value_counts().to_dict())
