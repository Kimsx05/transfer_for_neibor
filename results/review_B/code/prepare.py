from common import *
import re,subprocess,shutil

def approval():
 md=OUT/'00_approval/neighborhood_ring123_broad_codex_prompt_v7.md';gate=pd.read_csv(V6/'04_gate/candidate_review_TO_FILL.tsv',sep='\t',keep_default_na=False);rows=[]
 for line in md.read_text().splitlines():
  if not re.match(r'^\| niche[1-8] \|',line):continue
  cells=[x.strip() for x in line.strip('|').split('|')];n=cells[0]
  for tier,txt in zip(TIERS,cells[1:]):
   for c in txt.split('；'):
    c=c.strip()
    if c=='—':continue
    hit=gate[gate.scope.eq('ALL_DATASETS')&gate.niche.eq(n)&gate.celltype.eq(c)];assert len(hit)==1,(n,c);h=hit.iloc[0]
    rows.append(dict(master_candidate_id=h.candidate_id,pair_id=f'n{h.accepted_label}_c{TYPES.index(c):02d}',niche=n,accepted_label=int(h.accepted_label),celltype=c,cell_index=TYPES.index(c),review_tier=tier))
 p=pd.DataFrame(rows);assert len(p)==139 and p.master_candidate_id.is_unique and p.pair_id.is_unique and p.celltype.nunique()==29 and 'Epithelial' not in set(p.celltype);assert p.review_tier.value_counts().to_dict()==dict(zip(TIERS,[28,32,39,40]))
 tab('00_approval/approved_pairs_v7.tsv',p)
 expanded=gate.merge(p,on=['niche','accepted_label','celltype'],validate='many_to_one');assert len(expanded)==1112 and expanded.candidate_id.is_unique and expanded.groupby('master_candidate_id').scope.nunique().eq(8).all()
 expanded['decision']='GO';expanded['approved_scope']=expanded.scope;expanded['reason']='User conversation blanket approval of exact v7 appendix tiers; no manual row review';expanded['reviewer']='USER_CONVERSATION_APPROVAL';expanded['review_date']='2026-10-14';tab('00_approval/approved_gate_v7.tsv',expanded)
 receipt=dict(authorization='Explicit current user request authorizes all 139 combinations in all eight scopes',manual_row_review=False,pathology_review_implied=False,source='Only Markdown uploaded; reconstructed exact appendix per its final explicit fallback instruction; true IDs joined from original v6 Gate',n_pairs=139,n_expanded=1112,tiers=p.review_tier.value_counts().to_dict(),original_v6_gate_sha256=sha(V6/'04_gate/candidate_review_TO_FILL.tsv'),spec_sha256=sha(md),pairs_sha256=sha(OUT/'00_approval/approved_pairs_v7.tsv'),gate_sha256=sha(OUT/'00_approval/approved_gate_v7.tsv'),date='2026-10-14')
 js('00_approval/approval_receipt.json',receipt);js('00_approval/provenance_v7.json',receipt);log('139 exact approved pairs / 1112 true scoped IDs verified')

def prepare():
 if (OUT/'00_audit/input_manifest.tsv').exists():
  man=pd.read_csv(OUT/'00_audit/input_manifest.tsv',sep='\t');assert all(sha(r.path)==r.sha256 for r in man.itertuples()),'INPUT_CHANGED';log('M00 existing input hashes verified');return
 approval()
 for cmd,name in [(['df','-h',str(OUT)],'disk'),(['ps','-u','t070721','-o','pid,etime,pcpu,pmem,comm'],'processes')]: (OUT/'00_audit'/f'{name}_at_start.txt').write_text(subprocess.check_output(cmd,text=True))
 spec=dict(seed=123,configuration='SAMPLE_R12_LEIDEN',n_pairs=139,scopes=SCOPES,ring_min_valid=[4,8,12],ring_complete=[6,12,18],section_min_target=20,section_min_rest=20,association_min_centers=30,min_distinct=3,signed_min_units=5,signed_min_nonzero=5,signed_zero_tolerance=0.,signed_rule='n<50 no zero and no absolute ties exact, else continuity-corrected approx; Wilcox exclusion',composition_tolerance=1e-6,rank_variance_tolerance=1e-12,QR_relative_tolerance=1e-10,Fisher_clip=1e-7,TME_min_share_all_rings=.05,main_covariates=['intercept','rank(center_mean_Epi)','rank(log1p(center_nCount))'],center_type_sensitivity='add rank(center_mean_celltype)',block='floor(array_row/8), floor(array_col/8); stable SHA256(seed123|spot_uid) minimum per block; shared across all targets/rings/layers',families=dict(background_zero=417,ring_pairwise=417,association_per_method=1008),sets=SETS,main_plots='3 panels of 3 unconnected boxes: raw/excess/ringpairs; association 2 score layers x main/raw panels',only_sensitivity_tests='verified-patient subset separate BH; other sensitivities descriptive',approval='current user blanket authorization; original v6 gate read only')
 js('00_audit/spec_frozen.json',spec)
 paths=[(V6/'00_audit/input_manifest.tsv','v6 prior provenance'),(V6/'reports/METHODS.md','v6 definition'),(V6/'00_audit/analysis_unit_policy.tsv','unit policy'),(V6/'00_audit/spec_frozen.yaml','v6 thresholds'),(V6/'04_gate/candidate_review_TO_FILL.tsv','original gate'),(V6/'04_gate/all_niche_celltype_evidence.tsv','all 270 per scope evidence'),(V2/'03_features/common_anchors.parquet','original matrix UID axis'),(V2/'03_features/features.npz','frozen q1 q2'),(V2/'03_features/feature_columns.json','30 type axes'),(V2/'03_features/COMPLETE.json','producer hash contract'),(V2/'03_features/anchor_master.parquet','full tissue graph axes'),(V2/'code/01_prepare.py','physical graph producer'),(V2/'00_preflight/coverage63.tsv','63 sections'),(V3/'01_features/metadata.parquet','analysis units and common metadata'),(V3/'03_main/SAMPLE_R12_LEIDEN/labels.parquet','formal labels'),(V3/'02_R12_pilot/SAMPLE_WEAK_R12/LEIDEN/labels.parquet','v7 referenced alias verification'),(V4/'04_niche_comparison/anchor_MP10_all_layers.parquet','frozen scores'),(V5/'00_audit_and_spec/identity_and_independence.tsv','unit evidence'),(V5/'code/06_contrasts.py','prior continuous implementation'),(V5/'05_evidence_synthesis/evidence_matrix_all9.tsv','source limitations')]
 cov=pd.read_csv(V2/'00_preflight/coverage63.tsv',sep='\t')
 for r in cov.itertuples():paths.extend([(r.q05_path,'full tissue q05 '+r.section_id),(r.mean_path,'center posterior mean '+r.section_id)]+[(V2/'02_graphs'/r.dataset_id/(r.section_id+'__'+ring+'.npz'),'old '+ring+' graph') for ring in ['R1','R2']])
 paths += [(OUT/'00_approval/approved_pairs_v7.tsv','approved exact pairs'),(OUT/'00_approval/approved_gate_v7.tsv','approved scoped IDs')]
 tab('00_audit/input_manifest.tsv',[input_record(p,why) for p,why in paths]);shutil.copy2(V6/'00_audit/analysis_unit_policy.tsv',OUT/'00_audit/analysis_unit_policy.tsv');shutil.copy2(V6/'04_gate/all_niche_celltype_evidence.tsv',OUT/'00_approval/v6_all_evidence_readonly.tsv')
 a=pd.read_parquet(V3/'03_main/SAMPLE_R12_LEIDEN/labels.parquet').set_index('spot_uid');b=pd.read_parquet(V3/'02_R12_pilot/SAMPLE_WEAK_R12/LEIDEN/labels.parquet').set_index('spot_uid');assert a.index.is_unique and b.index.is_unique;assert np.array_equal(a.accepted_label,b.loc[a.index,'accepted_label']);tab('00_audit/label_alias_check.tsv',[dict(formal_path=str(V3/'03_main/SAMPLE_R12_LEIDEN/labels.parquet'),referenced_alias=str(V3/'02_R12_pilot/SAMPLE_WEAK_R12/LEIDEN/labels.parquet'),UID_label_identity=True)])
 (OUT/'00_audit/v5_covariate_difference.md').write_text('v5 used unranked standardized center Epi, log1p(model_UMI), and central corresponding type (only GC B/SPP1 TAM/Neutrophil), then regressed ranked outcomes. v7 explicit primary is intercept+rank(center mean Epi)+rank(log1p(center_nCount)); corresponding type rank is separate sensitivity for all 29 nonEpi types. No v5 coefficients or p-values are substituted. The array_row//8,array_col//8 blocks used in v5 direct-pair sensitivity provide the fixed physical 8x8 definition; v7 selects at most one scored eligible anchor per block with stable SHA256 seed123.\n')
 checkpoint('M00_APPROVAL_AND_INPUTS_VERIFIED')

def ring3(g1,g2):
 r=(g2.astype(np.int32)@g1.astype(np.int32));r.data[:]=1;r=r.tocsr();r=r-r.multiply(g1)-r.multiply(g2);r.setdiag(0);r.eliminate_zeros();assert (r.data==1).all();return r.astype(np.int8)

def rings():
 if done('M01_RING123_READY'):log('M01 checkpoint reused');return
 master=pd.read_parquet(V2/'03_features/anchor_master.parquet');m=pd.read_parquet(V3/'01_features/metadata.parquet');original=pd.read_parquet(V2/'03_features/common_anchors.parquet');old=np.load(V2/'03_features/features.npz');contract=json.loads((V2/'03_features/COMPLETE.json').read_text());assert sha(V2/'03_features/features.npz')==contract['sha_features'] and sha(V2/'03_features/common_anchors.parquet')==contract['sha_common'];lab=pd.read_parquet(V3/'03_main/SAMPLE_R12_LEIDEN/labels.parquet');mp=pd.read_parquet(V4/'04_niche_comparison/anchor_MP10_all_layers.parquet')
 for d in [m,original,lab,mp]:assert d.spot_uid.is_unique and set(d.spot_uid)==set(m.spot_uid)
 m=m.merge(lab[['spot_uid','accepted_label']],on='spot_uid',validate='one_to_one').merge(mp[['spot_uid']+LAYERS+['epi_score_valid','n_signature_available_model','n_signature_available_raw','center_nCount']],on='spot_uid',validate='one_to_one');q=np.full((len(m),3,30),np.nan);center=np.full((len(m),30),np.nan);validn=np.zeros((len(m),3),int);geom=np.zeros_like(validn);cov=pd.read_csv(V2/'00_preflight/coverage63.tsv',sep='\t');aud=[];coverage=[]
 oldq=np.stack([align(original.spot_uid,m.spot_uid,old['q1']),align(original.spot_uid,m.spot_uid,old['q2'])],axis=1)
 for sn,r in enumerate(cov.itertuples()):
  t=master[master.section_id.eq(r.section_id)].reset_index(drop=True);ix=m.index[m.section_id.eq(r.section_id)].to_numpy();local=m.loc[ix];take=pd.Index(t.spot_uid).get_indexer(local.spot_uid);assert (take>=0).all() and t.spot_uid.is_unique
  raw=pd.read_csv(r.q05_path,index_col=0).drop(columns=['section','dataset'],errors='ignore');assert raw.index.is_unique and set(raw.columns)==set(TYPES);raw=raw.reindex(index=t.barcode,columns=TYPES).to_numpy(float);ok=np.isfinite(raw).all(1)&(raw>=0).all(1)&(raw.sum(1)>0);assert np.array_equal(ok,t.composition_valid);p=np.zeros_like(raw);p[ok]=raw[ok]/raw[ok].sum(1,keepdims=True)
  g1=sparse.load_npz(V2/'02_graphs'/r.dataset_id/(r.section_id+'__R1.npz')).tocsr();g2=sparse.load_npz(V2/'02_graphs'/r.dataset_id/(r.section_id+'__R2.npz')).tocsr();g3=ring3(g1,g2)
  assert g1.shape==g2.shape==g3.shape==(len(t),len(t));assert not g1.multiply(g2).nnz and not g1.multiply(g3).nnz and not g2.multiply(g3).nnz;assert not g3.diagonal().any() and (g3-g3.T).nnz==0
  # Exact shortest-path BFS spotchecks in every section, independent of matrix product.
  for i in np.unique(np.r_[0,take[:2]]).astype(int):
   seen={i};front={i}
   for hop,gg in enumerate([g1,g2,g3],1):
    front={j for z in front for j in g1[z].indices}-seen;seen|=front;assert front==set(gg[i].indices)
  dest=OUT/'01_ring123/graphs'/r.dataset_id;dest.mkdir(parents=True,exist_ok=True)
  for rn,g in enumerate([g1,g2,g3]):
   n=np.asarray(g@ok.astype(int)).ravel();ng=np.asarray(g.sum(1)).ravel();qq=np.full((len(ix),30),np.nan);np.divide(g[take]@p,n[take,None],out=qq,where=n[take,None]>0);q[ix,rn]=qq;validn[ix,rn]=n[take];geom[ix,rn]=ng[take]
   if rn<2 and len(ix):assert np.nanmax(abs(qq-oldq[ix,rn]))<1e-6,(r.section_id,rn)
   sparse.save_npz(dest/f'{r.section_id}__{RINGS[rn]}.npz',g)
  t[['spot_uid','barcode','array_row','array_col','composition_valid']].to_parquet(dest/(r.section_id+'__node_axis.parquet'),index=False)
  mean=pd.read_csv(r.mean_path,index_col=0);mean.columns=[c.replace('meanscell_abundance_w_sf_','').replace('means_cell_abundance_w_sf_','') for c in mean.columns];assert set(TYPES).issubset(mean.columns) and mean.index.is_unique;mv=mean.reindex(index=local.barcode,columns=TYPES).to_numpy(float);valid=np.isfinite(mv).all(1)&(mv>=0).all(1)&(mv.sum(1)>0);center[ix[valid]]=mv[valid]/mv[valid].sum(1,keepdims=True)
  if len(ix):assert np.allclose(center[ix,EPI],local.epithelial_fraction,atol=1e-8)
  aud.append(dict(dataset_id=r.dataset_id,section_id=r.section_id,n_tissue_nodes=len(t),n_missing_composition=int((~ok).sum()),n_common=len(ix),q1_max_error=float(np.max(abs(q[ix,0]-oldq[ix,0]))) if len(ix) else 0,q2_max_error=float(np.max(abs(q[ix,1]-oldq[ix,1]))) if len(ix) else 0,R3_max_geometric=int(np.max(g3.sum(1))),rings_disjoint=True,BFS_exact=True))
  log(f'M01 section {sn+1}/63 {r.section_id}; anchors={len(ix)}')
 assert np.isfinite(q).all() and (q>=0).all() and np.allclose(q.sum(2),1,atol=1e-6)
 for j,r in enumerate(RINGS):m[r+'_valid_n']=validn[:,j];m[r+'_geometric_n']=geom[:,j]
 m['main_eligible']=m.accepted_label.ge(0)&(validn>=np.array([4,8,12])).all(1);m['complete123']=(validn==[6,12,18]).all(1)&(geom==[6,12,18]).all(1);m['score_eligible']=m.main_eligible&m.epi_score_valid&m.n_signature_available_model.eq(20)&m.n_signature_available_raw.eq(20)&np.isfinite(m[LAYERS]).all(1);m['TME_eligible']=m.main_eligible&(q[:,:,NON].sum(2)>=.05).all(1)
 for cr in cov.itertuples():
  g=m[m.section_id.eq(cr.section_id)]
  for k in [-1]+list(range(9)):
   z=g[g.accepted_label.eq(k)];coverage.append(dict(dataset_id=cr.dataset_id,section_id=cr.section_id,accepted_label=k,n_common=len(z),n_accepted=len(z) if k>=0 else 0,n_R123_main=int(z.main_eligible.sum()),n_complete123=int((z.main_eligible&z.complete123).sum()),n_score=int(z.score_eligible.sum()),n_TME=int(z.TME_eligible.sum()),n_R3_insufficient=int(z.R3_valid_n.lt(12).sum()),status='NO_COMMON_ANCHORS' if not len(g) else 'UNASSIGNED_EXCLUDED' if k==-1 else 'AVAILABLE'))
 m.to_parquet(OUT/'01_ring123/anchor_metadata.parquet',index=False);np.savez_compressed(OUT/'01_ring123/ring123_arrays.npz',spot_uid=m.spot_uid.to_numpy(dtype=str),celltypes=np.array(TYPES),rings=np.array(RINGS),q=q,center_mean=center)
 tab('01_ring123/coverage_all_steps.tsv',coverage);tab('01_ring123/graph_and_cache_audit.tsv',aud);tab('01_ring123/cache_hashes.tsv',[input_record(p,'new ring cache') for p in sorted((OUT/'01_ring123').rglob('*')) if p.is_file()]);checkpoint('M01_RING123_READY')
if __name__=='__main__':prepare();rings()
