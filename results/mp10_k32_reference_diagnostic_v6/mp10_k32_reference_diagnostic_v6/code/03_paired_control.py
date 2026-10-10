from common import *
import anndata as ad
from scipy import sparse
O=R/'04_paired_depth_control';old=S/'04_pseudospots/DEV';cfg=json.loads((R/'CONFIG.yaml').read_text());depths=json.loads((S/'04_pseudospots/SIMULATION_FROZEN.json').read_text())['depth_targets_full_RNA_UMI'];labels=json.loads((R/'00_manifest/labels_S31.json').read_text())
if (O/'READY.json').exists():
 st=json.loads((O/'READY.json').read_text());assert st['control_config_hash']==sha(R/'CONTROL_FROZEN.yaml')
 for rule,h in st['counts_hashes'].items():assert sha(O/rule/'counts.h5ad')==h
 say('VALIDATED_CONTROL_DONE');raise SystemExit
m=pd.read_csv(old/'manifest.tsv.gz',sep='\t');mm=pd.read_csv(old/'cell_membership_RNA.tsv.gz',sep='\t');source=pd.read_csv(S/'01_scope/actual_cell_split_manifest.tsv.gz',sep='\t',low_memory=False).set_index('cell_id')
def wanted(d):return ((d.scenario.eq('S1_regular')&d.N_epithelial.eq(10)&d.N_K32.isin([0,1,2,5,10]))|(d.scenario.eq('S3_pure_epithelial')&d.N_epithelial.eq(20)&d.N_K32.isin([0,10,20]))|(d.scenario.str.startswith('S2_')&d.N_epithelial.eq(10)&d.N_K32.eq(0)))
bags=m[wanted(m)&m['repeat'].isin([0,1])&m.depth_condition.eq('low')].copy().sort_values(['sample','scenario','N_K32','repeat']);assert bags.mixture_id.is_unique
members=[]
for z in bags.itertuples():
 q=mm[mm.pseudospot_id.eq(z.pseudospot_id)].sort_values('cell_id').copy();assert len(q)==20 and q.cell_id.is_unique
 orig=source.loc[q.cell_id];assert orig.split.eq('DEV').all() and orig['sample'].nunique()==1 and orig['sample'].iloc[0]==z.sample
 q=q[['cell_id','subtype','pre_full_UMI','pre_shared_UMI']];q['bag_id']=z.mixture_id;members.append(q)
bm=pd.concat(members,ignore_index=True);hs=bm.groupby('bag_id').cell_id.agg(lambda x:key(sorted(x)));bags['member_hash']=bags.mixture_id.map(hs);bags['bag_id']=bags.mixture_id;bags['old_anchor_spot']=bags.pseudospot_id;bags=bags.drop(columns=['pseudospot_id','depth_condition','target_full_UMI','thinning_probability','post_full_UMI','post_shared_UMI','depth_target_unattainable','any_cell_zero_shared_RNA','K32_shared_RNA_zero']);bags['unique_cell_set_first']=~bags.member_hash.duplicated();median=float(bags.drop_duplicates('member_hash').pre_full_UMI.median());probs={d:float(t/median) for d,t in depths.items()};wt(bags,O/'cell_bag_manifest.tsv');wt(bm,O/'cell_bag_members.tsv.gz');wt(bm.groupby('cell_id').size().rename('bag_reuse_count').reset_index(),O/'cell_reuse.tsv')
avail=pd.read_csv(old/'scenario_availability.tsv',sep='\t');avail=avail[wanted(avail)].copy();counts=bags.groupby(['sample','scenario','N_epithelial','N_K32']).size();avail['selected_bags']=[int(counts.get((x.sample,x.scenario,x.N_epithelial,x.N_K32),0)) for x in avail.itertuples()];assert (avail.loc[avail.available,'selected_bags']==2).all();assert (avail.loc[~avail.available,'selected_bags']==0).all();wt(avail,O/'scenario_availability.tsv')
plan={'N0':20,'selected_old_repeats':[0,1],'n_bags':len(bags),'n_unique_cell_sets':int(bags.member_hash.nunique()),'unique_definition_for_probability':'unique sorted cell-ID member_hash, not depth replicates','median_pre_full_UMI':median,'target_full_RNA_UMI':depths,'fixed_capture_probabilities':probs,'fixed_capture_feasible':{d:p<=1 for d,p in probs.items()},'full_RNA_definition':'same per-cell raw_UMI as v5; audited full-gene source, not G0','bag_manifest_sha256':sha(O/'cell_bag_manifest.tsv'),'bag_members_sha256':sha(O/'cell_bag_members.tsv.gz'),'source_split':'DEV only; unchanged conservative source groups','selection':'all feasible prescribed sample/scenario combinations, old repeats 0 and 1, exact same cell IDs across both rules and all depths','rng':'nested binomial over sorted descending distinct probabilities within bag; original G0 entries and non-G0 remainder sampled independently; correct Binomial(n,p) marginal for each condition','seeds':{b:seed('v6_control_thinning|'+b) for b in bags.bag_id},'not_old_predictions':'new counts and a smaller mapping collection; not claimed numeric replication of old predictions','mapping_pair_seed':seed('v6_control_mapping'),'mapping_posterior_pair_seed':seed('v6_control_posterior')}
js(R/'CONTROL_FROZEN.yaml',plan)
a=ad.read_h5ad(S/'02_inputs/DEV_counts.h5ad');ix=pd.Series(np.arange(a.n_obs),index=a.obs_names);raw=pd.read_csv(V/'00_audit/reference_counts_per_cell.tsv.gz',sep='\t').set_index('cell_id');out={r:{'rows':[],'X':[],'members':[],'truth':[],'RNA':[]} for r in ['TARGET_UMI','FIXED_CAPTURE']}
for z in bags.itertuples():
 mem=bm[bm.bag_id.eq(z.bag_id)].copy();cells=mem.cell_id.tolist();X=a.X[ix.loc[cells].to_numpy()].tocsr();full=mem.pre_full_UMI.to_numpy(np.int64);g0=np.asarray(X.sum(1)).ravel().astype(np.int64);assert np.array_equal(g0,mem.pre_shared_UMI) and np.array_equal(full,raw.loc[cells,'raw_UMI']);assert full.sum()==z.pre_full_UMI and g0.sum()==z.pre_shared_UMI
 cond=[]
 for d,t in depths.items():
  cond.append(('TARGET_UMI',d,min(1.,t/full.sum())))
  if probs[d]<=1:cond.append(('FIXED_CAPTURE',d,probs[d]))
 rng=np.random.default_rng(plan['seeds'][z.bag_id]);current=X.data.astype(np.int64).copy();remainder=full-g0;previous=1.;sampled={}
 for p in sorted({q[2] for q in cond},reverse=True):
  current=rng.binomial(current,p/previous);remainder=rng.binomial(remainder,p/previous);sampled[p]=(current.copy(),remainder.copy());previous=p
 truth=mem.subtype.value_counts().reindex(labels,fill_value=0)
 for rule,d,p in cond:
  vals,rem=sampled[p];Y=X.copy();Y.data=vals.astype(np.int32);Y.eliminate_zeros();kept=np.asarray(Y.sum(1)).ravel();fullkept=kept+rem;sid=rule+'__'+z.bag_id+'__'+d;t=out[rule];t['X'].append(sparse.csr_matrix(Y.sum(0),dtype=np.int32));info=z._asdict();info.pop('Index',None);info.update(pseudospot_id=sid,rule=rule,depth_condition=d,thinning_probability=p,target_full_UMI=depths[d],post_full_UMI=int(fullkept.sum()),post_shared_UMI=int(kept.sum()),depth_target_unattainable=bool(rule=='TARGET_UMI' and z.pre_full_UMI<depths[d]),any_cell_zero_shared_RNA=bool((kept==0).any()),K32_shared_RNA_zero=bool(z.N_K32>0 and kept[mem.subtype.eq('Epi_K32')].sum()==0));t['rows'].append(info);t['truth'].append({'pseudospot_id':sid,**truth.to_dict()})
  for j,cid in enumerate(cells):t['members'].append({'pseudospot_id':sid,'bag_id':z.bag_id,'cell_id':cid,'subtype':mem.subtype.iloc[j],'pre_full_UMI':int(full[j]),'pre_G0_UMI':int(g0[j]),'post_full_UMI':int(fullkept[j]),'post_G0_UMI':int(kept[j]),'zero_retained_G0_RNA':bool(kept[j]==0)})
  for label in labels:
   mask=mem.subtype.eq(label).to_numpy();t['RNA'].append({'pseudospot_id':sid,'subtype':label,'true_cells':int(truth[label]),'post_full_UMI':int(fullkept[mask].sum()),'post_G0_UMI':int(kept[mask].sum()),'true_positive_type_zero_RNA':bool(truth[label]>0 and kept[mask].sum()==0)})
checks=[];hashes={}
for rule,t in out.items():
 o=O/rule;o.mkdir(exist_ok=True);d=pd.DataFrame(t['rows']);mat=sparse.vstack(t['X'],format='csr');assert np.array_equal(np.asarray(mat.sum(1)).ravel(),d.post_shared_UMI);b=ad.AnnData(X=mat,obs=d[['pseudospot_id','sample']].set_index('pseudospot_id'),var=a.var.copy());b.write_h5ad(o/'counts.h5ad',compression='lzf');wt(d,o/'manifest.tsv.gz');wt(pd.DataFrame(t['members']),o/'cell_membership_RNA.tsv.gz');wt(pd.DataFrame(t['RNA']),o/'subtype_RNA_contribution.tsv.gz');wt(pd.DataFrame(t['truth']),o/'truth_counts.tsv.gz');hashes[rule]=sha(o/'counts.h5ad');checks.append({'rule':rule,'n_spots':len(d),'n_bags':d.bag_id.nunique(),'all_N20':bool(d.N_total.eq(20).all()),'only_DEV':bool(d.split.eq('DEV').all()),'member_hashes_match':bool(d.member_hash.eq(d.bag_id.map(hs)).all()),'no_predictive_truth_covariates':list(b.obs.columns)==['sample']});say('CONTROL_READY',rule,b.shape)
wt(pd.DataFrame(checks),O/'validation.tsv');js(O/'READY.json',{'control_config_hash':sha(R/'CONTROL_FROZEN.yaml'),'counts_hashes':hashes,'bag_manifest_sha256':sha(O/'cell_bag_manifest.tsv'),'n_bags':len(bags),'fit_not_started_during_creation':True})
