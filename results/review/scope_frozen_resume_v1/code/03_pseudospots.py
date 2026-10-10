from common import *
import anndata as ad
from scipy import sparse
cfg=json.loads((R/'CONFIG.yaml').read_text());freeze=json.loads((R/'02_inputs/INPUTS_READY.json').read_text());labels=json.loads((R/'02_inputs/labels.json').read_text());meta=pd.read_csv(R/'01_scope/actual_cell_split_manifest.tsv.gz',sep='\t',low_memory=False).set_index('cell_id');raw=pd.read_csv(V/'00_audit/reference_counts_per_cell.tsv.gz',sep='\t').set_index('cell_id');spqc=pd.read_csv(V/'00_audit/spatial_UMI_QC.tsv',sep='\t');depths=dict(zip(['low','medium','high'],[int(round(spqc[c].median())) for c in ['umi_p25','umi_median','umi_p75']]))
plan={'N0':20,'replicates':10,'depth_targets_full_RNA_UMI':depths,'thinning':'independent binomial per retained UMI; p=min(1,target/full_pre_UMI); no upscaling','within_source_only':True,'K32_truth':'frozen K32 membership, not MP10-high','TME':'uniform cell draws, requiring >=2 observed TME subtypes','scenario_counts':'S1 E5/10/16,K distinct 0/1/2/round(E/2)/E; S2 five strata same E; S3 E20 K0/1/2/10/20; S4 E10 K5 low/high median, identical other/background cells','source_split_hash':freeze['split_hash'],'member_label_hash':freeze['scope_member_label_hash'],'gene_order_hash':freeze['gene_order_hash'],'seeds':{s:seed('pseudospots_'+s) for s in ['DEV','TEST']}}
js(R/'04_pseudospots/SIMULATION_FROZEN.json',plan)
for part in ['DEV','TEST']:
 O=R/'04_pseudospots'/part
 if (O/'READY.json').exists():
  old_ready=json.loads((O/'READY.json').read_text());assert old_ready['plan_hash']==sha(R/'04_pseudospots/SIMULATION_FROZEN.json');assert sha(O/'counts.h5ad')==old_ready['counts_sha256'];assert sha(O/'truth_counts.tsv.gz')==old_ready['truth_sha256'];continue
 a=ad.read_h5ad(R/f'02_inputs/{part}_counts.h5ad');md=meta.loc[a.obs_names].copy();md['row']=np.arange(a.n_obs);rng=np.random.default_rng(seed('pseudospots_'+part));rows=[];members=[];truth=[];contrib=[];mat=[];avail=[];mixnum=0
 def choose(pool,n):return rng.choice(np.asarray(pool,dtype=int),n,replace=False) if n else np.array([],int)
 def tme(pool,n):
  if n==0:return np.array([],int)
  for attempt in range(100):
   ix=choose(pool,n)
   if md.iloc[ix].new_reference_label.nunique()>=2:return ix
  raise RuntimeError('TME rejection sampler failed; do not silently change composition')
 def emit(ids,info,rep,intensity='not_applicable',pair=''):
  global mixnum
  assert len(ids)==20 and len(set(ids))==20;z=md.iloc[ids];assert z['sample'].nunique()==1 and z.split.eq(part).all() and z.included_in_this_run.all();X=a.X[ids].tocsr();counts=z.new_reference_label.value_counts().reindex(labels,fill_value=0);full=raw.loc[z.index,'raw_UMI'].to_numpy().astype(np.int64);shared=np.asarray(X.sum(1)).ravel().astype(np.int64);assert np.all(full>=shared);mixnum+=1;mix=f'{part}_mix{mixnum:06d}'
  for depth,target in depths.items():
   p=min(1.,target/full.sum());Y=X.copy();Y.data=rng.binomial(Y.data.astype(np.int64),p).astype(np.int32);Y.eliminate_zeros();kept=np.asarray(Y.sum(1)).ravel();other=rng.binomial(full-shared,p);full_kept=kept+other;spot=f'{part}_spot{len(rows)+1:07d}';sumrow=sparse.csr_matrix(Y.sum(0),dtype=np.int32);mat.append(sumrow)
   d={**info,'pseudospot_id':spot,'mixture_id':mix,'repeat':rep,'depth_condition':depth,'target_full_UMI':target,'thinning_probability':p,'pre_full_UMI':int(full.sum()),'post_full_UMI':int(full_kept.sum()),'pre_shared_UMI':int(shared.sum()),'post_shared_UMI':int(kept.sum()),'depth_target_unattainable':bool(p==1 and full.sum()<target),'N_total':20,'N_K32':int(counts.Epi_K32),'N_other':int(counts.Epi_other),'N_epithelial':int(counts.Epi_K32+counts.Epi_other),'state_intensity':intensity,'intensity_pair_id':pair,'K32_mean_UCell':float(z.loc[z.new_reference_label=='Epi_K32','ucell10'].mean()),'K32_mean_original_UMI':float(full[z.new_reference_label.eq('Epi_K32').to_numpy()].mean()) if counts.Epi_K32 else np.nan,'cross_sample_mixture':False,'any_cell_zero_shared_RNA':bool((kept==0).any()),'K32_shared_RNA_zero':bool(counts.Epi_K32>0 and kept[z.new_reference_label.eq('Epi_K32').to_numpy()].sum()==0)};rows.append(d);truth.append({'pseudospot_id':spot,**counts.to_dict()})
   for j,cid in enumerate(z.index):members.append({'pseudospot_id':spot,'cell_id':cid,'subtype':z.new_reference_label.iloc[j],'pre_full_UMI':int(full[j]),'post_full_UMI':int(full_kept[j]),'pre_shared_UMI':int(shared[j]),'post_shared_UMI':int(kept[j])})
   for lab in labels:
    b=z.new_reference_label.eq(lab).to_numpy();contrib.append({'pseudospot_id':spot,'subtype':lab,'true_cells':int(counts[lab]),'post_full_UMI':int(full_kept[b].sum()),'post_shared_UMI':int(kept[b].sum()),'true_positive_type_zero_RNA':bool(counts[lab]>0 and kept[b].sum()==0)})
 for sample,z in md.groupby('sample',observed=True):
  kp=z.loc[z.new_reference_label=='Epi_K32','row'].to_numpy();op=z.loc[z.new_reference_label=='Epi_other','row'].to_numpy();bg=z.loc[~z.new_reference_label.isin(['Epi_K32','Epi_other']),'row'].to_numpy();base={'split':part,'dataset':z.dataset.iloc[0],'sample':sample,'source_group':z.source_group.iloc[0]}
  def condition(scene,E,K,pool):
   feasible=len(kp)>=K and len(pool)>=E-K and len(bg)>=20-E and (E==20 or md.iloc[bg].new_reference_label.nunique()>=2);avail.append({**base,'scenario':scene,'N_epithelial':E,'N_K32':K,'n_K32_pool':len(kp),'n_other_pool':len(pool),'n_TME_pool':len(bg),'available':feasible})
   if feasible:
    for rep in range(10):emit(np.concatenate([choose(kp,K),choose(pool,E-K),tme(bg,20-E)]),{**base,'scenario':scene},rep)
  for E in [5,10,16]:
   for K in sorted({0,1,2,round(E/2),E}):condition('S1_regular',E,K,op)
  for st,mask in [('C14',z.cluster.eq('C14')),('C7',z.cluster.eq('C7')),('C17',z.cluster.eq('C17')),('UCell_high_non_K32',z.UCell_high.eq(True)),('NMF_high_non_K32',z.NMF_high.eq(True))]:
   pool=z.loc[mask&z.new_reference_label.eq('Epi_other'),'row'].to_numpy()
   for E in [5,10,16]:condition('S2_'+st,E,0,pool)
  for K in [0,1,2,10,20]:condition('S3_pure_epithelial',20,K,op)
  zz=z[z.new_reference_label=='Epi_K32'];med=zz.ucell10.median();lo=zz.loc[zz.ucell10<med,'row'].to_numpy();hi=zz.loc[zz.ucell10>med,'row'].to_numpy();ok=min(len(lo),len(hi))>=5 and len(op)>=5 and len(bg)>=10 and md.iloc[bg].new_reference_label.nunique()>=2;avail.append({**base,'scenario':'S4_intensity','N_epithelial':10,'N_K32':5,'available':ok,'n_low_K32':len(lo),'n_high_K32':len(hi),'within_K32_UCell_median':med})
  if ok:
   for rep in range(10):
    fixed=np.concatenate([choose(op,5),tme(bg,10)]);pair=sample+'|'+str(rep)
    for st,pool in [('low',lo),('high',hi)]:emit(np.concatenate([choose(pool,5),fixed]),{**base,'scenario':'S4_intensity'},rep,st,pair)
  say('PSEUDOSPOTS',part,sample,len(rows))
 df=pd.DataFrame(rows);X=sparse.vstack(mat,format='csr');assert (np.asarray(X.sum(1)).ravel()==df.post_shared_UMI).all();wt(df,O/'manifest.tsv.gz');wt(pd.DataFrame(truth),O/'truth_counts.tsv.gz');wt(pd.DataFrame(members),O/'cell_membership_RNA.tsv.gz');wt(pd.DataFrame(contrib),O/'subtype_RNA_contribution.tsv.gz');wt(pd.DataFrame(avail),O/'scenario_availability.tsv')
 wt(pd.DataFrame(members).groupby('cell_id').size().rename('spot_reuse_count').reset_index(),O/'cell_reuse.tsv.gz')
 b=ad.AnnData(X=X,obs=df[['pseudospot_id','sample']].set_index('pseudospot_id'),var=a.var.copy());b.write_h5ad(O/'counts.h5ad',compression='lzf');js(O/'READY.json',{'n_spots':b.n_obs,'n_genes':b.n_vars,'counts_sha256':sha(O/'counts.h5ad'),'truth_sha256':sha(O/'truth_counts.tsv.gz'),'plan_hash':sha(R/'04_pseudospots/SIMULATION_FROZEN.json'),'member_source_check':'only scoped '+part+' cells; no cross sample; no repeated cell within spot'});say('PSEUDOSPOTS_READY',part,b.shape)
