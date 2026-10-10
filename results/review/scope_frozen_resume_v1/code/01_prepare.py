from common import *
import anndata as ad
from scipy import sparse
cfg={'version':'scope_frozen_resume_v1','authorization':'User explicitly excludes 18182 outside-frozen old Epithelial; retains every other old subtype member','master_seed':123,'labels':['Epi_K32','Epi_other'],'split_policy':'reuse original v5 unchanged; conservative sample groups, patients unverified','NB':{'initial_epochs':250,'extension_epochs':125,'max_epochs':1000,'batch_size':2500,'lr':.002,'train_size':1.,'convergence_relative_window_change':.01},'mapping':{'N_cells_per_location':20,'detection_alpha':20,'lr':.002,'batch_size':None,'initial_epochs':3000,'extension_epochs':2000,'max_epochs':30000,'convergence_relative_window_change':.005},'simulation':{'N0':20,'repeats':10,'epi_counts':[5,10,16],'depth':'median across 63 sections of p25/median/p75 UMI','cross_sample':False},'posterior_samples':1000,'ratio_min_epithelial_count':.1,'threshold_target_FPR':.05,'threshold_weighting':'maximum per-negative-stratum 95th percentile with equal source weights within stratum; strict > threshold','real_spatial_fit':False}
if (R/'CONFIG.yaml').exists():assert json.loads((R/'CONFIG.yaml').read_text())==cfg
else:js(R/'CONFIG.yaml',cfg)
(R/'SCOPE_REVISION.md').write_text('''# Authorized scope revision

The user's continuation explicitly replaces the incompatible old scope constraints. Keep all old non-Epithelial members and only the 123699 frozen epithelial-analysis members. Exclude 18182 outside-frozen old Epithelial cells from references and all main pseudospots, without assigning a negative state or QC-fail label. Epi_rest is renamed Epi_other, meaning other frozen epithelial clusters, not MP10-low. Preserve frozen Cycling-origin members by exact old-reference membership. The original BLOCKED audit remains unchanged and was not a biological failure. No additional scope confirmation is needed.
''')
old=pd.read_csv(V/'01_label_crosswalk/cell_crosswalk.tsv.gz',sep='\t',low_memory=False);split=pd.read_csv(V/'02_splits_and_genes/cell_split_manifest.tsv.gz',sep='\t');hist=pd.read_csv(V/'01_label_crosswalk/epithelial_historical_trace.tsv.gz',sep='\t',low_memory=False)
assert old.cell_id.is_unique and split.cell_id.is_unique
m=old.merge(split[['cell_id','group','split']],on='cell_id',validate='one_to_one');assert m[['group','split']].notna().all().all()
m['old_reference_label']=m.old_subtype;m['in_frozen_scope']=m.in_v4;m['frozen_cluster']=m.cluster;m['new_reference_label']=m.new_subtype.replace({'Epi_rest':'Epi_other'});m['included_in_this_run']=(m.old_subtype!='Epithelial')|m.in_v4;m['source_group']=m['group'];m['exclusion_reason']=np.where(m.included_in_this_run,'','outside_frozen_discovery_scope');hs=hist.set_index('cell_id').trace_status;m['exclusion_stage']=m.cell_id.map(hs).where(~m.included_in_this_run,'not_excluded');m['historical_specific_reason']=np.where(m.exclusion_stage=='OUTSIDE_CLEAN_OBJECT_ALREADY','UNVERIFIED',np.where(m.exclusion_stage=='REMOVED_AT_SMALL_CLUSTER_FILTER','documented_small_cluster_removal','not_applicable'))
wt(m,R/'01_scope/reference_scope_manifest.tsv.gz');wt(m[~m.included_in_this_run],R/'01_scope/scope_exclusions.tsv.gz');q=m[m.included_in_this_run].copy().reset_index(drop=True)
assert len(q)==281851 and q.new_reference_label.notna().all() and q.new_reference_label.nunique()==31 and q.K32.eq(True).sum()==3282 and q.new_reference_label.eq('Epi_other').sum()==120417
assert q.groupby('source_group').split.nunique().max()==1;assert q.loc[q['sample'].isin([P07,'CNP0000460|CNP0000460_P07N']),'split'].eq('TRAIN').all()
labels=sorted(q.new_reference_label.unique());js(R/'02_inputs/labels.json',labels)
wt(q,R/'01_scope/actual_cell_split_manifest.tsv.gz');wt(m.assign(new_reference_label=m.new_reference_label.fillna('EXCLUDED_SCOPE')).groupby(['old_reference_label','new_reference_label','included_in_this_run'],observed=True).size().rename('n').reset_index(),R/'01_scope/old_new_crosswalk.tsv')
wt(q.groupby(['split','dataset','sample','source_group','new_reference_label'],observed=True).size().rename('n_cells').reset_index(),R/'01_scope/source_support.tsv')
sp=pd.read_csv(V/'02_splits_and_genes/split_manifest.tsv',sep='\t');wt(sp,R/'01_scope/frozen_source_split_manifest.tsv');assert q[q.K32.eq(True)].groupby('split').size().to_dict()=={'DEV':634,'TEST':253,'TRAIN':2395}
raw=P/'subtype/v2/reference_training_v2/04_gene_harmonization/reference_v2_raw_counts_minimally_harmonized.h5ad';inp=pd.read_csv(V/'00_audit/INPUT_MANIFEST.tsv',sep='\t');expected=inp.set_index('path').loc[str(raw)];assert raw.stat().st_size==expected.bytes and raw.stat().st_mtime_ns==expected.mtime_ns
# Hash raw source once on resume; the previous integer audit is reused after exact identity verification.
say('VERIFY_RAW_HASH');assert sha(raw)==expected.sha256
frozen={'raw_counts_sha256':expected.sha256,'scope_member_label_hash':key(q.cell_id+'|'+q.new_reference_label),'split_hash':key(q.cell_id+'|'+q.split+'|'+q.source_group),'config_hash':sha(R/'CONFIG.yaml'),'original_split_sha256':sha(V/'02_splits_and_genes/cell_split_manifest.tsv.gz')}
if (R/'02_inputs/INPUTS_READY.json').exists():
 prior=json.loads((R/'02_inputs/INPUTS_READY.json').read_text());assert all(prior[k]==v for k,v in frozen.items());cache=json.loads((R/'02_inputs/COUNTS_CACHE_HASHES.json').read_text());assert all(sha(R/f'02_inputs/{part}_counts.h5ad')==value for part,value in cache.items());say('VALIDATED_INPUT_CACHE');raise SystemExit
say('READ_SPARSE_RAW');a=ad.read_h5ad(raw,backed='r');assert a.obs_names.equals(pd.Index(old.cell_id));prev=pd.read_csv(V/'02_splits_and_genes/gene_manifest.tsv',sep='\t',low_memory=False).set_index('gene').loc[a.var_names]
train=m.split.eq('TRAIN').to_numpy()&m.included_in_this_run.to_numpy();no=train&m['sample'].ne(P07).to_numpy();tot=np.zeros(a.n_vars,np.int64);nnz=np.zeros(a.n_vars,np.int64);ntot=np.zeros(a.n_vars,np.int64)
for start in range(0,a.n_obs,4000):
 end=min(a.n_obs,start+4000);X=a.X[start:end].tocsr();X.eliminate_zeros();z=X[train[start:end]];tot+=np.asarray(z.sum(0)).ravel().astype(np.int64);nnz+=z.getnnz(0);ntot+=np.asarray(X[no[start:end]].sum(0)).ravel().astype(np.int64)
mean=np.divide(tot,nnz,out=np.zeros(a.n_vars),where=nnz>0);passed=(nnz>train.sum()*.03)|((nnz>5)&(mean>1.12));keep=passed&prev.resolved_identity.to_numpy()&prev.shared_all_eligible_sections.to_numpy()&~prev.is_MT.to_numpy()&(ntot>0)
g=prev.reset_index(names='gene');g['TRAIN_UMI']=tot;g['TRAIN_n_cells_nonzero']=nnz;g['TRAIN_nonz_mean']=mean;g['TRAIN_NO_P07T_UMI']=ntot;g['passes_TRAIN_filter']=passed;g['selected']=keep;g['model_gene_status']='FROZEN_SCOPED_MAIN_TRAIN';g=g.drop(columns=['provisional_common_gene_keep']);g=annotate_MP_sources(g);wt(g,R/'02_inputs/gene_manifest.tsv');genes=a.var_names[keep].tolist();wt(pd.DataFrame({'gene':genes}),R/'02_inputs/genes.tsv');frozen['gene_order_hash']=key(genes);frozen['label_order_hash']=key(labels)
for part in ['TRAIN','DEV','TEST']:
 ix=np.flatnonzero((m.split==part)&m.included_in_this_run);blocks=[]
 for start in range(0,len(ix),4000):blocks.append(a.X[ix[start:start+4000]][:,keep].tocsr())
 X=sparse.vstack(blocks,format='csr');obs=m.iloc[ix][['cell_id','dataset','sample','source_group','split','new_reference_label']].set_index('cell_id').copy()
 b=ad.AnnData(X=X,obs=obs,var=pd.DataFrame(index=pd.Index(genes,name='gene')));b.uns['scope_input_hashes']=frozen;b.write_h5ad(R/f'02_inputs/{part}_counts.h5ad',compression='lzf');say('WRITE',part,b.shape)
a.file.close()
for branch in ['MAIN','NO_P07T']:
 z=q[q.split=='TRAIN'];z=z if branch=='MAIN' else z[z['sample']!=P07];support=z.groupby('new_reference_label',observed=True).size().reindex(labels,fill_value=0);wt(support.rename('n_cells').reset_index(),R/f'03_reference/{branch}/training_support.tsv');assert (support>0).all();wt(z[['cell_id','sample','source_group','new_reference_label']],R/f'03_reference/{branch}/training_members.tsv.gz')
js(R/'02_inputs/INPUTS_READY.json',{**frozen,'n_candidate_cells':len(q),'n_genes':len(genes),'n_TRAIN':int(train.sum()),'n_NO_P07T':int(no.sum())});say('INPUTS_READY',len(genes),int(train.sum()),int(no.sum()))

js(R/'02_inputs/COUNTS_CACHE_HASHES.json',{p:sha(R/f'02_inputs/{p}_counts.h5ad') for p in ['TRAIN','DEV','TEST']})
