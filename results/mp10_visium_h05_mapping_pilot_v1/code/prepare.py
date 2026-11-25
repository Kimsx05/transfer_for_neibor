from common import *
from scipy import sparse, io
import h5py, json
T=Path(os.environ['PILOT_TMP'])
seed=20261125
rng=np.random.default_rng(seed)
inputs=[]
def record(role,p,section='reference'):
 p=Path(p)
 inputs.append(dict(role=role,section=section,path=str(p),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns))
def strings(g):
 if isinstance(g,h5py.Dataset): return g.asstr()[:]
 if 'values' in g:
  assert not np.any(g['mask'][:])
  return g['values'].asstr()[:]
 cats=g['categories'].asstr()[:]; codes=g['codes'][:]; assert min(codes)>=0
 return cats[codes]
def pack(name,m,genes,ids):
 # CSR cells x genes is CSC genes x cells for R
 m=m.tocsr();m.sum_duplicates();m.eliminate_zeros();m.sort_indices()
 assert len(ids)==m.shape[0] and len(genes)==m.shape[1] and pd.Index(ids).is_unique and pd.Index(genes).is_unique
 assert (m.data>=0).all() and np.array_equal(m.data,m.data.astype(np.int32))
 p=T/name;p.mkdir(exist_ok=True)
 for k,x in [('i',m.indices),('p',m.indptr),('x',m.data)]:x.astype('<i4').tofile(p/k)
 for k,x in [('genes',genes),('ids',ids)]: (p/k).write_text('\n'.join(x)+'\n')
def slice_h5(f,idx):
 # Read small consecutive blocks, in physical order.
 x=f['X'];ptr=x['indptr'][:]; idx=np.asarray(idx)
 out=[];groups=np.split(np.sort(idx),np.where(np.diff(np.sort(idx))>24)[0]+1)
 for g in groups:
  a,b=ptr[g[0]],ptr[g[-1]+1]
  block=sparse.csr_matrix((x['data'][a:b],x['indices'][a:b],ptr[g[0]:g[-1]+2]-a),shape=(g[-1]-g[0]+1,x.attrs['shape'][1]))
  out.append(block[g-g[0]])
 return sparse.vstack(out,format='csr')[np.argsort(np.argsort(idx))]
a=read(A/'analysis_cell_metadata.tsv.gz')
members=read(A/'pseudobulk_membership.tsv.gz')
q=read(A/'specimen_exposure_eligibility.tsv')
for p in ['analysis_cell_metadata.tsv.gz','pseudobulk_membership.tsv.gz','specimen_exposure_eligibility.tsv']:record('frozen_cells',A/p)
record('harmonized_raw_counts',REF)
for c in ['barcode','sample_uid']: assert a[c].notna().all()
assert a.barcode.is_unique and q.sample_uid.is_unique and len(q)==40
epi=a.analysis_major_c2l_final.eq('Epithelial')
keys=lambda d:set(zip(d.sample_uid,d.barcode))
pool_path=R/'mp10_cellchat_high_definition_sensitivity_v1/tables/epi_full_scores_pools.tsv.gz'
record('authoritative_original_full_H05_flags',pool_path)
pool=read(pool_path)
assert keys(pool)==keys(a[epi]) and pool.barcode.is_unique
high=pool[pool.H05]
later_high=members[members.profile_type.eq('H05')]
boundary=a.loc[a.barcode.isin(set(later_high.barcode)^set(high.barcode)),['sample_uid','barcode','MP10_cell']].copy()
boundary['original_H05']=boundary.barcode.isin(high.barcode)
boundary['later_pseudobulk_H05']=boundary.barcode.isin(later_high.barcode)
save(boundary,'membership_version_difference.tsv')
assert keys(a[epi])==keys(members[members.profile_type.eq('Epi_all')])
# Serialized quantiles can round at the last digit. Membership is authoritative.
delta=a.MP10_cell-a.sample_uid.map(q.set_index('sample_uid').Q95)
mismatch=epi & ((delta>0)!=a.barcode.isin(high.barcode))
save(a.loc[mismatch,['barcode','sample_uid','MP10_cell']].assign(score_minus_serialized_Q95=delta[mismatch],frozen_high=a.loc[mismatch,'barcode'].isin(high.barcode)),'serialized_boundary_audit.tsv')
assert (delta[mismatch].abs() < 1e-12).all()
assert np.array_equal(high.groupby('sample_uid').size().reindex(q.sample_uid,fill_value=0),q.n_H05)
assert keys(high).issubset(keys(a[epi]))
a['label']=a.analysis_major_c2l_final
a.loc[epi,'label']='Epi_nonH05'; a.loc[a.barcode.isin(high.barcode),'label']='Epi_H05'
sp=read(K/'01_scope/frozen_source_split_manifest.tsv')
record('source_groups',K/'01_scope/frozen_source_split_manifest.tsv')
sp['sample_uid']=sp['sample'].str.replace('|','::',regex=False)
a=a.merge(sp[['sample_uid','group']],on='sample_uid',how='left',validate='many_to_one')
held=['CNP0000460::P01','InhouseData::T4463209']
hold=a['group'].isin(held)
assert set(a.loc[hold,'sample_uid'])=={'CNP0000460::CNP0000460_P01T','InhouseData::T4463209a','InhouseData::T4463209b'}
assert not set(a.loc[hold,'patient_uid_recorded']) & set(a.loc[~hold,'patient_uid_recorded'])
sample=a.loc[np.concatenate([rng.choice(g.index,min(len(g),35),replace=False) for _,g in a[~hold].groupby(['sample_uid','label'])])].copy()
sample=sample.sort_values(['sample_uid','label','barcode']).reset_index(drop=True)
with h5py.File(REF) as f:
 ids=pd.Index(strings(f['obs/_index'])); idx=ids.get_indexer(sample.barcode);assert (idx>=0).all() and ids.is_unique
 src=np.char.add(np.char.add(strings(f['obs/dataset'])[idx].astype(str),'::'),strings(f['obs/orig.ident'])[idx].astype(str))
 assert np.array_equal(src,sample.sample_uid)
 pack('reference',slice_h5(f,idx),strings(f['var/_index']),sample.barcode)
save(sample,'reference_sampling_ids.tsv.gz')
sizes=a.groupby(['sample_uid','label']).size().rename('eligible').to_frame().join(sample.groupby(['sample_uid','label']).size().rename('sampled')).fillna(0).reset_index()
sizes['sampled_share_reference']=sizes.sampled/len(sample)
sizes['within_specimen_share']=sizes.sampled/sizes.groupby('sample_uid').sampled.transform('sum').replace(0,np.nan)
save(sizes,'reference_sampling_proportions.tsv')
save(q[['sample_uid','Q95','n_epi','n_H05']],'frozen_specimen_thresholds.tsv')
record('spatial_manifest',F/'tables/input_manifest.tsv','four')
record('historical_query_counts',F/'full_sections.tsv','four')
join=[]
for sec,r in M.iterrows():
 p=R/'mp10_paper_distance_full_v1/spots'/f'{sec}.parquet'
 d=pd.read_parquet(p)
 for k in ['counts_path','genes_path','barcodes_path','image_path']:record(k,r[k],sec)
 record('full_v1_score_coordinates_and_c2l',p,sec)
 ids=pd.Index(Path(r.barcodes_path).read_text().splitlines())
 genes=pd.Index(Path(r.genes_path).read_text().splitlines())
 assert set(d.section_id)=={sec} and d.barcode.is_unique and ids.is_unique and set(ids)==set(d.barcode)
 d=d.set_index('barcode').loc[ids].rename_axis('barcode').reset_index()
 x=io.mmread(r.counts_path).tocsr()
 assert x.shape==(len(genes),len(ids))
 if not genes.is_unique:
  ix,genes=pd.factorize(genes)
  x=sparse.csr_matrix((np.ones(len(ix),dtype=int),(ix,np.arange(len(ix)))),shape=(len(genes),len(ix)))@x
 pack(sec,x.T,genes,ids)
 totals=np.asarray(x.sum(axis=0)).ravel()
 cols=['section_id','barcode','query','x','y','x_um','y_um','array_row','array_col','tissue_component_id','lesion_location','signed_distance_um','center_qc','analysis_eligible_legacy','EpiMP10_Top20_AUCell','Epithelial_fraction','nCount_RNA','epi_MP10_model_universe','spot_MP10_model_universe','epi_score_valid','source_context']
 b=d[cols].copy();b['raw_UMI']=totals
 broad=[k+'_fraction' for k in ['Epithelial','B','Endothelial','Fibroblast','Mast','Mural','Myeloid','Plasma','T_NK']]
 b['dominant_c2l_major']=d[broad].idxmax(axis=1).str.replace('_fraction','',regex=False)
 assert int(d['query'].sum())==S[sec] and d[['x','y']].notna().all().all()
 for gene in TOP+'EPCAM KRT8 KRT18 KRT19 PTPRC COL1A1 PECAM1'.split():
  b['UMI_'+gene]=x[genes.get_loc(gene)].toarray().ravel() if gene in genes else np.nan
 save(b,f'spots/{sec}.input.tsv.gz')
 join.append(dict(section_id=sec,n_expression=len(ids),n_frozen=len(d),query_observed=d['query'].sum(),query_expected=S[sec],max_abs_depth_difference=np.max(np.abs(totals-d.nCount_RNA)),unmatched=0))
save(join,'barcode_audit.tsv')
save(inputs,'input_manifest.tsv')
# Old TEST selects two whole source groups, up to 80 medium mixtures/specimen, plus their paired low realizations.
def prepare_pseudo():
 base=K/'04_pseudospots/TEST'
 for name in ['counts.h5ad','manifest.tsv.gz','cell_membership_RNA.tsv.gz']:record('reused_TEST',base/name,'pseudo')
 m=read(base/'manifest.tsv.gz'); mem=read(base/'cell_membership_RNA.tsv.gz')
 m=m[m.source_group.isin(held)&m.depth_condition.isin(['medium','low'])].copy()
 m['sample_uid']=m['sample'].str.replace('|','::',regex=False)
 m=m[m.sample_uid.isin(a.loc[hold,'sample_uid'])]
 mem=mem[mem.pseudospot_id.isin(m.pseudospot_id)].merge(a[['barcode','sample_uid','label']],left_on='cell_id',right_on='barcode',how='left',validate='many_to_one')
 mem=mem.merge(m[['pseudospot_id','sample_uid']],on='pseudospot_id',suffixes=('','_target'),validate='many_to_one')
 mem['known']=mem.sample_uid.eq(mem.sample_uid_target);mem['hit']=mem.label.eq('Epi_H05')
 g=mem.groupby('pseudospot_id');m=m.join(g.known.all().rename('known'),on='pseudospot_id').join(g.hit.sum().rename('n_H05'),on='pseudospot_id').join(g.size().rename('n_members'),on='pseudospot_id')
 m=m[m.known.eq(True)&m.n_members.eq(m.N_total)]
 pairs=m.groupby('mixture_id').filter(lambda g:len(g)==2 and g.n_H05.nunique()==1)
 medium=pairs[pairs.depth_condition=='medium'];rand=np.random.default_rng(seed+1)
 chosen=np.concatenate([rand.choice(g.mixture_id,min(80,len(g)),replace=False) for _,g in medium.groupby('sample_uid')])
 d=pairs[pairs.mixture_id.isin(chosen)].sort_values(['sample_uid','mixture_id','depth_condition']).copy()
 assert len(d)==2*len(chosen) and len(d)>0
 d['H05_fraction']=d.n_H05/d.N_total
 # Validate exact cell multisets across depth realizations.
 memberlists=mem.groupby('pseudospot_id').cell_id.agg(lambda x:tuple(sorted(x)))
 assert (d.assign(members=d.pseudospot_id.map(memberlists)).groupby('mixture_id').members.nunique()==1).all()
 with h5py.File(base/'counts.h5ad') as f:
  idx=pd.Index(strings(f['obs/pseudospot_id'])).get_indexer(d.pseudospot_id);assert min(idx)>=0
  for depth,grp in d.groupby('depth_condition'):
   pos=pd.Index(d.pseudospot_id).get_indexer(grp.pseudospot_id)
   pack('pseudo_'+depth,slice_h5(f,idx[pos]),strings(f['var/gene']),grp.pseudospot_id)
 save(d[['pseudospot_id','mixture_id','sample_uid','source_group','depth_condition','scenario','post_shared_UMI','N_total','n_H05','H05_fraction']],'pseudo_truth.tsv')
 save(mem.loc[mem.pseudospot_id.isin(d.pseudospot_id),['pseudospot_id','cell_id','sample_uid','label']],'pseudo_members.tsv.gz')
try: prepare_pseudo()
except Exception:
 import traceback
 (O/'logs/pseudo_prepare_failure.txt').write_text(traceback.format_exc())
save(inputs,'input_manifest.tsv')
print('prepared',len(sample),join,flush=True)

