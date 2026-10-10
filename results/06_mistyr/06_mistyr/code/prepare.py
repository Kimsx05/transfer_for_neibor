import bootstrap
from bootstrap import ROOT, RUN
import hashlib, json, os, shutil, platform
from pathlib import Path
import pandas as pd
import numpy as np
import h5py
from scipy import sparse

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    return h.hexdigest()
def save(d,p): d.to_csv(p,sep='\t',index=False,na_rep='NA')
def dump(d,p): Path(p).write_text(json.dumps(d,indent=2,ensure_ascii=False)+'\n')

contract=json.loads((RUN/'00_inputs/project_contract.json').read_text())
h0=json.loads((RUN/'00_inputs/handoff.json').read_text())
h1=json.loads((RUN/'01_mp_display/handoff.json').read_text())
assert h0['status'].startswith('PASS') and h1['status'].startswith('PASS')
mem=int(next(l.split()[1] for l in Path('/proc/meminfo').read_text().splitlines() if l.startswith('MemAvailable:')))*1024
disk=shutil.disk_usage(ROOT).free
dump(dict(available_memory_bytes=mem,available_disk_bytes=disk,minimum_memory_bytes=4*1024**3,minimum_disk_bytes=5*1024**3,cpu_budget=8,jobs_parallel=1,blas_threads=1,gpu=False),ROOT/'resource_preflight.json')
assert mem>4*1024**3 and disk>5*1024**3
manifest=pd.read_csv(contract['eligible_manifest'],sep='\t')
c2l=pd.read_csv(contract['c2l_manifest'],sep='\t').set_index('section_id')
mapping=pd.read_csv(contract['major_map'],sep='\t')
assert not mapping.subtype.duplicated().any() and set(mapping.subtype)==set(contract['factor_order'])
spots=pd.read_parquet(contract['frozen_spot_table'])
scores=pd.read_parquet(h1['spot_table'],filters=[('MP','==','EpiMP10')])
assert not spots.duplicated(['section_id','barcode']).any()
assert not scores.duplicated(['section_id','barcode']).any()
merged=spots.merge(scores[['section_id','barcode','mp_score_raw','analysis_eligible_legacy']],on=['section_id','barcode'],suffixes=('','_01'),validate='one_to_one')
assert len(merged)==len(spots)
assert np.array_equal(merged.mp_score_raw,merged.EpiMP10_Top20_AUCell,equal_nan=True)
assert (merged.analysis_eligible_legacy==merged.analysis_eligible_legacy_01).all()
majors=sorted(mapping.major.unique()); factors=contract['factor_order']
dictionary=[]
for i,f in enumerate(factors): dictionary.append(dict(variable=f'F{i+1:02}',source_label=f,role='B_TME_predictor' if f!='Epithelial' else 'excluded_epithelial_factor',units='posterior_mean_expected_cells_per_spot',transform='none'))
for i,m in enumerate(majors): dictionary.append(dict(variable=f'M{i+1:02}',source_label=m,role='A_target_and_other_major_predictor',units='posterior_mean_expected_cells_per_spot',transform='sum approved subtypes; no log or closure'))
dictionary += [dict(variable='MP10',source_label='EpiMP10 raw Top20 AUCell',role='B_target_only',units='raw_AUCell',transform='none'),dict(variable='epi_fraction_gate',source_label='historical c2l mean fraction',role='B_technical_predictor',units='fraction',transform='none'),dict(variable='log1p_nCount',source_label='raw X row sum over all frozen features',role='B_technical_predictor',units='log1p_UMI_count',transform='log1p; fixed, not fitted')]
save(pd.DataFrame(dictionary),ROOT/'variable_dictionary.tsv')
sources=[];audit=[]; summaries=[]
def register(p,role):
    p=Path(p); sources.append(dict(path=str(p.resolve()),role=role,bytes=p.stat().st_size,sha256=sha(p)))
for p in [RUN/'00_inputs/project_contract.json',RUN/'00_inputs/handoff.json',RUN/'01_mp_display/handoff.json',Path(contract['eligible_manifest']),Path(contract['c2l_manifest']),Path(contract['major_map']),Path(contract['frozen_spot_table']),Path(h1['spot_table'])]: register(p,'frozen_upstream')
# Verify hashes of directly consumed frozen tables against their handoffs.
for hand in [h0,h1]:
    expected={o['path']:o['sha256'] for o in hand['outputs']}
    for s in sources:
        if s['path'] in expected:
            ok=s['sha256']==expected[s['path']]
            audit.append(dict(path=s['path'],expected_sha256=expected[s['path']],observed_sha256=s['sha256'],status='PASS' if ok else 'BLOCKED'))
assert all(x['status']=='PASS' for x in audit)
save(pd.DataFrame(audit),ROOT/'upstream_hash_validation.tsv')
for row in manifest.itertuples():
    sec=row.section_id; folder=ROOT/'by_dataset'/row.dataset/sec;folder.mkdir(parents=True,exist_ok=True)
    d=merged[merged.section_id==sec].copy()
    ab=pd.read_csv(c2l.loc[sec,'mean_path'],index_col=0)
    assert ab.index.is_unique and set(d.barcode)<=set(ab.index)
    assert list(ab.columns)==factors
    ab=ab.loc[d.barcode]; assert np.isfinite(ab.values).all() and (ab.values>=0).all()
    base=d[['section_id','barcode','dataset','cohort','patient_uid','specimen','array_row','array_col','x','y','epi_fraction_gate','analysis_eligible_legacy','mp_score_raw']].rename(columns={'mp_score_raw':'MP10'}).reset_index(drop=True)
    for i,f in enumerate(factors): base[f'F{i+1:02}']=ab[f].to_numpy()
    for i,m in enumerate(majors): base[f'M{i+1:02}']=ab[mapping.loc[mapping.major==m,'subtype']].sum(axis=1).to_numpy()
    # Only read counts; never load or rewrite Seurat/AnnData objects.
    with h5py.File(row.raw_input_path,'r') as f:
        ix=f['obs'][f['obs'].attrs['_index']]
        if isinstance(ix,h5py.Dataset): idx=ix.asstr()[:]
        elif 'values' in ix:
            assert not ix['mask'][:].any()
            idx=ix['values'].asstr()[:]
        else: idx=ix['categories'].asstr()[:][ix['codes'][:]]
        x=f['X']; enc=x.attrs['encoding-type']; shape=tuple(x.attrs['shape'])
        val=x['data'][:]; assert np.isfinite(val).all() and (val>=0).all()
        assert np.allclose(val,np.rint(val),rtol=0,atol=1e-6)
        klass={'csr_matrix':sparse.csr_matrix,'csc_matrix':sparse.csc_matrix}[enc]
        mat=klass((val,x['indices'][:],x['indptr'][:]),shape=shape)
        nc=pd.Series(np.asarray(mat.sum(axis=1)).ravel(),index=idx)
        assert nc.index.is_unique and set(base.barcode)<=set(nc.index)
        base['nCount']=nc.loc[base.barcode].to_numpy()
        del val,mat
    base['log1p_nCount']=np.log1p(base.nCount)
    # Deterministic geometry-only five contiguous bands along widest true-lattice axis.
    xr=base.array_col.to_numpy()/2; yr=base.array_row.to_numpy()*np.sqrt(3)/2
    axis=xr if np.ptp(xr)>=np.ptp(yr) else yr
    fold_scope=base.analysis_eligible_legacy & base[['MP10','epi_fraction_gate','log1p_nCount']].notna().all(axis=1)
    bins=np.unique(np.quantile(axis[fold_scope],[.2,.4,.6,.8])) if fold_scope.any() else np.array([])
    base['spatial_fold']=np.searchsorted(bins,axis,side='right')+1
    save(base,folder/'model_input.tsv')
    register(c2l.loc[sec,'mean_path'],'c2l_posterior_mean_export')
    # Hash actual counts input bytes for reproducible nCount, even though only row sums enter models.
    register(row.raw_input_path,'raw_counts_for_nCount')
    register(folder/'model_input.tsv','derived_model_input')
    n=int(base.analysis_eligible_legacy.sum())
    summaries.append(dict(section_id=sec,dataset=row.dataset,cohort=row.cohort,patient_uid=row.patient_uid,specimen=row.specimen,n_tissue=len(base),n_legacy_eligible=n,n_mp10_available=int(base.MP10.notna().sum()),n_complete_B=int((base.analysis_eligible_legacy & base[['MP10','epi_fraction_gate','log1p_nCount']].notna().all(axis=1)).sum()),major_count=len(majors),TME_factor_count=len(factors)-1,fold_axis='array_col/2' if np.ptp(xr)>=np.ptp(yr) else 'sqrt(3)*array_row/2',fold_boundaries=json.dumps(bins.tolist()),status='PASS'))
    print('PREPARED',sec,len(base),n,flush=True)
save(pd.DataFrame(summaries),ROOT/'section_input_summary.tsv')
save(pd.DataFrame(sources),ROOT/'input_manifest.tsv')
params=dict(module='06',pilot_section='GSE319536_BC15',cpu_budget=8,jobs_parallel=1,future_plan='sequential',ranger_threads=1,blas_threads=1,gpu=False,seed=6102026,num_trees=100,importance='impurity',min_node_size=5,mtry='floor(sqrt(p))',internal_meta_folds=10,spatial_folds=5,spatial_design='eligible complete-case geometry-only contiguous quantile bands along widest tissue lattice axis; no buffer',min_spots=100,min_train_spots=50,min_test_spots=20,major_variables=[f'M{i+1:02}' for i in range(len(majors))],tme_variables=[f'F{i+1:02}' for i,f in enumerate(factors) if f!='Epithelial'],technical_variables=['epi_fraction_gate','log1p_nCount'],target='MP10',transform_A='raw mean abundance',transform_B='raw mean abundance; frozen raw AUCell target; fixed log1p nCount only',importance_interpretability_A='unclipped internal meta R2 > 0 and OOB R2 > 0; internal only',importance_interpretability_B='spatial pooled full R2 > 0 AND full R2 - baseline R2 > 0',aggregation='equal-section median of per-target standardized importance; no patient inference',allowed_scope_override='Current user explicitly requests Module 06; 00 contract allowed_modules records previous session scope only')
dump(params,ROOT/'parameters.json')
dump(dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,h5py=h5py.__version__,base_environment_modified=False),ROOT/'python_versions.json')
