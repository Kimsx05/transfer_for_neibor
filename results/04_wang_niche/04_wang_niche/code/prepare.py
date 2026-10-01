import sys, json, hashlib, shutil, time
from pathlib import Path
ROOT = Path('/home/data/t070721/codex_workspace/Bladder Metabolism/spatial_ecology_top20_v2')
sys.path.insert(0, str(ROOT/'code'))
from common import *
OUT = ROOT/'04_wang_niche'
SEED = 20260930
ATTACHMENT = Path('/home/data/t070721/.codex/attachments/3c23f9ae-e6fe-4f14-8ba3-ad7b8cafbeec/04_wang_composition_niche.md')

def prepare():
    OUT.mkdir(exist_ok=True)
    contract_path = ROOT/'00_inputs/project_contract.json'
    handoff_path = ROOT/'00_inputs/handoff.json'
    contract = json.loads(contract_path.read_text())
    upstream = json.loads(handoff_path.read_text())
    assert upstream['status'] in ('PASS', 'PASS_WITH_WARNING')
    assert sha(ROOT/'00_inputs/input_provenance.tsv') == upstream['input_signature']
    assert not contract['gpu_allowed'] and contract['jobs_parallel'] == 1
    mem = dict(line.split(':',1) for line in Path('/proc/meminfo').read_text().splitlines())
    mem_gb = int(mem['MemAvailable'].split()[0])/1024**2
    disk_gb = shutil.disk_usage(OUT).free/1024**3
    assert mem_gb >= 32 and disk_gb >= 2, 'Insufficient declared resource budget'
    savejson(OUT/'resource_preflight.json', dict(available_memory_GiB=mem_gb, available_disk_GiB=disk_gb,
        required_memory_GiB=32, required_disk_GiB=2, cpu_budget_total=8, simultaneous_datasets=1,
        candidate_workers=4,BLAS_OpenMP_threads=1,graph_threads=1,GPU=False,
        checked_UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
    factors = contract['factor_order']
    assert len(factors) == 30 and len(set(factors)) == 30 and 'Epithelial' in factors
    sm = readcsv(contract['eligible_manifest'],sep='\t')
    cm = readcsv(contract['c2l_manifest'],sep='\t').set_index('section_id')
    provenance_df = readcsv(ROOT/'00_inputs/input_provenance.tsv',sep='\t')
    # Explicit column projection prevents MP/gate/CNV values from entering this module.
    metadata_cols = ['section_id','barcode','dataset','cohort','patient_uid','specimen','x','y']
    base = pd.read_parquet(contract['frozen_spot_table'],columns=metadata_cols)
    assert not base.duplicated(['section_id','barcode']).any()
    assert len(base)==int(sm.n_spots.sum())
    ident_path = P/'cell2location_final_validation_gate_v1/09_final_status/DOWNSTREAM_USAGE_TABLE.csv'
    ident = readcsv(ident_path)
    assert set(ident.factor)==set(factors)
    tsv(OUT/'factor_interpretation.tsv',ident)
    parameter = dict(module='04', authorization='User explicitly requested Module 04 with attached task specification',
        authorization_document=str(ATTACHMENT), authorization_document_sha256=sha(ATTACHMENT),
        upstream_contract=str(contract_path), upstream_contract_sha256=sha(contract_path),
        upstream_handoff=str(handoff_path), upstream_handoff_sha256=sha(handoff_path),
        upstream_input_signature=upstream['input_signature'],
        prior_contract_allowed_modules='Historical 00/01 scope retained unchanged; current user request additionally authorizes 04',
        factor_order=factors, posterior_summary='original accepted decimal CSV export of posterior mean abundance',
        abundance_units='expected_cells_per_spot', transform='row sum normalization over all 30 factors',
        invalid_row_rule='any nonfinite or negative abundance, nonfinite/zero total -> NA proportions; retain full spot record',
        features='all 30 subtype proportions ONLY', integration=None, log=None, zscore=None, ILR=None, PCA=None,
        pooling='all accepted sections within each frozen dataset; no cross-dataset label pooling',
        graph='scran::buildSNNGraph', input_orientation='spots x factors; transposed=TRUE', d='NA_integer_',
        SNN_weighting='rank', neighbor_search='BiocNeighbors::KmknnParam(distance="Euclidean"); exact',
        k=[10,20,30], resolutions=[round(v/10,1) for v in range(1,16)],
        clustering='igraph::cluster_louvain; undirected, finite strictly positive weights',
        seed=SEED, silhouette_distance='Euclidean in unchanged 30-factor proportions',
        silhouette_evaluation='all valid spots if N<=5000; otherwise fixed proportional section-stratified 5000 spots',
        silhouette_selection='max evaluation-SPOT mean; invalid/all-one-class/all-singleton candidates excluded',
        tie_tolerance=1e-8, tie_break_order=['lower resolution','lower k'],
        BH_family='all niche x factor two-sided tests within each dataset',
        wilcoxon='Mann-Whitney rank sum, asymptotic with ties and continuity correction; no paired samples',
        FDR_cutoff=0.05, annotation_role='defining composition description, not independent validation',
        niche_numbering='dataset-scoped; decreasing full cluster size, then smallest composite spot key',
        label_policy='top three mean subtype fractions with numbers; no forced biological names or MP labels',
        low_identifiability='retain all factors as features; flag original status in interpretation',
        centroid_redundancy_diagnostic_threshold=0.03, diagnostic_only=True, cpu_budget_total=8,
        jobs_parallel=1,candidate_workers=4,BLAS_OpenMP_threads=1,graph_threads=1,gpu_allowed=False,
        pdf=True, png_dpi=600, no_MP_niche_tests=True, stop_after_module='04')
    savejson(OUT/'parameters.json',parameter)
    savejson(OUT/'software_versions_python.json',versions())
    code_files = sorted((OUT/'code').glob('*.*')) + [ROOT/'code/common.py']
    savejson(OUT/'code_checksums.json',{str(p):sha(p) for p in code_files if p.is_file()})
    savejson(OUT/'upstream_freeze_checksums.json',{str(p):sha(p) for p in [contract_path,handoff_path,
        ROOT/'00_inputs/input_provenance.tsv',Path(contract['eligible_manifest']),Path(contract['c2l_manifest']),
        Path(contract['major_map']),Path(contract['histology_manifest']),Path(contract['frozen_spot_table'])]})
    rows=[]; signature_inputs=[]; dataset_dirs=[]
    # Smallest dataset is an actual four-section complete mechanical smoke.
    order = sm.groupby('dataset').n_spots.sum().sort_values().index.tolist()
    for ds in order:
        folder=OUT/'by_dataset'/ds;folder.mkdir(parents=True,exist_ok=True)
        metadata=[]; compositions=[]
        for _,sec in sm.loc[sm.dataset.eq(ds)].iterrows():
            c=cm.loc[sec.section_id]; src=Path(c.mean_path)
            expected=provenance_df.loc[provenance_df.path.eq(str(src)),'sha256'].dropna().unique()
            assert len(expected)==1, f'No unique upstream mean checksum: {src}'
            checksum=sha(src);assert checksum==expected[0],f'Mean source changed: {src}'
            d=readcsv(src);assert d.columns[1:].tolist()==factors
            assert json.loads(c.factor_names)==factors
            d=d.rename(columns={d.columns[0]:'barcode'})
            assert not d.barcode.duplicated().any()
            m=base.loc[base.section_id.eq(sec.section_id)].copy().reset_index(drop=True)
            assert len(m)==int(sec.n_spots) and set(m.barcode)==set(d.barcode)
            assert m.dataset.eq(ds).all()
            a=d.set_index('barcode').loc[m.barcode,factors].to_numpy(dtype=np.float64)
            total=a.sum(axis=1)
            valid=np.isfinite(a).all(axis=1)&(a>=0).all(axis=1)&np.isfinite(total)&(total>0)
            reasons=np.select([~np.isfinite(a).all(axis=1),(a<0).any(axis=1),~np.isfinite(total),total<=0],
                ['NONFINITE_ABUNDANCE','NEGATIVE_ABUNDANCE','NONFINITE_TOTAL','ZERO_TOTAL'],default='VALID')
            prop=np.full_like(a,np.nan);prop[valid]=a[valid]/total[valid,None]
            assert np.allclose(prop[valid].sum(axis=1),1,rtol=0,atol=1e-12)
            m['spot_key']=m.section_id+'::'+m.barcode
            m['composition_valid']=valid;m['composition_reason']=reasons;m['total_mean_abundance']=total
            metadata.append(m);compositions.append(pd.DataFrame(prop,columns=factors))
            r=dict(section_id=sec.section_id,dataset=ds,cohort=sec.cohort,patient_uid=sec.patient_uid,
                specimen=sec.specimen,n_frozen_spots=len(m),n_valid_compositions=int(valid.sum()),
                n_invalid_compositions=int((~valid).sum()),n_factors=len(factors),mean_source=str(src),
                mean_source_sha256=checksum,posterior_mean_key=c.mean_key,
                transformation='posterior mean / total of all 30 factors; no other transform',
                abundance_units='expected_cells_per_spot',composition_units='fraction; row sum=1',
                maximum_row_sum_error=float(np.max(np.abs(prop[valid].sum(axis=1)-1))) if valid.any() else np.nan,
                discovery_gate='all valid compositions; no epithelial or CNV gate',status='PASS' if valid.all() else 'PASS_WITH_WARNING')
            rows.append(r);signature_inputs.append(dict(path=str(src),sha256=checksum))
        m=pd.concat(metadata,ignore_index=True);c=pd.concat(compositions,ignore_index=True)
        full=pd.concat([m,c],axis=1)
        full.to_parquet(folder/'spot_composition_inputs.parquet',index=False,compression='zstd')
        valid=full.loc[full.composition_valid].reset_index(drop=True)
        assert valid.spot_key.is_unique
        valid[['spot_key']+factors].to_csv(folder/'composition_features.csv.gz',index=False,
            float_format='%.17g',compression={'method':'gzip','mtime':0})
        n=len(valid);target=min(5000,n);rng=np.random.default_rng(SEED)
        sizes=valid.groupby('section_id',sort=True).size()
        if n<=target: selected=np.arange(n)
        else:
            quota=sizes.to_numpy()*target/n
            allocation=np.maximum(1,np.floor(quota).astype(int))
            while allocation.sum()<target:
                j=int(np.argmax(quota-allocation));allocation[j]+=1
            while allocation.sum()>target:
                candidates=np.where(allocation>1,allocation-quota,-np.inf)
                allocation[int(np.argmax(candidates))]-=1
            picks=[]
            for sec,n_take in zip(sizes.index,allocation):
                ids=np.flatnonzero(valid.section_id.eq(sec).values)
                picks.extend(rng.choice(ids,size=int(n_take),replace=False).tolist())
            selected=np.sort(np.array(picks,dtype=int))
        evaluation=valid.loc[selected,['spot_key','section_id','barcode']].copy()
        evaluation['input_row_1based']=selected+1
        tsv(folder/'silhouette_evaluation_spots.tsv',evaluation)
        alloc=evaluation.groupby('section_id').size().rename('n_evaluation').to_frame().join(sizes.rename('n_valid'))
        tsv(folder/'silhouette_section_allocation.tsv',alloc.reset_index())
        signature=dict(dataset=ds,features_sha256=sha(folder/'composition_features.csv.gz'),
            metadata_sha256=sha(folder/'spot_composition_inputs.parquet'),
            evaluation_sha256=sha(folder/'silhouette_evaluation_spots.tsv'),
            parameters_sha256=sha(OUT/'parameters.json'),code_checksums_sha256=sha(OUT/'code_checksums.json'),
            upstream_input_signature=upstream['input_signature'],n_valid=n,n_evaluation=target,
            evaluation_mode='all_valid_spots' if n<=5000 else 'fixed_section_stratified_sample_approximation')
        signature['signature']=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()
        savejson(folder/'input_signature.json',signature);dataset_dirs.append(str(folder))
        print(f'PREPARED {ds}: {len(m)} frozen / {n} valid / {target} evaluation spots',flush=True)
    tsv(OUT/'input_composition_manifest.tsv',rows)
    savejson(OUT/'input_checksums.json',signature_inputs+[
        dict(path=str(P/'ST/code/03_spatial_QC_summary.R'),sha256=sha(P/'ST/code/03_spatial_QC_summary.R')),
        dict(path=str(P/'ST/code/16_run_full_RCTD_one.R'),sha256=sha(P/'ST/code/16_run_full_RCTD_one.R')),
        dict(path=str(ident_path),sha256=sha(ident_path)),dict(path=str(ATTACHMENT),sha256=sha(ATTACHMENT))])
    savejson(OUT/'dataset_execution_order.json',dict(smoke_dataset=order[0],dataset_order=order,dataset_dirs=dataset_dirs))
    return dataset_dirs

def refresh_code_signatures():
    # Used only before smoke when input/parameter checks remain unchanged.
    for rec in json.loads((OUT/'input_checksums.json').read_text()):
        assert sha(rec['path'])==rec['sha256']
    code_files=sorted((OUT/'code').glob('*.*'))+[ROOT/'code/common.py']
    savejson(OUT/'code_checksums.json',{str(p):sha(p) for p in code_files if p.is_file()})
    order=json.loads((OUT/'dataset_execution_order.json').read_text())
    for folder in order['dataset_dirs']:
        folder=Path(folder);sig=json.loads((folder/'input_signature.json').read_text());sig.pop('signature')
        for filename,key in [('composition_features.csv.gz','features_sha256'),
            ('spot_composition_inputs.parquet','metadata_sha256'),('silhouette_evaluation_spots.tsv','evaluation_sha256')]:
            assert sha(folder/filename)==sig[key]
        assert sha(OUT/'parameters.json')==sig['parameters_sha256']
        sig['code_checksums_sha256']=sha(OUT/'code_checksums.json')
        sig['signature']=hashlib.sha256(json.dumps(sig,sort_keys=True).encode()).hexdigest()
        savejson(folder/'input_signature.json',sig)

if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='--refresh-code-signatures': refresh_code_signatures()
    else: prepare()
