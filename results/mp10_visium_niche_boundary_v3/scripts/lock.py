from common3 import *
C=OUT/'classification';labels=C/'section_classification.tsv';meta=read(labels)
files=[O/'cache/immune_metadata.tsv',O/'scripts/cluster.py',O/'scripts/features.py',O/'cache/immune_features.tsv',O/'cache/immune_windows.parquet',O/'cache/immune_links.tsv.gz',C/'feature_selection.tsv',C/'X_pam.tsv',C/'X_z.tsv',C/'bootstrap_runs.tsv.gz',C/'sensitivity_models.tsv',labels]
lock=dict(seed=seed,parameters=dict(grid=GRID,axis=AXIS,pool=list('ABCDEFG'),drop_abs_spearman_above=.8,family_weight='1/sqrt(m)',PAM='deterministic BUILD then best single medoid SWAP until no improvement; Euclidean K 3 (2,4 comparisons)',GMM='same X_z, K 1–4, diagonal primary and tied sensitivity, regularizer .001 (.0001,.01), 20 starts, 500 iterations also for bootstrap',bootstrap=100),hashes={str(f):sha(f) for f in files},labels=dict(zip(meta.section_id,meta.exploratory_cluster)),interpretation='C1–C3 retained. Amount/enrichment/heterogeneity profiles overlap; neither medoid nor local dispersion establishes distinct coexisting entry and exclusion states. Published TLS or cohort findings do not validate these three labels at section level.')
path=OUT/'classification_lock.json'
if path.exists():
    old=json.loads(path.read_text())
    if old['labels']!=lock['labels'] or old['seed']!=seed:raise RuntimeError('Different classification: use a new V3_MODEL_OUT')
path.write_text(json.dumps(lock,ensure_ascii=False,indent=2))

