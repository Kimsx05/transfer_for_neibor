from prepare import *
import html, zipfile
from PIL import Image

PAPER='https://link.springer.com/article/10.1186/s12943-025-02377-9'

def markdown_table(frame):
    cols=frame.columns.tolist()
    lines=['| '+' | '.join(cols)+' |','| '+' | '.join(['---']*len(cols))+' |']
    lines += ['| '+' | '.join(str(v).replace('|','\\|') for v in row)+' |' for row in frame.itertuples(index=False,name=None)]
    return '\n'.join(lines)

def methods():
    items=[
        ('cellular input','cell2location cellular compositions','accepted c2l outputs; all 30 factors','PAPER_CORE','model itself is not retrained'),
        ('posterior summary','NOT_REPORTED','original accepted posterior mean decimal export','PROJECT_ADAPTATION','not q05 or recovered major means'),
        ('transformation','NOT_REPORTED','all-factor row proportions; invalid rows NA','PROJECT_ADAPTATION','no log/z/ILR/PCA/Harmony'),
        ('feature space','cellular compositions','all frozen subtype factors including Epithelial','PROJECT_ADAPTATION','no TME-only subset, no neighborhood averages'),
        ('pooling','similar-composition spots across tissue slices','pool sections within each frozen dataset','PROJECT_ADAPTATION','shared dictionary within dataset; separate labels across datasets'),
        ('discovery universe','NOT_REPORTED as an epithelial/CNV filter','all frozen valid composition spots','PROJECT_ADAPTATION','no Epi gate, MP/CNV mask, histology, coordinate or patient feature'),
        ('SNN graph','scran::buildSNNGraph','scran::buildSNNGraph','PAPER_CORE','composition neighbors, not physical lattice neighbors'),
        ('k','10,20,30','10,20,30','PAPER_CORE','three composition-space graphs'),
        ('SNN weights','NOT_REPORTED','type=rank (explicit installed default)','PROJECT_ADAPTATION','finite positive undirected graph checked'),
        ('matrix orientation','NOT_REPORTED','spots x 30 factors; transposed=TRUE','PROJECT_ADAPTATION','d=NA explicitly disables PCA'),
        ('nearest-neighbor backend','NOT_REPORTED','KmknnParam exact Euclidean; serial','PROJECT_ADAPTATION','fixed seed for graph construction'),
        ('clustering algorithm','Louvain','igraph::cluster_louvain','PAPER_CORE','not Leiden; same graph reused for its 15 resolutions'),
        ('resolution grid','NOT_REPORTED','0.1,0.2,...,1.5','PROJECT_ADAPTATION','45 total candidates per dataset'),
        ('seed','NOT_REPORTED','20260930','PROJECT_ADAPTATION','same seed reset before every graph and candidate'),
        ('compute scheduling','NOT_REPORTED','one dataset at a time; four candidate workers; graph/BLAS/OpenMP each 1 thread','PROJECT_ADAPTATION','CPU budget 8; serial-versus-parallel 45-candidate membership and silhouette parity verified on first two datasets'),
        ('parameter selection','maximize mean silhouette; full aggregation not specified','evaluation-spot weighted mean silhouette','PROJECT_ADAPTATION','per-cluster means reported separately'),
        ('silhouette distance','NOT_REPORTED','Euclidean in original 30 proportion features','PROJECT_ADAPTATION','no UMAP, PCA or physical-coordinate distances'),
        ('evaluation set','NOT_REPORTED','all if <=5000, else fixed proportional section-stratified 5000','PROJECT_ADAPTATION','six datasets use sample-based approximation; no full dataset pair matrix'),
        ('silhouette ties','NOT_REPORTED','tolerance 1e-8; lower resolution then lower k','PROJECT_ADAPTATION','declared before fitting'),
        ('degenerate candidates','NOT_REPORTED','exclude one class / all singleton / invalid silhouette','PROJECT_ADAPTATION','excluded reasons retained'),
        ('annotation','niche vs rest Wilcoxon, FDR<0.05','two-sided unpaired rank sum on proportions; BH<0.05','PAPER_CORE','defining annotation, not independent biological validation'),
        ('Wilcoxon numerical options','NOT_REPORTED','asymptotic, tie-adjusted variance, continuity correction','PROJECT_ADAPTATION','optimized calculation independently checked against scipy'),
        ('FDR family','NOT_REPORTED','all niche x factor tests within each dataset','PROJECT_ADAPTATION','more transparent than assuming author correction scope'),
        ('naming and cross-dataset IDs','NOT_REPORTED','dataset-scoped IDs and numerical top-three factor descriptions','PROJECT_ADAPTATION','no rare-factor rich labels from small p-values; no cross-dataset label merging'),
        ('low-identifiability factors','NOT_REPORTED for bladder reference','retain factors; attach accepted interpretation status','PROJECT_ADAPTATION','no automatic factor merging'),
        ('scran version','1.32.0','1.38.1 installed','PROJECT_ADAPTATION','installed interface recorded; no base-environment changes'),
        ('R version','4.3.0','4.5.3 installed','PROJECT_ADAPTATION','versions of all graph packages recorded'),
        ('upstream SCT/Harmony','paper upstream expression preprocessing','existing accepted c2l inputs reused; no rerun','OUT_OF_SCOPE','do not apply expression integration to abundance'),
        ('spatial NMF','separate paper analysis','not performed','OUT_OF_SCOPE','no NMF or alternative method competition'),
        ('author analysis repository','NOT_RECOVERED','no public author implementation located','NOT_REPORTED','publisher article/data availability and DOI/author searches; supplement link retrieval failed'),
    ]
    tsv(OUT/'methods_parity.tsv',[dict(component=a,paper_or_author_source=b,project_implementation=c,
        evidence_status=d,note=e,source_url=PAPER) for a,b,c,d,e in items])
    savejson(OUT/'public_source_search.json',dict(article_url=PAPER,doi='10.1186/s12943-025-02377-9',
        search_queries=['"s12943-025-02377-9" code github','"FABP1" "PLG" "PLAT" "buildSNNGraph"',
        '"s12943-025-02377-9" ("code" OR "GitHub") Wang','"Yiqiu Wang" "FABP1" code'],
        findings='Publisher methods and data-availability inspected; no GitHub link or author code repository recovered. Public searches returned the article and method extracts, no implementation.',
        supplementary_status='Publisher DOCX link exists; web link retrieval returned Internal Error, content not reviewed.',
        author_posterior_summary='NOT_REPORTED',author_transform='NOT_REPORTED',author_resolution_grid='NOT_REPORTED',
        author_silhouette_aggregation='NOT_REPORTED in sufficient implementation detail',
        limitation='Search result is not proof that no public code exists; no claim of exact parameter reproduction.'))

def finalize():
    order=json.loads((OUT/'dataset_execution_order.json').read_text());folders=[Path(p) for p in order['dataset_dirs']]
    methods()
    required=['selected_parameters','silhouette_grid','niche_centroids','niche_annotation_wilcoxon','niche_sample_coverage',
        'niche_descriptions','niche_major_summary','niche_diagnostics','per_cluster_silhouette','graph_diagnostics',
        'centroid_redundancy','figure_index','wilcoxon_numerical_validation','plot_geometry_validation']
    for stem in required:
        tsv(OUT/f'{stem}.tsv',pd.concat([readcsv(f/f'{stem}.tsv',sep='\t') for f in folders],ignore_index=True))
    assign=pd.concat([pd.read_parquet(f/'spot_niche_assignments.parquet') for f in folders],ignore_index=True)
    assert not assign.duplicated(['section_id','barcode']).any()
    assign.to_parquet(OUT/'spot_niche_assignments.parquet',index=False,compression='zstd')
    original=pd.read_parquet(ROOT/'00_inputs/frozen_spot_inputs.parquet',columns=['section_id','barcode','dataset','patient_uid','specimen','cohort','x','y'])
    assert set(zip(assign.section_id,assign.barcode))==set(zip(original.section_id,original.barcode))
    left=assign.set_index(['section_id','barcode']).sort_index();right=original.set_index(['section_id','barcode']).sort_index()
    for col in right:
        assert left[col].fillna('__MISSING__').equals(right[col].fillna('__MISSING__')),f'Metadata changed: {col}'
    assert assign.loc[assign.composition_valid,'niche_id'].notna().all()
    assert assign.loc[~assign.composition_valid,'niche_id'].isna().all()
    checks=[dict(check='full frozen spot key universe retained',n_checked=len(assign),max_error=0,status='PASS'),
        dict(check='dataset, patient, specimen, cohort and image coordinates unchanged',n_checked=len(assign),max_error=0,status='PASS')]
    selected=readcsv(OUT/'selected_parameters.tsv',sep='\t');grid=readcsv(OUT/'silhouette_grid.tsv',sep='\t')
    assert len(grid)==45*len(folders)
    for _,r in selected.iterrows():
        g=grid.loc[grid.dataset.eq(r.dataset)&grid.eligible].copy()
        best=g.mean_silhouette.max();candidate=g.loc[g.mean_silhouette.ge(best-1e-8)].sort_values(['resolution','k']).iloc[0]
        assert candidate.candidate==r.candidate
        f=OUT/'by_dataset'/r.dataset
        features=readcsv(f/'composition_features.csv.gz');ev=readcsv(f/'silhouette_evaluation_spots.tsv',sep='\t')
        assert features.columns[1:].tolist()==json.loads((OUT/'parameters.json').read_text())['factor_order']
        assert features.spot_key.iloc[ev.input_row_1based.to_numpy()-1].tolist()==ev.spot_key.tolist()
        e=float(np.max(np.abs(features.iloc[:,1:].sum(axis=1)-1)))
        assert e<1e-12
        checks.append(dict(check=f'{r.dataset}: 45-candidate selection, shared evaluation keys, 30-feature row sums',
            n_checked=len(features),max_error=e,status='PASS'))
    for rec in json.loads((OUT/'input_checksums.json').read_text()):
        assert sha(rec['path'])==rec['sha256'],f'Input changed: {rec["path"]}'
    for p,h in json.loads((OUT/'upstream_freeze_checksums.json').read_text()).items():
        assert sha(p)==h,f'Frozen upstream changed: {p}'
    for p,h in json.loads((OUT/'code_checksums.json').read_text()).items():
        assert sha(p)==h,f'Code changed during execution: {p}'
    checks.append(dict(check='all mean input SHA and upstream freeze signatures unchanged',n_checked=63,max_error=0,status='PASS'))
    figure=readcsv(OUT/'figure_index.tsv',sep='\t');image_checks=[]
    for _,r in figure.iterrows():
        assert Path(r.pdf).read_bytes()[:4]==b'%PDF'
        with Image.open(r.png) as im:
            im.verify()
        with Image.open(r.png) as im:
            dpi=im.info.get('dpi');assert dpi and abs(dpi[0]-600)<1
            image_checks.append(dict(dataset=r.dataset,section_id=r.section_id,view=r.view,png=r.png,
                width=im.width,height=im.height,dpi_x=dpi[0],dpi_y=dpi[1],status='PASS'))
    tsv(OUT/'figure_file_validation.tsv',image_checks);tsv(OUT/'acceptance_validation.tsv',checks)
    visual_path=OUT/'visual_review.tsv'
    assert visual_path.is_file(),'Human-readable agent visual review required before final handoff'
    visual=readcsv(visual_path,sep='\t');assert visual.status.eq('PASS').all()
    assert set(visual.dataset)==set(selected.dataset)
    summary=[]
    for f in folders:
        discovery=json.loads((f/'discovery_status.json').read_text());annotation=json.loads((f/'annotation_status.json').read_text())
        completion=json.loads((f/'completion.json').read_text());assert completion['status']=='PASS_WITH_WARNING'
        for p,h in completion['output_checksums'].items():assert sha(p)==h,f'Dataset cache changed: {p}'
        summary.append(dict(dataset=discovery['dataset'],n_spots=annotation['n_assignments'],n_valid=discovery['n_valid'],
            n_niches=discovery['n_niches'],selected_candidate=discovery['candidate'],n_evaluation=discovery['n_evaluation'],
            evaluation_mode=discovery['evaluation_mode'],discovery_seconds=discovery['elapsed_seconds'],
            peak_parent_RSS_GiB=discovery['peak_RSS_GiB'],max_candidate_worker_RSS_GiB=discovery['max_candidate_worker_RSS_GiB'],
            conservative_sum_peak_RSS_GiB=discovery['conservative_sum_peak_RSS_GiB'],
            candidate_workers=discovery['candidate_workers'],patient_coverage_status=annotation['patient_coverage_status'],status='PASS_WITH_WARNING'))
    tsv(OUT/'dataset_completion.tsv',summary)
    warning=['Posterior mean proportions, dataset pools, resolution grid, SNN rank weighting, evaluation sampling and detailed selection rules are transparent project adaptations; author code not recovered.',
        'Six datasets use a fixed 5000-spot silhouette approximation; GSE171351 uses all valid spots.',
        'Frozen patient_uid absent; patient coverage is SKIPPED_NOT_ESTIMABLE; section coverage is available.',
        'Low-identifiability subtype factors retained and flagged; defining Wilcoxon/BH annotations are not independent validation.',
        'CanDis, Inhousedata and NPJ dataset identities retain frozen UNVERIFIED labels.',
        'Installed scran 1.38.1 differs from paper 1.32.0; actual interfaces and versions recorded.']
    lines=['# Module 04 — PASS_WITH_WARNING','',
        f'完成 {len(folders)} 个冻结 dataset、{assign.section_id.nunique()} 张切片、{len(assign):,} 个 spot；有效组成 {int(assign.composition_valid.sum()):,}，无效组成 {int((~assign.composition_valid).sum()):,}。',
        f'每个 dataset 完成 3 个 k × 15 个 resolution = 45 个候选，共 {len(grid)} 个候选；选出 dataset-scoped niches 共 {int(sum(r["n_niches"] for r in summary))} 个。',
        '',markdown_table(pd.DataFrame(summary)[['dataset','n_spots','n_niches','selected_candidate','n_evaluation','evaluation_mode','status']]),'',
        '输入来自 Module 00 冻结的原始完整 30-factor posterior mean abundance CSV；全 factor 行和归一化为比例。按真实 (section_id, barcode) 对齐，原始 abundance 文件 SHA 与上游冻结 SHA 均未改变。所有有效组成进入发现，无 Epi/CNV 预筛；MP、H&E、空间坐标和患者标签不进入 feature matrix。未做 log、z、ILR、PCA、Harmony 或 c2l 重训。',
        '',f'文献核心来自 [Wang et al. 2025]({PAPER})：cell2location compositions、scran SNN 的 k=10/20/30、Louvain、silhouette 选参和 niche-vs-rest Wilcoxon/FDR。',
        '公开作者实现未恢复；posterior mean proportions、dataset 分池、resolution=0.1–1.5、rank weighting、精确 Euclidean NN、spot-weighted silhouette、固定 section 分层评估集和 BH family 都是透明项目适配，见 methods_parity.tsv。并非逐参数原样复刻。',
        '', 'silhouette 只在不超过 5,000 个固定评估 spot 的原始 30 维组成空间计算。N≤5,000 的 GSE171351 全量评估，其余六组是 sample-based approximation；同 dataset 的所有候选使用完全相同的评估集。预先设定 1e-8 tie tolerance，先低 resolution 再小 k。cluster size、section dominance、centroid redundancy 仅记录诊断，未据此改解或粗化。',
        '', 'CPU 预算8：dataset 逐组串行，组内4个 candidate workers，graph/BLAS/OpenMP均1线程，无嵌套并行或GPU。前两组的单 worker 和4-worker全部候选 memberships、silhouette逐值一致，见 execution_parallelism_parity.json。内存为逐进程 VmHWM；父进程+4倍最大worker峰值是保守上界，会重复计入fork共享内存，不宣称是真实合计RSS。',
        '', 'Wilcoxon 为独立两组的双侧 rank-sum，使用 ties 与 continuity correction；BH 在各 dataset 的所有 niche×factor tests 内校正。输出 U、rank-biserial effect、组内/组外 median fraction、差值、mean fraction、绝对 mean abundance 和 FDR。优化实现与 scipy 独立算法核对通过。注释描述组成特征，不是独立验证，也不能将空间 spot 当独立生物重复。',
        '', 'niche 描述使用实际的前三个 mean fractions；ID 仅在 dataset 内定义，跨 dataset 同编号没有对应意义。未自动合并标签；major 聚合只用于解释图。低可辨识性原始 subtype 标记随表提供，所有 T/NK factors 保留。',
        '', f'空间图 {assign.section_id.nunique()} 张，组成/major/选参/coverage 图 {4*len(folders)} 张，共 {len(figure)} views，均 PDF + 600 dpi PNG，另有缩略图。优先复用原 R save_both 函数，调整 PNG 至600 dpi；H&E alpha=0.25，spot alpha=1。零值黑色，NA 灰色，使用冻结 lowres 坐标，无再次缩放。视觉抽查覆盖七个 dataset。',
        '', '警告表示完成的分析存在解释或复刻边界，不表示运行报错：']
    lines += ['- '+w for w in warning]
    lines += ['', 'patient_uid 全缺失，因此患者层覆盖标 SKIPPED_NOT_ESTIMABLE；切片层覆盖正常输出。未推断患者编号。',
        '', '交付：input_composition_manifest.tsv、selected_parameters.tsv、silhouette_grid.tsv、spot_niche_assignments.parquet、niche_centroids.tsv、niche_annotation_wilcoxon.tsv、niche_sample_coverage.tsv、methods_parity.tsv；另有图、代码、版本、校验、逐 dataset 日志、output_manifest.tsv、handoff.json 和 review_bundle.zip。',
        '', '**停止位置：Module 04 完成。未运行 MP–niche 检验、CNV 关联或下一模块。原始数据、历史 MP 图与 Module 00/01 输出均未覆盖。**','']
    (OUT/'RUN_SUMMARY.md').write_text('\n'.join(lines))
    desc=readcsv(OUT/'niche_descriptions.tsv',sep='\t')
    page=['<!doctype html><meta charset="utf-8"><title>Module 04 review</title><style>body{font:15px system-ui;margin:24px;background:#fafafa}nav{display:flex;gap:16px;flex-wrap:wrap}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:16px}.card{background:white;padding:12px;border:1px solid #ddd}img{width:100%}table{border-collapse:collapse}td,th{padding:6px;border:1px solid #ddd}</style>',
        '<h1>Module 04 · Wang-style composition niches</h1><p>PASS_WITH_WARNING · Dataset-scoped labels · No MP association tests</p>',
        '<p>Bundle contains tables, provenance, thumbnails and representative PDFs. Full 600 dpi PNGs and all PDFs are listed in figure_index.tsv and remain in the module output directory.</p>',
        '<nav>'+''.join(f'<a href="#{html.escape(ds)}">{html.escape(ds)}</a>' for ds in selected.dataset)+'</nav>']
    for ds in selected.dataset:
        page.append(f'<h2 id="{html.escape(ds)}">{html.escape(ds)}</h2>')
        page.append(desc.loc[desc.dataset.eq(ds),['niche_id','n_spots','composition_description','low_identifiability_in_top_factors']].to_html(index=False,escape=True))
        page.append('<div class="grid">')
        for _,r in figure.loc[figure.dataset.eq(ds)].iterrows():
            rel=Path(r.thumbnail).relative_to(OUT).as_posix()
            page.append(f'<div class="card"><a href="{html.escape(rel)}"><img loading="lazy" src="{html.escape(rel)}"></a><p>{html.escape(r.stem)}</p></div>')
        page.append('</div>')
    (OUT/'review.html').write_text('\n'.join(page))
    # Compact standalone review: no raw abundance/features or full spot matrices.
    archive=OUT/'review_bundle.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(OUT.iterdir()):
            if p.is_file() and p.suffix in ('.tsv','.json','.md','.html') and p.name not in ('handoff.json','output_manifest.tsv'):
                z.write(p,p.relative_to(OUT))
        for p in sorted((OUT/'code').glob('*.*')):z.write(p,p.relative_to(OUT))
        for f in folders:
            for p in sorted((f/'thumbnails').glob('*.png')):z.write(p,p.relative_to(OUT))
            for name in ('discovery.log','plots.log','installed_interfaces.txt','software_versions_R.json','silhouette_section_allocation.tsv'):
                p=f/name;z.write(p,p.relative_to(OUT))
            for stem in ('niche_subtype_composition','niche_major_composition'):
                p=f/'figures'/f'{stem}.pdf';z.write(p,p.relative_to(OUT))
            r=figure.loc[figure.dataset.eq(f.name)&figure.view.eq('spatial_niche_labels')].iloc[0]
            p=Path(r.pdf);z.write(p,p.relative_to(OUT))
    manifest=output_manifest(OUT)
    savejson(OUT/'handoff.json',dict(module='04',status='PASS_WITH_WARNING',
        completed_scope='7 dataset-pooled 30-factor composition niche discovery, annotation, section coverage and spatial/composition figures',
        upstream_module='00',upstream_handoff=str(ROOT/'00_inputs/handoff.json'),
        dependencies=['00'],copykat_required=False,MP_display_required=False,
        authorization_document=str(ATTACHMENT),authorization_document_sha256=sha(ATTACHMENT),
        input_signature=sha(OUT/'input_checksums.json'),parameters_signature=sha(OUT/'parameters.json'),
        code_signature=sha(OUT/'code_checksums.json'),spot_table=str(OUT/'spot_niche_assignments.parquet'),
        spot_key_columns=['section_id','barcode'],n_frozen_spots=len(assign),n_valid_spots=int(assign.composition_valid.sum()),
        n_datasets=len(folders),n_sections=int(assign.section_id.nunique()),
        n_niches=sum(r['n_niches'] for r in summary),niche_id_scope='dataset only; never align numeric suffixes across datasets',
        primary_features='all 30 posterior mean row proportions; no additional transformation',
        excluded_composition_rows='NA assignment with composition_reason, never zero-filled',
        abundance_units='expected_cells_per_spot',composition_units='fraction of total 30-factor posterior mean',
        silhouette_units='dimensionless Euclidean-feature silhouette; evaluation mode per dataset',
        patient_coverage_status='SKIPPED_NOT_ESTIMABLE',warnings=warning,outputs=manifest,
        review_bundle=str(archive),review_html=str(OUT/'review.html'),MP_niche_tests_performed=False,
        stop_after_module='04',next_module_authorized=False))
    savejson(OUT/'completion.json',dict(module='04',status='PASS_WITH_WARNING',n_datasets=len(folders),
        n_sections=int(assign.section_id.nunique()),n_spots=len(assign),n_niches=sum(r['n_niches'] for r in summary),
        handoff=str(OUT/'handoff.json'),review_bundle=str(archive),stop_after_module='04'))
    # Refresh manifest once after completion.json; the manifest/handoff intentionally exclude themselves.
    manifest=output_manifest(OUT)
    handoff=json.loads((OUT/'handoff.json').read_text());handoff['outputs']=manifest;savejson(OUT/'handoff.json',handoff)
    print('MODULE04_COMPLETE',json.dumps(summary,ensure_ascii=False),flush=True)

if __name__=='__main__':finalize()
