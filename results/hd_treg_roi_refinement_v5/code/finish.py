"""Verify saved outputs and package the review materials, preserving recorded visual judgements."""
from common5 import *
from load_result import load_result
from zipfile import ZipFile,ZIP_DEFLATED
tables=R/'tables'
visual=rd(tables/'geometry_review.tsv')
for phase in ('round1','round2'):
    im=rd(tables/(phase+'_images.tsv'))
    visual=visual.merge(im[['section','roi','path']].rename(columns={'path':phase+'_file'}),on=['section','roi'],validate='one_to_one')
visual['score_panel_note']=[
    '主体内散在可评分上皮。焦点附近原始红格与多个核相邻。',
    '孤立点周围缺少可评分上皮和连续组织。',
    '局部深染碎片。两个圆内缺少可解释的细胞邻域。',
    '孤立染色点。未见成片可评分上皮。',
    '组织小突起进入圆内，较大邻域跨空白，局部上皮机会受限。',
    '焦点在边缘，部分圆在组织之外，边界处红格归属未定。',
    '主体核密集，上皮分数的可见范围依组织位置变化。十个辅助中心另附独立邻域图。',
    '离散组织片段，局部稀疏蓝格。焦点附近缺少可辨识的上皮分数。']
visual['action_after_MP10']='身份和几何冻结，无评分后修订'
put(visual,tables/'roi_visual_review.tsv')
pages=[]; quality=[]; cover=[]; dispo=[]; check=[]
for s in sections:
    o=R/s
    res=load_result(s,counts=True)
    d=res['cells_scored'];c=res['candidate_ledger'];x=res['counts']
    a=rd(o/'counts_audit.tsv.gz')
    assert c.unit.is_unique and c.roi.gt(0).all() and c.disposition.notna().all()
    assert d.cell_id.is_unique and np.array_equal(d.id,np.arange(1,len(d)+1))
    assert np.array_equal(np.asarray(x.sum(axis=1)).ravel(),d.nCount)
    assert np.array_equal(np.asarray(x.sum(axis=0)).ravel(),a.assigned)
    assert np.array_equal(a.input,a.assigned+a.excluded+a.unassigned)
    ep=d[d.score_ok];sx=x[ep.index].tocsr()
    assert ep.MP10_cell_score.notna().all() and d.loc[~d.score_ok,'MP10_cell_score'].isna().all()
    avail=[]; tied=[]
    for i in range(sx.shape[0]):
        vals=sx.data[sx.indptr[i]:sx.indptr[i+1]]
        _,freq=np.unique(vals,return_counts=True)
        freq=freq[::-1];mean_rank=np.cumsum(freq)-(freq-1)/2
        avail.append(freq[mean_rank<905].sum());tied.append(freq[freq>1].sum()/len(vals))
    tv=ep.MP10_cell_score
    symbols=res['genes'].symbol.to_numpy(); assert all(np.count_nonzero(symbols==g)==1 for g in top)
    top_cols=[np.flatnonzero(symbols==g)[0] for g in top]
    assert min(top_cols)>=0
    gene_width=sx[:,top_cols].getnnz(axis=1)
    for name,arr in [('UMI',ep.nCount),('detected_genes',ep.nFeature),('top20_detected',gene_width),('positive_genes_avg_rank_lt_905',avail),('fraction_detected_genes_in_ties',tied),('MP10',tv)]:
        q=np.percentile(arr,[0,25,50,75,100])
        quality.append([s,name,len(arr),*q,len(np.unique(arr)),int(np.sum(np.asarray(arr)==0))])
    cc=pd.read_parquet(o/'candidate_coverage.parquet')
    assert np.array_equal(c.unit,cc.unit)
    for roi,g in cc.groupby('roi'):
        for rad in [50,200]:
            den=g[f'disk_{rad}um_sites'];f=g[f'supplied_{rad}um_sites']/den
            cover.append([s,roi,rad,len(g),f.median(),sum(f>=.5)])
    for name,n in c.disposition.value_counts().items():dispo.append([s,name,n])
    z=pd.read_parquet(o/'zoom_index.parquet')
    assert not z.duplicated(['file','slot']).any()
    for j,(file,g) in enumerate(z.groupby('file',sort=True)):
        direct=(s=='Cervilla' and j<44) or (s=='D1_Johnson' and g.roi.iloc[0]==2)
        p=Path(file) if direct else R/'figures'/'review_sheets'/f'{s}_{((j-44 if s=="Cervilla" else j)//4)+1:02d}.jpg'
        assert p.is_file() and Path(file).is_file()
        pages.append([s,g.roi.iloc[0],file,str(p),len(g),'AI辅助拼版目视复核，逐对象身份按冻结的RNA、参考及几何规则记录'])
    chk=rd(o/'historical_score_check.tsv')
    assert np.allclose(chk.MP10,chk.recomputed,atol=1e-7)
    check.append([s,len(c),len(d),x.shape[1],int(x.sum()),len(chk),max(abs(chk.MP10-chk.recomputed)),len(z),z.file.nunique()])
    d.loc[d.group.eq('RNA_candidate')].to_csv(o/'auxiliary_candidates.tsv',sep='\t',index=False)
def table(data,cols,name):
    put(pd.DataFrame(data,columns=cols.split()),tables/name)
table(pages,'section roi source_sheet viewed_file fields review_mode','detail_review.tsv')
table(quality,'section metric n min q25 median q75 max distinct zeros','mp10_measurement.tsv')
table(cover,'section roi radius_um screening_entries median_supplied_fraction n_at_least_half_supplied','coverage_summary.tsv')
table(dispo,'section linked_cell_disposition_combinations screening_entries','source_dispositions.tsv')
table(check,'section unique_screen_entries cells genes assigned_UMI historical_comparisons max_score_error indexed_zoom_fields reviewed_source_pages','verification.tsv')
table([['screen_and_geometry','complete','8 merged ROIs; all observed candidates have ROI'],['proseg','not_run','Rust toolchain mismatch; H&E and UMI gradient assignment executed in both sections'],['nuclear_detection','reused','registered StarDist nuclear masks; no new nucleus detector run'],['primary_contrast','not_estimable','0 supported Treg in either section'],['auxiliary_contrast','descriptive_only','D1 200 um only; two matched blocks'],['AI_image_review','complete','8 ROIs and 150 indexed detail sheets; molecular identity remains uncertain']],'branch state reason','run_status.tsv')
with ZipFile(R/'review.zip','w',ZIP_DEFLATED,compresslevel=5) as out:
    files=[*tables.glob('*.tsv'),*(R/'code').glob('*.py'),*(R/'figures').glob('*section.png'),*(R/'figures').glob('*.pdf'),R/'figures/D1_Johnson_RNA_candidates.jpg']
    files.extend(R/name for name in ['REPORT.md','METHODS.md','frozen_config.yaml','input_manifest.tsv','rerun.sh'])
    for s in sections:files.extend(R/s/file for file in ['auxiliary_candidates.tsv','identity_funnel.tsv','historical_score_check.tsv'])
    for f in sorted(set(files)):out.write(f,f.relative_to(R))
with ZipFile(R/'review.zip') as out:
    assert out.testzip() is None
print(pd.read_csv(tables/'verification.tsv',sep='\t').to_string(index=False))
