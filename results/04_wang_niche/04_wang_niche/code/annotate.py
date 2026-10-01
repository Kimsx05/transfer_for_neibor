from prepare import *
from scipy.stats import rankdata, norm, mannwhitneyu
from scipy.spatial.distance import pdist, squareform

def bh(p):
    p=np.asarray(p,float);order=np.argsort(p,kind='stable');n=len(p)
    q=np.minimum.accumulate((p[order]*n/np.arange(1,n+1))[::-1])[::-1]
    out=np.empty(n);out[order]=np.minimum(q,1);return out

def annotate(folder):
    folder=Path(folder);status=json.loads((folder/'discovery_status.json').read_text())
    ds=status['dataset'];params=json.loads((OUT/'parameters.json').read_text());factors=params['factor_order']
    full=pd.read_parquet(folder/'spot_composition_inputs.parquet')
    assign=readcsv(folder/'selected_assignments.tsv',sep='\t')
    assert assign.spot_key.is_unique and set(assign.spot_key)==set(full.loc[full.composition_valid,'spot_key'])
    cols=['spot_key','niche_id','raw_louvain_cluster']
    full=full.merge(assign[cols],on='spot_key',how='left',validate='one_to_one',sort=False)
    full['assignment_reason']=np.where(full.composition_valid,'ASSIGNED',full.composition_reason)
    output_cols=['section_id','barcode','spot_key','dataset','cohort','patient_uid','specimen',
        'composition_valid','composition_reason','total_mean_abundance','niche_id','raw_louvain_cluster',
        'assignment_reason','x','y']
    full[output_cols].to_parquet(folder/'spot_niche_assignments.parquet',index=False,compression='zstd')
    valid=full.loc[full.composition_valid].reset_index(drop=True)
    x=valid[factors].to_numpy(float);labels=valid.niche_id.to_numpy();n=len(valid)
    niches=sorted(assign.niche_id.unique())
    assert n==status['n_valid'] and len(niches)==status['n_niches']
    ident=readcsv(OUT/'factor_interpretation.tsv',sep='\t').set_index('factor')
    major_map=readcsv(ROOT/'00_inputs/major_map.tsv',sep='\t').set_index('subtype').major.to_dict()
    # Re-read accepted sources for absolute-abundance annotation only; never change discovery features.
    cm=readcsv(ROOT/'00_inputs/c2l_manifest.tsv',sep='\t').set_index('section_id')
    abundance=np.empty_like(x)
    for sec,indices in valid.groupby('section_id',sort=False).groups.items():
        d=readcsv(cm.loc[sec,'mean_path']);d=d.rename(columns={d.columns[0]:'barcode'}).set_index('barcode')
        abundance[indices]=d.loc[valid.loc[indices,'barcode'],factors].to_numpy(float)
    assert np.allclose(abundance/abundance.sum(axis=1)[:,None],x,rtol=0,atol=1e-15)
    # One exact rank vector per factor; identical to asymptotic two-sided Wilcoxon rank sum.
    ranks=np.column_stack([rankdata(x[:,j],method='average') for j in range(len(factors))])
    tie_term=[]
    for j in range(len(factors)):
        cnt=np.unique(x[:,j],return_counts=True)[1].astype(float)
        tie_term.append(float(np.sum(cnt**3-cnt)))
    annotations=[];centroids=[];summaries=[];major_rows=[];coverage=[];diagnostics=[]
    numerical_checks=[];means=[]
    for ni,niche in enumerate(niches):
        mask=labels==niche;nx=int(mask.sum());ny=n-nx
        mean=x[mask].mean(axis=0);means.append(mean)
        median=np.median(x[mask],axis=0);restmedian=np.median(x[~mask],axis=0)
        absmean=abundance[mask].mean(axis=0);absmedian=np.median(abundance[mask],axis=0)
        inside=valid.loc[mask];sec_counts=inside.section_id.value_counts()
        patients=inside.patient_uid.dropna().unique()
        patient_known=len(patients)>0
        top=np.argsort(-mean,kind='stable')[:3]
        descriptor='; '.join(f'{factors[j]} {mean[j]:.1%}' for j in top)
        low_top=[factors[j] for j in top if ident.loc[factors[j],'status']=='LOW_IDENTIFIABILITY']
        summaries.append(dict(dataset=ds,niche_id=niche,n_spots=nx,fraction_dataset=nx/n,
            n_sections=inside.section_id.nunique(),n_patients_known=len(patients) if patient_known else np.nan,
            n_spots_patient_unknown=int(inside.patient_uid.isna().sum()),composition_description=descriptor,
            low_identifiability_in_top_factors='; '.join(low_top),interpretation_status='PASS_WITH_WARNING',
            name_source='top 3 mean subtype fractions, dataset-scoped; descriptive only'))
        diagnostics.append(dict(dataset=ds,niche_id=niche,n_spots=nx,fraction_dataset=nx/n,
            dominant_section_id=sec_counts.index[0],dominant_section_fraction=int(sec_counts.iloc[0])/nx,
            n_sections=inside.section_id.nunique(),n_patients_known=len(patients) if patient_known else np.nan,
            interpretation='section dominance and cluster size are diagnostics only; no coarsening'))
        for j,factor in enumerate(factors):
            u=float(ranks[mask,j].sum()-nx*(nx+1)/2)
            variance=nx*ny/12*((n+1)-tie_term[j]/(n*(n-1)))
            p=min(1.,float(2*norm.sf((abs(u-nx*ny/2)-0.5)/np.sqrt(variance)))) if variance>0 else 1.
            effect=2*u/(nx*ny)-1
            delta=float(median[j]-restmedian[j])
            annotations.append(dict(dataset=ds,niche_id=niche,factor=factor,major=major_map[factor],
                n_in=nx,n_rest=ny,U_statistic=u,pvalue=p,rank_biserial_effect=effect,
                median_fraction_in=float(median[j]),median_fraction_rest=float(restmedian[j]),
                median_fraction_difference=delta,mean_fraction_in=float(mean[j]),
                mean_fraction_rest=float(x[~mask,j].mean()),mean_abundance_in=float(absmean[j]),
                median_abundance_in=float(absmedian[j]),mean_abundance_rest=float(abundance[~mask,j].mean()),
                effect_direction='higher_in_niche' if effect>0 else 'lower_in_niche' if effect<0 else 'equal_ranks',
                factor_interpretation_status=ident.loc[factor,'status'],
                original_downstream_usage=ident.loc[factor,'downstream_usage'],
                interpretation_reason=ident.loc[factor,'reason'],test='two-sided Wilcoxon rank sum; ties+continuity',
                BH_family='all niche x 30 factor tests within dataset',role='defining annotation; not independent validation'))
            centroids.append(dict(dataset=ds,niche_id=niche,factor=factor,major=major_map[factor],n_spots=nx,
                mean_fraction=float(mean[j]),median_fraction=float(median[j]),
                mean_abundance=float(absmean[j]),median_abundance=float(absmedian[j]),
                factor_interpretation_status=ident.loc[factor,'status']))
            # Independent library implementation verifies the optimized rank calculations.
            if ni==0 or j in (8,12):
                ref=mannwhitneyu(x[mask,j],x[~mask,j],alternative='two-sided',method='asymptotic',use_continuity=True)
                numerical_checks.append(dict(dataset=ds,niche_id=niche,factor=factor,
                    optimized_pvalue=p,scipy_pvalue=float(ref.pvalue),absolute_difference=abs(p-float(ref.pvalue)),
                    U_difference=abs(u-float(ref.statistic)),status='PASS'))
                assert abs(p-float(ref.pvalue))<1e-12 and u==float(ref.statistic)
        for major in sorted(set(major_map.values())):
            cols=[j for j,f in enumerate(factors) if major_map[f]==major]
            major_rows.append(dict(dataset=ds,niche_id=niche,major=major,n_spots=nx,
                mean_fraction=float(x[mask][:,cols].sum(axis=1).mean()),
                median_fraction=float(np.median(x[mask][:,cols].sum(axis=1))),
                mean_abundance=float(abundance[mask][:,cols].sum(axis=1).mean()),
                role='sum of original subtype means; interpretation only, not clustering input'))
        for sec,g in full.groupby('section_id',sort=False):
            hit=g.niche_id.eq(niche);den=int(g.composition_valid.sum())
            coverage.append(dict(dataset=ds,niche_id=niche,coverage_scope='section',section_id=sec,
                patient_uid=g.patient_uid.iloc[0],specimen=g.specimen.iloc[0],cohort=g.cohort.iloc[0],
                n_spots=int(hit.sum()),n_frozen_spots=len(g),n_valid_compositions=den,
                fraction_valid_spots=float(hit.sum()/den) if den else np.nan,status='PASS'))
        if full.patient_uid.notna().any():
            for patient,g in full.loc[full.patient_uid.notna()].groupby('patient_uid',sort=False):
                den=int(g.composition_valid.sum())
                coverage.append(dict(dataset=ds,niche_id=niche,coverage_scope='patient',patient_uid=patient,
                    n_spots=int(g.niche_id.eq(niche).sum()),n_frozen_spots=len(g),n_valid_compositions=den,
                    fraction_valid_spots=float(g.niche_id.eq(niche).sum()/den) if den else np.nan,status='PASS'))
        else:
            coverage.append(dict(dataset=ds,niche_id=niche,coverage_scope='patient_unknown',patient_uid=np.nan,
                n_spots=nx,n_frozen_spots=len(full),n_valid_compositions=n,fraction_valid_spots=np.nan,
                status='SKIPPED_NOT_ESTIMABLE',reason='No frozen patient identifiers; never infer patient from section/specimen'))
    annotation_df=pd.DataFrame(annotations);annotation_df['FDR']=bh(annotation_df.pvalue)
    annotation_df['significant_FDR_lt_0_05']=annotation_df.FDR.lt(.05)
    annotation_df['overrepresented_for_description']=annotation_df.significant_FDR_lt_0_05&annotation_df.rank_biserial_effect.gt(0)
    # Do not label a rare factor "rich" from significance alone.
    annotation_df['notable_fraction_for_description']=annotation_df.mean_fraction_in.ge(.05)
    tsv(folder/'niche_annotation_wilcoxon.tsv',annotation_df)
    for name,rows in [('niche_centroids',centroids),('niche_descriptions',summaries),('niche_major_summary',major_rows),
        ('niche_sample_coverage',coverage),('niche_diagnostics',diagnostics),('wilcoxon_numerical_validation',numerical_checks)]:
        tsv(folder/f'{name}.tsv',rows)
    distances=squareform(pdist(np.asarray(means),metric='euclidean'))
    pairs=[]
    for i in range(len(niches)):
        for j in range(i+1,len(niches)):
            pairs.append(dict(dataset=ds,niche_id_a=niches[i],niche_id_b=niches[j],
                centroid_euclidean_distance=float(distances[i,j]),
                diagnostic_close=bool(distances[i,j]<params['centroid_redundancy_diagnostic_threshold']),
                action='none; diagnostic only; labels remain selected'))
    tsv(folder/'centroid_redundancy.tsv',pd.DataFrame(pairs,columns=['dataset','niche_id_a','niche_id_b',
        'centroid_euclidean_distance','diagnostic_close','action']))
    hist=readcsv(ROOT/'00_inputs/histology_manifest.tsv',sep='\t')
    tsv(folder/'spatial_plot_metadata.tsv',hist.loc[hist.dataset.eq(ds)])
    tsv(folder/'spatial_plot_spots.tsv',full[['section_id','barcode','x','y','niche_id','composition_valid']])
    nerr=max([r['absolute_difference'] for r in numerical_checks],default=0)
    savejson(folder/'annotation_status.json',dict(dataset=ds,status='PASS_WITH_WARNING',
        n_assignments=len(full),n_valid=n,n_invalid=len(full)-n,n_niches=len(niches),
        n_wilcoxon_tests=len(annotation_df),n_significant=int(annotation_df.significant_FDR_lt_0_05.sum()),
        max_wilcoxon_validation_difference=nerr,patient_coverage_status='SKIPPED_NOT_ESTIMABLE' if full.patient_uid.isna().all() else 'PASS',
        warnings=['Annotation is derived from clustering features, not independent validation.',
            'Low-identifiability factors retained in features and explicitly flagged.']))
    print(f'ANNOTATED {ds}: {len(niches)} niches / {len(annotation_df)} defining tests',flush=True)

if __name__=='__main__': annotate(sys.argv[1])
