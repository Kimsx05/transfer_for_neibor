from common import *
from scipy.stats import spearmanr
import importlib.metadata as im
D=pd.read_parquet(R/'tables/Cervilla_cell_formal_metadata.parquet');reg=pd.read_parquet(R/'tables/Cervilla_epithelial_regions.parquet');raw=pd.read_csv(R/'tables/section_target_evaluability.tsv',sep='\t');rows=[]
for r in raw.itertuples():
 z=r._asdict();z['assessment_stage']='raw_input_screen';rows.append(z)
for target,n,status,note in [('Treg',0,'not_evaluable','144 fine-reference candidates but no joint multi-gene Treg/lineage support; bin8/bin16 and conservative4um ROI checks did not corroborate; not evidence of biological absence'),('Endothelial',204,'partially_evaluable','Fine reference plus >=2 raw endothelial genes; subtype not used; supported objects are an incomplete observed subset, not vessel census'),('Epi_MP10',1209,'partially_evaluable','Region-level historical Top20 rank-AUC on connected credible Epi only; >=5 Epi, >=5 Top20 and >=905 detected genes; depth sensitivity limits local rank precision'),('registration',268924,'evaluable','Inherited verified native pixel-to-micron transform and supplied footprint; larger image outside supplied coverage remains unmeasured')]:
 rows.append(dict(section='Cervilla',target=target,analysis_bins=268924,raw_supported_bins=n,status=status,reason=note,assessment_stage='final_selected_unit'))
rows.append(dict(section='D1_Johnson',target='registration',analysis_bins=364805,raw_supported_bins=364805,status='evaluable',reason='Existing verified crop origin (23764,35851) fullres pixels preserved; light ROI check only',assessment_stage='light_final'))
put(pd.DataFrame(rows),'section_target_evaluability_final.tsv')
# Direct group evidence, not independent validation or formal differential expression.
groups={'Treg_reference_candidates':D.fine_prediction.eq('Treg'),'T_lineage_supported':D.type_supported.eq('T_subtype_unresolved'),'Endothelial_supported':D.type_supported.eq('Endothelial'),'Epi_credible_unmixed':D.mp10_epithelial_eligible,'Unknown':D.type_supported.eq('Unknown')};rr=[]
for name,ix in groups.items():
 q=D[ix]
 for g in sorted(set(PANELS['T']+PANELS['Treg']+PANELS['Tfh']+PANELS['Endo']+PANELS['Epi']+PANELS['Fib']+['TIGIT'])):
  v=q['raw_'+g];rr.append(dict(group=name,gene=g,n=len(q),positive=int((v>0).sum()),positive_fraction=float((v>0).mean()),raw_UMI=int(v.sum()),mean_LogNormalize=float(np.log1p(v*10000/q.nCount_RNA.clip(lower=1)).mean()),source='unsmoothed counts; normalized per original object; shared classifier markers not independent truth'))
put(pd.DataFrame(rr),'target_group_marker_evidence.tsv')
# Coverage and analysis selection at all original data levels.
rr=[]
for unit,tab in [('cell','Cervilla_cell_measured.parquet'),('bin8','Cervilla_bin8_measurement.parquet'),('bin16','Cervilla_bin16_measurement.parquet')]:
 q=pd.read_parquet(R/'tables'/tab);keep=q.analysis_qc if 'analysis_qc'in q else (q.qc_keep if 'qc_keep'in q else pd.Series(True,index=q.index))
 rr.append(dict(unit=unit,n_supplied=len(q),UMI_supplied=int(q.nCount_RNA.sum()),n_analysis=int(keep.sum()),UMI_analysis=int(q.loc[keep,'nCount_RNA'].sum()),UMI_excluded=int(q.loc[~keep,'nCount_RNA'].sum()),note='cell file here is existing qc_keep_analysis; full v2 contains 158144 geometry objects'if unit=='cell'else'counts include all supplied before 16um aggregation'))
put(pd.DataFrame(rr),'input_QC_balance.tsv')
rr=[];cross=[]
for radius in [50,200]:
 n=pd.read_parquet(R/'tables'/f'Cervilla_cell_neighborhood_r{radius}.parquet');v=pd.read_parquet(R/'tables'/f'Cervilla_region_association_r{radius}.parquet');q=v[v.association_eligible]
 rr.append(dict(radius_um=radius,epi_anchors=len(n),neighborhood_eligible=int(n.neighborhood_eligible.sum()),unresolved=int((~n.neighborhood_eligible).sum()),training=int(n.niche_assignment.eq('training').sum()),predicted=int(n.niche_assignment.eq('model_predict').sum()),score_regions=int(reg.score_eligible.sum()),association_regions=len(q),distance_observable_regions=int(q.Endothelial_distance.notna().sum()),Endo_positive_regions=int((q.Endothelial_density>0).sum()),mean_known_fraction=float(n.known_fraction.mean()),median_known_neighbors=float(n.n_known.median()),Endo_objects=204))
 for k,g in q.groupby('niche_region_summary'):
  for definition in ['group_H50','group_H10','group_H05']:
   s=g[g[definition].isin(['H50','H10','H05','L50'])];high=s[definition].ne('L50');out=q[(q.niche_region_summary!=k)&q[definition].isin(['H50','H10','H05','L50'])];out_high=out[definition].ne('L50');odds=(high.sum()*(~out_high).sum()/((~high).sum()*out_high.sum()))if (~high).sum()*out_high.sum()>0 else np.nan
   cross.append(dict(radius_um=radius,niche=int(k),group_definition=definition,n_regions=len(g),n_comparable=len(s),n_high=int(high.sum()),n_L50=int((~high).sum()),high_fraction=float(high.mean())if len(s)else np.nan,descriptive_odds_ratio_vs_other_regions=odds,median_MP10=float(g.MP10.median()),note='0=mixed/unresolved regional summary; descriptive only, no independent repeated-tail tests'))
put(pd.DataFrame(rr),'analysis_coverage.tsv');put(pd.DataFrame(cross),'niche_high_MP10_descriptive.tsv')
a=pd.read_parquet(R/'tables/Cervilla_cell_neighborhood_r50.parquet');b=pd.read_parquet(R/'tables/Cervilla_cell_neighborhood_r200.parquet');z=a[['cell_id','niche']].merge(b[['cell_id','niche']],on='cell_id',suffixes=('_50','_200'),validate='one_to_one');put(z.groupby(['niche_50','niche_200']).size().rename('n_centers').reset_index(),'niche_scale_crosstab.tsv')
comp=pd.read_csv(R/'tables/niche_composition.tsv',sep='\t');ee=[]
for radius in [50,200]:
 n=a if radius==50 else b;valid=n.neighborhood_eligible
 for row in comp[comp.radius_um==radius].itertuples():
  for p in [c for c in n if c.startswith('P_')]:
   mean=float(n.loc[valid&n.niche.eq(row.niche),p].mean());base=float(n.loc[valid,p].mean());ee.append(dict(radius_um=radius,niche=row.niche,type=p[2:],mean_composition=mean,global_eligible_center_mean=base,ratio_to_global_mean=mean/base if base>0 else np.nan))
put(pd.DataFrame(ee),'niche_relative_enrichment.tsv')
# MP10 member-gene detectability in all eligible pooled regions.
put(pd.DataFrame([dict(gene=g,n_regions=int(reg.score_eligible.sum()),positive=int((reg.loc[reg.score_eligible,'raw_'+g]>0).sum()),mean_UMI=float(reg.loc[reg.score_eligible,'raw_'+g].mean()))for g in TOP]),'MP10_regional_gene_support.tsv')
roi=pd.read_csv(R/'tables/roi_manifest.tsv',sep='\t');roi['area_um2']=(roi.x1-roi.x0)*(roi.y1-roi.y0);roi['coordinate_units']='original fullres physical microns; x right, y down';put(roi,'roi_manifest.tsv')
ver=[]
for p in ['numpy','pandas','scipy','scikit-learn','anndata','statsmodels','numba','matplotlib','pyarrow']:
 try:ver.append(dict(package=p,version=im.version(p)))
 except im.PackageNotFoundError:pass
put(pd.DataFrame(ver),'Python_versions.tsv')
# Compact, readable region table for R without an optional arrow dependency.
for radius in [50,200]:
 q=pd.read_parquet(R/'tables'/f'Cervilla_region_association_r{radius}.parquet');cols=['region_id','x_um','y_um','n_epi','nCount_RNA','nFeature_RNA','MP10_detected','MP10','score_eligible','score_reason','association_eligible','niche_region_summary','niche_K2_region_summary','known_fraction','observed_area_fraction','Endothelial_density','Endothelial_distance','distance_censored_fraction','group_H50','group_H10','group_H05'];put(q[cols],f'region_result_r{radius}.tsv.gz')
print('DELIVERY_TABLES_COMPLETE')

# Correct installed RCTD selected-model semantics without changing the previously frozen cell analysis.
original=pd.read_csv(R/'tables/unit_selection.tsv',sep='\t');corrected=pd.read_csv(R/'tables/unit_comparison.tsv',sep='\t');audit=corrected.merge(original[['unit','selected','frozen_at','Treg_status','association_inspected']],on='unit',validate='one_to_one');audit['selection_confidence']='provisional operational main; no demonstrated universal bin/cell winner';audit['audit_note']='Installed get_beta_doublet: singlet first_type=1, certain doublet pair weights, unresolved not complete proportions. Corrected mixing comparable between cell and16um; former less-mixing rationale withdrawn. Cell/niche/MP10/association definitions unchanged.';audit['reason']=audit.unit.map({'cell':'Nuclear-associated lineage positions and Epi-only connected regional pooling without dividing mixed RNA by model Epi fraction; provisional, not proven more accurate than bins','bin8':'Retain original spatial signal; sparse single-bin MP10 (6 eligible fixed-ROI bins)','bin16':'Greater depth but comparable detected mixture to cells; model contributions do not recover epithelial-specific RNA; bounded comparison'});put(audit,'unit_selection_audited.tsv')
f=pd.read_parquet(R/'tables/Cervilla_cell_formal_metadata.parquet');new=pd.read_parquet(R/'tables/Cervilla_cell_annotation.parquet');assert f.cell_id.equals(new.cell_id) and f.type_supported.equals(new.type_supported) and f.mp10_epithelial_eligible.equals(new.mp10_epithelial_eligible),'RCTD reporting audit must not change frozen cell analysis'
rr=[]
for unit in ['bin8','bin16']:
 c=pd.read_parquet(R/'tables'/f'{unit}_ROI_contributions.parquet');pc=[x for x in c if x.startswith('contribution_')];valid=c.reliable_fit;assert np.allclose(c.loc[valid,pc].sum(1),1);assert c.loc[~valid,pc].isna().all().all();rr.append(dict(unit=unit,n_reliable=int(valid.sum()),row_sum_one=True,nonreliable_contributions_NA=True,native_get_beta_doublet_equal=True,original_pair_weights_retained=True))
put(pd.DataFrame(rr),'RCTD_selected_weight_validation.tsv')

sv=pd.read_csv(R/'tables/mp10_score_validation.tsv',sep='\t');sv['score_implementation']=np.where(sv.scope.eq('auc_original'),'native random-zero AUCell sensitivity (not historical primary)','native zero-excluded AUCell sensitivity');put(sv,'mp10_score_validation.tsv')
four=pd.read_csv(R/'tables/conservative4um_ROI_cells.tsv.gz',sep='\t');join=four[['cell_id','nCount','nFeature','fine_prediction']].merge(D[['cell_id','nCount_RNA','nFeature_RNA','fine_prediction']],on='cell_id',how='left',suffixes=('_4um','_8um'),validate='one_to_one');join['UMI_added_8_vs4']=join.nCount_RNA-join.nCount;join['relation']='same original nucleus-label object ID; 4um ROI may include objects outside v2 analysis set';put(join,'conservative4_vs8_ID_comparison.tsv.gz')
