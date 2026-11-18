from fu_common import *
from analysis_stats import fit_effect
if (R/'frozen_config.yaml').exists():print('ALREADY FROZEN');raise SystemExit
checks=[]
for path in ['code/08b_common_cell_comparison.py','frozen_config.yaml','tables/association_primary.tsv']:
 a=(V/path).read_text().strip();b=(R/'sources'/('baseline_'+path.replace('/','__'))).read_text().strip();checks.append(dict(path=path,baseline_commit='e16e986fc776329cbb6f77b686131d1aecc98c42',exact_text_matches=a==b));assert a==b
put(pd.DataFrame(checks),'baseline_version_comparison.tsv')
base=pd.read_csv(V/'tables/input_manifest.tsv',sep='\t');put(base,'input_manifest.tsv');put(pd.read_csv(V/'tables/association_primary.tsv',sep='\t'),'v1_five_section_main.tsv')
b=pd.read_parquet(V/'Cervilla/T32_neighborhoods.parquet');c=pd.read_parquet(V/'Cervilla/cell_region_neighborhood_comparison.parquet');rows=[];design=[]
for rad in [50,200]:
 bb=b[b.radius_um.eq(rad)&b.source_rule.eq('exclude_shared')&b.association_eligible].set_index('T_region_id');cc=c[c.radius_um.eq(rad)&c.association_eligible].set_index('T_region_id');ids=bb.index.intersection(cc.index);d=bb.loc[ids].copy();d['MP10_bin']=d.Epi_MP10;d['MP10_cell']=cc.loc[ids,'Epi_MP10'];d['depth_bin']=d.Epi_nCount;d['depth_cell']=cc.loc[ids,'Epi_nCount'];d['n_Epi_cell']=cc.loc[ids,'Epi_n'];d['n_Epi_bin']=d.Epi_n;d['radius_um']=rad
 for f in ['Treg_main','T_lineage','T_nCount','x_um','y_um','component','observed_circle_fraction','T_background_balance','Tfh','Activated_T']:assert np.allclose(d[f],cc.loc[ids,f],equal_nan=True),f
 d.reset_index().to_parquet(R/'tables'/f'Cervilla_common_centers_r{rad}.parquet',index=False)
 for name,v,cov in [('A_bin_main_full',bb,None),('B_common_bin_full',bb.loc[ids],None),('v1_common_bin_named',bb.loc[ids],['T_lineage','T_nCount','Epi_nCount','observed_circle_fraction','T_background_balance']),('v1_common_cell_named',cc.loc[ids],['T_lineage','T_nCount','Epi_nCount','observed_circle_fraction','T_background_balance'])]:
  q=fit_effect(v,'Treg_main',block_um=400,covars=cov);rows.append(dict(comparison=name,radius_um=rad,**q))
 # Freeze common matrix before any new effect is inspected. Rank each source depth once; same rows used for both outcomes.
 for version,extra in [('common',[]),('competition',['Tfh','Activated_T'])]:
  cov=['T_lineage','T_nCount','observed_circle_fraction','T_background_balance','depth_bin','depth_cell']+extra
  q=d.dropna(subset=cov+['Treg_main','MP10_bin','MP10_cell']);C=[np.ones(len(q))];names=['intercept'];decisions=[]
  for f in cov:
   if q[f].nunique()<3:decisions.append(dict(column=f,retained=False,reason='constant_or_less_than3_values'));continue
   z=zrank(q[f]);X=np.column_stack(C);rem=z-X@np.linalg.lstsq(X,z,rcond=None)[0];keep=np.var(rem)>.1;decisions.append(dict(column=f,retained=keep,residual_variance=float(np.var(rem)),reason='fixed_order_VIF_less10'if keep else 'collinear_VIF_ge10'))
   if keep:C.append(z);names.append(f)
  X=np.column_stack(C+[zrank(q.Treg_main)]);names+=['Treg_main'];out=pd.DataFrame(X,columns=names,index=q.index).reset_index();out.to_parquet(R/'tables'/f'design_r{rad}_{version}.parquet',index=False)
  design.extend([dict(radius_um=rad,version=version,**z)for z in decisions])
put(pd.DataFrame(rows),'baseline_reproduction.tsv');put(pd.DataFrame(design),'common_design_decisions.tsv')
old=pd.read_csv(V/'tables/cell_bin_common_effects.tsv',sep='\t');new=pd.DataFrame(rows)
for rad in [50,200]:
 for unit in ['bin','cell']:
  v=new[(new.radius_um==rad)&new.comparison.eq('v1_common_'+unit+'_named')].iloc[0];w=old[(old.radius_um==rad)&old.block_um.eq(400)&old.unit.str.startswith(unit)].iloc[0];assert v.n==w.n and abs(v.effect-w.effect)<1e-12
roi=pd.read_csv(V/'tables/roi_manifest.tsv',sep='\t');assert len(roi)==24 and roi.roi.str.startswith('candidate').sum()==14;put(roi,'roi_manifest_original.tsv');roi['center_x']=(roi.x0+roi.x1)/2;roi['center_y']=(roi.y0+roi.y1)/2;roi['display_x0']=roi.center_x-256;roi['display_x1']=roi.center_x+256;roi['display_y0']=roi.center_y-256;roi['display_y1']=roi.center_y+256;put(roi,'roi_manifest_display512.tsv')
F={'version':'hd_treg_mp10_spatial_followup_v2','frozen_at':stamp(),'status':'post_v1_exploratory_followup_not_independent_validation','sections':SECTIONS,'seed':123,'threads_total':4,'new_primary_test':'Cervilla200_common_design_delta_cell_minus_bin_zero','radii_um':[50,200],'block_um':[400,800],'bootstrap_reps':999,'rank_transform':'average ranks, population standardization fixed before resampling','paired_bootstrap':'same block draws/indices both outcomes','null':'restricted common-design model of zrank(cell)-zrank(bin); HC2 block Rademacher wild bootstrap','primary_covariates':['T_lineage','T_nCount','observed_circle_fraction','T_background_balance','depth_bin','depth_cell'],'competition_extra':['Tfh','Activated_T'],'missing':'complete case jointly for both outcomes and common covariates','collinearity':'one frozen common design, ordered nuisance residual variance>.1; no effect-driven deletion','area':'not comparable; omitted both; no unknown-to-zero','hotspot_exclusion':'T centroids inside union of original Cervilla candidate256um ROI rectangles; keep remaining neighbor outcomes fixed','no_new_5section_continuous_associations':True,'signature':SIG,'maxRank':1300,'old_threshold':cfg['Treg_like_threshold_training_competitor_q95'],'mixture':{'partitions':'original train versus heldout sample; never mix reference cell between partitions','backgrounds':[{'name':'Epi40','Epi':.4,'other':None,'other_fraction':0.},{'name':'Epi70','Epi':.7,'other':None,'other_fraction':0.},{'name':'Epi90','Epi':.9,'other':None,'other_fraction':0.},{'name':'Epi40_Fib30','Epi':.4,'other':'Fibroblast','other_fraction':.3},{'name':'Epi40_Myeloid30','Epi':.4,'other':'Myeloid','other_fraction':.3}],'competitors':['CD4 Helper/Conventional','Tfh/CXCL13-like','CD8 Exhausted'],'Treg_T_RNA_fractions':[0,.1,1/3,1.],'depths':[500,2000,5000],'replicates':60,'cells_per_profile':5,'threshold':'one pooled common threshold; maximum of background/competitor/depth train-negative95th quantiles among multigene/core eligible negatives; old joint conditions retained','target_train_FPR_max_each_condition':.05,'heldout_acceptance':'all adequately sampled negative conditions empirical FPR<=.05 and Wilson upper<=.15; gain sensitivity>=.05 overall positives, no condition loss>.02; no refit to holdout','applicability':'500-5000 total UMI; nFeature>=325; other old measurement and joint conditions unchanged; out-of-range uncalibrated','uncertainty':'Wilson intervals describe simulated mixtures conditional on finite sample pools; not patients'}}
p=R/'frozen_config.yaml'
if p.exists():
 oldF=yaml.safe_load(p.read_text());F['frozen_at']=oldF['frozen_at'];assert F==oldF,'Definition changed: create an explicit new version'
else:p.write_text(yaml.safe_dump(F,sort_keys=False))
print('FROZEN',F['frozen_at'])
