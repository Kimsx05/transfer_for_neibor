from common import *
import html,zipfile,shutil,subprocess

def synthesis():
 p=pairs();g=pd.read_csv(OUT/'00_approval/approved_gate_v7.tsv',sep='\t');bg=pd.read_csv(OUT/'03_ring_statistics/ring123_background_tests.tsv',sep='\t');rp=pd.read_csv(OUT/'03_ring_statistics/ring123_pairwise_tests.tsv',sep='\t');mp=pd.read_csv(OUT/'04_MP10_association/MP10_scope_associations.tsv',sep='\t');bg=bg[bg.unit_subset.eq('ALL_POLICY_UNITS')];rp=rp[rp.unit_subset.eq('ALL_POLICY_UNITS')];mp=mp[mp.unit_subset.eq('ALL_POLICY_UNITS')&mp.method.eq('ADJUSTED_MAIN')];b=bg.set_index(['scope','pair_id','ring']);r=rp.set_index(['scope','pair_id','ring_pair']);a=mp.set_index(['scope','query_id','ring','score_layer']);sens=pd.read_csv(OUT/'05_sensitivity/TME_conditional_sensitivity.tsv.gz',sep='\t').set_index(['analysis_set','scope','pair_id','ring']);rows=[]
 for scope in SCOPES:
  for pr in p.itertuples():
   old=g[g.scope.eq(scope)&g.master_candidate_id.eq(pr.master_candidate_id)].iloc[0];row=dict(scope=scope,master_candidate_id=pr.master_candidate_id,candidate_id=old.candidate_id,pair_id=pr.pair_id,niche=pr.niche,celltype=pr.celltype,review_tier=pr.review_tier,v6_unit_effect_pp=old.unit_effect_pp,v6_unit_q=old.unit_q_BH,v6_paper_q=old.paper_q_BH,v6_source_limitation=old.source_limitation,composition_figure='../06_figures/'+scope+'/composition/'+pr.pair_id+'_excess.png',within_MP10_figure='../06_figures/'+scope+'/MP10/'+pr.pair_id+'_association.png',ALL_MP10_figure='../06_figures/'+scope+f'/MP10/all_c{pr.cell_index:02d}_association.png');tags=['COMPOSITION_CHARACTERISTIC'];denom=False
   for ring in RINGS:
    z=b.loc[(scope,pr.pair_id,ring)];row.update({ring+'_excess_pp':100*z.median_effect,ring+'_excess_mean_pp':100*z.mean_effect,ring+'_background_q':z.background_zero_q,ring+'_n_units':z.n_units,ring+'_n_verified_units':z.n_verified_units,ring+'_n_other_units':z.n_other_units,ring+'_n_positive':z.n_positive,ring+'_n_negative':z.n_negative,ring+'_status':z.status})
    if z.background_zero_q<.05:tags.append('EXCESS_ABOVE_SAMPLE_BACKGROUND' if z.median_effect>0 else 'DEPLETION_BELOW_SAMPLE_BACKGROUND')
    if z.n_units<5:tags.append('INSUFFICIENT_COVERAGE')
    sm=sens.loc[('TME_MATCHED_RAW',scope,pr.pair_id,ring)];st=sens.loc[('TME_CONDITIONAL',scope,pr.pair_id,ring)];row[ring+'_TME_matched_raw_excess']=sm.median_excess;row[ring+'_TME_conditional_excess']=st.median_excess;row[ring+'_TME_n_units']=st.n_units
    if np.isfinite(sm.median_excess) and np.isfinite(st.median_excess) and np.sign(sm.median_excess)!=np.sign(st.median_excess):denom=True
    for range_name,query in [('ALL',f'all_c{pr.cell_index:02d}'),('WITHIN',pr.pair_id)]:
     vals=[];qs=[]
     for layer,short in zip(LAYERS,['spot','model_epi']):
      z=a.loc[(scope,query,ring,layer)];pre=range_name+'_'+ring+'_'+short;row.update({pre+'_rho':z.median_effect,pre+'_mean_rho':z.mean_effect,pre+'_IQR':z.IQR,pre+'_q':z.association_q,pre+'_n_units':z.n_units,pre+'_n_positive':z.n_positive,pre+'_n_negative':z.n_negative,pre+'_status':z.status});vals.append(z.median_effect);qs.append(z.association_q)
      if z.association_q<.05:tags.append('MP10_POSITIVE_ASSOCIATION_EXPLORATORY' if z.median_effect>0 else 'MP10_NEGATIVE_ASSOCIATION_EXPLORATORY')
     finite=np.isfinite(vals).all();concordant=bool(finite and np.sign(vals[0])==np.sign(vals[1]));row[range_name+'_'+ring+'_two_layer_same_sign']=concordant if finite else None;row[range_name+'_'+ring+'_both_layers_q_lt05']=bool(np.all(np.array(qs)<.05));row[range_name+'_'+ring+'_one_layer_q_lt05']=bool(np.sum(np.array(qs)<.05)==1)
     if concordant:tags.append('TWO_LAYER_DIRECTION_CONCORDANT')
     if finite and np.sign(vals[0])*np.sign(vals[1])<0:tags.append('OPPOSITE_SCORE_LAYER_DIRECTIONS')
     if np.sum(np.array(qs)<.05)==1:tags.append('ONE_LAYER_STATISTICAL_SUPPORT_ONLY')
   for pair in ['R1-R2','R1-R3','R2-R3']:
    z=r.loc[(scope,pr.pair_id,pair)];row[pair+'_excess_difference_pp']=100*z.median_effect;row[pair+'_q']=z.ring_pairwise_q;row[pair+'_n_units']=z.n_units
   if denom:tags.append('TME_DENOMINATOR_SENSITIVE')
   ds=bg[bg.scope.ne('ALL_DATASETS')&bg.pair_id.eq(pr.pair_id)&bg.n_units.gt(0)]
   heterogeneous=any((dd.median_effect.gt(0).any() and dd.median_effect.lt(0).any()) for _,dd in ds.groupby('ring'))
   if heterogeneous:tags.append('DATASET_HETEROGENEOUS')
   row['interpretation_tags']=';'.join(dict.fromkeys(tags));row['TME_direction_change_same_center_set']=denom;row['dataset_direction_heterogeneous']=heterogeneous;rows.append(row)
 e=pd.DataFrame(rows)
 for col in ['composition_figure','within_MP10_figure','ALL_MP10_figure']:
  e[col+'_in_dataset_package']=e[col]
  non=e.scope.ne('ALL_DATASETS');e.loc[non,col]=['../reports/SCOPE_'+sc+'.html#'+pid for sc,pid in zip(e.loc[non,'scope'],e.loc[non,'pair_id'])]
 tab('07_synthesis/candidate_evidence_139.tsv',e[e.scope.eq('ALL_DATASETS')]);tab('07_synthesis/candidate_dataset_evidence.tsv',e[e.scope.ne('ALL_DATASETS')]);tab('07_synthesis/all_scope_evidence.tsv',e)
 # Numeric dataset support for composition and both MP10 ranges retained independently.
 tab('07_synthesis/composition_dataset_directions.tsv',bg);tab('07_synthesis/MP10_dataset_directions.tsv',mp)
 fam=[]
 fullbg=pd.read_csv(OUT/'03_ring_statistics/ring123_background_tests.tsv',sep='\t');fullrp=pd.read_csv(OUT/'03_ring_statistics/ring123_pairwise_tests.tsv',sep='\t');fullmp=pd.read_csv(OUT/'04_MP10_association/MP10_scope_associations.tsv',sep='\t')
 for subset in ['ALL_POLICY_UNITS','VERIFIED_PATIENTS_ONLY']:
  for scope in SCOPES:
   for famname,df,n,qcol in [('BACKGROUND_ZERO',fullbg,417,'background_zero_q'),('RING_PAIRWISE',fullrp,417,'ring_pairwise_q')]:
    z=df[df.scope.eq(scope)&df.unit_subset.eq(subset)];fam.append(dict(unit_subset=subset,scope=scope,family=famname,planned_n=n,n_p_available=int(z.p.notna().sum()),q_field=qcol,seven_dataset_n=n*7,eight_scope_n=n*8))
   for method in ['ADJUSTED_MAIN','RAW','ADJUSTED_CENTER_TYPE']:
    z=fullmp[fullmp.scope.eq(scope)&fullmp.method.eq(method)&fullmp.unit_subset.eq(subset)];fam.append(dict(unit_subset=subset,scope=scope,family='MP10_'+method,planned_n=1008,n_p_available=z.p.notna().sum(),q_field='association_q',seven_dataset_n=7056,eight_scope_n=8064))
 tab('00_audit/family_manifest.tsv',fam)
 m,q,center=load();cov=pd.read_csv(V2/'00_preflight/coverage63.tsv',sep='\t');summary=[]
 for scope in SCOPES:
  z=m if scope=='ALL_DATASETS' else m[m.dataset_id.eq(scope)];zmain=z[z.main_eligible];zs=z[z.score_eligible];bt=bg[bg.scope.eq(scope)];mt=mp[mp.scope.eq(scope)];counts={}
  for layer,short in zip(LAYERS,['spot','model_epi']):
   for ran,rs in [('ALL_ELIGIBLE_EPI','ALL'),('WITHIN_SELECTED_NICHE','WITHIN')]:
    t=mt[mt.score_layer.eq(layer)&mt.range.eq(ran)];counts[rs+'_'+short+'_positive_q05']=int((t.association_q.lt(.05)&t.median_effect.gt(0)).sum());counts[rs+'_'+short+'_negative_q05']=int((t.association_q.lt(.05)&t.median_effect.lt(0)).sum());counts[rs+'_'+short+'_NA']=int(t.p.isna().sum())
  summary.append(dict(scope=scope,manifest_sections=len(cov) if scope=='ALL_DATASETS' else int(cov.dataset_id.eq(scope).sum()),main_sections=zmain.section_id.nunique(),accepted_anchors=int(z.accepted_label.ge(0).sum()),main_R123_anchors=len(zmain),complete_R123_anchors=int((z.main_eligible&z.complete123).sum()),TME_sensitivity_anchors=int(z.TME_eligible.sum()),scored_anchors=len(zs),main_units=zmain.analysis_unit_uid.nunique(),verified_units=zmain.loc[zmain.source_patient_identity_known,'analysis_unit_uid'].nunique(),other_units=zmain.loc[~zmain.source_patient_identity_known,'analysis_unit_uid'].nunique(),background_positive_q05=int((bt.background_zero_q.lt(.05)&bt.median_effect.gt(0)).sum()),background_negative_q05=int((bt.background_zero_q.lt(.05)&bt.median_effect.lt(0)).sum()),background_NA=int(bt.p.isna().sum()),ring_pair_significant=int(rp[rp.scope.eq(scope)].ring_pairwise_q.lt(.05).sum()),**counts))
 summary=pd.DataFrame(summary);tab('07_synthesis/scope_summary.tsv',summary);tab('07_synthesis/tier_completion.tsv',[dict(scope=sc,review_tier=t,n_pairs=int(p.review_tier.eq(t).sum()),planned_background_tests=int(p.review_tier.eq(t).sum()*3),status='ALL_NUMERIC_RESULTS_PRESENT_INCLUDING_NA') for sc in SCOPES for t in TIERS])
 js('00_audit/software_versions.json',dict(python=sys.version,numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,R=subprocess.check_output(['Rscript','--version'],text=True,stderr=subprocess.STDOUT).strip()))
 checkpoint('M05_SENSITIVITY_DONE');return e,summary

METHODS='''# v7 Methods
Frozen SAMPLE_R12_LEIDEN accepted 0–8; mean-Epithelial/all30 gate >=0.25 inherited. Labels, scores and v2–v6 outputs are read-only. Matrix identity is secured by original COMPLETE hashes and explicit spot_uid joins. The alias in the v7 attachment and formal v3 path have identical UID/accepted labels.
Full tissue physical R1/R2 graphs retain missing-composition nodes. Exact R3 is reachability via R2×R1, excluding R1/R2/self; sparse, no all-pairs dense distances. Every section has independent BFS checks. Recomputed full q1/q2 match original float32 caches within 1e-6. Per-spot q05 normalization -> equal valid-neighbor ring means; missing rows excluded from denominator, never imputed as zero. R12 remains equal-ring average.
S_s: accepted common anchors with valid R1/R2/R3 neighbors >=4/8/12. All three rings use identical centers. Full 6/12/18 geometric and valid counts form the separate complete-rings sensitivity. Background is mean q over S_s within section, including target and niche0. e=target-background; rest excludes target; e=(1-target_weight)*(target-rest), checked numerically. Section target/rest each >=20 for inferential summaries. Raw/excess/ringpair boxes use the identical inferential section/unit set. Separate descriptive summaries retain smaller targets.
Same analysis_unit_uid has equal-section averages, not anchor weighting. Global units collapse across datasets. Source verified patients and existing anonymous/user-assumed policy units are counted separately. Background tests use unit e vs zero, ringpairs use paired e_R1-e_R2 etc. Two-sided signed-rank, Wilcox zero removal; >=5 units and >=5 nonzero values; all zero with >=5 units p=1. Exact if n<50 and no zeros/absolute ties, otherwise tie and continuity-corrected approximation. Numerical zero tolerance is 0, matching v6; no jitter added to values. R/Python fixed cases cross-checked. Effect output includes median/mean/IQR, not a claim of a mean test.
BH: per scope 139×3=417 background tests and separate 417 ringpairs. Fixed planned n even if some P missing. Across only seven datasets n=2919; all eight scopes n=3336. All tiers in same family. Main q used on figures. Four stars <0.0001, three <0.001, two <0.01, one <0.05, otherwise ns; unavailable N/A. Raw composition panel has no stars. Ringpair brackets only in their own panel.
S_score: S_s with complete Top20, epi_score_valid and both frozen model-universe spot/Epi scores finite. B5 remains composition-eligible but partial-signature scores are excluded. Same centers across two layers. Queries: 29 unique ALL_ELIGIBLE_EPI types (includes niche0) plus139 WITHIN_SELECTED_NICHE pairs. Each section/range >=30 complete centers, at least3 distinct score and composition values. Average-tie ranks. Main partial correlation: residualize rank(q),rank(score) against intercept,rank(center mean Epi),rank(log1p(center_nCount)); correlate residuals. Constant covariates removed and noted; pivoted QR rank audited, collinear design recorded NA. Main QR tolerance1e-10 relative; residual numerical variance tolerance1e-12. Corresponding central mean celltype rank is extra sensitivity only. V5 used unranked covariates and corresponding type in its limited diagnostic; v7 follows its new explicit design.
Score layers use identical usable sections for each range/ring/method; invalid in either layer removes that section in both. Unit correlation=tanh(mean(arctanh(clip(section_rho,-1+1e-7,1-1e-7)))); raw correlation/clipping retained. Unit signed-rank vs0. Per scope and method planned family(29+139)×3×2=1008; raw and center-type are independent auxiliary families; seven-dataset n7056, eight-scope n8064. No anchor P used. Verified patient subset uses same planned families with separate q.
Sensitivity: complete rings; verified patients; fixed TME share>=0.05 in all3 rings; normalize ring means by their nonEpi sum (not normalize individual neighbors). Recompute matched-center raw-q background in parallel with TME conditional background; both composition and association retain full planned rows. TME and complete/block effects are descriptive, no substituted main stars. Fixed 8×8 physical blocks reuse v5 array coordinate definition; select SHA256(seed123|spot_uid) minimum per block among main score-eligible centers, then restrict each query; common draw across types/rings/layers. Blocks are not independent patients or a full permutation test.
Figures: each composition file contains three separate adjacent-box panels raw/excess/ringpairs; each MP10 file contains two score-layer rows and primary/raw columns. Each panel has three unconnected R1/R2/R3 boxes. Median/Q1/Q3, whiskers1.5IQR; showfliers=False with all actual unit points plotted exactly once; verified circle, other triangle. N is test-unit count. PNG300dpi transparent and vector PDF. This layout consolidates files without omitting requested views.
Reference boundary: attachment-provided Khaliq et al. 2024 supports box appearance; its Fig3d unpaired t-test is not our paired signed-rank procedure. Rings, thresholds,139 selections and families are project adaptations. No new literature claims, model training, score recomputation or communication analysis.
'''
UNCERTAIN='''# What remains uncertain
- All139 pairs were chosen after inspecting composition-defined niches; new tests remain exploratory with selection effects. R1/R2 were clustering inputs; R3 remains spatially correlated, not an independent cohort.
- Overlapping neighborhoods and section-local dependence remain. Fixed sparse block direction checks do not eliminate spatial correlation. Unit independence follows inherited evidence tiers, including user-assumed or anonymous policies; these are not all confirmed patients.
- Relative q05 compositions have closure constraints; lower proportions do not prove lower absolute counts or active exclusion. TME-conditional changes use a different denominator and are descriptive sensitivity.
- Model-attributed Epi MP10 is not measured pure epithelial truth. v5 reported allocation error, low counts/rank ties, reference overlap, unpropagated model uncertainty and synthetic Epi RNA underallocation. Full spot/model Epi share the same centers and gene background but are not independent validation.
- Pathology remains unreviewed; the original Epi>=0.25 gate does not establish pure epithelium or malignancy. No new pathology approval was implied or required.
- ALL associations include between-niche environment variation; within-niche analyses restrict range and condition on composition. A nonsignificant within-niche result does not prove no relationship. ns is not equivalence or return to background.
- Same-sign score-layer point estimates, one layer with q<0.05, and both layers with q<0.05 are separate fields. No effect was required to match its old review tier.
- Sparse scopes and absent target/rest remain NA with actual coverage. Unit boxes summarize available inference-eligible sections, while separate descriptive tables retain smaller groups.
'''

def reports(e,summary):
 (OUT/'reports/METHODS.md').write_text(METHODS);(OUT/'reports/WHAT_IS_STILL_UNCERTAIN.md').write_text(UNCERTAIN)
 cols=['scope','main_sections','main_R123_anchors','complete_R123_anchors','scored_anchors','main_units','verified_units','other_units','background_positive_q05','background_negative_q05','background_NA','ring_pair_significant']
 text='# 阶段 B 实际结果\n\n状态：全部已授权139组合、总体及七dataset的数值结果与预定敏感性已完成；图像和包完整性见 validation_results.tsv。\n\n'+summary[cols].to_markdown(index=False)+'\n\n正向/负向与环间显著计数均按“组合×圈”或“组合×环对”计，不是独立患者或机制数。下面另列主调整连续 MP10 结果（按范围×类型/组合×圈计）。\n\n'+summary[['scope']+[c for c in summary if c.endswith(('positive_q05','negative_q05','_NA')) and not c.startswith('background')]].to_markdown(index=False)
 text+='\n\n28核心正向、32宽口径正向、39探索、40负向对照全量运行。29个全体范围关联只计算一次；139个群内分支单列。所有结果包括NA，不按显著性删除。\n\n新增附件两份审阅表已完整读取，139项与v7附表分层一致；两份approved TSV由附表及原v6真实ID生成，记录用户整体授权，不声称逐行人工手填。\n\n三个证据层独立：组成高/低于同片背景；圈间差异；连续MP10关联。任何一层不能替代其他层，更不能推出主动招募、排斥或机制。\n\n图像组织：每个组合有raw/excess/ringpairs三个面板，每面板三个并列箱体；每个关联图有两评分层×主调整/raw四个面板。无距离曲线，无单位点连线。\n'
 (OUT/'reports/RUN_SUMMARY.md').write_text(text)
 cmd='bash "/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring123_broad_v7/code/run_ring123_broad_v7.sh"'
 (OUT/'reports/RESUME.md').write_text('# 可用续跑及验收\n\n```bash\n'+cmd+' --stage all --approval "/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring123_broad_v7/00_approval/approved_gate_v7.tsv"\n'+cmd+' --verify-only\n```\n\n检查点绑定输入manifest和冻结spec哈希。完整性验证通过后跳过一致阶段；绘图只补缺失文件。命令实际执行全量三环与关联，不再调用v6 Gate-only接口。原始输入变化会报冲突，不能静默复用缓存。\n')
 # Portable main review index includes ALL figures + all scope overviews, not broken dataset full-image links.
 index=['# v7 阶段 B 审阅入口','', '[运行摘要](RUN_SUMMARY.md) · [方法](METHODS.md) · [限制](WHAT_IS_STILL_UNCERTAIN.md) · [续跑](RESUME.md)','', '[总体139行证据](../07_synthesis/candidate_evidence_139.tsv) · [七dataset 973行证据](../07_synthesis/candidate_dataset_evidence.tsv) · [覆盖与显著计数](../07_synthesis/scope_summary.tsv)','', '[全部总体箱线图索引](INDEX_ALL_DATASETS.html)','', '先看组成raw/excess/ringpairs，再看MP10主调整/原始两层关联；两类星号不能互借。按既定tier排序，不按本轮显著性筛图。','']
 for scope in SCOPES:
  index.append(f'- {scope}：[缩略总览](../06_figures/{scope}/overview.png) · [Epi/TME结构性背景](../06_figures/{scope}/structural_Epi_TME.png)')
 index+=['','七个dataset完整图包分别为 figures_<dataset>.zip。每包独立解压，打开其 reports/INDEX_<dataset>.html。主review包仅带总体全图及所有scope总览，避免重复千张图。','']
 (OUT/'reports/REVIEW_INDEX.md').write_text('\n'.join(index))
 for scope in SCOPES[1:]:
  h=['<!doctype html><meta charset="utf-8"><h1>'+scope+'</h1><p>Full boxplots are in figures_'+scope+'.zip. Extract that package separately and open reports/INDEX_'+scope+'.html; the exact filenames below identify each candidate. This landing page keeps review_B links portable.</p>', '<p><a href="../06_figures/'+scope+'/overview.png">Overview</a> | <a href="../06_figures/'+scope+'/structural_Epi_TME.png">Epi/TME structure</a></p>']
  for pr in pairs().itertuples():h.append('<h3 id="'+pr.pair_id+'">'+html.escape(pr.niche+' / '+pr.celltype)+'</h3><p>Composition: '+pr.pair_id+'_excess.png; within MP10: '+pr.pair_id+'_association.png; ALL MP10: all_c'+f'{pr.cell_index:02d}'+'_association.png</p>')
  (OUT/'reports'/('SCOPE_'+scope+'.html')).write_text('\n'.join(h))
 for scope in SCOPES:
  mapping=pd.read_csv(OUT/'06_figures'/scope/'display_name_map.tsv',sep='\t');h=['<!doctype html><html><head><meta charset="utf-8"><title>'+scope+' v7 review</title><style>body{font-family:Arial;max-width:1200px;margin:30px auto}td,th{padding:6px;border-bottom:1px solid #ddd}table{border-collapse:collapse}a{color:#175b87}</style></head><body>',f'<h1>{scope} — all authorized results</h1>','<p>Composition: separate raw/background-zero/ring-pair panels. MP10: primary/raw × whole-spot/model-attributed Epi panels. N/A retained. No connected distance curves.</p>',f'<p><a href="../06_figures/{scope}/overview.png">Overview</a> | <a href="../06_figures/{scope}/structural_Epi_TME.png">Epi/TME structure</a></p>','<table><tr><th>Tier</th><th>Range / niche</th><th>Celltype</th><th>Kind</th><th>PNG</th><th>PDF</th></tr>']
  order={t:i for i,t in enumerate(TIERS+['SHARED_ALL_RANGE'])};mapping['order']=mapping.review_tier.map(order)
  for r in mapping.sort_values(['order','niche','celltype','kind']).itertuples():h.append('<tr><td>'+html.escape(r.review_tier)+'</td><td>'+html.escape(r.niche)+'</td><td>'+html.escape(r.celltype)+'</td><td>'+r.kind+'</td><td><a href="../'+r.png+'">PNG</a></td><td><a href="../'+r.pdf+'">PDF</a></td></tr>')
  h.append('</table></body></html>');(OUT/'reports'/('INDEX_'+scope+'.html')).write_text('\n'.join(h))


def package():
 main=OUT/'reports/review_B.zip'
 with zipfile.ZipFile(main,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for folder in ['00_approval','00_audit','02_background_composition','03_ring_statistics','04_MP10_association','05_sensitivity','07_synthesis','code']:
   for p in sorted((OUT/folder).rglob('*')):
    if p.is_file() and 'section_cache' not in str(p) and '__pycache__' not in str(p) and p.suffix!='.xlsx':z.write(p,p.relative_to(OUT))
  for p in [OUT/'01_ring123/coverage_all_steps.tsv',OUT/'01_ring123/graph_and_cache_audit.tsv',OUT/'01_ring123/cache_hashes.tsv']:z.write(p,p.relative_to(OUT))
  for p in sorted((OUT/'reports').iterdir()):
   if p.is_file() and p.name!='delivery_manifest.tsv' and (p.suffix in ['.md','.json','.tsv'] or p.name=='INDEX_ALL_DATASETS.html' or p.name.startswith('SCOPE_')):z.write(p,p.relative_to(OUT))
  for scope in SCOPES:
   for p in sorted((OUT/'06_figures'/scope).rglob('*')):
    if p.is_file() and (scope=='ALL_DATASETS' or p.stem in ['overview','structural_Epi_TME']):z.write(p,p.relative_to(OUT))
 for scope in SCOPES[1:]:
  with zipfile.ZipFile(OUT/'reports'/('figures_'+scope+'.zip'),'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
   for p in sorted((OUT/'06_figures'/scope).rglob('*')):
    if p.is_file():z.write(p,p.relative_to(OUT))
   p=OUT/'reports'/('INDEX_'+scope+'.html');z.write(p,p.relative_to(OUT))
  log('packaged '+scope)

if __name__=='__main__':
 e,s=synthesis();reports(e,s);package()
