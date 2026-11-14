import pathlib,json,hashlib,datetime,re,html,sys
import pandas as pd
P=pathlib.Path.cwd();B=P/'results/mp10_hypothesis_followup_v1';O=P/'results/mp10_immune_endothelial_mechanisms_v1';R=P/'resources/mp10_regulatory_priors'
def put(d,n): pd.DataFrame(d).to_csv(B/'tables'/n,sep='\t',index=False)
def old(n): return pd.read_csv(O/'tables'/n,sep='\t')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
 return h.hexdigest()
if (B/'config/M0_FROZEN.sha256').exists():
 assert sha(B/'config/analysis_plan.yaml')==(B/'config/M0_FROZEN.sha256').read_text().strip(), 'Frozen M0 config identity conflict'
 print('Existing frozen M0 preserved; use resume.py, not a new freeze.')
 sys.exit(0)
M20=(P/'neighborhood_mp10_calibration_v5/01_score_fidelity/top20.txt').read_text().split()
sets={};sources={}
h=pd.read_csv(R/'hallmark.tsv.gz',sep='\t');react=pd.read_csv(R/'reactome.tsv.gz',sep='\t')
for name in ['MTORC1_SIGNALING','OXIDATIVE_PHOSPHORYLATION','G2M_CHECKPOINT','MYC_TARGETS_V1','TNFA_SIGNALING_VIA_NFKB','ANGIOGENESIS']:
 sets[name]=sorted(set(h.loc[h.gs_name=='HALLMARK_'+name,'gene_symbol']));sources[name]='MSigDB Hallmark 2026.1.Hs'
sets['FA_SYNTH_DESAT']='ACACA ACACB FASN SCD SCD5 FADS1 FADS2 ELOVL6'.split();sources['FA_SYNTH_DESAT']='UniProt reviewed human enzyme functions; CAN-12-1329; predefined synthesis/desaturation/elongation panel'
sets['STEROL_SYNTH']=sorted(set(react.loc[react.gs_name=='REACTOME_CHOLESTEROL_BIOSYNTHESIS','gene_symbol']));sources['STEROL_SYNTH']='Reactome CHOLESTEROL_BIOSYNTHESIS / MSigDB 2026.1.Hs; includes shared mevalonate upstream'
sets['LIPID_UPTAKE_TRANSPORT']=sorted(set(react.loc[react.gs_name.isin(['REACTOME_TRANSPORT_OF_FATTY_ACIDS','REACTOME_PLASMA_LIPOPROTEIN_CLEARANCE']),'gene_symbol'])|{'CD36'});sources['LIPID_UPTAKE_TRANSPORT']='Reactome transport of fatty acids + plasma lipoprotein clearance 2026.1.Hs; CD36 from doi:10.1038/s41590-019-0589-5; not lipid export'
sets['PRENYLATION']='FNTA FNTB PGGT1B RABGGTA RABGGTB CHM CHML GGPS1 FDPS ICMT RCE1'.split();sources['PRENYLATION']='UniProt reviewed enzyme functions; doi:10.1038/s41586-021-03235-6; precursor/prenylation/postprocessing panel'
sets['TREG_IMMUNOREG']='CTLA4 ICOS IL2RA PDCD1 NT5E TIGIT FGL2'.split();sources['TREG_IMMUNOREG']='Literature-defined transcription panel: nature12297 (CTLA4 ICOS IL2RA context); s41586-018-0846-z (PDCD1 NT5E TIGIT FGL2); not validated suppression assay'
sets['MP10']=M20;sources['MP10']='Frozen project Top20; score not recomputed'
rows=[]
for name,genes in sets.items():
 for mask in ['Full','M20']:
  gg=sorted(set(genes)-(set(M20) if mask=='M20' else set()))
  for g in gg: rows.append(dict(signature=name,mask=mask,gene=g,n_expected=len(gg),source=sources[name],shared_M20=g in M20,shared_mevalonate=g in 'ACAT2 HMGCS1 HMGCR MVK PMVK MVD IDI1 IDI2 FDPS GGPS1'.split(),direction='higher transcript program; not metabolic flux',score_min_genes=5,score_min_fraction=.5))
put(rows,'signature_registry.tsv')
# Registry is built before any newly estimated association.
X=['MP10_Q50','SREBF1_M20','SREBF2_M20','MTORC1_SIGNALING_M20','NECTIN2','FA_SYNTH_DESAT_M20','STEROL_SYNTH_M20','LIPID_UPTAKE_TRANSPORT_M20','PRENYLATION_M20']
models=[]
def add(fam,x,y,ct='Treg',method='HC3',cohort='own',role='primary'):
 models.append(dict(FDR_family=fam,x=x,y=y,celltype=ct,method=method,cohort=cohort,role=role,setting='main',min_cells=50,patient='main',omit_dataset='none',covariate='none'))
for y in X[1:]: add('epi_continuous_primary','MP10_Q50',y,'Epi_all')
for x in X[1:3]:add('epi_continuous_primary',x,'NECTIN2','Epi_all')
for x in X:
 add('bridge_primary',x,'OXIDATIVE_PHOSPHORYLATION_Full');add('bridge_primary',x,'cycling_fraction','Treg','quasibinomial')
 add('bridge_noncycling',x,'OXIDATIVE_PHOSPHORYLATION_Full','Treg_NC')
 add('bridge_matched_all',x,'OXIDATIVE_PHOSPHORYLATION_Full','Treg','HC3','NC_matched','control')
 for y in ['FA_SYNTH_DESAT','STEROL_SYNTH','LIPID_UPTAKE_TRANSPORT','PRENYLATION','TREG_IMMUNOREG']:add('treg_secondary',x,y+'_Full',role='secondary')
for ct in ['Treg','Treg_NC','Treg_CY']:
 for y in ['OXIDATIVE_PHOSPHORYLATION_Full','G2M_CHECKPOINT_Full']:add('treg_decomposition','MP10_Q50',y,ct,role='secondary')
add('treg_decomposition','MP10_Q50','G2M_CHECKPOINT_Full','Treg','HC3','NC_matched','control')
for ct in ['Tfh/CXCL13-like','Capillary EC','Arterial EC']:
 for y in ['MTORC1_SIGNALING_Full','MYC_TARGETS_V1_Full','TNFA_SIGNALING_VIA_NFKB_Full','OXIDATIVE_PHOSPHORYLATION_Full','ANGIOGENESIS_Full','G2M_CHECKPOINT_Full']:add('parallel_'+ct.replace('/','_').replace(' ','_'),'MP10_Q50',y,ct,role='secondary')
base=list(models);datasets=sorted(old('specimen_exposure_eligibility.tsv').dataset.unique())
for m in base:
 settings=['min30','replace_b','epi_RNA','receiver_RNA']+['LODO_'+d for d in datasets]
 if m['x']=='MP10_Q50':settings+=['Q95']
 if m['celltype']=='Treg_NC':settings+=['receiver_MT']
 if m['celltype']=='Treg' and m['method']=='HC3': settings+=['cycling_adjust']
 if m['celltype']=='Epi_all':settings+=['epi_major_fraction']
 if m['x'].endswith('_M20') or (m['celltype']=='Epi_all' and m['y'].endswith('_M20')):settings+=['Full']
 if m['x'] in ['SREBF1_M20','SREBF2_M20'] or m['y'] in ['SREBF1_M20','SREBF2_M20']:settings+=['DoRothEA_AB','TF_unique','remove_NECTIN2']
 for setting in settings:
  z=m.copy();z.update(setting=setting,role='sensitivity',FDR_family=m['FDR_family']+'__'+setting)
  if setting=='min30':z['min_cells']=30
  elif setting=='replace_b':z['patient']='replacement'
  elif setting.startswith('LODO_'):z['omit_dataset']=setting[5:]
  elif setting in ['epi_RNA','receiver_RNA','receiver_MT','cycling_adjust','epi_major_fraction']:z['covariate']=setting
  elif setting=='Q95':z['x']='MP10_Q95'
  models.append(z)
reg=pd.DataFrame(models);reg.insert(0,'model_id',['F%05d'%i for i in range(1,len(reg)+1)]);reg['planned_n']=reg.groupby('FDR_family').model_id.transform('size');reg['min_patient_units']=4;reg['min_residual_df']=3;reg['test']='two_sided';put(reg,'model_registry.tsv')
plan=dict(task='mp10_hypothesis_followup_v1',frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),nature='retrospective directed internal follow-up; not prospective or independent validation',seed=123,threads_total=4,threads_per_process=1,scoring=dict(method='GSVA ssGSEA',version='2.4.9',alpha=.25,normalize=False,min_genes=5,min_fraction=.5,universe='edgeR filterByExpr on original lineage profiles >=30, reference fixed',mask='remove Top20 from downstream target AND ranking/scoring universe; no modification of counts or full-RNA normalization',TMM='full RNA library sizes retained for new scores; fixed reference profile per lineage; no re-normalization after masking',centering='gene means from main-a >=50 reference profiles; fixed for new TF ULM and all sensitivity',scale='All-Treg main-a >=50 scores set receiver endpoint scale reused for NC/CY; Epi_all main-a reference SD; original MP10 SD retained'),cycle=dict(method='Seurat CellCycleScoring / cc.genes.updated.2019',pool='all detected genes in all original Treg normalized RNA',seed=123,nbin=24,ctrl='min(number of mapped S, mapped G2M)',validity='>=20 and >=50% mapped genes in each list; finite scores; Phase S/G2M operational cycling; G1 non-cycling; no G0 claim'),models=dict(continuous='OLS HC3 t(df residual) for P and CI',proportion='quasibinomial glm; Pearson dispersion; t(df residual) P and CI',eligibility='>=4 recorded patient units, full rank, residual df>=3; n_cells>=50 main >=30 sensitivity',failed_ratio='nonconvergence, nonfinite covariance, constant outcome, separation or coefficient SE>10 => NE',FDR='BH fixed planned family n including NE as unobserved; NE P/q NA',clinical='only reliable fields with crossed dataset coverage; otherwise NE; no invented staging',interaction='descriptive cross-celltype comparison; no celltype specificity claim without comparable lineage score universe'),gsea=dict(minSize=15,maxSize=500,scoreType='std',eps=1e-10,nPermSimple=[10000,100000],escalation='all pathways of every affected model/library/mask; second budget if NA pval or nonfinite log2err with pval>eps; never based on candidate significance',family='preserve old within model + original joint primary Full family; repaired values separate'),source_result=str(O),source_pb=str(O/'objects/pseudobulk.qs'),source_pb_sha256=sha(O/'objects/pseudobulk.qs'),registry_sha256=sha(B/'tables/model_registry.tsv'),signatures_sha256=sha(B/'tables/signature_registry.tsv'))
# JSON is a strict YAML 1.2 subset, deliberately avoids an extra YAML dependency.
(B/'config/analysis_plan.yaml').write_text(json.dumps(plan,ensure_ascii=False,indent=2))
put([dict(context='existing v1 TF',mask='M20-masked',definition='remove Top20 network targets and centered expression rows after original lineage filter/TMM; original gene centering preserved',normalization='v1 keep.lib.sizes=FALSE; preserved as historical method'),dict(context='new Epi scores/ULM',mask='M20',definition='remove Top20 targets and scoring-universe rows after full RNA TMM/logCPM; never recount libraries',normalization='new frozen Epi_all reference; not claimed identical to v1 activities'),dict(context='existing/repaired GSEA',mask='M20-masked',definition='remove Top20 from DE rank vector and pathway membership; DE unchanged'),dict(context='Treg endpoint',mask='Full',definition='full frozen All-Treg expression universe used for All/NC/CY; no M20 removal')],'masking_definition.tsv')
# Broad old evidence extracts: retain nulls and sensitivity, not a cherry-picked list.
a=old('activity_associations.tsv.gz');a=a[a.source.isin(['SREBF1','SREBF2','MTOR','HIF1A','MYC'])].copy();a['evidence_type']='old_regulon_association';a['variable']=a.source;a['effect']=a.activity_effect;a['P']=a.sample_model_p;a['q']=a.sample_model_q;a['source_path']=str(O/'tables/activity_associations.tsv.gz');a['FDR_family']='old model/resource/mask TF family; preserved'
g=old('gsea_results.tsv.gz');g=g[g.pathway.isin(['HALLMARK_'+n for n in ['MTORC1_SIGNALING','OXIDATIVE_PHOSPHORYLATION','G2M_CHECKPOINT','MYC_TARGETS_V1','TNFA_SIGNALING_VIA_NFKB','ANGIOGENESIS']])].copy();g['evidence_type']='old_GSEA';g['variable']=g.pathway;g['effect']=g.NES;g['P']=g.pval;g['q']=g.padj;g['source_path']=str(O/'tables/gsea_results.tsv.gz');g['FDR_family']='old model/library/mask; joint_main_q separate'
e=pd.concat([a,g],ignore_index=True);meta=old('model_status.tsv');e=e.merge(meta[['model_id','kind','target','n_specimens','n_recorded_patient_units','n_datasets','status','residual_df']],on='model_id',suffixes=('','_model'));e['statistical_unit']='paired specimen or dataset-adjusted recorded patient unit; provisional identity';e['evidence_nature']='same-data association/enrichment';e['shared_data']='all v1 outputs share scRNA profiles; sensitivities are not independent replicates';e['gap']='no kinase/protein/flux/function or directional causality established';e['CI_low']=float('nan');e['CI_high']=float('nan')
extra=[]
for name in ['abundance_associations.tsv','scenic_all_associations.tsv.gz']:
 d=old(name);d['source_path']=str(O/'tables'/name);d['evidence_type']=name;d['gap']='historical full table retained; no new inference';extra.append(d)
e=pd.concat([e,*extra],ignore_index=True);e.insert(0,'evidence_id',['E%06d'%i for i in range(1,len(e)+1)]);put(e,'evidence_matrix.tsv')
notes=[('Epi SREBF2 SCENIC','not_evaluable','no motif-pruned regulon; never substitute SREBF1'),('forward NicheNet','not_evaluable','receiver positive targets only 0-2; no relaxation'),('Treg abundance','limited_support','existing main abundance q not significant'),('Capillary ANGIOGENESIS','direction_not_supportive','main NES negative and not significant; Q95 negative'),('Treg own mTORC1/TF','limited_support','no primary sample-model q<0.05; distinguish pathway GSEA'),('HD local contact/metabolite transfer','pending','not run in this task')];put([dict(relation=a,status=b,reason=c) for a,b,c in notes],'evidence_limitations.tsv')
records=json.loads((B/'logs/literature_source_records.json').read_text());context=[('cell-intrinsic models','mTORC1/S6K-SREBP metabolic transcription','not proof of BLCA Epi-to-Treg signaling'),('human BLCA cells/xenograft','FGFR3 knockdown/SCD perturbation; SREBP cleavage and growth','does not establish MP10-specific FGFR3 or immune mechanism'),('mouse Treg-specific Rptor deletion','mevalonate/lipid programs, CTLA4/ICOS and suppressive function','Treg-intrinsic; not cancer-cell lipid supply'),('gastric cancer human/mouse RHOA','tumor PI3K-AKT-mTOR/free fatty acids and Treg advantage','different cancer/genotype; no BLCA transfer demonstrated'),('tumor mouse Treg SCAP/FASN deletion','tumor growth, IFNg, PD1, prenylation and Treg function','Treg-intrinsic; not Epi SREBF causal bridge'),('tumor Treg CD36 genetic perturbation','mitochondrial adaptation, survival and antitumor responses','does not make OXPHOS RNA a flux/suppression assay'),('murine tumors/human tumor observations','Treg lipid synthesis/glycolysis/oxidation and expansion','mixed contexts; not MP10 exposure association'),('mouse Treg complex III deletion','suppressive function, methylation/metabolites; stable proliferation','respiration requirement; cannot infer low RNA means stronger suppression')]
lit=[]
for i,j in enumerate(records):
 c=context[i];lit.append(dict(title=html.unescape(re.sub('<[^>]+>','',j['title'])),year=j.get('pubYear'),DOI=j.get('doi'),PMID=j.get('id'),PMCID=j.get('pmcid'),URL='https://pubmed.ncbi.nlm.nih.gov/'+j['id']+'/',context=c[0],actual_readouts=c[1],cannot_extrapolate=c[2],verified_level='Europe PMC primary-paper bibliographic record and abstract; selected PMC full text',retrieved_utc=plan['frozen_utc']))
put(lit,'literature_evidence.tsv')
(B/'EVIDENCE_REVIEW.md').write_text('''# 既有证据整理\n\n本任务是看过发现性结果后的定向内部检验。完整正、负和敏感性结果见 tables/evidence_matrix.tsv；旧 P/q 未修改。NECTIN2 的上皮配对上调、SREBF 程序、Treg OXPHOS 负向/G2M 正向是起点，不是已成立的跨细胞机制。\n\n主 Epi 配对覆盖39/33/24，Treg/Tfh/Capillary EC 为32/23/12个记录患者单位；身份仍有 provisional 限制。Arterial EC 残差自由度不足。Treg总丰度、Treg自身调控活动、正向 NicheNet targets 和内皮 ANGIOGENESIS 的限制保留，见 evidence_limitations.tsv。SREBF2 的 Epi SCENIC 为 NE。\n\n文献表区分癌细胞外在环境作用与 Treg 内在依赖。胃癌 RHOA 结果不能直接套用到本项目；复合体 III 的功能必要性尤其不支持把 OXPHOS RNA 下降直接写为免疫抑制增强。空间既有结果仅作背景：mTORC1 去共享基因后减弱，未作为本轮RNA暴露或新验证。\n''')
(B/'HYPOTHESIS_v1.md').write_text('''# 冻结工作假说\n\n项目已有相关：MP10较高 Epi 的脂质/SREBF程序及NECTIN2配对变化；跨标本Treg OXPHOS负向/G2M正向。一般文献机制：mTORC1–SREBP调控脂质程序，Treg具有内在代谢依赖。待验证连接：Epi接触信号或局部代谢改变→Treg状态；平行检查毛细血管内皮。\n\n预测：去Top20后上皮程序仍有支持；若NC共同样本中OXPHOS保持负向，提高状态内变化优先级；若只见cycling比例变化，优先考虑组成；若NECTIN2或特定脂质分支在HC3和预设敏感性中方向稳定，可提出后续验证。反向/弱效应降低相应桥接假说；覆盖不足为NE。不得用显著/不显著差异代替效应差异检验。\n\n替代解释：共享基因、Epi亚群与RNA含量、Treg cycling组成、共同环境、反向作用、数据集依赖。新连续分数与历史NES为不同统计量。\n\n图：Epi MP10 ↔ mTORC1/SREBF/脂质程序（项目相关）；mTORC1 → SREBP（文献一般机制）；Epi NECTIN2 / 代谢环境 ⇢ Treg cycling / OXPHOS（项目待检验）；⇢ Capillary EC（平行支线）。虚线不表示已证实因果。\n''')
put([dict(path=str(p),bytes=p.stat().st_size,mtime=p.stat().st_mtime,sha256=sha(p),reuse='read_only') for p in [O/'objects/pseudobulk.qs',O/'tables/analysis_cell_metadata.tsv.gz',O/'tables/pseudobulk_membership.tsv.gz',O/'tables/specimen_exposure_eligibility.tsv',O/'tables/patient_relationship_audit.tsv',R/'hallmark.tsv.gz',R/'reactome.tsv.gz',R/'collectri.tsv',R/'dorothea_AB.tsv']],'input_reuse_audit.tsv')
(B/'config/M0_FROZEN.sha256').write_text(sha(B/'config/analysis_plan.yaml')+'\n')
print('M0 FROZEN',plan['frozen_utc'],'models',len(reg),'evidence',len(e))
