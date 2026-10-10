import bootstrap
from bootstrap import ROOT, RUN
import pandas as pd
import numpy as np
import json,hashlib,zipfile,html,platform
from pathlib import Path
from PIL import Image
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
 return h.hexdigest()
def dump(x,p):Path(p).write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n')
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA')
man=pd.read_csv(ROOT/'section_input_summary.tsv',sep='\t');st=pd.read_csv(ROOT/'section_status.tsv',sep='\t');perf=pd.read_csv(ROOT/'mp10_baseline_vs_tme.tsv',sep='\t');dic=pd.read_csv(ROOT/'variable_dictionary.tsv',sep='\t');labels=dic.set_index('variable').source_label.to_dict()
checks=[]
def check(name,ok,detail):
 checks.append(dict(check=name,status='PASS' if ok else 'BLOCKED',detail=detail))
check('63 frozen sections and 7 datasets',len(man)==63 and man.dataset.nunique()==7,f'{len(man)}, {man.dataset.nunique()}')
check('both branches have section status',len(st)==2*len(man),'A and B separately retained')
check('no technical block',(st.status!='BLOCKED').all(),str(st.status.value_counts().to_dict()))
check('pilot first and passed',json.loads((ROOT/'pilot_validation.json').read_text())['status']=='PASS','BC15')
check('audited learner equals installed native learner',json.loads((ROOT/'native_learner_parity.json').read_text())['status']=='PASS','exact OOB predictions and impurity importances on real pilot')
ca=pd.read_csv(ROOT/'native_coefficient_serialization_audit.tsv',sep='\t')
check('native coefficient serialization defect audited and corrected',ca.coefficient_p_verified.all() and ca.native_importance_reproduced.all(),f'{len(ca)} target-model records; native outputs retained')
check('only permitted predictors',True,'checked below against exact allowlists, every target and fold')
params=json.loads((ROOT/'parameters.json').read_text())
allowA=set(params['major_variables']);allowB=set(params['technical_variables']+params['tme_variables'])
au=pd.read_csv(ROOT/'predictor_exclusion_audit.tsv',sep='\t');ok=True
for row in au.itertuples():
 allowed=allowA-{row.target} if row.model=='wang_style_cellcell_intra' else set(params['technical_variables']) if row.model=='mp10_technical_baseline' else allowB
 ok=ok and row.predictor in allowed and row.predictor!=row.target
checks[-1]['status']='PASS' if ok else 'BLOCKED'
for row in man.itertuples():
 folder=ROOT/'by_dataset'/row.dataset/row.section_id
 f=folder/'mp10_spatial_predictions.tsv'
 if not f.exists() or f.stat().st_size<10:continue
 d=pd.read_csv(f,sep='\t');keys=['section_id','barcode','spatial_fold'];a=d[d.model=='baseline'];b=d[d.model=='tme']
 same=set(map(tuple,a[keys].to_numpy()))==set(map(tuple,b[keys].to_numpy()))
 check(row.section_id+' same evaluation spots',same and not d.duplicated(keys+['model']).any(),f'{len(a)} baseline; {len(b)} full')
 pp=perf[perf.section_id==row.section_id]
 for tag,g in d.groupby('model'):
  stored=pp[pp.model==tag].iloc[0]
  if not stored.spatial_cv_complete:continue
  r2=1-np.sum((g.MP10-g.prediction)**2)/np.sum((g.MP10-g.MP10.mean())**2)
  check(row.section_id+' '+tag+' recomputed R2',np.isclose(r2,stored.spatial_r2,atol=1e-10,rtol=1e-10),str(r2))
 # Verify no training/test key overlap from the emitted fold audit.
 ff=pd.read_csv(folder/'spatial_fold_audit.tsv',sep='\t')
 check(row.section_id+' spatial fold audit',ff.fold_overlap.eq(0).all() and ff.target_excluded.all(),'training-only OOB calibration; fixed folds')
 # Independently reproduce native importance standardization from raw importances and coefficient P values.
 for tag in ['wang_style_cellcell_intra','mp10_target_extension']:
  fn=folder/('cellcell_importance.tsv' if tag.startswith('wang') else 'mp10_importance.tsv')
  if not fn.exists():continue
  imp=pd.read_csv(fn,sep='\t');coef=pd.read_csv(folder/tag/'coefficients.txt',sep=r'\s+',skiprows=1,names=['target','intercept','intra','extra_NA','p.intercept','p.intra'])
  matches=[]
  for target,g in imp.groupby('Target'):
   raw=pd.read_csv(folder/tag/f'importances_{target}_intra.txt').set_index('target').imp
   z=(raw-raw.mean())/raw.std(ddof=1)*(1-float(coef.loc[coef.target==target,'p.intra'].iloc[0]))
   valid=g[g.Predictor.isin(z.index)]
   matches.append(np.allclose(valid.standardized_importance,z.loc[valid.Predictor],atol=1e-8,rtol=1e-8,equal_nan=True))
  check(row.section_id+' '+tag+' importance formula',all(matches),'sample-SD standardization times (1-p_intra)')
figs=pd.read_csv(ROOT/'figure_index.tsv',sep='\t');fcheck=[]
for row in figs.itertuples():
 with Image.open(row.png) as im:
  dpi=im.info.get('dpi',(0,0));fcheck.append(dict(path=row.png,width=im.width,height=im.height,dpi_x=dpi[0],dpi_y=dpi[1],status='PASS' if abs(dpi[0]-600)<1 else 'BLOCKED'))
 check(row.stem+' PDF',Path(row.pdf).read_bytes()[:4]==b'%PDF',row.pdf)
save(pd.DataFrame(fcheck),ROOT/'figure_validation.tsv')
check('600 dpi PNGs',all(x['status']=='PASS' for x in fcheck),str(len(fcheck)))
save(pd.DataFrame(checks),ROOT/'acceptance_validation.tsv')
ds=pd.read_csv(ROOT/'dataset_performance_summary.tsv',sep='\t');b=perf[perf.model=='tme'];n_complete=int(b.spatial_cv_complete.sum());n_report=int(b.importance_interpretable.sum());n_inc=int((b.delta_r2_tme_minus_baseline>0).sum())
status='PASS_WITH_WARNING' if all(x['status']=='PASS' for x in checks) else 'BLOCKED'
lines=['# Module 06 已完成：intra-spot 预测依赖', '',f'状态：**{status}**。覆盖 {len(man)} 张冻结切片、{man.dataset.nunique()} 个数据集、{man.n_tissue.sum():,} 个 tissue spots。仅运行 Module 06；未启动后续模块。','',f'A：按批准映射汇总 9 个 major posterior mean abundances，逐切片训练；状态分布：{st[st.branch=="A"].status.value_counts().to_dict()}。',f'B：冻结 MP10 raw AUCell 为目标，技术 baseline = 原 Epi fraction + log1p(nCount)；完整模型增加全部 29 个非上皮 TME factors。同一 complete-case 集合与相同空间 folds。{n_complete}/{len(man)} 张切片完成五折空间评估；{n_inc} 张 delta R² > 0；{n_report} 张同时达到 full R² > 0、delta R² > 0 的 importance 解释条件。','','\n'.join(['| '+' | '.join(ds.columns)+' |','| '+' | '.join(['---']*len(ds.columns))+' |']+['| '+' | '.join(str(round(v,5)) if isinstance(v,float) else str(v) for v in row)+' |' for row in ds.itertuples(index=False,name=None)]),'','以上是切片层面的描述性预测结果，不是患者独立重复的推断、显著性检验、通讯或因果证据。负提升保留，不计为技术失败。','', '## 主要交付','', '- `cellcell_importance.tsv`：逐切片 major 网络 raw 与 standardized importance。','- `mp10_baseline_vs_tme.tsv`：相同 spots/folds 的未截断空间预测表现与增量。','- `cellcell_dataset_median_standardized_importance.tsv`、`mp10_dataset_median_standardized_importance.tsv`：按 dataset 的中位数及有效切片数，同时保留全量与性能合格汇总。','- `abundance_mp_technical_associations.tsv`：与 importance 分开的方向关联表。','- `by_dataset/`：模型、完整预测、输入、性能、审计与 PDF + 600 dpi PNG。','- `review.html`、`review_bundle.zip`：紧凑审阅入口。','- `variable_dictionary.tsv`、`EVALUATION_METHODS.md`、`parameters.json`、`input_manifest.tsv`、`acceptance_validation.tsv`、`handoff.json`。','', '## 限制与来源','', 'MISTyR 1.18.0（原文 1.10.0）与 distances 已在模块私有 R 库安装。官方集成输出会截断负 R²，本项目另存未截断 OOB/内部 meta-fold 指标；空间评估完全独立重训训练 fold 的森林与 OOB 校准器。A 仅内部验证；B 为切片内连续区域外推，无 buffer，不代表新患者验证。', '', '标准化 importance 的负值表示低于相对平均贡献，不表示抑制。丰度与 counts 同变、变量相关性和混合 spots 均限制生物解释。patient_uid 和未核实 dataset 标签沿用冻结输入，不以 section 代替患者。','', '方法与解释详见 [EVALUATION_METHODS.md](EVALUATION_METHODS.md)。已在预测依赖/增量预测信息层面停止。']
lines += ['', '原生输出缺陷修正：MISTyR 1.18.0 的系数表为 5 列表头/6 列数据，原生 importance 误用了截距 P 值。已逐 target 用保存的 OOB 数据重建校准回归，验证真实 slope P 值并修正汇总权重；保留 native_standardized_importance、原生文件与完整缺陷审计。此修正不改变预测性能。']
(ROOT/'RUN_SUMMARY.md').write_text('\n'.join(lines)+'\n')
body=['<!doctype html><meta charset="utf-8"><title>Module 06 review</title><style>body{font:16px system-ui;max-width:1200px;margin:30px auto;line-height:1.5}img{max-width:100%}table{border-collapse:collapse;font-size:13px}td,th{padding:7px;border:1px solid #ddd}</style><h1>Module 06 · intra-spot prediction</h1>',f'<p>{len(man)} sections / {man.dataset.nunique()} datasets. MP10: {n_complete} completed spatial CV; {n_inc} positive increments; {n_report} reportable importance models.</p><p>Predictive dependence only. Negative standardized importance does not mean inhibition. Grey = unavailable/not reportable.</p>',ds.to_html(index=False),'<p><a href="RUN_SUMMARY.md">运行总结</a> · <a href="EVALUATION_METHODS.md">方法与限制</a> · <a href="mp10_baseline_vs_tme.tsv">完整性能表</a></p>']
for row in figs.itertuples():
 rel=str(Path(row.thumbnail).relative_to(ROOT));pdf=str(Path(row.pdf).relative_to(ROOT));body.append(f'<h2>{html.escape(str(Path(row.pdf).parent.parent.relative_to(ROOT)))} / {html.escape(row.stem)}</h2><a href="{html.escape(pdf)}">PDF</a><img src="{html.escape(rel)}">')
(ROOT/'review.html').write_text('\n'.join(body))
source=dict(mistyR_release_url='https://bioconductor.org/packages/3.22/bioc/src/contrib/mistyR_1.18.0.tar.gz',archive=str(ROOT/'sources/mistyR_complete.tar.gz'),archive_sha256=sha(ROOT/'sources/mistyR_complete.tar.gz'),git_last_commit='1401996',installed_version='1.18.0',paper_version='1.10.0',paper_doi='10.1186/s12943-025-02377-9',retrieval_date='2026-10-05',installation_library=str(ROOT/'runtime/Rlib'),base_environment_modified=False)
dump(source,ROOT/'source_registry.json')
code={str(f.relative_to(ROOT)):sha(f) for f in (ROOT/'code').glob('*') if f.is_file()};dump(code,ROOT/'code_checksums.json')
inputsig=sha(ROOT/'input_manifest.tsv');codesig=sha(ROOT/'code_checksums.json')
warnings=['MISTyR 1.18.0 differs from paper 1.10.0.','Native coefficient serialization has an extra NA field: importance weighting corrected with OOB-verified slope P values; unmodified native values preserved.','A internal validation only; B held-out regions within section, no buffer or independent-patient validation.','Frozen patient_uid absent; section is not a patient. Frozen __UNVERIFIED dataset labels retained.','Shared counts and correlated abundances can drive predictive dependencies; no causal or communication claims.','Native negative R2 clipping is bypassed for reporting via retained unclipped metrics.','Negative standardized importance means below-average contribution, not inhibition.','Performance-qualified median is conditional; all-estimable medians also supplied.']
handoff=dict(module='06',status=status,completed_scope=dict(sections=len(man),datasets=man.dataset.nunique(),branches=['wang_style_cellcell_intra','mp10_target_extension'],views=['intra'],spatial_cv_complete_sections=n_complete,positive_increment_sections=n_inc,mp10_reportable_sections=n_report),upstream_handoffs=[str(RUN/'00_inputs/handoff.json'),str(RUN/'01_mp_display/handoff.json')],input_signature=inputsig,code_signature=codesig,parameter_signature=sha(ROOT/'parameters.json'),source_registry=str(ROOT/'source_registry.json'),field_units=dict(raw_importance='ranger impurity importance; target scale dependent',standardized_importance='per-target z of impurity importance times (1-p_intra); no effect sign',spatial_r2='dimensionless fraction, not percent; may be negative',spatial_rmse='raw AUCell',mean_abundance='expected cells per spot',mp_score='raw frozen AUCell',nCount='raw UMI counts across all input features'),spot_key_columns=['section_id','barcode'],warnings=warnings,stop_after_module='06',next_module_started=False,outputs={f:str(ROOT/f) for f in ['RUN_SUMMARY.md','review.html','review_bundle.zip','output_manifest.tsv','cellcell_importance.tsv','mp10_baseline_vs_tme.tsv','cellcell_dataset_median_standardized_importance.tsv','mp10_dataset_median_standardized_importance.tsv','variable_dictionary.tsv','EVALUATION_METHODS.md','acceptance_validation.tsv']})
dump(handoff,ROOT/'handoff.json')
with zipfile.ZipFile(ROOT/'review_bundle.zip','w',zipfile.ZIP_DEFLATED) as z:
 for f in sorted(ROOT.rglob('*')):
  rel=f.relative_to(ROOT)
  if not f.is_file():continue
  include=(len(rel.parts)==1 and f.suffix in ['.tsv','.md','.json','.html']) or 'thumbnails' in rel.parts or ('figures' in rel.parts and f.suffix=='.pdf') or rel.parts[0]=='code'
  if include and '__pycache__' not in rel.parts and f.name!='output_manifest.tsv':z.write(f,str(rel))
outs=[]
for f in sorted(ROOT.rglob('*')):
 rel=f.relative_to(ROOT)
 if not f.is_file() or rel.parts[0] in ['runtime','sources'] or '__pycache__' in rel.parts or f.name=='output_manifest.tsv':continue
 role='model' if f.suffix=='.rds' else 'figure' if f.suffix in ['.pdf','.png'] else 'code' if rel.parts[0]=='code' else 'table' if f.suffix=='.tsv' else 'documentation_or_log'
 outs.append(dict(path=str(f),bytes=f.stat().st_size,sha256=sha(f),role=role))
save(pd.DataFrame(outs),ROOT/'output_manifest.tsv')
print(status,'checks',len(checks),'failed',sum(x['status']!='PASS' for x in checks),'outputs',len(outs),'CV',n_complete,'positive increment',n_inc,'reportable',n_report)
assert status!='BLOCKED','Final validation failed; see acceptance_validation.tsv'
