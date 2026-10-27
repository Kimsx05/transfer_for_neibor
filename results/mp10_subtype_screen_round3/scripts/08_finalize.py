from common import *
import shutil,zipfile,subprocess
subprocess.run(['/usr/bin/python3',str(O/'scripts/07_report.py')],check=True)
obj=read('object_index.tsv');cov=read('node_coverage.tsv');edge=read('total_node_edges.tsv.gz');g=read('global_incoming_outgoing.tsv');c=read('epi_allTME_total_contrasts.tsv');order=read('node_order.tsv');lr=read('treg_LR_contribution_decomposition.tsv.gz');par=read('aggregation_parity.tsv')
# Enrich edge coverage using the same object/node membership without changing numerical W/C/R.
lookup=cov.set_index(['object_id','node'])
for side in ['source','target']:
 for field in ['n_sampled_cells','n_contributing_recorded_IDs','n_matched_specimens','n_matched_recorded_IDs']:
  idx=pd.MultiIndex.from_arrays([edge.object_id,edge[side]]);edge[side+'_'+field]=lookup[field].reindex(idx).to_numpy()
save(edge,'total_node_edges.tsv.gz')
g['n_computable_other_partners']=np.where(g.computed,g.n_computable_partners-1,0);save(g,'global_incoming_outgoing.tsv')
checks=[]
def check(name,ok,details=''):
 checks.append(dict(check=name,passed=bool(ok),details=str(details)));assert ok,(name,details)
check('all primary objects present',len(obj)==12 and obj.groupby('network_level').scope.nunique().eq(6).all())
check('only shared within-specimen objects',obj.comparison_type.eq('epi_state_within_specimen').all()&obj.condition.eq('shared').all())
check('native aggregation parity',par.filter(regex='^max_').to_numpy().max()<1e-10,par.filter(regex='^max_').max().to_dict())
check('NE edges remain missing',edge.loc[edge.status.ne('computed'),['W','C','R']].isna().all().all())
check('NE node global sums missing',g.loc[~g.computed,['incoming_all_W','outgoing_all_W','self_W']].isna().all().all())
# Compare standard incl-self IO with v2. There is exactly one Treg row per fine object.
v2io=pd.read_csv(V/'tables/network_incoming_outgoing.tsv',sep='\t');v2io=v2io[v2io.comparison_type.eq('epi_state_within_specimen')];io=g.merge(v2io,left_on=['object_id','node'],right_on=['network_id','node'],how='inner')
check('incoming agrees with v2',np.allclose(io.incoming_all_W,io.incoming_supported,equal_nan=True));check('outgoing agrees with v2',np.allclose(io.outgoing_all_W,io.outgoing_supported,equal_nan=True))
check('single shared Treg per fine object',g[g.node.eq('Treg')&g.network_level.eq('fine')].groupby('object_id').size().eq(1).all())
for o in obj.itertuples():
 a=edge[edge.object_id.eq(o.object_id)];gg=g[g.object_id.eq(o.object_id)].set_index('node');nn=order[order.network_level.eq(o.network_level)].node.tolist();check(o.object_id+' full node pairs',len(a)==len(nn)**2)
 for met in ['W','C','R']:
  inc=a.groupby('target')[met].sum(min_count=1).reindex(gg.index);out=a.groupby('source')[met].sum(min_count=1).reindex(gg.index)
  check(o.object_id+' '+met+' row/column directions',np.allclose(inc,gg['incoming_all_'+met],equal_nan=True) and np.allclose(out,gg['outgoing_all_'+met],equal_nan=True))
  check(o.object_id+' '+met+' self decomposition',np.allclose(gg['incoming_all_'+met]-gg['self_'+met],gg['incoming_other_'+met],equal_nan=True) and np.allclose(gg['outgoing_all_'+met]-gg['self_'+met],gg['outgoing_other_'+met],equal_nan=True))
 for dr in ['epi_to_TME','TME_to_epi']:
  z=c[c.object_id.eq(o.object_id)&c.direction.eq(dr)]
  for state in ['low','high']:
   shares=z[state+'_share'];total=z[state+'_W'].sum(min_count=1);check(o.object_id+' '+dr+' '+state+' share normalization',np.isclose(shares.sum(min_count=1),1) if total>0 else shares.isna().all())
  for met in ['W','C','R']:check(o.object_id+' '+dr+' '+met+' difference',np.allclose(z['delta_'+met],z['high_'+met]-z['low_'+met],equal_nan=True))
check('Treg full LR decomposition',read('treg_contribution_parity.tsv')[['W_absolute_error','R_absolute_error']].to_numpy().max()<1e-10)
# Related reverse DB record remains transparent, never counted as independent validation.
sup=read('supplementary_NECTIN2_TIGIT.tsv');fwd=sup[sup.LR_id.eq('NECTIN2_TIGIT')];rev=sup[sup.LR_id.eq('TIGIT_NECTIN2')];pair=fwd.merge(rev,on=['scope','epi_state'],suffixes=('_forward','_reverse'));pair['raw_difference']=pair.raw_probability_forward-pair.raw_probability_reverse;save(pair,'supplementary_reverse_record_comparison.tsv');check('reverse DB probabilities compared',len(pair)==12)
# All positive matrix edges must be represented in the A circles; B/C use directed masks only.
cs=read('circle_scales.tsv')
for o in obj.itertuples():
 ee=edge[edge.object_id.eq(o.object_id)]
 for measure in ['W','C']:
  sc=cs[cs.object_id.eq(o.object_id)&cs.measure.eq(measure)];check(o.object_id+' '+measure+' shared circle scale',sc.edge_weight_max.nunique()==1)
  for view in ['A_full','B_epi_to_TME','C_TME_to_epi']:
   mask=np.ones(len(ee),bool) if view=='A_full' else ee.source.isin(EPI)&~ee.target.isin(EPI) if view=='B_epi_to_TME' else ~ee.source.isin(EPI)&ee.target.isin(EPI);n=int((ee[measure].gt(0)&mask).sum());check(o.object_id+' '+measure+' '+view+' no edge filtering',sc.loc[sc.view.eq(view),'n_edges_drawn'].iloc[0]==n)
save(checks,'numerical_checks.tsv')
figs=pd.concat([read('figure_index_'+x+'.tsv') for x in ['circles','global','treg','supplement']],ignore_index=True);check('figure paths exist',all((O/p).is_file() for p in figs.path));save(figs,'figure_index.tsv');save(checks,'numerical_checks.tsv')
lines=['# Round3 逐队列图索引','', '主入口是全LR聚合circle与完整热图。W为支持总强度，C为支持LR记录数量；NE不是0。major/fine来自独立对象，不比较其绝对W。','']
for ds in ['ALL']+sorted(obj[obj.scope.ne('ALL')].scope.unique()):
 lines+=['## '+ds,'']
 for lev in ['major','fine']:
  lines+=['### '+lev,'']
  for row in figs[figs.scope.eq(ds)&figs.network_level.eq(lev)&~figs.path.str.contains('supplementary')].itertuples():lines.append(f'- [{Path(row.path).stem}]({row.path}) — {row.question}')
  lines.append('')
lines+=['## NECTIN2–TIGIT 补充','']
for row in figs[figs.path.str.contains('supplementary')].itertuples():lines.append(f'- [{Path(row.path).stem}]({row.path})')
(O/'FIGURE_INDEX.md').write_text('\n'.join(lines)+'\n')
review=O/'review';review.mkdir(exist_ok=True)
for name in ['README.md','methods.md','round3_report.md']:shutil.copy2(O/name,review/name)
# Core PDFs: all six circle views per resolution, full heatmaps, global IO/ranks/roles,
# epithelial all-group comparisons and Treg panels, plus the small supplementary axis.
keep=figs[~figs.path.str.contains('self_fraction_|global_IO_.*_other_|roles_.*_other')].copy()
for row in keep.itertuples():
 dest=review/row.path;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(O/row.path,dest)
# Include every small/medium result table, but no old qs, expression object or PNG.
(review/'tables').mkdir(exist_ok=True)
for fp in sorted(T.glob('*')):
 if fp.is_file():shutil.copy2(fp,review/'tables'/fp.name)
ri=['# 轻量包逐队列图索引','', '固定选图规则覆盖全部五队列和ALL；不只选阳性图。含自环与去自环完整数值均在表中，完整输出目录另有去自环及自环比例图。','']
for ds in ['ALL']+sorted(obj[obj.scope.ne('ALL')].scope.unique()):
 ri+=['## '+ds,'']
 for row in keep[keep.scope.eq(ds)].itertuples():ri.append(f'- [{Path(row.path).stem}]({row.path}) — {row.network_level}; {row.question}')
 ri.append('')
(review/'FIGURE_INDEX.md').write_text('\n'.join(ri)+'\n')
(review/'README.md').write_text('# Round3 阅读包\n\n先读[中文报告](round3_report.md)，再看[逐队列图索引](FIGURE_INDEX.md)和[方法](methods.md)。本包只包含核心PDF、主要表及文字文档；全部TME与NE/反向结果均保留。完整透明PNG、小型聚合矩阵与重现脚本在round3目录。\n\n主图是全LR聚合网络；NECTIN2–TIGIT只在补充部分。Treg全局总量只有一套，不是high/low两组Treg。\n')
with zipfile.ZipFile(O/'review.zip','w',zipfile.ZIP_DEFLATED,compresslevel=6) as zz:
 for fp in sorted(review.rglob('*')):
  if fp.is_file():zz.write(fp,fp.relative_to(O))
with zipfile.ZipFile(O/'review.zip') as zz:
 assert zz.testzip() is None;assert not any(n.endswith(('.qs','.png','.parquet','.mtx','.bin')) for n in zz.namelist())
print('PASS',len(checks),'numerical/plot checks; figures',len(figs),'review PDFs',len(keep),'zip MB',round((O/'review.zip').stat().st_size/1024**2,2))
(O/'logs/COMPLETED.txt').write_text('Round3 modules 0–4 completed using 12 read-only v2 shared-state networks; no inference rerun or external upload.\n')
