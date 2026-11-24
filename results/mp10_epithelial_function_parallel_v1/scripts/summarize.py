from public_common import *
import platform,sys,openpyxl,scipy
from importlib.metadata import version
def fmt(v):return 'NE' if pd.isna(v) else f'{v:.3g}'
m=read('tables/project_statistics.tsv');profiles=read('plot_data/project_profiles.tsv.gz')
gen=read('tables/project_gene_drivers.tsv.gz').query('contrast=="H05"').set_index('gene')
summary=[]
for row in profiles.query('contrast=="H05"').drop_duplicates(['module','mask']).itertuples():
 key=row.module;mask=row.mask
 r=m[(m.contrast=='H05')&(m.module==key)&(m['mask']==mask)]
 main=r[r.model=='paired'];exclude=r[r.model.str.startswith('omit_')]
 g=[a for a in sets.get(key,[]) if a in gen.index and (mask=='Full' or a not in top)]
 drivers=gen.reindex(g).delta.dropna().sort_values(key=abs,ascending=False)
 summary.append(dict(module=key,mask=mask,covered=row.covered,expected=row.expected,covered_exact=' '.join(g) if key in sets else 'PROGENy_coverage.tsv',drivers='; '.join(f'{k} {v:+.2f}' for k,v in drivers.head(3).items()),effect=main.effect.iloc[0] if len(main) else np.nan,lo=main.lo.iloc[0] if len(main) else np.nan,hi=main.hi.iloc[0] if len(main) else np.nan,p=main.p.iloc[0] if len(main) else np.nan,q=main.q.iloc[0] if len(main) else np.nan,omit_min=exclude.effect.min(),omit_max=exclude.effect.max(),status='computed' if len(main) else 'NE: no independent members'))
summary=pd.DataFrame(summary)
con=read('tables/PROGENy_H05_exact_contributions.tsv')
for (module,msk),g in con.groupby(['module','mask']):
 ix=(summary.module==module)&(summary['mask']==msk)
 summary.loc[ix,'drivers']='; '.join(f'{r.gene} {r.mean_contribution:+.3f}' for r in g.head(3).itertuples())
 summary.loc[ix,'covered_exact']=' '.join(g.gene)
output(summary,'tables/project_modules_H05.tsv')
d=read('tables/candidate_profile_detection.tsv.gz');d=d[(d.profile_type=='H05')&d.main_patient&d.eligible_50&d.gene.isin(['HMGCS1','SREBF2','CHMP4B','GPX4','SCAP'])]
gate=d.groupby('gene').det_fraction.agg(['count','min','median','max']).reset_index()
gate['matched_independent_Epi_RNA_target_perturbation']='未取得靶点匹配的独立 BLCA Epi 扰动转录组'
gate['network_patient_holdout']='现有网络未证明按患者留出后的重构稳定性';gate['decision']='NE';output(gate,'tables/virtual_KO_gate.tsv')
effects=[];influence=[]
for src in ['CTRP','DepMap']:
 a=read('tables/'+src+'_associations.tsv');a['platform']=src;effects.append(a)
 units=read('plot_data/'+src+'_units.tsv');units=units[units.eligible]
 inf=read('tables/'+src+'_influence.tsv')
 for k,g in units.groupby('endpoint'):
  # Descriptive influential points: highest absolute residual from the linear trend in ranks.
  xx=stats.rankdata(g.score); yy=stats.rankdata(g.response);res=yy-np.polyval(np.polyfit(xx,yy,1),xx)
  z=g[['endpoint','models','score','response']].copy();z['rank_residual']=res;z['abs_residual_rank']=stats.rankdata(-abs(res))
  influence.append(z.merge(inf[inf.endpoint.eq(k)],left_on=['endpoint','models'],right_on=['endpoint','omit']))
output(pd.concat(influence),'tables/model_influence_and_residuals.tsv')
output(pd.concat(effects),'tables/pharmacology_dependency.tsv')
phy=read('tables/pharmacology_dependency.tsv');mask_summary=summary[(summary['mask']=='M20-masked')&~summary.module.isin(['MP10','cycling','stress'])]
def tab(d,cols):
 lines=['| '+' | '.join(cols.values())+' |','|'+'|'.join(['---']*len(cols))+'|']
 for _,r in d.iterrows():lines.append('| '+' | '.join(fmt(r[k]) if isinstance(r[k],(float,np.floating)) else str(r[k]) for k in cols)+' |')
 return '\n'.join(lines)
h=tab(mask_summary,{'module':'模块','covered':'去重叠后覆盖','effect':'High−Low','lo':'95%下限','hi':'上限','q':'BH q'})
e=tab(phy[phy.cohort.eq('BLCA')],{'endpoint':'实测终点','n':'独立记录单位','effect':'ρ','lo':'95%下限','hi':'上限','p':'p','q':'q'})
a=read('tables/CTRP_intersections.tsv');i=tab(a,{'endpoint':'药物','eligible':'合格模型','no_expression':'缺表达','no_auc':'缺AUC','models':'交集模型','units':'合并后单位','auc_min':'最低AUC','auc_max':'最高AUC'})
tmpl=(O/'scripts/results_template.md').read_text()
(O/'RESULTS.md').write_text(tmpl.replace('{EXTERNAL_TABLE}',e).replace('{COVERAGE_TABLE}',i).replace('{MODULE_TABLE}',h))
env=[('python',sys.version.replace('\n',' ')),('platform',platform.platform())]+[(k,version(k)) for k in ['numpy','pandas','scipy','openpyxl']]
output(pd.DataFrame(env,columns=['component','version']),'logs/python_versions.tsv')

