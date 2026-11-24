from public_common import *
from PIL import Image
import subprocess
checks=[]
def check(label,ok):
 checks.append((label,bool(ok)))
 if not ok:raise ValueError(label)
def equal(a,b):return np.allclose(a,b,atol=1e-7,rtol=0)
rows=read('tables/project_statistics.tsv')
for state,n in [('H05',24),('H10',33),('H50',39)]:
 check(state+' 配对数',set(rows.query('contrast==@state and model=="paired"').n)=={n})
profiles=read('plot_data/project_profiles.tsv.gz')
for r in rows.query('contrast=="H05" and model=="paired"').itertuples():
 sub=profiles[(profiles.contrast=='H05') & (profiles.module==r.module) & (profiles['mask']==r.mask)].pivot(index='sample_uid',columns='profile_type',values='score')
 check('配对差 '+r.module+' '+r.mask,equal(r.effect,(sub.H05-sub.L50).mean()))
 check('区间包含点估计 '+r.module+' '+r.mask,r.lo<=r.effect<=r.hi)
con=read('tables/PROGENy_H05_exact_contributions.tsv')
for (k,msk),g in con.groupby(['module','mask']):
 est=rows[(rows.contrast=='H05')&(rows.module==k)&(rows['mask']==msk)&(rows.model=='paired')].effect.item()
 check('原足迹和贡献 '+k+' '+msk,equal(g.mean_contribution.sum(),est))
for platform in ('CTRP','DepMap'):
 u=read('plot_data/'+platform+'_units.tsv');a=read('tables/'+platform+'_associations.tsv')
 for r in a[a.cohort=='BLCA'].itertuples():
  g=u[(u.endpoint==r.endpoint)&u.eligible]
  check(platform+' 有效单位和相关 '+r.endpoint,len(g)==g.unit.nunique()==r.n and equal(stats.spearmanr(g.score,g.response).statistic,r.effect))
check('28 个独立定义的条件对比',len(read('tables/GSE149933_contrasts.tsv'))==28)
check('18 个急性 RNA GSM 一一对应',read('tables/GSE182638_samples.tsv').GSM.nunique()==18)
check('12 个 CHMP2A RNA GSM 一一对应',read('tables/GSE175863_samples.tsv').geo_accession.nunique()==12)
check('Brunello 映射',read('tables/GSE126696_mapping.tsv').forward_unique_match.item()==77441)
check('官方缓存',read('tables/resources_final.tsv').final_status.str.contains('由完整 family SOFT 替代').sum()==1)
for pdf in sorted((O/'figures').glob('*.pdf')):
 png=pdf.with_suffix('.png')
 with Image.open(png) as im:
  check(png.name+' 可读',im.width>500 and im.height>500);im.verify()
 check(pdf.name+' 单页可读','Pages: 1' in subprocess.check_output(['mutool','info',str(pdf)],text=True).splitlines())
check('11 组 PDF 和 PNG',len(list((O/'figures').glob('*.pdf')))==len(list((O/'figures').glob('*.png')))==11)
for n in ['RESULTS.md','METHODS.md','FIGURE_INDEX.md']:
 check(n+' 内容完整',(O/n).stat().st_size>1000 and '{MODULE_TABLE}' not in (O/n).read_text())
output(pd.DataFrame(checks,columns=['check','passed']),'logs/verification.tsv')
print(f'{len(checks)} checks passed')

