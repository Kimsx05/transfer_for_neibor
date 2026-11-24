from public_common import *
import zipfile
z=zipfile.ZipFile(next((C/'CTRPv2.1').glob('*.zip')))
def tab(n,**kw):return pd.read_csv(z.open('v21.'+n+'.txt'),sep='\t',**kw)
feat=tab('meta.gex_features');meta=tab('meta.per_cell_line'); drug=tab('meta.per_compound')
ids=[63578,154846,609060,609110];drug=drug[drug.master_cpd_id.isin(ids)].copy()
assert list(drug.cpd_name)==['erastin','ML162','1S,3R-RSL-3','ML210']
drug['bottom_umol_nominal']=drug.top_test_conc_umol/(2**15)
output(drug,'tables/CTRP_compounds.tsv')
# Only unambiguously annotated features enter the platform background.
f=feat[feat.gene_symbol.eq(1)&feat.gene_primary_name.notna()].set_index('idx_gene_feature')
output(f[f.gene_primary_name.isin(top)].reset_index(),'tables/CTRP_MP10_features.tsv')
g=tab('data.gex_avg_log2',dtype={'master_ccl_id':'int32','idx_gene_feature':'int32','mrna_expression_avg_log2':'float32'})
mat=g.pivot(index='idx_gene_feature',columns='master_ccl_id',values='mrna_expression_avg_log2');del g
mat=mat.loc[mat.index.intersection(f.index)]; mat.index=f.loc[mat.index,'gene_primary_name'];mat=mat.groupby(level=0).mean()
r=mat.rank(pct=True);s=r.reindex(top).mean()
cov=r.reindex(top).notna().sum()
output(pd.DataFrame(dict(master_ccl_id=s.index,score=s.values,n_top=cov.values,background=mat.notna().sum().values)),'tables/CTRP_scores.tsv')
m=meta.rename(columns={'ccl_name':'name','ccle_primary_site':'lineage'}); m['unit']=m.name.replace({'253JBV':'253J'})
m['eligible']=(m.lineage=='urinary_tract')&(m.ccle_hist_subtype_1=='transitional_cell_carcinoma')&~m.name.isin(['SCABER'])
from identity import annotate
m=annotate(m)
m=m.merge(s.rename('score'),how='left',left_on='master_ccl_id',right_index=True)
a=tab('data.auc_sensitivities');a=a[a.master_cpd_id.isin(ids)]
a=a.groupby(['master_ccl_id','master_cpd_id']).area_under_curve.agg(['mean','size']).rename(columns={'mean':'response','size':'experiments'}).reset_index()
joined=m.merge(drug[['master_cpd_id','cpd_name']].rename(columns={'cpd_name':'endpoint'}),how='cross').merge(a,on=['master_ccl_id','master_cpd_id'],how='left')
u=joined[joined.lineage=='urinary_tract'];output(u,'tables/CTRP_urinary_crosswalk.tsv')
audit=[]
for drug,x in u.groupby('endpoint'):
 y=x[x.eligible];t=y.dropna(subset=['response','score']);audit.append(dict(endpoint=drug,urinary=len(x),eligible=len(y),no_expression=y.score.isna().sum(),no_auc=y.response.isna().sum(),both_missing=(y.response.isna()&y.score.isna()).sum(),models=len(t),units=t.unit.nunique(),score_min=t.score.min(),score_max=t.score.max(),auc_min=t.response.min(),auc_max=t.response.max()))
output(pd.DataFrame(audit),'tables/CTRP_intersections.tsv')
associations(joined,'CTRP')

