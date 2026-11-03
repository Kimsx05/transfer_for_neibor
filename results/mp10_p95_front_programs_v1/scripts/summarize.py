from analyze import *
from scipy.stats import ttest_1samp

def bh(p):
 p=np.asarray(p,float);q=np.full(len(p),np.nan);ids=np.flatnonzero(np.isfinite(p));ix=ids[np.argsort(p[ids])]
 q[ix]=np.minimum.accumulate((p[ix]*len(ix)/np.arange(1,len(ix)+1))[::-1])[::-1].clip(0,1);return q

def main():
 for name in ['coverage','effects','bands','niche','matching']:
  frames=[]
  for s in M.section_id:
   p=O/'tables'/(f'{s}_{name}.tsv'+('' if name=='coverage' else '.gz'))
   p=Path(str(p))
   if not p.exists():continue
   try:frames.append(pd.read_csv(p,sep='\t'))
   except pd.errors.EmptyDataError:pass
  a=pd.concat(frames,ignore_index=True)
  if name=='effects':
   a['panel']=np.where(a.program.str.contains('__NO_MP10|TOP100'),'sensitivity','main')
   a['test_p']=np.where(a.quadratic_beta.notna(),a.overall_p,a.p)
   a['q']=a.groupby(['section_id','family','panel','adjustment']).test_p.transform(bh)
  save(a,O/'tables'/f'{name}.tsv.gz')
 e=pd.read_csv(O/'tables/effects.tsv.gz',sep='\t').merge(M[['section_id','dataset_id','source_context','analysis_unit_uid','patient_uid','eligible_for_patient_balanced_analysis']],on='section_id',validate='m:1')
 # Preserve quadratic effects and dataset coverage, rather than selecting slope signs alone.
 e2=e[e.quadratic_beta.notna()].copy();e2['effect']=e2.quadratic_beta;e2['analysis']=e2.analysis+'_quadratic';e=pd.concat([e,e2],ignore_index=True)
 keys=['family','panel','analysis','adjustment','program'];out=[];units=[];datasets=[]
 for source,g in e.groupby('source_context'):
  for key,v in g.groupby(keys):
   for pop,b,idcol in [('sections',v,'section_id'),('units',v[v.analysis_unit_uid.notna()],'analysis_unit_uid'),('verified',v[v.eligible_for_patient_balanced_analysis.eq(True)&v.patient_uid.notna()],'patient_uid')]:
    b=b.dropna(subset=['effect']);b=b[np.isfinite(b.effect)]
    if b.empty:continue
    u=b.groupby(['dataset_id',idcol]).effect.mean().reset_index();ds=u.groupby('dataset_id').effect.mean();common=dict(source=source,population=pop,**dict(zip(keys,key)))
    for rr in u.itertuples(index=False):units.append(dict(common,dataset=rr[0],unit=rr[1],effect=rr[2]))
    for dataset,w in u.groupby('dataset_id'):
     vals=w.effect.to_numpy();datasets.append(dict(common,dataset=dataset,effect=vals.mean(),n_units=len(vals),n_sections=b[b.dataset_id==dataset].section_id.nunique(),positive_units=(vals>0).sum(),negative_units=(vals<0).sum()))
    rng=np.random.default_rng(20261103);boots=[]
    for _,w in u.groupby('dataset_id'):
     vals=w.effect.to_numpy();boots.append(rng.choice(vals,(1999,len(vals)),replace=True).mean(1))
    boots=np.array(boots);bt=boots[rng.integers(0,len(ds),(len(ds),1999)),np.arange(1999)].mean(0)
    loo=[ds.drop(i).mean() for i in ds.index] if len(ds)>1 else [np.nan]
    out.append(dict(common,effect=ds.mean(),low=np.quantile(bt,.025),high=np.quantile(bt,.975),p=ttest_1samp(ds,0).pvalue if len(ds)>=3 and ds.std()>0 else np.nan,n_sections=b.section_id.nunique(),n_units=len(u),n_verified_patients=b.loc[b.eligible_for_patient_balanced_analysis.eq(True),'patient_uid'].nunique(),n_datasets=len(ds),positive_units=(u.effect>0).sum(),negative_units=(u.effect<0).sum(),positive_datasets=(ds>0).sum(),negative_datasets=(ds<0).sum(),loo_min=np.nanmin(loo),loo_max=np.nanmax(loo)))
 a=pd.DataFrame(out);a['q']=a.groupby(['source','population','family','panel','adjustment']).p.transform(bh)
 save(a,O/'tables/summary_effects.tsv');save(units,O/'tables/unit_effects.tsv.gz');save(datasets,O/'tables/dataset_effects.tsv.gz')
 coverage=pd.read_csv(O/'tables/coverage.tsv.gz',sep='\t').merge(M[['section_id','dataset_id','source_context','analysis_unit_uid','patient_uid','eligible_for_patient_balanced_analysis']],on='section_id');save(coverage,O/'full_sections.tsv')
 # Gene mapping/coverage and matching diagnostics.
 a=pd.concat([pd.read_csv(p,sep='\t') for p in (O/'scores').glob('*_match.tsv')]);save(a,O/'tables/new_gene_set_coverage.tsv')
 a=pd.concat([pd.read_csv(p,sep='\t') for p in (O/'scores').glob('*_gene_detection.tsv')]);save(a,O/'tables/gene_detection.tsv.gz')
 save([dict(program=n,n_original=g.original_symbol.nunique(),n_mapped=g.gene_symbol.nunique(),n_overlap=g.loc[g.MP10_overlap,'gene_symbol'].nunique(),overlap='|'.join(g.loc[g.MP10_overlap,'gene_symbol']),source=g.source.iloc[0],panel=g.panel.iloc[0]) for n,g in G.groupby('gs_name')],O/'tables/signature_index.tsv')
 print('summary complete',len(e),len(a),flush=True)
if __name__=='__main__':main()
