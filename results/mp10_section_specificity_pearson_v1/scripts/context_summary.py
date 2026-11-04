from core import *
from aggregate import adjust
# Same dataset/unit hierarchical implementation as the accepted previous analysis.
def main():
 e=pd.read_parquet(O/'tables/context85_effects.parquet');extra=e[e.quadratic_beta.notna()].copy();extra['effect']=extra.quadratic_beta;extra['analysis']=extra.analysis+'_quadratic';e=pd.concat([e,extra],ignore_index=True);keys=['family','program_family','analysis','adjustment','program'];out=[];units=[];datasets=[]
 for source,g in e.groupby('source_context'):
  for key,v in g.groupby(keys):
   for pop,b,idcol in [('sections',v,'section_id'),('units',v[v.analysis_unit_uid.notna()],'analysis_unit_uid'),('verified',v[v.eligible_for_patient_balanced_analysis.eq(True)&v.patient_uid.notna()],'patient_uid')]:
    b=b[np.isfinite(b.effect)]
    if b.empty:continue
    u=b.groupby(['dataset_id',idcol]).effect.mean().reset_index();ds=u.groupby('dataset_id').effect.mean();common=dict(source=source,population=pop,**dict(zip(keys,key)))
    for rr in u.itertuples(index=False):units.append(dict(common,dataset=rr[0],unit=rr[1],effect=rr[2]))
    for dataset,w in u.groupby('dataset_id'):
     vals=w.effect.to_numpy();datasets.append(dict(common,dataset=dataset,effect=vals.mean(),n_units=len(vals),n_sections=b[b.dataset_id==dataset].section_id.nunique(),positive_units=(vals>0).sum(),negative_units=(vals<0).sum()))
    rng=np.random.default_rng(20261103);boots=[]
    for _,w in u.groupby('dataset_id'):
     vals=w.effect.to_numpy();boots.append(rng.choice(vals,(1999,len(vals)),replace=True).mean(1))
    boots=np.array(boots);bt=boots[rng.integers(0,len(ds),(len(ds),1999)),np.arange(1999)].mean(0);loo=[ds.drop(i).mean() for i in ds.index] if len(ds)>1 else [np.nan]
    out.append(dict(common,effect=ds.mean(),low=np.quantile(bt,.025),high=np.quantile(bt,.975),p=stats.ttest_1samp(ds,0).pvalue if len(ds)>=3 and ds.std()>0 else np.nan,n_sections=b.section_id.nunique(),n_units=len(u),n_verified_patients=b.loc[b.eligible_for_patient_balanced_analysis.eq(True),'patient_uid'].nunique(),n_datasets=len(ds),positive_units=(u.effect>0).sum(),negative_units=(u.effect<0).sum(),positive_datasets=(ds>0).sum(),negative_datasets=(ds<0).sum(),loo_min=np.nanmin(loo),loo_max=np.nanmax(loo)))
 a=pd.DataFrame(out);adjust(a,'p','q_context85',['source','population','family','program_family','adjustment']);save(a,O/'tables/context85_hierarchical_summary.tsv');save(units,O/'tables/context85_units.tsv.gz');save(datasets,O/'tables/context85_datasets.tsv.gz');print('Context hierarchy complete',flush=True)
if __name__=='__main__':main()
