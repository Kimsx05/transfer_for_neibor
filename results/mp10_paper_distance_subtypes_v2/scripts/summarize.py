from run import *
from scipy.stats import ttest_1samp
def bh(p):
 p=np.asarray(p);q=np.full(len(p),np.nan);idx=np.where(np.isfinite(p))[0];order=idx[np.argsort(p[idx])]
 q[order]=np.minimum.accumulate((p[order]*len(order)/np.arange(1,len(order)+1))[::-1])[::-1].clip(0,1);return q
def main():
 for name in ['effects','targets','interface','niche','strata']:
  frames=[]
  for sec in M.section_id:
   try:frames.append(pd.read_csv(O/'tables'/f'{sec}_{name}.tsv.gz',sep='\t'))
   except pd.errors.EmptyDataError:pass
  a=pd.concat(frames,ignore_index=True)
  if name=='effects':
   a['family']=np.where(a.kind=='interaction','interaction',np.where(a.axis=='MP10','MP10_program',np.where(a.program==RAW,'MP10_distance','program_distance')))
   a.loc[a.axis.isin(['local_interface','observed_outer','strict_outer']),'family']='interface_'+a.kind
   a['reference']=~a.axis.str.contains('_P95$')|a.axis.str.startswith(('CAF_','TAM_'))
   a['q']=a.groupby(['section_id','family','reference'],group_keys=False).p.transform(lambda p:bh(p))
  write(a,O/'tables'/f'{name}.tsv.gz')
 e=pd.read_csv(O/'tables/effects.tsv.gz',sep='\t').merge(M[['section_id','dataset_id','source_context','analysis_unit_uid','patient_uid','eligible_for_patient_balanced_analysis']],on='section_id')
 group=['kind','axis','program']
 e['value']=np.where(e.kind=='interaction',e.beta,e.adj)
 e.loc[e.kind=='matched','value']=e.loc[e.kind=='matched','r']
 out=[];units=[]
 for source,g in e.groupby('source_context'):
  for key,z in g.groupby(group):
   for population,base,idcol in [('sections',z,'section_id'),('units',z[z.analysis_unit_uid.notna()],'analysis_unit_uid'),('verified',z[z.eligible_for_patient_balanced_analysis.eq(True)&z.patient_uid.notna()],'patient_uid')]:
    b=base.dropna(subset=['value']).copy()
    if not len(b):continue
    isr=key[0]!='interaction';b['z']=np.arctanh(b.value.clip(-.999999,.999999)) if isr else b.value
    unit=b.groupby(['dataset_id',idcol]).z.mean().reset_index()
    for row in unit.itertuples(index=False):units.append(dict(source=source,population=population,kind=key[0],axis=key[1],program=key[2],dataset=row[0],unit=row[1],value=np.tanh(row[2]) if isr else row[2]))
    ds=unit.groupby('dataset_id').z.mean();v=ds.mean()
    rng=np.random.default_rng(20261103);draws=[]
    for dsid,sub in unit.groupby('dataset_id'):
     vals=sub.z.to_numpy();draws.append(rng.choice(vals,size=(399,len(vals)),replace=True).mean(1))
    draws=np.array(draws);boots=draws[rng.integers(0,len(ds),size=(len(ds),399)),np.arange(399)].mean(0)
    lo,hi=np.quantile(boots,[.025,.975])
    out.append(dict(source=source,population=population,kind=key[0],axis=key[1],program=key[2],n_sections=b.section_id.nunique(),n_units=len(unit),n_datasets=len(ds),effect=np.tanh(v) if isr else v,low=np.tanh(lo) if isr else lo,high=np.tanh(hi) if isr else hi,positive_datasets=(ds>0).sum(),negative_datasets=(ds<0).sum(),p=ttest_1samp(ds,0).pvalue if len(ds)>=3 else np.nan))
 a=pd.DataFrame(out);a['family']=np.where(a.kind=='interaction','interaction',np.where(a.axis=='MP10','MP10_program',np.where(a.program==RAW,'MP10_distance','program_distance')))
 a.loc[a.axis.isin(['local_interface','observed_outer','strict_outer']),'family']='interface_'+a.kind
 a['reference']=~a.axis.str.endswith('_P95')|a.axis.str.startswith(('CAF_','TAM_'))
 a['q']=a.groupby(['source','population','family','reference']).p.transform(lambda v:bh(v))
 write(a,O/'tables/summary.tsv');write(units,O/'tables/unit_effects.tsv.gz')
 genes=set('HMGCS1 MSMO1 INSIG1 FDFT1 SQLE FDPS IDI1 ACSS2 ELOVL5 LDLR MVD SCD HMGCR ACAT2 SLC20A1 ACSL3 ACSL1 ADIPOR2 MXD1 NPC1'.split())
 write([dict(program=n,n_genes=g.gene_symbol.nunique(),n_overlap=len(genes&set(g.gene_symbol)),overlap='|'.join(sorted(genes&set(g.gene_symbol)))) for n,g in SETS.groupby('gs_name')],O/'tables/gene_overlap.tsv')
if __name__=='__main__':main()
