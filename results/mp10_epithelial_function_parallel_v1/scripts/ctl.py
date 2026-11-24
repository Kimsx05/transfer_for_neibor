from membrane import diff,targets,f
from public_common import *
import gzip,re
a='GSE149933';samples=[]
for line in gzip.open(f(a,'*soft*'),'rt'):
 l=line.strip()
 if l.startswith('^SAMPLE ='):samples.append(dict(GSM=l.split(' = ')[1]))
 if l.startswith('!Sample_title ='):
  title=l.split(' = ')[1];r=re.match(r'(MUS_\d+)_(\S+)',title);samples[-1].update(title=title,col=r[1],model=r[2])
 if l.startswith('!Sample_characteristics_ch1 ='):
  key,value=l.split(' = ',1)[1].split(': ',1);samples[-1][key]=value
meta=pd.DataFrame(samples);assert len(meta)==175 and meta.col.is_unique
output(meta,'tables/GSE149933_samples.tsv')
orth=read('tables/mouse_orthology.tsv');valid=orth[orth.eligible==True].drop_duplicates(['symbol','human_symbol']).set_index('symbol').human_symbol
results=[];coverage=[];contrasts=[];guide_qc=[]
for library,tag in [('mTKO','mTKO'),('mTKO-validation','mVal')]:
 lib=pd.read_csv(f(a,'*'+tag+'_library*'),sep='\t');raw=pd.read_csv(f(a,'*'+tag+'_readCounts*'),sep='\t')
 x=raw.merge(lib[['GUIDE_ID','SEQUENCE','GENE']].drop_duplicates(),left_on=['Guide.ID','Guide.Seq','Gene.symbol'],right_on=['GUIDE_ID','SEQUENCE','GENE'],how='left')
 assert len(x)==len(raw)
 guide_qc.append([tag,len(raw),len(lib),x.GUIDE_ID.notna().sum(),raw['Guide.ID'].duplicated().sum()])
 m=meta[meta['crispr library']==library];assert set(m.col)<=set(raw.columns)
 extra=sorted(set(raw.filter(regex='^MUS_').columns)-set(m.col));output(pd.DataFrame({'unmapped_column':extra}),'tables/GSE149933_'+tag+'_unmapped.tsv')
 counts=x[m.col];norm=np.log2(counts.div(counts.sum())*1e6+.5)
 human=x.GENE.map(valid);gok=human.notna();norm=norm[gok];norm.index=human[gok];by=norm.groupby(level=0).median()
 for g in targets:
  coverage.append([tag,g,','.join(orth.loc[orth.human_symbol==g,'symbol'].dropna().astype(str)),(human==g).sum()])
 for model,mm in m.groupby('model'):
  arms=mm[~mm.treatment.isin(['Control','none'])]
  ctrl=mm[mm.treatment=='Control']
  if not len(arms):continue
  def times(z):return sorted(z['time point'].unique(),key=lambda t:int(t[1:]))
  for tr,tt in arms.groupby('treatment'):
   assert len(times(tt))==len(times(ctrl))==2
   for stage,t1,t0 in zip(['mid','end'],times(tt),times(ctrl)):
    test=tt[tt['time point']==t1].sort_values('title');base=ctrl[ctrl['time point']==t0].sort_values('title')
    assert len(test)==len(base)==3
    start=mm[mm.treatment=='none'];contrasts.append([tag,model,stage,tr,t1,t0,','.join(test.col),','.join(base.col),','.join(start.col)])
    # Median guide log abundance per replicate. Welch intervals model replicate variability; guides are not replicates.
    for gene,row in by.iterrows():
     d=diff(row[test.col],row[base.col])
     results.append(dict(library=tag,model=model,stage=stage,treatment=tr,test_day=t1,control_day=t0,gene=gene,n_test=len(test),n_control=len(base),control_minus_T0=row[base.col].mean()-row[start.col].mean(),**d))
output(pd.DataFrame(guide_qc,columns=['library','count_rows','library_rows','matched_id_seq_gene','duplicate_ids']),'tables/GSE149933_guide_QC.tsv')
output(pd.DataFrame(contrasts,columns=['library','model','stage','arm','test_day','ctrl_day','test','ctrl','T0']),'tables/GSE149933_contrasts.tsv')
output(pd.DataFrame(coverage,columns=['library','gene','mapped_mouse','guides']),'tables/GSE149933_coverage.tsv')
d=adjust(pd.DataFrame(results),['library','model','stage','treatment']);output(d,'tables/GSE149933_candidates.tsv')

