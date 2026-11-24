from public_common import *
import csv,gzip
def file(a,q):return next((C/a).glob(q))
a='GSE182638'
md={}
with gzip.open(file(a,'*series_matrix*'),'rt') as f:
 for row in csv.reader(f,delimiter='\t'):
  if not row:continue
  if row[0]=="!Sample_description":md["raw_column"]=row[1:]
  if row[0]=='!Sample_geo_accession':md['GSM']=row[1:]
  elif row[0]=='!Sample_title':md['title']=row[1:]
  elif row[0]=='!Sample_characteristics_ch1':
   key=row[1].split(': ')[0];md[key]=[x.split(': ',1)[1] for x in row[1:]]
p=pd.DataFrame(md);p['arm']=np.select([p.treatment.str.contains('Fer-1'),p.treatment.str.contains('RSL3')],['FR','RSL3'],'UT')
p['sample']=p['cell line']+'-T3-'+p['cell passage']+'-'+p.arm
assert (p.raw_column==p["sample"]).all()
assert set(p['sample'])==set(pd.read_csv(file(a,'*raw*'),sep='\t',nrows=1).columns[1:]);output(p,'tables/GSE182638_samples.tsv')
a='GSE126696';raw=pd.read_csv(file(a,'*786O*'),sep='\t');ref=pd.read_csv(file(a,'*reference*'),header=None,names=['bc','construct','seq'])
lib=pd.read_csv(file('Brunello_2016','*.txt'),sep='\t')
seq='sgRNA Target Sequence'; gene='Target Gene Symbol';raw=raw.rename(columns={'Construct Barcode':'bc','Construct IDs':'construct'})
assert not raw.bc.duplicated().any() and not ref.bc.duplicated().any()
raw=raw.merge(ref,on=['bc','construct'],how='left',validate='one_to_one')
unique=lib.drop_duplicates([seq,gene]);unique=unique[~unique[seq].duplicated(keep=False)]
raw['gene']=raw.seq.map(unique.set_index(seq)[gene]);rev=raw.seq.str.translate(str.maketrans('ACGT','TGCA')).str[::-1]
output(pd.DataFrame([dict(count_guides=len(raw),unique_barcodes=raw.bc.nunique(),reference_matched=raw.seq.notna().sum(),barcode_eq_sequence=(raw.bc==raw.seq).sum(),library_rows=len(lib),unique_sequences=lib[seq].nunique(),ambiguous_dropped=len(lib)-len(unique),forward_unique_match=raw.gene.notna().sum(),reverse_match=rev.isin(unique[seq]).sum())]),'tables/GSE126696_mapping.tsv')
cols=[x for x in raw if x.startswith('786-O')];counts=raw[cols];rpm=counts.div(counts.sum())*1e6;log=np.log2(rpm+.5)
targets=set(sum(sets.values(),[])+['CHMP2A','POR','SREBF2','SCAP'])
diff=log.iloc[:,1:].sub(log.iloc[:,0],axis=0);diff['gene']=raw.gene
u=diff.groupby('gene').median();perc=u.rank(pct=True)
arr=[]
for g in sorted(targets):
 for d in u.columns:
  b=diff.loc[diff.gene.eq(g),d]
  arr.append(dict(gene=g,sample=d,nguides=len(b),median_log2_RPM_ratio=b.median(),positive_guides=(b>0).sum(),percentile=perc.at[g,d] if g in u.index else np.nan,biological_replicates=1))
output(pd.DataFrame(arr),'tables/GSE126696_conditional.tsv')
output(pd.concat([raw[['bc','construct','gene']+cols],diff.drop(columns='gene').add_suffix('_delta')],axis=1).loc[raw.gene.isin(targets)],'tables/GSE126696_candidate_guides.tsv')

