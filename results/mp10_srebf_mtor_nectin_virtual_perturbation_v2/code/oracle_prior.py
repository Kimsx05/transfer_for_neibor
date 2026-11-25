from common import *
import pyarrow.parquet as pq
from collections import defaultdict
f=Path('cache/promoter.parquet')
assert hashlib.sha1(b'blob '+str(f.stat().st_size).encode()+b'\0'+f.read_bytes()).hexdigest()==C['celloracle']['prior_blob']
h=rd(prior/'hgnc_complete_set_20261118.txt',low_memory=False).query("status=='Approved'").set_index('symbol')
alias=defaultdict(set)
for col in ['prev_symbol','alias_symbol','ensembl_gene_id']:
 for g,s in h[col].dropna().items():
  for a in s.split('|'):alias[a].add(g)
def canonical(g):
 if g in h.index:return g
 u=alias.get(g,());return next(iter(u)) if len(u)==1 else None
genes=set(rd('tables/panel_members.tsv').query('P_clean').gene);pf=pq.ParquetFile(f)
tf={s:canonical(s) for s in pf.schema.names if s not in ['peak_id','gene_short_name']}
used=[s for s,g in tf.items() if g in genes]
t=pd.read_parquet(f,columns=['gene_short_name',*used]);t['target']=t.gene_short_name.map(canonical)
t=t[t.target.isin(genes)].groupby('target')[used].max()
rows=[(tf[s],g) for s in used for g in t.index[t[s].values>0] if tf[s]!=g]
edges=pd.DataFrame(sorted(set(rows)),columns=['source','target']);save(edges,'cache/base_edges.tsv.gz')
save(pd.DataFrame([dict(species='human',genome='hg38',database='CellOracle hg38_gimmemotifsv5_fpr2_threshold_10_20210630',blob_sha1=C['celloracle']['prior_blob'],sha256=sha(f),peak_rows=pf.metadata.num_rows,raw_regulator_columns=len(tf),matched_TF_columns=len(used),mapped_base_edges=len(edges),sources=edges.source.nunique(),targets=edges.target.nunique())]),'qc/promoter_coverage.tsv')

