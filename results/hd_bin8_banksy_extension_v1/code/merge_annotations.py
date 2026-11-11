from pathlib import Path
import pandas as pd,sys,json
r=Path(__file__).resolve().parents[1];s=sys.argv[1];o=r/s
p=pd.read_csv(o/'predictions.tsv.gz',sep='\t');d=pd.read_csv(o/'bin_metadata.tsv.gz',sep='\t');d=d[d.qc_keep].copy();assert d.bin_id.reset_index(drop=True).equals(p.bin_id);assert p.bin_id.is_unique;assert {'singler_hpca_marker_count','epi_fib_ge2','marker_conflict'}.issubset(p.columns)
d=d.merge(p,on='bin_id',sort=False,validate='one_to_one');summary=json.loads((o/'annotation_summary.json').read_text());assert len(d)==summary['n'] and int(d.nCount_RNA.sum())==summary['UMI'];tmp=o/'analysis_metadata.partial.tsv.gz';d.to_csv(tmp,sep='\t',index=False);tmp.replace(o/'analysis_metadata.tsv.gz');print('METADATA_MERGE_COMPLETE',s,len(d))
