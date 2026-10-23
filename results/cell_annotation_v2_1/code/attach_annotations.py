"""Read original v2 object, align v2.1 proposals by cell_id, optionally save a NEW object.
Usage: python attach_annotations.py /path/to/v2/SECTION/cells.h5ad --output /new/path.h5ad
Without --output, only validates and prints the aligned annotation summary.
"""
from pathlib import Path
import argparse,pandas as pd,anndata as ad,numpy as np
p=argparse.ArgumentParser();p.add_argument('input');p.add_argument('--labels',default=str(Path(__file__).resolve().parents[1]/'cell_annotations.parquet'));p.add_argument('--output');args=p.parse_args();source=Path(args.input).resolve();table=pd.read_parquet(args.labels);assert table.cell_id.is_unique
a=ad.read_h5ad(source,backed='r' if not args.output else None);ids=a.obs['cell_id'].astype(str) if 'cell_id' in a.obs else pd.Series(a.obs_names,index=a.obs_names);q=table.set_index('cell_id').reindex(ids.to_numpy());assert q.proposed_major.notna().all() and len(q)==a.n_obs
cols=['proposed_major','proposed_subtype','display_label','state','annotation_support','evidence_source','unknown_reason','quality_flag','artifact_review_flag','global_cluster_v21','lineage_cluster_v21','review_status','review_entry']
print(q.proposed_major.value_counts().to_string());print('Aligned by cell_id:',len(q))
if args.output:
 out=Path(args.output).resolve();assert out!=source and not out.exists(),'Choose a new output path; inputs and existing files cannot be overwritten.'
 for c in cols:a.obs['v21_'+c]=q[c].to_numpy()
 a.write_h5ad(out,compression='gzip');print(out)
else:a.file.close()
