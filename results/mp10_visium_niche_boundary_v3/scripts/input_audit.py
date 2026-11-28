from common3 import *
import importlib.metadata as im
import pyarrow.parquet as pq
remote=read(O/'audit/reviewed_git_blobs.tsv').set_index('path');local=read(V/'file_manifest.tsv').set_index('path')
inputs=['README.md','METHODS.md','RESULTS.md','config.json','file_manifest.tsv','classification/section_classification.tsv','classification/immune_input.tsv.gz','classification/immune_windows.tsv.gz','classification/window_links.tsv.gz','audit/input_manifest.tsv','audit/section_clinical.tsv','audit/source_studies.tsv','audit/R12_dictionary.tsv','audit/score_coverage_audit.tsv','tables/section_bins.tsv.gz','tables/section_effects.tsv.gz','tables/paired_boundaries.tsv.gz','tables/region_niche.tsv.gz','tables/all_section_levels.tsv','tables/unit_exposure_groups.tsv','scripts/metrics.py','scripts/discover.py','scripts/base.py','scripts/measure.py','scripts/classify.py']+[str(Path(k)/(s+'.parquet')) for k in ['boundary','spots'] for s in M.index]
rows=[]
for n in inputs:
    p=V/n;data=p.read_bytes();gh=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\x00'+data).hexdigest();actual=hashlib.sha256(data).hexdigest()
    if n in local.index:assert actual==local.at[n,'sha256']
    expected=remote.git_sha1.get(n,None)
    rows.append(dict(path=n,bytes=len(data),local_sha256=actual,git_blob=gh,reviewed_git_blob=expected,status='matches_review_commit' if expected==gh else 'in_local_only_or_remote_unlisted' if expected is None else 'changed'))
save(rows,O/'audit/actual_inputs_vs_commit.tsv')
idx=read(O/'cache/point_index.tsv').merge(M.reset_index(),on='section_id',validate='1:1');save(idx,O/'cache/sections.tsv')
schema=pq.read_schema(V/'spots'/(M.index[0]+'.parquet'))
save([dict(field=f.name,dtype=str(f.type)) for f in schema],O/'cache/point_schema.tsv')
save([dict(package=k,version=im.version(k)) for k in ['pandas','numpy','scipy','scikit-learn','matplotlib','pyarrow','shapely'] ],O/'audit/software.tsv')
print(pd.DataFrame(rows).status.value_counts().to_dict())

