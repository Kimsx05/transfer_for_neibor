from pathlib import Path
import csv,json,platform
D=Path(__file__).resolve().parents[1];src=D/'logs/resource_use.tsv'
values={}
with src.open() as f:
 for v in csv.DictReader(f,delimiter='\t'):
  if v['wall_s']:values[v['stage']]=v
res={'timing_scope':'sequential warm-cache replay, latest successful invocation per stage; excludes downloads and superseded attempts','wall_seconds':round(sum(float(v['wall_s']) for v in values.values()),2),'peak_RSS_GiB':round(max(int(v['peak_RSS_KiB']) for v in values.values())/(1024**2),3),'stages_measured':list(values),'workers_per_analysis':1,'disk_MiB':round(sum(p.stat().st_size for p in D.rglob('*') if p.is_file())/1024**2,1)}
(D/'logs/resource_summary.json').write_text(json.dumps(res,indent=2))

