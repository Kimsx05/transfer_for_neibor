"""Sample proportional set size of this workflow and descendants, keeping fork-shared pages proportional."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import os
out=Path(__file__).resolve().parents[1];tree={};workers=set()
for p in Path('/proc').glob('[0-9]*'):
 try:
  if p.stat().st_uid!=os.getuid():continue
  info=dict(line.split(':',1) for line in (p/'status').read_text().splitlines() if ':' in line)
  cmd=(p/'cmdline').read_bytes().split(b'\0')
  pid=int(p.name);tree[pid]=int(info['PPid'])
  if info['Name'].strip() in ('R','python','python3','python3.12') and any(str(out/'code').encode() in arg.replace(b'~+~',b' ') for arg in cmd[1:]) and pid!=os.getpid():workers.add(pid)
 except OSError:pass
for i in range(7):workers.update(k for k,v in tree.items() if v in workers)
pss=rss=0
for pid in workers:
 try:
  vals=dict((x.split(':')[0],int(x.split()[1])) for x in (Path('/proc')/str(pid)/'smaps_rollup').read_text().splitlines() if x.startswith(('Pss:','Rss:')))
  pss+=vals['Pss'];rss+=vals['Rss']
 except (OSError,KeyError): pass
avail=next(int(x.split()[1]) for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:'))
p=out/'logs/resources.tsv';new=not p.exists()
with p.open('a') as f:
 if new:f.write('timestamp\tprocesses\tPSS_GiB\tRSS_GiB\thost_available_GiB\n')
 f.write(f'{datetime.now(ZoneInfo("Asia/Shanghai")).isoformat()}\t{len(workers)}\t{pss/1048576:.2f}\t{rss/1048576:.2f}\t{avail/1048576:.2f}\n')

