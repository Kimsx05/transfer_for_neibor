import os,sys,time,signal,json
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
pid=int(sys.argv[1]);out=Path(__file__).resolve().parents[1]/'logs/runtime.json';paused=False
while True:
 try:os.kill(pid,0)
 except ProcessLookupError:break
 now=datetime.now(ZoneInfo('Asia/Shanghai'));allowed=now.hour>=23 or now.hour<7
 if not allowed and not paused:os.kill(pid,signal.SIGSTOP);paused=True
 if allowed and paused:os.kill(pid,signal.SIGCONT);paused=False
 out.write_text(json.dumps(dict(time=now.isoformat(),pid=pid,state='paused_outside_window' if paused else 'running_not_complete')))
 time.sleep(30)
out.write_text(json.dumps(dict(time=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),state='process_exited_check_outputs')))
