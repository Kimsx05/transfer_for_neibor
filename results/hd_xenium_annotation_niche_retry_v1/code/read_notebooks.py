from pathlib import Path
import json,urllib.request,subprocess,time
O=Path(__file__).resolve().parents[1]
ar=[]
repo='/tmp/blca-crohn-nov27'
commit=subprocess.check_output(['git','-C',repo,'rev-parse','HEAD']).decode().strip()
files=subprocess.check_output(['git','-C',repo,'ls-tree','-r','--name-only','HEAD','xenium']).decode().splitlines()
for f in files:
 if not f.endswith('.ipynb') or not Path(f).name.startswith(('2_','3a','3b','4a','4b','4e','5_')):continue
 for retry in range(3):
  try:
   url='https://raw.githubusercontent.com/Teichlab/crohns-fibroblast-atlas/'+commit+'/'+f
   obj=json.load(urllib.request.urlopen(url,timeout=80))
   ar.append(f)
   for c in obj['cells']:
    if c['cell_type']=='code':
     ar.extend(l.strip() for l in ''.join(c['source']).splitlines() if any(w in l.lower() for w in ['celltypist','scvi.','aggregate_neighbors','cellcharter','cc.tl.','cc.gr.','sc.tl.leiden','manual','read_csv']))
   break
  except Exception as e:print(f,retry,repr(e),flush=True)
 (O/'shared/crohn_notebook_settings.txt').write_text(commit+'\n'+'\n'.join(ar))

