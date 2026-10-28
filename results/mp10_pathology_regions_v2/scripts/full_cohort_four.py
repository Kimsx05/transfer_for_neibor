import os,sys,time,subprocess,concurrent.futures,csv,json
from pathlib import Path
R=Path(__file__).resolve().parents[2];O=R/'mp10_pathology_regions_v2';rows=list(csv.DictReader(open(O/'tables/section_manifest.tsv'),delimiter='\t'));rows.sort(key=lambda r:-int(r['n_spots']));env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
adopt={v['cmd'].split('--args ')[1].split()[0]:v['pid'] for v in json.load(open(O/'logs/concurrency_four.json'))['kept']}
rows.sort(key=lambda r:(r['section_id'] not in adopt,-int(r['n_spots'])))
def job(r):
 sec=r['section_id'];p=Path(r['output_dir']);start=time.time()
 if sec in adopt:
  pid=adopt[sec]
  while True:
   try:
    stat=Path(f'/proc/{pid}/stat').read_text().rsplit(')',1)[1].split()
    if stat[0]=='Z':break
   except FileNotFoundError:break
   time.sleep(5)
 if not ((p/'spacet.done').exists() and (p/'banksy.done').exists()):
  with open(O/'logs'/f'{sec}_full.log','w') as log:subprocess.run(['/usr/bin/time','-v','/usr/local/bin/Rscript','--vanilla',str(O/'scripts/run_section_two.R'),sec,'6'],cwd=R,env=env,stdout=log,stderr=subprocess.STDOUT)
 with open(O/'logs'/f'{sec}_annotation.log','w') as log:subprocess.run(['/usr/bin/python3.12',str(O/'scripts/annotate_plot.py'),sec],cwd=R,env=env,stdout=log,stderr=subprocess.STDOUT)
 print(time.strftime('%H:%M:%S'),sec,'SpaCET',int((p/'spacet.done').exists()),'BANKSY',int((p/'banksy.done').exists()),'plots',int((p/'annotation.done').exists()),'seconds',round(time.time()-start),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:list(ex.map(job,rows))
