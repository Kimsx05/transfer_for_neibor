import os,sys,json,time,subprocess,concurrent.futures
from pathlib import Path
import csv
R=Path(__file__).resolve().parents[2];O=R/'mp10_pathology_regions_v2';rows=list(csv.DictReader(open(O/'tables/section_manifest.tsv'),delimiter='\t'));pilot=json.loads((O/'tables/pilot_sections.json').read_text());phase=sys.argv[1] if len(sys.argv)>1 else 'pilot';secs=pilot if phase=='pilot' else [r['section_id'] for r in rows];workers=int(sys.argv[2]) if len(sys.argv)>2 else 2;lookup={r['section_id']:r for r in rows}
env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
def run(sec):
 dest=Path(lookup[sec]['output_dir'])
 if (dest/'spacet.done').exists() and (dest/'banksy.done').exists():return sec,'cached'
 start=time.time()
 with open(O/'logs'/f'{sec}_algorithms.log','w') as log:
  p=subprocess.run(['/usr/bin/time','-v','/usr/local/bin/Rscript','--vanilla',str(O/'scripts/run_section_two.R'),sec,'6'],cwd=R,env=env,stdout=log,stderr=subprocess.STDOUT)
 status='complete' if (dest/'spacet.done').exists() and (dest/'banksy.done').exists() else 'partial_or_failed'
 print(time.strftime('%H:%M:%S'),sec,status,round(time.time()-start,1),flush=True);return sec,status
with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:list(ex.map(run,secs))
