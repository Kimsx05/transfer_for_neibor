from common import *
import time,subprocess
expected=['S30_FULL_DEV','S31_TARGET_UMI','S31_FIXED_CAPTURE','S30_TARGET_UMI','S30_FIXED_CAPTURE']
while True:
 completed=[t for t in expected if (R/f'05_deconvolution/{t}/MAPPING_DONE.json').exists()]
 js(R/'logs/finish_analysis_progress.json',{'stage':'WAITING_MAPPING','completed':completed,'pending':[t for t in expected if t not in completed],'time':datetime.datetime.now().isoformat()})
 if len(completed)==len(expected):break
 errors=list((R/'logs').glob('ERROR_*.json'))
 if errors:raise RuntimeError('Mapping error records require audit: '+str(errors))
 time.sleep(15)
for script in ['06_evaluate.py','07_reference_and_history_audit.py','09_additional_evidence.py','08_figures.py','10_extra_figures.py','11_validate_delivery.py']:
 js(R/'logs/finish_analysis_progress.json',{'stage':script,'time':datetime.datetime.now().isoformat()});say('START',script)
 with open(R/'logs'/('final_'+script.replace('.py','.log')),'w') as f:p=subprocess.run(['bash',str(R/'code/python.sh'),str(R/'code'/script)],stdout=f,stderr=subprocess.STDOUT)
 if p.returncode:raise RuntimeError(script+' failed; see final log')
js(R/'logs/finish_analysis_progress.json',{'stage':'ANALYSIS_READY_FOR_REPORT','time':datetime.datetime.now().isoformat()});say('ANALYSIS_READY_FOR_REPORT')
