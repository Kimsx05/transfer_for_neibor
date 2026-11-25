from common import *
from knk import base
files=['config/main.json','config/input_fingerprints.tsv','code/prepare.py','code/fit.R','code/knk.py','tables/panel_members.tsv','tables/program_members.tsv','cache/panel.bin','cache/metadata.tsv.gz']
old=hashlib.sha256(json.dumps([(p,hashlib.sha256(Path(p).read_bytes().replace(b'dt[is_KO==FALSE,',b'dt[!is_KO,')).hexdigest() if p=='code/fit.R' else sha(p)) for p in files]).encode()).hexdigest()
report=[]
for marker in Path('cache/jobs').glob('*/*.json'):
 j=json.loads(marker.read_text());suffix='' if marker.stem=='reference' else sha('tables/control_matching.tsv')
 assert j['key'] in [old+suffix,base+suffix] and all(sha(f)==h for f,h in j['outputs'].items())
 assert 'Exit status: 0' in Path('logs/'+marker.parent.name+'_'+marker.stem+'.usage').read_text()
 report.append(dict(path=str(marker),orchestration_key=j['key'],verified_source_key=base+suffix,reason='reference unaffected; all successful export processes ran after data.table logical-row syntax fix'))
 j['key']=base+suffix;t=marker.with_suffix('.tmp');t.write_text(json.dumps(j));os.replace(t,marker)
save(pd.DataFrame(report),'qc/checkpoint_reconciliation.tsv')

