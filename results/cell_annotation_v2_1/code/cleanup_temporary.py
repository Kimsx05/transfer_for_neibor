from common import *
import pandas as pd
checks=pd.read_csv(ROOT/'metadata/delivery_checks.tsv',sep='\t');assert checks.passed.all()
removed=[]
for name in ['reference_shared_counts.h5','Cervilla_shared_counts.h5','D1_Johnson_shared_counts.h5']:
 p=ROOT/'temporary'/name
 if p.exists():removed.append({'file':str(p),'bytes':p.stat().st_size});p.unlink()
js(ROOT/'metadata/temporary_cleanup.json',{'removed_only_regenerable_v21_bridges':removed,'recreate_with':'code/prepare_inputs.py','original_v1_v2_inputs_touched':False})
print('Removed regenerable bridge bytes:',sum(x['bytes'] for x in removed))
