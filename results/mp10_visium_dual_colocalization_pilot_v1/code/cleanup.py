from common import *
for folder in ['Epi','TME',*S]:
 p=O/'.work'/folder
 for f in ['genes','ids','x','i','p']:(p/f).unlink(missing_ok=True)
 if p.exists() and not any(p.iterdir()):p.rmdir()
(O/'.work/helper_fragment').unlink(missing_ok=True)

