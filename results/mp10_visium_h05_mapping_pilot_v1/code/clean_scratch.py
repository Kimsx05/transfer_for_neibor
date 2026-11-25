from common import *
# Delete only this workflow's known sparse interchange files after verification.
root=Path(os.environ['PILOT_TMP'])
for name in ['reference',*S,'pseudo_medium','pseudo_low']:
 folder=root/name
 if not folder.is_dir():continue
 for k in ['genes','ids','x','p','i']:(folder/k).unlink(missing_ok=True)
 if not any(folder.iterdir()):folder.rmdir()
if root.is_dir() and not any(root.iterdir()):root.rmdir()

