#!/bin/bash
set -eu
cd "$(dirname "$0")/.."
py="../../highres_VisiumHD/analysis/hd_epi_niche_multimethod_v2/.venv/bin/python"
"$py" -m pip install --no-deps --target shared/deps tacco==0.5.0 POT==0.9.7.post1 sparse==0.19.3 pyamg==5.4.1 fast-histogram==0.14
"$py" - <<'PY'
from pathlib import Path
p=Path('shared/deps/tacco/preprocessing/_qc.py');s=p.read_text().replace('profiles = profiles.to_numpy()','profiles = profiles.to_numpy(copy=True)');p.write_text(s)
p=Path('shared/deps/tacco/utils/_utils.py');s=p.read_text()
if '    cost = cost.copy()\n    cost /= cost.max()' not in s:s=s.replace('    cost /= cost.max()','    cost = cost.copy()\n    cost /= cost.max()')
p.write_text(s)
PY

