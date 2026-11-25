#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
source code/env.sh
export XDG_CONFIG_HOME="$PWD/tmp/xdgconfig" XDG_CACHE_HOME="$PWD/tmp/xdgcache" MPLCONFIGDIR="$PWD/tmp/mpl" NUMBA_CACHE_DIR="$PWD/tmp/numba"
# Preserve verified prerequisite files and all previous results. An input mismatch stops reuse.
if [ -f config/preparation_checkpoint.json ]; then
 "$VP_PY" code/check_inputs.py
else
 "$VP_PY" code/prepare.py
 "$VP_PY" code/oracle_prior.py
 "$VP_PY" code/check_inputs.py --seal
fi
"$VP_PY" code/knk.py
if [ ! -x cache/co310/bin/python ]; then bash code/bootstrap_oracle.sh; fi
"$VP_PY" code/oracle_qc_compat.py
"$VP_PY" code/run_oracle.py
"$VP_PY" code/knk_tables.py
"$VP_PY" code/module.py
"$VP_PY" code/diagnostics.py
"$VP_PY" code/extra_summaries.py
"$VP_PY" code/oracle_tables.py
"$VP_PY" code/finish.py
"$VP_PY" code/plot_inputs.py
Rscript --vanilla code/replot.R
"$VP_PY" code/verify.py

