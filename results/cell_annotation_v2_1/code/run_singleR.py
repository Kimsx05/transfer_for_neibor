from common import *
import subprocess
# Only dependency availability is checked here; the completed repair is recorded in logs/sessionInfo.
ready=subprocess.run(['Rscript','-e','quit(status=if(requireNamespace("scrapper",quietly=TRUE))0 else 1)'],capture_output=True)
if ready.returncode:
 raise RuntimeError('SingleR wilcox/aggregation requires scrapper. Restore the documented existing R library or install this missing dependency; do not rebuild the project environment.')
raise SystemExit(subprocess.call(['Rscript',str(ROOT/'code/run_singler.R'),str(ROOT)]))
