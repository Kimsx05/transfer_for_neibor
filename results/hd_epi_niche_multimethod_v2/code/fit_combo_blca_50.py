from pathlib import Path
import fcntl
root=Path(__file__).resolve().parents[1]
with (root/'logs/fit_combo_blca_50.lock').open('w') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX)
 from fit_composition import main
 main('blca',50)
