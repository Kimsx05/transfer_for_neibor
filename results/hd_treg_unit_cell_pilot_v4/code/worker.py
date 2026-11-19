from base import *
import importlib.util,sys
spec=importlib.util.spec_from_file_location('night_supervise',A/'hd_epi_niche_multimethod_v2/code/night_supervise.py')
w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w);w.R=R;w.py=sys.executable
w.run('model')
