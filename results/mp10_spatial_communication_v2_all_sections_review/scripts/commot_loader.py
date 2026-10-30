"""Load official COMMOT inference/database modules without optional downstream tools.
Source files are unmodified. Namespace package avoids importing tradeSeq/pysal/karateclub.
"""
from runtime import *
import types,importlib,importlib.util
root=V/'vendor/COMMOT/commot'
for name,path in [('commot',root),('commot.tools',root/'tools'),('commot.preprocessing',root/'preprocessing')]:
 m=types.ModuleType(name);m.__path__=[str(path)];m.__file__=str(path/'__init__.py');sys.modules[name]=m
spatial=importlib.import_module('commot.tools._spatial_communication')
database=importlib.import_module('commot.preprocessing._ligand_receptor_database')
# Timing only: wrap the four official solver stages without changing arguments/results.
import time,resource
cotmod=importlib.import_module('commot._optimal_transport._cot')
def _timed_stage(fn,name):
 def call(*a,**kw):
  t=time.time();print('COT_STAGE_START',name,flush=True)
  out=fn(*a,**kw)
  print('COT_STAGE_DONE',name,'seconds',round(time.time()-t,2),'peak_rss_GiB',round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,3),flush=True)
  return out
 return call
for name in ['cot_sparse','cot_row_sparse','cot_col_sparse','cot_blk_sparse']:
 setattr(cotmod,name,_timed_stage(getattr(cotmod,name),name))
