from common import *
import argparse
p=argparse.ArgumentParser();p.add_argument('--resume',action='store_true');p.add_argument('--verify-only',action='store_true');args=p.parse_args()
if args.verify_only:
 from validate import validate
 validate(final=True)
 from finish import verify_packages
 verify_packages();log('VERIFY_ONLY_PASSED: no permutations rerun')
else:
 from prepare import prepare,geometry
 from hotspots import pilot,full as hotspots_full,materialize_pilot_point_audit
 from relations import full as relations_full
 from figures import full as figures_full
 from validate import validate
 from finish import finish
 from sensitivity import enrich
 prepare();geometry();pilot()
 if not done('M03_ALL_HOTSPOTS_COMPLETE'):hotspots_full()
 materialize_pilot_point_audit()
 if not done('M06_RELATIONS_STATISTICS_COMPLETE'):relations_full()
 if not done('M07_FIGURES_COMPLETE'):figures_full()
 enrich();validate(final=True);finish()
