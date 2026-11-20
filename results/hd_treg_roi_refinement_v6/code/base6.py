from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
V5=ROOT.parent/'hd_treg_roi_refinement_v5'
sys.path.insert(1,str(V5/'code'))
from common5 import *
R=ROOT if os.environ.get("HD6_BRANCH","conservative")=="conservative" else ROOT/"proseg_observed"
direct=T[:4]
def tab(df,name):put(df,R/'tables'/name)
