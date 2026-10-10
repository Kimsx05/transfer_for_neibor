"""Read compatible existing libraries without importing system site packages."""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT.parent
PROJECT = RUN.parent
sys.path[:0] = [str(RUN/'runtime'), str(PROJECT/'cell2location_setup_v1/cell2location_env/lib/python3.12/site-packages')]
