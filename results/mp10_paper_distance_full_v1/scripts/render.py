from program_figures import *
from concurrent.futures import ProcessPoolExecutor
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("phase",choices=["core","program"])
parser.add_argument("--workers",type=int,default=8)
parser.add_argument("--sections",nargs="*")
parser.add_argument("--resume",action="store_true")
args=parser.parse_args()
cols=([raw,"epi_MP10_model_universe"]+[t+"_fraction" for t in targets]+[t+"_P80_all_distance_um" for t in targets]) if args.phase=="core" else focus+["Treg_fraction","Treg_abundance"]
frames=[]
for s in secs:frames.append(pd.read_parquet(out/"spots"/(s+".parquet"),columns=cols))
allvals=pd.concat(frames,ignore_index=True)
scales.update({c:max(v,1e-9) if np.isfinite(v) else 1 for c,v in allvals.quantile(.99).items()})
save([dict(feature=k,lower=0,upper=v) for k,v in scales.items()],args.phase+"_display_scales.tsv")
del frames,allvals
def render(s):
    folder=out/"figures"/section_ds[s]/s
    names=["01_context","02_targets","03_distance","04_regions"] if args.phase=="core" else ["05_program_maps","06_program_maps","07_program_maps","08_program_trends","09_program_effects"]
    if args.resume and all((folder/(k+ext)).is_file() for k in names for ext in [".png",".pdf","_legend_only.pdf","_panel_only.pdf"]):return s
    D[s]=pd.read_parquet(out/"spots"/(s+".parquet"))
    if args.phase=="core":core(s)
    else:maps(s);trends(s)
    del D[s]
    return s
if __name__=="__main__":
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for sec in pool.map(render,args.sections or secs):print(sec,flush=True)
