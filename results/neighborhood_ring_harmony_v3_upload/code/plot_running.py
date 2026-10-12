from base import *
import matplotlib.pyplot as plt
F=OUT/'04_comparison/figures';F.mkdir(exist_ok=True)
for od in sorted((OUT/'02_R12_pilot').iterdir()):
 if not od.is_dir() or not (od/'objective.tsv').exists():continue
 a=pd.read_csv(od/'objective.tsv',sep='\t');fig,ax=plt.subplots(figsize=(6,4));ax.plot(a.iteration,a.objective,marker='o',markersize=3);ax.set(xlabel='Harmony iteration',ylabel='Objective',title=od.name);fig.tight_layout();fig.savefig(F/(od.name+'__objective.png'),dpi=160,transparent=True);fig.savefig(F/(od.name+'__objective.pdf'),transparent=True);plt.close(fig)
