from pathlib import Path
import pandas as pd,numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
font_manager.fontManager.addfont('/home/data/t070721/.local/share/fonts/corefonts-arial/Arial.TTF')
R=Path(__file__).resolve().parents[1];composition=pd.read_csv(R/'tables/niche_composition.tsv',sep='\t',dtype={'niche':str});E=R.parent/'hd_bin8_banksy_extension_v1';F=R/'figures';F.mkdir(exist_ok=True);plt.rcParams.update({'font.family':'Arial','font.size':12,'pdf.fonttype':42});rows=[]
for s in ['Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2']:
 source=E/s/'fixed_rois.tsv'
 if not source.exists():source=R.parent/'hd_bin8_comparison_v1'/s/'fixed_rois.tsv'
 rois=pd.read_csv(source,sep='\t');nodes=pd.read_parquet(R/s/'nodes.parquet')
 for _,roi in rois.iterrows():
  name=roi['roi'];x,y,X,Y=[float(roi[c]) for c in ['x','y','X','Y']];a=nodes[nodes.epi_anchor&nodes.x_um.between(x,X)&nodes.y_um.between(y,Y)];bg=nodes[nodes.x_um.between(x,X)&nodes.y_um.between(y,Y)]
  if not len(a):rows.append(dict(section=s,roi=name,anchors=0,x=x,y=y,X=X,Y=Y,note='No shared anchors; not substituted by a new favorable ROI'));continue
  fig,axs=plt.subplots(2,4,figsize=(16,8))
  for i,ref in enumerate(['blca','hpca']):
   for j,(method,radius) in enumerate([('A',50),('B',50),('A',200),('B',200)]):
    lab=pd.read_parquet(R/f'{ref}_r{radius}'/'labels.parquet');allids=sorted(set(lab[method+'_niche'])-{'unresolved'});lab=lab[lab.section==s].set_index('bin_id');z=lab.loc[a.bin_id];pal={k:plt.get_cmap('hsv')((n*.61803398875+(0 if method=='A' else .23+radius/1000+(0.11 if ref=='hpca' else 0)))%1) for n,k in enumerate(allids)};
    if method=='A' and len(allids)==2:
     cc=composition[(composition.reference==ref)&(composition.radius==radius)&(composition.method=='A')&(composition.niche!='unresolved')];high=cc.loc[cc.Epithelial.idxmax(),'niche'];pal={k:('#f5b264' if k==high else '#8752a3') for k in allids}
    pal['unresolved']='#444444';ax=axs[i,j];ax.scatter(bg.x_um,bg.y_um,c='#e8e0d0',s=6);ax.scatter(a.x_um,a.y_um,c=[pal[k] for k in z[method+'_niche']],s=9);ax.set_xlim(x,X);ax.set_ylim(Y,y);ax.set_aspect('equal');ax.set_xlabel(f'{ref} {method} {radius} μm');ax.set_ylabel('y (μm)');ax.plot([x+10,x+74],[Y-10,Y-10],c='black',lw=1.2);ax.text(x+10,Y-18,'64 μm',fontsize=9)
  fig.savefig(F/f'{s}_ROI_{name}.png',dpi=300,bbox_inches='tight');fig.savefig(F/f'{s}_ROI_{name}.pdf',bbox_inches='tight');plt.close(fig)
  raw=E/s/'figures'/f'ROI_single_genes_{name}.png';rows.append(dict(section=s,roi=name,anchors=len(a),x=x,y=y,X=X,Y=Y,original_raw_marker_figure=str(raw.relative_to(R.parent)) if raw.exists() else 'See upstream fixed-ROI raw expression figures',note='Matched existing ROI; center labels only; no anatomy assertion for MPBC'))
pd.DataFrame(rows).to_csv(R/'tables/ROI_review_index.tsv',sep='\t',index=False)
