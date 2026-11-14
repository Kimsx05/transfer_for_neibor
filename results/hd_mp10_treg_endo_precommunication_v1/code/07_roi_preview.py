from common import *
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
Image.MAX_IMAGE_PIXELS=None
import matplotlib.font_manager as fm
fm.fontManager.addfont('/home/data/t070721/.local/share/fonts/corefonts-arial/Arial.TTF');plt.rcParams.update({'font.family':'Arial','font.size':12,'axes.linewidth':1,'pdf.fonttype':42,'figure.facecolor':'white','savefig.facecolor':'white'})
d=pd.read_parquet(R/'tables/Cervilla_bin8_measurement.parquet');ro=pd.read_csv(R/'tables/roi_manifest.tsv',sep='\t');ass=json.loads((A/'segmentation_assignment_v2/Cervilla/metadata/assets.json').read_text());im=Image.open(H/ass['image_project_relative']);mpp=.8636122748145504
for r in ro.itertuples():
 q=d[d.x_um.between(r.x0,r.x1)&d.y_um.between(r.y0,r.y1)];fig,axes=plt.subplots(1,5,figsize=(15,3),layout='constrained');crop=im.crop((int(r.x0/mpp),int(r.y0/mpp),int(r.x1/mpp),int(r.y1/mpp)));axes[0].imshow(crop,extent=[r.x0,r.x1,r.y1,r.y0]);axes[0].set_xlabel('H&E')
 for ax,key in zip(axes[1:],['Epi','T','Treg','Endo']):
  val=np.log1p(q[key+'_UMI']);v=ax.scatter(q.x_um,q.y_um,c=val,s=12,marker='s',cmap='viridis',vmin=0,rasterized=True);fig.colorbar(v,ax=ax,shrink=.7,label='log1p raw UMI');ax.set_xlabel(key+' panel')
 for ax in axes:ax.set_xlim(r.x0,r.x1);ax.set_ylim(r.y1,r.y0);ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([])
 axes[0].plot([r.x0+10,r.x0+60],[r.y1-12]*2,'k-',lw=2);axes[0].text(r.x0+35,r.y1-17,'50 μm',ha='center',fontsize=10)
 fig.savefig(R/'figures'/f'ROI_raw_preview_{r.roi}.png',dpi=300);plt.close(fig)
print(ro.to_string(index=False))
