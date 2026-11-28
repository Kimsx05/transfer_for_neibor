from common3 import *
import matplotlib.pyplot as plt
from PIL import Image
theme();w=pd.read_parquet(O/'cache/immune_windows.parquet').query("grid==@GRID and eligible");cl=read(OUT/'classification/section_classification.tsv').set_index('section_id')
qa=[];fmax=w['frac_'+AXIS].max()
for s,r in cl.iterrows():
    d=spot(s,['x_um','y_um','block','all_E','region_main','window500','T_f']);a=w[w.section_id==s].set_index('window')
    width=1000 if GRID=='1000' else 500;dx=250 if GRID in ['x250','xy250'] else 0;dy=250 if GRID in ['y250','xy250'] else 0
    grid_id=d.block.astype(str)+'|'+np.floor((d.x_um-dx)/width).astype(int).astype(str)+'|'+np.floor((d.y_um-dy)/width).astype(int).astype(str)
    f,ax=plt.subplots(1,4,figsize=(13.5,3.2),layout='constrained')
    if Path(r.image_path).is_file():ax[0].imshow(Image.open(r.image_path))
    ax[0].set_title('H&E · native image',fontsize=9)
    for i,key,lbl in [(1,'frac_T','T fraction inside E'),(2,'neighbor_T','R1 T fraction around E'),(3,'penetration_T','Entry / periphery (log2)')]:
        vals=grid_id.map(a[key.replace('_T','_'+AXIS)]).where(d.all_E)
        ax[i].scatter(d.x_um,d.y_um,c='#dfdfdc',s=3,lw=0)
        ok=vals.notna();art=ax[i].scatter(d.x_um[ok],d.y_um[ok],c=vals[ok],s=4,lw=0,cmap='viridis' if i<3 else 'coolwarm',vmin=0 if i<3 else None,vmax=fmax if i==1 else None);f.colorbar(art,ax=ax[i],fraction=.04);ax[i].set_title(lbl.replace('T fraction',AXIS+' fraction'),fontsize=9);ax[i].invert_yaxis()
    for aa in ax:aa.set_aspect('equal');aa.axis('off')
    f.suptitle(s+' · '+r.exploratory_cluster+' · assignment frequency '+str(round(r.stability,2)),fontsize=11);f.savefig(OUT/'review'/(s+'.jpg'),dpi=125);plt.close(f)
    qa.append(dict(section_id=s,group=r.exploratory_cluster,medoid=r.medoid,stability=r.stability,eligible_windows=len(a),fraction_min=a['frac_'+AXIS].min(),fraction_max=a['frac_'+AXIS].max(),file='review/'+s+'.jpg',he_available=Path(r.image_path).is_file()))
save(qa,OUT/'review/index.tsv')

