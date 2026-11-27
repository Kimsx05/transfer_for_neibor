from lib import *
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Patch
from PIL import Image
import argparse
p=argparse.ArgumentParser();p.add_argument('--version',default='primary');arg=p.parse_args()
out=O if arg.version=='primary' else O/'versions'/arg.version
(out/'review').mkdir(exist_ok=True)
state=read(out/'classification/window_states.tsv.gz').query("level=='window500'")
cl=read(out/'classification/section_classification.tsv').set_index('section_id')
colors=dict(ENN='#cb645c',ENS='#db994d',EIN='#b96eae',EIS='#875a92',NIN='#377f99',NIS='#56a59e',NNS='#94ae68',NNN='#bbb8ab',Unknown='#e6e6e6')
rows=[]
for sec,r in M.iterrows():
    d=load(f'spots/{sec}.parquet');edge=load(f'boundary/{sec}.parquet');mm=load(f'metrics/{sec}_local_moran.parquet')
    st=state[state.section_id==sec].set_index('local_id')
    fig,aa=plt.subplots(2,5,figsize=(17,7.8));ax=aa.ravel()
    def draw(i,title,c,**kw):
        ax[i].set_title(title,fontsize=9)
        return ax[i].scatter(d.x_um,d.y_um,c=c,s=2.7,linewidths=0,**kw)
    img=Path(r.image_path)
    if img.is_file():ax[0].imshow(Image.open(img))
    ax[0].set_title('H&E | native image',fontsize=9)
    draw(1,'E/I/S regions',d.region_main.map(colors).astype(str));ax[1].legend(handles=[Patch(color=v,label=k) for k,v in colors.items() if k in d.region_main.unique()],fontsize=6,loc='upper center',bbox_to_anchor=(.5,0),ncol=3,frameon=False)
    for i,key in [(2,'T_f'),(3,'I_f'),(4,'local_Moran_T')]:
        art=draw(i,{'T_f':'Total T fraction','I_f':'All immune fraction','local_Moran_T':'Local T Moran I'}[key],mm[key] if i==4 else d[key],cmap='coolwarm' if i==4 else 'viridis');fig.colorbar(art,ax=ax[i],fraction=.035)
    label=d.window500.map(st.local_state);ok=d.window500.map(st.neighbor_support).eq(True)&d.all_E
    pal={'infiltrated-like':'#297bb1','non-infiltrated-like':'#d07d44','infiltrated':'#297bb1','non_infiltrated':'#d07d44'}
    cs=label.map(pal).fillna('#bcbcbc').where(ok,'#d8d8d8')
    draw(5,'500 μm local pattern | supported E',cs);ax[5].legend(handles=[Patch(color='#297bb1',label='entry'),Patch(color='#d07d44',label='limited entry'),Patch(color='#d8d8d8',label='other / unresolved')],fontsize=6,frameon=False,loc='upper center',bbox_to_anchor=(.5,.01))
    for i,b in enumerate(B,6):
        draw(i,'All E border' if i==6 else 'ENN + ENS border',np.where(d[b],'#b8d2d9','#e0dfdb'))
        if len(edge):
            e=edge[edge.boundary==b];lines=e[['x1','y1','x2','y2']].to_numpy().reshape(-1,2,2);ax[i].add_collection(LineCollection(lines,colors=np.where(e.contact_type=='to_EIN_EIS','#9c397d','#2c4967'),linewidths=.52))
        ax[i].set_xlabel('Plum = EIN/EIS contact',fontsize=7)
    draw(8,'MP10 score | E domain','#dedede');v=d[d.new_epi_domain];art=ax[8].scatter(v.x_um,v.y_um,c=v.MP10,s=3,lw=0,cmap='magma');fig.colorbar(art,ax=ax[8],fraction=.035)
    draw(9,'MP10 coverage',np.where(d.new_epi_domain,np.where(d.legacy_query,'#518f8d','#d07c3f'),'#ddd'));ax[9].legend(handles=[Patch(color='#518f8d',label='E + old query'),Patch(color='#d07c3f',label='E now visible')],fontsize=6,frameon=False,loc='upper center',bbox_to_anchor=(.5,.01))
    for i,a in enumerate(ax):
        a.set_aspect('equal');a.set_xticks([]);a.set_yticks([])
        if i:a.invert_yaxis()
        for s in a.spines.values():s.set_visible(False)
    fig.suptitle(sec+'   |   '+cl.at[sec,'provisional_architecture_group']+'   |   '+r.source_context.split('_source')[0],fontsize=13)
    fig.tight_layout(rect=[0,.03,1,.95]);fig.savefig(out/'review'/f'{sec}.jpg',dpi=118);plt.close(fig)
    rows.append(dict(section_id=sec,image='review/'+sec+'.jpg',n_spots=len(d),E_scored=len(v),n_segments=len(edge),he_available=img.is_file(),group=cl.at[sec,'provisional_architecture_group']))
save(rows,out/'review/index.tsv')

