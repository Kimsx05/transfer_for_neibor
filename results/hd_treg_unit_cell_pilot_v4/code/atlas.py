from base import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib import font_manager
from source import kinds
import subprocess
font_manager.fontManager.addfont(subprocess.check_output(['fc-match','-f','%{file}','Arial'],text=True))
plt.rcParams.update({'font.family':'Arial','font.size':12,'pdf.fonttype':42,'figure.facecolor':'white'})
def crop(d,q):return d[d.x_um.between(q.x0,q.x1)&d.y_um.between(q.y0,q.y1)]
def points(ax,d,color,size=5,marker='o',**kw):
 if len(d):ax.scatter(d.x_um,d.y_um,c=color,s=size,marker=marker,lw=0,rasterized=True,**kw)
def atlas():
 with PdfPages(R/'figures/ROI_atlas.pdf')as pages:
  for s in sections:
   root=R/s;src=pd.read_parquet(root/'sources.parquet');src['T']=src[T].sum(axis=1)
   tr=js(N/s/'metadata/transform.json');scale=np.array(tr['local_to_model_scale_xy']);off=np.array(tr['fullres_crop_xyxy'][:2]);mpp=tr['native_mpp']
   img=np.load(N/s/'inputs/model_image.npy',mmap_mode='r');mask=sparse.load_npz(S/s/'inputs/nuclei.npz').tocsr()
   objs={k:pd.read_parquet(root/(k+'_evidence.parquet')) for k in kinds}
   historic=pd.read_parquet(Z/s/'T_regions.parquet');ep=pd.read_parquet(V/s/'Epi64_regions.parquet',columns=['region_id','measurement_eligible','MP10','H05']);em=pd.read_parquet(V/s/'Epi64_members.parquet')
   native=pd.read_parquet(V/s/'bin_evidence.parquet',columns=['bin_id','x_um','y_um']);native=native.merge(em,on='bin_id',how='left').merge(ep,on='region_id',how='left')
   stats=rd(R/'cell_assignment_audit.tsv').query('section==@s')
   for q in rd(R/'roi_manifest.tsv').query('section==@s').itertuples():
    fig,axs=plt.subplots(3,4,figsize=(13.6,10.8));axs=axs.ravel();srcq=crop(src,q);nat=crop(native,q)
    l,b=np.floor(((np.array([q.x0,q.y0])/mpp-off)+.5)*scale).astype(int)
    h,t=np.ceil(((np.array([q.x1,q.y1])/mpp-off)+.5)*scale).astype(int)
    ext=[(l/scale[0]-.5+off[0])*mpp,(h/scale[0]-.5+off[0])*mpp,(t/scale[1]-.5+off[1])*mpp,(b/scale[1]-.5+off[1])*mpp]
    rgb=img[b:t,l:h];lab=mask[b:t,l:h].toarray();edge=(lab!=np.roll(lab,1,axis=0))|(lab!=np.roll(lab,1,axis=1));edge[0,:]=False;edge[:,0]=False
    axs[0].imshow(rgb,extent=ext);over=np.zeros((*lab.shape,4));over[edge]=[0,.45,.4,.8];axs[0].imshow(over,extent=ext)
    for i,(gena,genb)in enumerate([('T','FOXP3'),('IL2RA','CTLA4')],1):
     ax=axs[i];points(ax,srcq,'#eeeeee',1,'s');points(ax,srcq[srcq[gena]>0],'#2676ac',9);points(ax,srcq[srcq[genb]>0],'#c33459',12,'D')
     if i==1:
      outline=over.copy();outline[edge]=[.17,.28,.28,.44];ax.imshow(outline,extent=ext,zorder=1.1)
     points(ax,srcq[(srcq[gena]>0)&(srcq[genb]>0)],'#59328c',13)
    d=crop(objs['T8_from2'],q);points(axs[3],d,'#e6e6e6',7,'s');points(axs[3],d[d.T_UMI>0],'#2676ac',13,'s');points(axs[3],d[d.FOXP3>0],'#c33459',14,'D');points(axs[3],d[(d.T_UMI>0)&(d.FOXP3>0)],'#59328c',15)
    for i,k in enumerate(kinds,4):
     d=crop(objs[k],q);ax=axs[i];points(ax,d,'#d7dddd',5);points(ax,d[d.nCount==0],'#555555',7,'s');points(ax,d[d.FOXP3_raw_screen],'#dc911a',20,'D');points(ax,d[d.model_candidate],'#8e51a2',25);points(ax,d[d.supported_tier],'#1a6852',30,'^')
     f=d[d.FOXP3_raw_screen&d.model_candidate];ax.scatter(f.x_um,f.y_um,facecolors='none',edgecolors='#dc911a',s=50,linewidths=1)
    d=crop(historic,q);points(axs[9],d,'#e0e0e0',11,'s');points(axs[9],d[d.T_eligible],'#88baca',22,'s');points(axs[9],d[d.FOXP3_supported],'#dc911a',30,'D');points(axs[9],d[d.strict],'#8e51a2',35,'^')
    d=crop(stats,q);points(axs[10],d,'#dedede',6);old=d.fine_prediction.eq('Treg')if s=='Cervilla' else d.singler_blca_labels.eq('T/NK');points(axs[10],d[old],'#247daa',20)
    points(axs[11],nat,'#e2e2e2',9,'s');d=nat[nat.measurement_eligible.eq(True)];sc=axs[11].scatter(d.x_um,d.y_um,c=d.MP10,vmin=0,vmax=ep.MP10.max(),cmap='viridis',s=10,marker='s',linewidths=0,rasterized=True)
    for ax,title in zip(axs,['H&E and nuclear outlines','2 µm · T / FOXP3','2 µm · IL2RA / CTLA4','8 µm · T / FOXP3','8 µm · classifier','16 µm · classifier','Nucleus · classifier','4 µm extension · classifier','8 µm extension · classifier','Historical T32 screening','Old fine Treg' if s=='Cervilla'else 'Old BLCA T/NK','Frozen Epi64 MP10']):
     ax.set(xlim=(q.x0,q.x1),ylim=(q.y1,q.y0),aspect='equal');ax.set_xticks([]);ax.set_yticks([]);ax.set_title(title,fontsize=12)
    axs[0].plot([q.x0+10,q.x0+60],[q.y1-12]*2,c='black',lw=2);axs[0].text(q.x0+12,q.y1-17,'50 µm',fontsize=10)
    fig.subplots_adjust(.025,.12,.985,.97,.06,.19)
    handles=[Line2D([],[],marker=m,color=c,lw=0,label=t)for m,c,t in [('D','#dc911a','FOXP3 screen'),('o','#8e51a2','Model Treg candidate'),('^','#1a6852','Reference-supported tier'),('o','#2676ac','First raw signal'),('D','#c33459','Second raw signal'),('o','#59328c','Both in same grid')]]
    legend=fig.legend(handles=handles,loc='lower center',ncol=3,frameon=False,bbox_to_anchor=(.4,.005))
    cax=fig.add_axes([.77,.045,.18,.013]);fig.colorbar(sc,cax=cax,orientation='horizontal',label='Frozen Epi64 MP10')
    path=R/'figures'/f'{s}_{q.roi}';fig.savefig(str(path)+'.png',dpi=300);pages.savefig(fig);fig.savefig(str(path)+'.pdf');legend.remove();cax.remove();fig.savefig(str(path)+'_panel_only.pdf');plt.close(fig)
    lf=plt.figure(figsize=(9,1.9));lf.legend(handles=handles,loc='upper center',ncol=3,frameon=False);cx=lf.add_axes([.33,.23,.34,.095]);lf.colorbar(matplotlib.cm.ScalarMappable(norm=sc.norm,cmap=sc.cmap),cax=cx,orientation='horizontal',label='Frozen Epi64 MP10');lf.savefig(str(path)+'_legend_only.pdf');plt.close(lf)
def confusion():
 d=rd(R/'confusion.tsv');print(d.columns)
 for depth in [100]:
  z=d[d.depth.astype(str)==str(depth)].pivot(index='label',columns='predicted',values='n').fillna(0);idx=z.index;z=z.reindex(columns=idx,fill_value=0);fig,ax=plt.subplots(figsize=(10,9));im=ax.imshow(z.div(z.sum(axis=1),axis=0),cmap='viridis',vmin=0,vmax=1);ax.set_xticks(range(len(z)),idx,rotation=90,fontsize=8);ax.set_yticks(range(len(z)),idx,fontsize=8);ax.set(xlabel='Predicted class',ylabel='Held-out reference class');fig.colorbar(im,ax=ax,label='Fraction within reference class');fig.tight_layout();fig.savefig(R/'figures/reference_confusion.png',dpi=300);fig.savefig(R/'figures/reference_confusion.pdf');plt.close(fig)
if __name__=='__main__':atlas();confusion()
