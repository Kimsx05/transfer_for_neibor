from base6 import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
plt.rcParams.update({'font.size':12,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','lines.linewidth':1.2})
def save(fig,name):
 fig.savefig(R/'figures'/(name+'.png'),dpi=300,bbox_inches='tight');fig.savefig(R/'figures'/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
def run(s):
 o=R/s;d=pd.read_parquet(o/'cells_scored.parquet');q=pd.read_parquet(o/'sources.parquet');e=d[d.score_ok];geo=np.load(o/'D.npz');reg=geo['roi'];ori=geo['origin'];step=8
 fig,ax=plt.subplots(1,2,figsize=(12,6),layout='constrained');ax[0].imshow(np.ma.masked_equal(reg,0),extent=(ori[0]-4,ori[0]+step*(reg.shape[1]-.5),ori[1]+step*(reg.shape[0]-.5),ori[1]-4),cmap='Pastel1',interpolation='nearest',alpha=.8)
 ax[0].scatter(q.x,q.y,s=.05,c='#38455c',alpha=.25,rasterized=True);ax[0].set_title('Screening sources and 200 µm union')
 ax[1].scatter(d.x[d.in_D],d.y[d.in_D],s=.05,c='#ddd',rasterized=True);ax[1].scatter(e.x,e.y,s=.15,c='#f0b080',rasterized=True)
 for k,col,marker in [('control','#2979af','.'),('A','#881233','o'),('B','#008174','s')]:
  z=d[d[k]];ax[1].scatter(z.x,z.y,s=8 if k=='control' else 23,c=col,marker=marker,label=f'{k} ({len(z)})',rasterized=False)
 ax[1].legend(fontsize=9);ax[1].set_title('Scorable epithelium and T centres')
 for a in ax:a.set_aspect('equal');a.set_xlabel('x (µm)');a.set_ylabel('y (µm)');a.set_ylim(sorted(a.get_ylim(),reverse=True))
 fig.suptitle(s);save(fig,s+'_coverage')
 # One display per predefined domain, anchored on a reviewed ID when possible, including negative domains.
 nb=pd.read_parquet(o/'neighborhoods.parquet');rm=rd(o/'roi_manifest.tsv');selection=[]
 for roi,focus in [(r,None) for r in rm.roi]+([(int(d.loc[d.id.eq(i),'roi'].iloc[0]),i) for i in (3391,62224,32362,138645)] if s=='D1_Johnson' else []):
  pool=d[(d.roi==roi)&d.review_action.ne('pending')].copy()
  if len(pool):
   # Prefer pre-existing traced ID, then reviewed main, then lowest ID. MP10 is unavailable to this choice.
   pool['old_priority']=pool.v5_old_id.isin([106665 if s=='Cervilla' else 49654])
   order=pool.sort_values(['old_priority','main','id'],ascending=[False,False,True]);z=order.iloc[0]
  else:
   z=d[d.roi.eq(roi)].sort_values('id').iloc[0] if d.roi.eq(roi).any() else q[q.roi.eq(roi)].iloc[0]
  if focus is not None:z=d.loc[d.id.eq(focus)].iloc[0]
  selection.append(dict(roi=roi,id=z.get('id'),x=z.x,y=z.y));cx=z.x;cy=z.y
  box=d[(d.x-cx).abs().le(225)&(d.y-cy).abs().le(225)];fig,ax=plt.subplots(1,2,figsize=(9.2,4.7),layout='constrained')
  color=plt.get_cmap('viridis');upper=np.quantile(e.MP10_cell_score,.99)
  for a in ax:
   a.scatter(box.x,box.y,s=3,color='#e2e2e2',label='Other nuclei',rasterized=True)
   for r in [50,200]:a.add_patch(Circle((cx,cy),r,fill=False,ls='--' if r==50 else '-',lw=1,color='#41464e'))
   a.set(xlim=(cx-225,cx+225),ylim=(cy+225,cy-225),aspect='equal',xlabel='x (µm)',ylabel='y (µm)')
  v=box[box.score_ok];im=ax[0].scatter(v.x,v.y,c=v.MP10_cell_score,vmin=0,vmax=max(upper,.001),cmap=color,s=12,rasterized=True);fig.colorbar(im,ax=ax[0],shrink=.72,label='MP10 (clipped at Q99)')
  for k,c,marker in [('H05','#cb7114','^'),('A','#9b103d','o'),('B','#168075','s'),('control','#2363ad','.')]:
   g=box[box[k].fillna(False)];ax[1].scatter(g.x,g.y,color=c,s=21 if k!='control' else 11,marker=marker,label=f'{k} ({len(g)})')
  ax[1].plot(cx,cy,'k+',ms=9);ax[1].legend(fontsize=8,loc='upper right');ax[0].set_title('Scorable epithelial cells');ax[1].set_title('H05 and T centres')
  fig.suptitle(f'{s}   domain {roi}   50 / 200 µm');save(fig,f'{s}_ROI{roi}'+(f'_cell{focus}' if focus else '')+'_MP10')
 put(pd.DataFrame(selection),o/'display_windows.tsv')
def effects():
 df=pd.concat([rd(R/s/'effects.tsv') for s in sections]);d=df[df.controls.eq('other_T')&df.metric.eq('H05_fraction')];fig,axs=plt.subplots(1,2,figsize=(10,4.6),sharey=True,layout='constrained')
 for ax,s in zip(axs,sections):
  a=d[d.section.eq(s)];ax.axvline(0,color='#bbb');ax.set_title(s)
  for j,g in enumerate(['A','B','main','strict','extended']):
   for j0,(r,col,m) in enumerate([(50,'#117733','o'),(200,'#882255','s')]):
    z=a[(a.group==g)&(a.radius==r)].iloc[0];y=j+(j0-.5)*.24
    if pd.notna(z.difference):ax.scatter(z.difference*100,y,c=col,marker=m);ax.annotate(str(z.matched),(z.difference*100,y),xytext=(5,3),textcoords='offset points',fontsize=8)
    else:ax.text(.04,y,'NA',transform=ax.get_yaxis_transform(),color=col,fontsize=8,va='center')
  ax.set_yticks(range(5),['A','B','A ∪ B','Strict','Extended']);ax.set_xlabel('H05 fraction difference (percentage points)');ax.legend(handles=[Line2D([],[],color=c,marker=m,ls='',label=f'{r} µm') for r,c,m in [(50,'#117733','o'),(200,'#882255','s')]],fontsize=9)
 fig.suptitle('Matched cases minus non-Treg T controls (numbers = matched cases)');save(fig,'matched_H05')
if __name__=='__main__':
 for s in sections:run(s)
 effects()
