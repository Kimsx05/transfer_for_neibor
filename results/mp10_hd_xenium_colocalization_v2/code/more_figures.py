from common import *
from figures import plt,save,colors
from graph import polygons
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon
from scipy.spatial import cKDTree
def run(s):
 v,d=folders(s);out=OUT/'figures'/s;out.mkdir(exist_ok=True);z=labels(s);xy=z[['x_um','y_um']].values;idt=pd.read_parquet(d/'identity.parquet')
 mid=xy[np.argmin(((xy-np.median(xy,axis=0))**2).sum(1))];sel=np.flatnonzero((abs(xy-mid)<130).all(1));g=sparse.load_npz(v/'graph_50.npz')[sel][:,sel].tocoo();ix=g.row<g.col
 fig,ax=plt.subplots(figsize=(5.7,4.8));ax.add_collection(LineCollection(np.stack([xy[sel[g.row[ix]]],xy[sel[g.col[ix]]]],axis=1),lw=.15,color='#a4a6ae',alpha=.5))
 for key in z.major.iloc[sel].unique():
  a=sel[z.major.iloc[sel].eq(key).values];ax.scatter(*xy[a].T,s=7,c=colors.get(key,'#787878'),label=key,lw=0)
 ax.set_aspect('equal');ax.invert_yaxis();ax.set(xlabel='x (µm), right →',ylabel='y (µm), down');ax.legend(loc='upper left',bbox_to_anchor=(1,.98),frameon=False,fontsize=7);ax.plot([mid[0]-120,mid[0]-70],[mid[1]+120]*2,color='k',lw=2);ax.text(mid[0]-120,mid[1]+115,'50 µm');save(fig,out/'fixed_graph_ROI')
 if s in XE:
  po=polygons(s,z.index);bd=pd.read_parquet(v/'boundary_geometry.parquet').reindex(z.index);cx=xy[int(np.argmax(bd.distance_centroid_to_polygon_um.values))];ii=np.flatnonzero((abs(xy-cx)<60).all(1));fig,ax=plt.subplots(figsize=(4.4,4.2))
  for j in ii:
   ax.add_patch(Polygon(np.array(po[j].exterior.coords),fc='#d39196' if bd.review_boundary.iloc[j] else '#d8e3e7',ec='#657383',lw=.5));ax.plot(*xy[j],marker='.',color='#536172',ms=1)
  ax.set(xlim=(cx[0]-50,cx[0]+50),ylim=(cx[1]+50,cx[1]-50),xlabel='x (µm)',ylabel='y (µm)');ax.set_aspect('equal');ax.plot([cx[0]-44,cx[0]-24],[cx[1]+45]*2,color='black',lw=2);ax.text(cx[0]-44,cx[1]+43,'20 µm',fontsize=7);save(fig,out/'polygon_QC_flagged')
 n=OUT/'niche'/s;df=pd.read_csv(n/'niche_overlay.tsv',sep='\t');df=df[(df.method=='MENDER')&(df.endpoint=='density')]
 mat=df.pivot(index='niche',columns='state',values='high_fraction');fig,ax=plt.subplots(figsize=(4.6,3.5));im=ax.imshow(mat,cmap='YlOrRd',vmin=0);ax.set(yticks=range(len(mat)),yticklabels=mat.index,xticks=range(len(mat.columns)),xticklabels=mat.columns,ylabel='MENDER niche');plt.setp(ax.get_xticklabels(),rotation=25,ha='right');fig.colorbar(im,ax=ax,label='High / eligible Epi');save(fig,out/'MENDER_high_fraction')
 for metric in ['density','Treg/T']:
  t=pd.read_csv(n/'within_niche_effects.tsv',sep='\t');t=t[t.endpoint==metric].copy();fig,ax=plt.subplots(figsize=(5.2,max(3,len(t)*.25+1)));yy=np.arange(len(t));ax.hlines(yy,t.CI_low,t.CI_high,lw=1,color='#677f91');ax.scatter(t.effect,yy,s=14,color='#903454');ax.axvline(0,color='gray',ls='--',lw=.7);ax.set(yticks=yy,yticklabels=['N'+str(nn)+' · '+m for nn,m in zip(t.niche,t.state_definition)],xlabel='High − low, '+metric);save(fig,out/('MENDER_effect_'+metric.replace('/','_')))
if __name__=='__main__':
 for s in sys.argv[1:] or HD+XE:run(s)
