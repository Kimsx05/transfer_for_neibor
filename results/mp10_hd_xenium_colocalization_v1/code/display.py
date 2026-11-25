from base import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
matplotlib.rcParams.update({'font.family':'Arial','pdf.fonttype':42,'font.size':10,'axes.spines.right':False,'axes.spines.top':False})
colors=json.loads((O/'shared/colors.json').read_text())
def save(fig,s,k):
    d=O/'figures'/s;d.mkdir(exist_ok=True);fig.savefig(d/(k+'.png'),dpi=300,bbox_inches='tight',facecolor='white');fig.savefig(d/(k+'.pdf'),bbox_inches='tight',facecolor='white');plt.close(fig)
def spatial(s,d,k,values,unit=None,limits=None):
    fig,ax=plt.subplots(figsize=(6.7,5.4));x=d.x_um.to_numpy();y=d.y_um.to_numpy();sz=1.35 if s in HD else 2.3
    if unit is None:
        v=pd.Series(np.asarray(values)).fillna('Unknown');cats=sorted(v.unique());handles=[]
        for c in sorted(cats,key=lambda c:(c not in ['Unknown','Mixed'],c)):
            j=v.eq(c);color=colors.get(c,'#888888')
            ax.scatter(x[j],y[j],s=sz,c=color,linewidths=0)
            handles.append(Line2D([],[],marker='o',color=color,linestyle='',markersize=4,label=c))
        ax.legend(handles=handles,loc='center left',bbox_to_anchor=(1,.5),frameon=False,fontsize=7)
    else:
        v=np.asarray(values,float);keep=np.isfinite(v);ax.scatter(x[~keep],y[~keep],s=sz,color='#dddddd',linewidths=0)
        z=ax.scatter(x[keep],y[keep],s=sz,c=v[keep],cmap='viridis',linewidths=0,vmin=None if limits is None else limits[0],vmax=None if limits is None else limits[1])
        fig.colorbar(z,ax=ax,shrink=.55,label=unit)
    ax.set_aspect('equal');ax.invert_yaxis();ax.set_xlabel('x (µm)');ax.set_ylabel('y (µm)');save(fig,s,k)
def scatter(s,a,b,x,y,k):
    v=pd.DataFrame({'x':a,'y':b}).dropna()
    fig,ax=plt.subplots(figsize=(4.5,3.8));ax.scatter(v.x,v.y,s=2,alpha=.15,color='#4477aa',linewidths=0)
    if v.x.nunique()>2:
        from statsmodels.nonparametric.smoothers_lowess import lowess
        z=lowess(v.y,v.x,frac=.28,it=0,delta=(v.x.max()-v.x.min())/100)
        ax.plot(z[:,0],z[:,1],color='#9a352a',lw=1.6)
    ax.set(xlabel=x,ylabel=y);save(fig,s,k)

