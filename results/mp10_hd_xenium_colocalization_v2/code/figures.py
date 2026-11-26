from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.collections import LineCollection
from scipy.ndimage import gaussian_filter1d
from matplotlib.lines import Line2D
from scipy.stats import binned_statistic
colors=json.loads((V/'shared/colors.json').read_text())
plt.rcParams.update({'font.family':'Arial','font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
def save(fig,p):
 fig.text(.02,.012,'Provisional · Codex-reviewed · human review pending',fontsize=7,color='#666666');fig.tight_layout(rect=(0,.044,1,1));fig.savefig(p.with_suffix('.png'),dpi=210);fig.savefig(p.with_suffix('.pdf'));plt.close(fig)
def mapfig(z,values,p,label,categorical=False):
 fig,ax=plt.subplots(figsize=(6,5));x=z.x_um;y=z.y_um;ax.scatter(x,y,s=.5,c='#dedede',rasterized=True,linewidth=0)
 if categorical:
  palette={'high':'#9d214a','low':'#77a5b0','supported':'#632448','possible':'#d48340','unresolved':'#8c80a3','unsupported':'#b2b2b2','outside':'#ededed'}
  cats=sorted(values.dropna().unique(),key=lambda c:(c=='high' or c=='supported',str(c)))
  for j,k in enumerate(cats):
   m=(values==k).fillna(False);ax.scatter(x[m],y[m],s=3 if k in ['high','supported'] else 1,c=[palette.get(str(k),plt.get_cmap('tab10')(j%10))],rasterized=True,linewidth=0,label=str(k))
  ax.legend(markerscale=4,bbox_to_anchor=(1,.95),loc='upper left',frameon=False,title=label)
 else:
  m=values.notna();im=ax.scatter(x[m],y[m],s=2,c=values[m],cmap='viridis',linewidth=0,rasterized=True);fig.colorbar(im,ax=ax,shrink=.65,label=label)
 ax.set_aspect('equal');ax.invert_yaxis();ax.set(xlabel='x (µm), right →',ylabel='y (µm), down');lo,hi=ax.get_xlim();bottom,top=ax.get_ylim();sx=lo+.05*(hi-lo);sy=bottom-.03*(bottom-top);ax.plot([sx,sx+200],[sy,sy],color='black',lw=2);ax.text(sx,sy-.02*(bottom-top),'200 µm',fontsize=7);save(fig,p)
def field(s):
 hd=s in HD;v,d=folders(s);z=labels(s);i=pd.read_parquet(d/'identity.parquet');out=OUT/'figures'/s;out.mkdir(exist_ok=True);methods=['MP10'] if hd else STATES
 mapfig(z,i.grade,out/'Treg_evidence','Treg evidence',True)
 for res in ([16,32] if hd else ['cell']):
  cs=pd.read_parquet(d/f'states_{res}.parquet');e=pd.read_parquet(d/f'centers_{res}.parquet');nn=pd.read_parquet(d/f'neighborhood_{res}_50.parquet').set_index('id').loc[e.index];base=cs if hd else z
  for method in methods:
   ylabel={'MP10':'MP10 score','MP10_covered':'Panel MP10 score','Chol_biosynthesis_8':'Cholesterol biosynthesis score','H05_transfer':'Transfer support'}[method]
   mapfig(base,cs[method],out/f'{res}_{method}_score',ylabel);h=e[method+'_high'].astype(bool);mapfig(base,pd.Series(np.where(h,'high','low'),index=e.index).reindex(base.index),out/f'{res}_{method}_group','Epi group',True)
   for end,y in [('density',nn.supported*1e6/e.area_50),('Treg_T',nn.supported/nn.reliable_T.replace(0,np.nan))]:
    valid=nn.valid&y.notna();x=e.loc[valid,method].values;yy=y[valid].values;fig,ax=plt.subplots(figsize=(4.7,3.6))
    if len(x):
     ax.scatter(x,yy,s=4,c='#384e65',alpha=.15,lw=0,rasterized=True);bins=np.unique(np.quantile(x,np.linspace(0,1,17)))
     if len(bins)>2:
      j=np.minimum(np.searchsorted(bins[1:],x,side='right'),len(bins)-2);trend=pd.DataFrame({'x':x,'y':yy,'j':j}).groupby('j').mean();ax.plot(trend.x,gaussian_filter1d(trend.y.to_numpy(),1),color='#ac3252',lw=2)
    ax.set(xlabel=ylabel,ylabel=('Treg component bins/mm²' if hd else 'Treg cells/mm²') if end=='density' else 'Treg / reliable T',ylim=(0,None));save(fig,out/f'{res}_{method}_{end}_continuous')
   env=pd.read_csv(d/'distance_envelope.tsv',sep='\t');ev=env[(env.res.astype(str)==str(res))&(env.state==method)&(env.definition=='supported')&(env.block==400)].copy()
   if len(ev):
    rr=ev.distance.str.replace('CDF_','').astype(float);fig,ax=plt.subplots(figsize=(4.7,3.6));ax.fill_between(rr,ev.global_lower,ev.global_upper,color='#c5d1db');ax.plot(rr,ev.observed,'o-',color='#82395a');ax.axhline(0,color='gray',lw=.6);ax.set(xlabel='Distance to Treg (µm)',ylabel='CDF difference, high − low');save(fig,out/f'{res}_{method}_distance')
  if not hd:
   mask=e[[m+'_high' for m in methods]].astype(bool).to_numpy();pat=[' + '.join(np.array(['MP10','Chol8','H05'])[row]) if row.any() else 'none' for row in mask]
   mapfig(base,pd.Series(pat,index=e.index).reindex(base.index),out/'three_state_high_overlap','High status',True)
 eff=pd.read_csv(d/'effects.tsv.gz',sep='\t')
 for end in ['density','Treg/T','Treg/all','any contact','restricted distance um']:
  sub=eff[(eff.endpoint==end)&(eff.block_um==400)&(eff.identity_definition=='supported')].copy()
  if sub.empty:continue
  labs=sub.state_definition+' · '+sub.unit_resolution.astype(str)+' · '+sub.radius_or_contact_um.astype(str)+' µm'
  fig,ax=plt.subplots(figsize=(5.5,max(2.5,len(sub)*.38+1.3)));ax.hlines(np.arange(len(sub)),sub.CI_low,sub.CI_high,color='#547384');ax.scatter(sub.effect,np.arange(len(sub)),s=18,color='#842d50');ax.axvline(0,ls='--',color='#999',lw=.7);ax.set(yticks=np.arange(len(sub)),yticklabels=labs,xlabel='High − low, '+end);save(fig,out/('effect_'+end.replace('/','_').replace(' ','_')))
 p=OUT/'niche'/s
 for f in list(p.glob('*.parquet'))+list(p.glob('BuildNicheAssay.tsv.gz')):
  t=pd.read_csv(f,sep='\t').set_index('id') if str(f).endswith('gz') else pd.read_parquet(f)
  if 'niche' in t:mapfig(z,t.niche.reindex(z.index).astype(str),out/('niche_'+f.name.split('.')[0]),'Niche',True)
 for name,f in [('major',v/'major_50_matrix.tsv.gz'),('subtype',d/'squidpy_50.parquet')]:
  if not f.exists():continue
  t=pd.read_parquet(f).reset_index() if name=='subtype' else pd.read_csv(f,sep='\t');m=t.pivot(index='from',columns='to',values='Z' if name=='subtype' else 'z_score');a=m.values.copy();np.fill_diagonal(a,np.nan);mx=max(1,np.nanquantile(abs(a),.97));fig,ax=plt.subplots(figsize=(max(5,len(m)*.27+2),max(4,len(m)*.255+1)));im=ax.imshow(a,cmap='RdBu_r',vmin=-mx,vmax=mx);ax.set(xticks=np.arange(len(m)),yticks=np.arange(len(m)),yticklabels=m.index,xticklabels=m.columns);ax.tick_params(labelsize=6);plt.setp(ax.get_xticklabels(),rotation=90);fig.colorbar(im,ax=ax,label='Neighborhood enrichment Z-score',shrink=.6);save(fig,out/('matrix_'+name))
if __name__=='__main__':
 for s in sys.argv[1:] or HD+XE:field(s)
