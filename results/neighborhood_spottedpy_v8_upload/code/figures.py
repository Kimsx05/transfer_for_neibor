from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm
from PIL import Image
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
NICHE_COLORS=plt.cm.hsv(np.arange(9)/9*.94);HOT='#D55E00';COLD='#0072B2';BG='#d0d0d0'
def save(fig,path):
 path=OUT/'06_figures'/path;path.parent.mkdir(parents=True,exist_ok=True);fig.savefig(str(path)+'.png',dpi=300,transparent=True,bbox_inches='tight');fig.savefig(str(path)+'.pdf',transparent=True,bbox_inches='tight');plt.close(fig)
def figure_done(png):
 pdf=png.with_suffix('.pdf')
 if not png.exists() or not pdf.exists():return False
 with open(pdf,'rb') as f:
  f.seek(max(0,pdf.stat().st_size-1024));return b'%%EOF' in f.read()
def spatial(ax,xy,title):
 ax.set_aspect('equal');ax.set_xlim(xy[:,0].min()-1,xy[:,0].max()+1);ax.set_ylim(xy[:,1].max()+1,xy[:,1].min()-1);ax.set_xticks([]);ax.set_yticks([]);ax.set_title(title,fontsize=8)
def source_map(sec):
 t,a,g=load_section(sec);meta=json.loads((OUT/'01_geometry'/sec/'metadata.json').read_text());xy=a['xy_pitch'];base=Path(meta['dataset_id'])/'sections'/sec;dest=OUT/'06_figures'/base/'MP10_overview.png'
 if figure_done(dest):return
 fig,axs=plt.subplots(3,4,figsize=(15,11));fig.suptitle(sec+' | SOURCE_G, k=10 | true physical coordinates (pitch)\nHot/cold = local statistic of neighbors; center score need not be high',fontsize=12)
 for row,(var,layer) in enumerate(zip(VARIABLES[:2],LAYERS)):
  d=pd.read_parquet(OUT/'02_hotspots'/sec/f'{var[0]}_k10.parquet');ids=d.node_index.to_numpy();ax=axs[row,0];ax.scatter(*xy.T,s=2,c=BG,rasterized=False);ee=np.flatnonzero(t.E0);vals=t[layer].to_numpy(float)[ee];positive=vals>0
  if positive.any():im=ax.scatter(*xy[ee[positive]].T,s=4,c=vals[positive],cmap='plasma',vmin=0,vmax=vals.max(),rasterized=False);fig.colorbar(im,ax=ax,fraction=.035,pad=.01)
  ax.scatter(*xy[ee[~positive]].T,s=4,c='black',rasterized=False);spatial(ax,xy,var[0]+f' frozen | E0={len(ee)}; 0=black')
  for col,track in enumerate(TRACKS,1):
   ax=axs[row,col];ax.scatter(*xy.T,s=2,c='#e5e5e5',rasterized=False);ax.scatter(*xy[ids].T,s=3,c='#aaa',rasterized=False)
   for kind,color in [('cold',COLD),('hot',HOT)]:
    mask=d[f'{track}_{kind}'].to_numpy();ax.scatter(*xy[ids[mask]].T,s=5,c=color,label=f'{kind}: {mask.sum()}; C={d.loc[mask,f"{track}_{kind}_component"].nunique()}',rasterized=False)
   spatial(ax,xy,track+f' | E={int(d.Z.notna().sum())}; R1 size>=5');ax.legend(loc='upper right',fontsize=6,framealpha=.6,markerscale=2)
  ax=axs[row,3];ax.scatter(*xy.T,s=2,c=BG,rasterized=False)
  for j in range(9):sel=t.accepted_label.to_numpy()==j;ax.scatter(*xy[sel].T,s=4,c=[NICHE_COLORS[j]],rasterized=False)
  hh=ids[d.SOURCE_P05_hot.to_numpy()];ax.scatter(*xy[hh].T,s=9,facecolors='none',edgecolors='black',linewidths=.3,rasterized=False);spatial(ax,xy,'Frozen niche0–8 + SOURCE_P05 hot outlines')
  ax=axs[2,row];v=d.value.to_numpy();bins=np.linspace(v.min(),v.max(),31) if len(v) and v.max()>v.min() else 10
  if len(v):ax.hist(v,bins=bins,density=True,histtype='step',color='#666',label='E supported');
  for track,color in zip(TRACKS,[HOT,COLD]):
   h=d.loc[d[f'{track}_hot'],'value'];
   if len(h):ax.hist(h,bins=bins,density=True,histtype='step',color=color,label=f'{track} hot ({len(h)})')
  ax.set_title(var[0]+' center values (no extra score cutoff)');ax.set_xlabel('Frozen MP10');ax.legend(fontsize=6)
 ax=axs[2,2];inventory=pd.read_csv(V3/'05_user_review/image_availability.tsv',sep='\t');ir=inventory[inventory.section_id.eq(sec)].iloc[0]
 if isinstance(ir.HE,str) and Path(ir.HE).exists():
  img=Image.open(ir.HE);ax.imshow(img);ax.scatter(t.imagecol*ir.scale_used,t.imagerow*ir.scale_used,s=.5,c='#00A6D6',alpha=.3,rasterized=False);ax.set_aspect('equal');ax.set_title('Existing H&E + spot centers\nRegistration inherited; no pathology validation');ax.axis('off')
 else:ax.text(.5,.5,'H&E unavailable',ha='center');ax.axis('off')
 axs[2,3].axis('off');handles=[Line2D([],[],marker='o',linestyle='',color=NICHE_COLORS[j],label=f'niche{j}') for j in range(9)];axs[2,3].legend(handles=handles,ncol=3,loc='upper left',fontsize=8);axs[2,3].text(0,.5,'Orange = hot; blue = cold\nGrey = unavailable / tissue background\nSOURCE_P05: native P < .05\nLOCAL_BH05: two-sided local BH q < .05\nLocal cutoffs are not unit-test stars.\nDistances never cross tissue components.',va='top',fontsize=8);fig.tight_layout(rect=[0,0,1,.95]);save(fig,base/'MP10_overview')

def cell_maps(sec):
 sources=[pd.read_parquet(OUT/'02_hotspots'/sec/f'{var[0]}_k10.parquet') for var in VARIABLES[:2]]
 t,a,g=load_section(sec);meta=json.loads((OUT/'01_geometry'/sec/'metadata.json').read_text());xy=a['xy_pitch'];base=Path(meta['dataset_id'])/'sections'/sec
 for page in range(3):
  if figure_done(OUT/'06_figures'/base/f'cell_targets_{page+1}.png'):continue
  types=NON[page*10:(page+1)*10];fig,axs=plt.subplots(5,4,figsize=(14,16));axs=axs.ravel();fig.suptitle(sec+f' | full-tissue composition hotspots | k10 | page {page+1}/3\nSOURCE_P05 (left), LOCAL_BH05 (right); orange/blue: target hot/cold\nBlack circle outlines: MP10_SPOT hot; green square outlines: MP10_EPI hot (same track)',fontsize=11)
  for i,j in enumerate(types):
   d=pd.read_parquet(OUT/'02_hotspots'/sec/f'CELL_{j:02d}_k10.parquet');ids=d.node_index.to_numpy()
   for ti,track in enumerate(TRACKS):
    ax=axs[i*2+ti];ax.scatter(*xy.T,s=1,c='#ddd',rasterized=False)
    for direction,col in [('cold',COLD),('hot',HOT)]:mask=d[f'{track}_{direction}'].to_numpy();ax.scatter(*xy[ids[mask]].T,s=3,c=col,rasterized=False)
    for sd,col,marker in zip(sources,['black','#009E73'],['o','s']):
     sh=sd.loc[sd[f'{track}_hot'],'node_index'].to_numpy();ax.scatter(*xy[sh].T,s=7,facecolors='none',edgecolors=col,marker=marker,linewidths=.25,rasterized=False)
    h=d[d[f'{track}_hot']];txt=f'{TYPES[j]} | {track}\nH={len(h)}; target mean={h.value.mean():.3g}' if len(h) else f'{TYPES[j]} | {track}\nNO_TARGET_HOTSPOT';spatial(ax,xy,txt)
  for ax in axs[2*len(types):]:ax.axis('off')
  fig.tight_layout(rect=[0,0,1,.96]);save(fig,base/f'cell_targets_{page+1}')

def scope_figures(scope,tests,units):
 base=Path(scope)/'summaries';sub=tests[tests.scope.eq(scope)&tests.k.eq(10)&tests.analysis_domain.eq('MAIN')&tests.subset.eq('ALL_POLICY_UNITS')]
 for track in TRACKS:
  for li,layer in enumerate(LAYERS):
   layername=['MP10_SPOT','MP10_EPI'][li]
   for kind,targets in [('NICHE',[f'niche{i}' for i in range(9)]),('CELL',[TYPES[j] for j in NON])]:
    prefix=f'{kind}_{layername}_{track}';dd=sub[sub.track.eq(track)&sub.layer.eq(layer)&sub.target_kind.eq(kind)];uu=units[units.scope.eq(scope)&units.layer.eq(layer)&units.track.eq(track)&units.target_kind.eq(kind)&units.analysis_domain.eq('MAIN')]
    if not figure_done(OUT/'06_figures'/base/(prefix+'_heatmap.png')):
     fig,axs=plt.subplots(1,3,figsize=(10,4 if kind=='NICHE' else 10));fig.suptitle(scope+' | '+prefix+'\nUnit mean effect; relation-test q_scope stars; NA grey',fontsize=10)
     for mi,metric in enumerate(['O','D','D_out']):
      ax=axs[mi];z=dd[dd.metric.eq(metric)].set_index('target_id').reindex(targets);values=z.mean_effect.to_numpy()* (100 if metric=='O' else 1);mx=max(np.nanmax(abs(values)) if np.isfinite(values).any() else 1,1e-8);cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#bdbdbd');im=ax.imshow(np.ma.masked_invalid(values[:,None]),cmap=cmap,vmin=-mx,vmax=mx,aspect='auto');ax.set_yticks(np.arange(len(targets)),targets if mi==0 else []);ax.set_xticks([]);ax.set_title({'O':'Overlap excess (pp)\npositive = enriched','D':'Distance excess (pitch)\nnegative = nearer','D_out':'Outside distance (pitch)\nnegative = nearer'}[metric],fontsize=8)
      for i,(_,r) in enumerate(z.iterrows()):ax.text(0,i,f'{values[i]:.2g} {stars(r.q_scope)}\nN={int(r.n_units)}' if np.isfinite(values[i]) else f'NA\nN={int(r.n_units)}',ha='center',va='center',fontsize=6,color='white' if np.isfinite(values[i]) and abs(values[i])>.6*mx else 'black')
      fig.colorbar(im,ax=ax,fraction=.04,pad=.03)
     fig.tight_layout(rect=[0,0,1,.93]);save(fig,base/(prefix+'_heatmap'))
    if not figure_done(OUT/'06_figures'/base/(prefix+'_unit_boxes.png')):
     fig,axs=plt.subplots(1,3,figsize=(13,5 if kind=='NICHE' else 13));fig.suptitle(scope+' | '+prefix+' | main k10\nSame eligible unit rows as tests; 1.5 IQR boxes; circles verified, triangles other units; q_scope',fontsize=10)
     for mi,metric in enumerate(['O','D','D_out']):
      ax=axs[mi];fac=100 if metric=='O' else 1
      for i,target in enumerate(targets):
       rows=uu[uu.k.eq(10)&uu.metric.eq(metric)&uu.target_id.eq(target)];vals=rows.effect.to_numpy()*fac
       if len(vals):
        ax.boxplot([vals],positions=[i],orientation='horizontal',widths=.55,showfliers=False,whis=1.5,manage_ticks=False,boxprops={'color':'#666'},medianprops={'color':'black'})
        for known,marker in [(True,'o'),(False,'^')]:r=rows[rows.known_patient.eq(known)];ax.scatter(r.effect*fac,np.full(len(r),i),marker=marker,s=10,c=COLD if known else HOT,alpha=.65,zorder=3)
       stat=dd[dd.metric.eq(metric)&dd.target_id.eq(target)].iloc[0];ax.text(1.01,i,f'N={len(vals)} {stars(stat.q_scope)}',transform=ax.get_yaxis_transform(),fontsize=6,va='center')
      ax.axvline(0,color='black',ls=':',lw=.8);ax.set_yticks(range(len(targets)),targets if mi==0 else []);ax.invert_yaxis();ax.set_title(metric+(' (pp)' if metric=='O' else ' (pitch)'));ax.grid(axis='x',alpha=.2)
     fig.tight_layout(rect=[0,0,1,.94]);save(fig,base/(prefix+'_unit_boxes'))
    if not figure_done(OUT/'06_figures'/base/(prefix+'_k_sensitivity.png')):
     fig,axs=plt.subplots(1,3,figsize=(13,6 if kind=='NICHE' else 14));fig.suptitle(scope+' | '+prefix+' | fixed scales, no sensitivity P selection\nAdjacent unit boxes: k6 blue / k10 black / k18 orange; coverage differs',fontsize=10)
     for mi,metric in enumerate(['O','D','D_out']):
      ax=axs[mi];fac=100 if metric=='O' else 1
      for i,target in enumerate(targets):
       for offset,k,col in [(-.22,6,COLD),(0,10,'black'),(.22,18,HOT)]:
        rows=uu[uu.k.eq(k)&uu.metric.eq(metric)&uu.target_id.eq(target)];vals=rows.effect.to_numpy()*fac
        if len(vals):ax.boxplot([vals],positions=[i+offset],orientation='horizontal',widths=.18,showfliers=False,whis=1.5,manage_ticks=False,boxprops={'color':col},medianprops={'color':col},whiskerprops={'color':col},capprops={'color':col})
      ax.axvline(0,color='black',ls=':',lw=.8);ax.set_yticks(range(len(targets)),targets if mi==0 else []);ax.invert_yaxis();ax.set_title(metric+(' (pp)' if metric=='O' else ' (pitch)'))
     fig.tight_layout(rect=[0,0,1,.95]);save(fig,base/(prefix+'_k_sensitivity'))
 log('scope figures '+scope)

def spatial_worker(sec):
 source_map(sec);cell_maps(sec);log('spatial figures '+sec)
def scope_worker(scope):
 tests=pd.read_parquet(OUT/'07_synthesis/all_scope_tests.parquet');units=pd.read_parquet(OUT/'07_synthesis/all_unit_effects.parquet');scope_figures(scope,tests,units)
def full():
 assert done('M06_RELATIONS_STATISTICS_COMPLETE')
 from concurrent.futures import ProcessPoolExecutor
 with ProcessPoolExecutor(max_workers=2) as pool:list(pool.map(spatial_worker,COV.section_id))
 with ProcessPoolExecutor(max_workers=2) as pool:list(pool.map(scope_worker,SCOPES))
 checkpoint('M07_FIGURES_COMPLETE')
if __name__=='__main__':full()
