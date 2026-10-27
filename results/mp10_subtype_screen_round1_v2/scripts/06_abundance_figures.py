from pathlib import Path
import os,sys,re
P=Path('/home/data/t070721/codex_workspace/Bladder Metabolism');O=P/'mp10_subtype_screen_round1_v2';T=O/'tables'
sys.path=[q for q in sys.path if 'dist-packages' not in q]
sys.path[:0]=[str(P/'spatial_ecology_top20_v2/runtime'),str(P/'cell2location_setup_v1/cell2location_env/lib/python3.12/site-packages')]
os.environ['MPLCONFIGDIR']=str(O/'logs/mplconfig');os.environ['MPLBACKEND']='Agg'
import pandas as pd,numpy as np,matplotlib.pyplot as plt
from matplotlib.colors import to_hex
plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'svg.fonttype':'none','axes.titleweight':'normal','axes.spines.top':False,'axes.spines.right':False})
s=pd.read_csv(T/'samples_v2.tsv',sep='\t');p=pd.read_csv(T/'proportions_v2.tsv.gz',sep='\t');p=p.merge(s[['sample_uid','tissue_status','T_category','sample_MP10_group_allEpi','MP10_sample_allEpi']],on='sample_uid',validate='many_to_one')
labels=sorted(set(p.subtype)|{'Malignant_like_epithelial','Normal_epithelial','Other_epithelial_unresolved'})
existing=dict(re.findall(r'"([^"]+)"="(#[A-Fa-f0-9]+)"',(P/'subtype/v1/final_annotation_preview_v1/code/00_preview_config.R').read_text()))
pool=[to_hex(plt.get_cmap(c)(i)) for c in ['tab20','tab20b','tab20c'] for i in range(20)]
colors={n:existing[n] for n in labels if n in existing}
for n in labels:
 if n not in colors:
  colors[n]=next(c for c in pool if c not in colors.values())
colors.update({'Epi_MP10_low':'#4575b4','Epi_MP10_high':'#d73027','B':'#984ea3','iCAF':'#e6ab02','Matrix_CAF':'#d95f02','rCAF':'#a6761d','vCAF':'#e7298a','SPP1 TAM':'#e31a1c','C1QC TAM':'#6a3d9a','FOLR2 TAM':'#1b9e77','Classical Mono':'#ff7f00','T/NK':'#377eb8','Epithelial':'#d95f02','Malignant_like_epithelial':'#d73027','Normal_epithelial':'#4575b4','Other_epithelial_unresolved':'#aaaaaa'})
pd.DataFrame({'cell_type':list(colors),'color':list(colors.values())}).to_csv(T/'cell_colors.tsv',sep='\t',index=False)
def safe(n):return re.sub(r'[^a-zA-Z0-9_-]+','_',n)
def dest(ds):return O/'figures/overall' if ds=='COMMON_DATASETS' or ds=='ALL' else O/'figures/by_dataset'/ds
def save(fig,ds,name):
 d=dest(ds);d.mkdir(parents=True,exist_ok=True)
 for ext in ['pdf','png']:fig.savefig(d/(name+'.'+ext),dpi=300,transparent=True,bbox_inches='tight')
 plt.close(fig)
comps={'Normal_Tumor':('tissue_status',['Normal','Tumor']),'T1_T2':('T_category',['T1','T2']),'sample_MP10_group_allEpi':('sample_MP10_group_allEpi',['low','high']),'T_overview':('T_category',['Ta','Tis','T1','T2','T3','T4'])}
families=[('major','all_cells','Atlas'),('major','TME','TME')]+[('subtype','parent',pa) for pa in sorted(p.loc[p.level.eq('subtype'),'major'].unique())]
rows=[];boundaries=[];skips=[]
for comp,(col,grps) in comps.items():
 for level,den,parent in families:
  z=p[p.level.eq(level)&p.denominator_type.eq(den)&p[col].isin(grps)&p.proportion.notna()].copy()
  if level=='subtype':z=z[z.major.eq(parent)]
  if comp in ['T1_T2','T_overview']:z=z[z.tissue_status.eq('Tumor')]
  if z.empty:continue
  order=sorted(z.subtype.unique());z['group']=z[col]
  b=z.groupby(['dataset','group','subtype'],observed=True).agg(mean_proportion=('proportion','mean'),n_specimens=('sample_uid','nunique'),n_patients_recorded=('patient_uid','nunique')).reset_index()
  b['n_datasets']=1;b['dataset_set']=b.dataset;b['summary_method']='equal_specimen_mean'
  # Common-dataset aggregation is computed separately per denominator/parent.
  common=[ds for ds,a in z.groupby('dataset') if set(grps).issubset(set(a.group))]
  if common:
   c=b[b.dataset.isin(common)].groupby(['group','subtype']).mean_proportion.mean().reset_index();c['dataset']='COMMON_DATASETS';c['n_datasets']=len(common);c['dataset_set']=';'.join(common);c['summary_method']='common_datasets_equal_weight'
   for g in grps:
    zz=z[z.dataset.isin(common)&z.group.eq(g)];c.loc[c.group.eq(g),'n_specimens']=zz.sample_uid.nunique();c.loc[c.group.eq(g),'n_patients_recorded']=zz.patient_uid.nunique()
   b=pd.concat([b,c],ignore_index=True)
  else:skips.append({'comparison':comp,'parent':parent,'reason':'no_common_datasets_for_all_requested_groups'})
  b['comparison']=comp;b['parent']=parent;b['level']=level;b['denominator_type']=den
  rows.append(b)
  for ds,a in b.groupby('dataset'):
   groups=[g for g in grps if g in a.group.unique()]
   mat=a.pivot(index='subtype',columns='group',values='mean_proportion').reindex(index=order,columns=groups).fillna(0)
   assert np.allclose(mat.sum(axis=0),1,atol=1e-8)
   fig,ax=plt.subplots(figsize=(8.3,max(4.8,len(order)*.24)))
   for j,g in enumerate(groups):
    bottom=0
    for typ,v in mat[g].items():
     ax.bar(j,v,bottom=bottom,width=.48,color=colors[typ],edgecolor='white',linewidth=.35,label=typ if j==0 else None);bottom+=v
   for j in range(len(groups)-1):
    left=mat[groups[j]].cumsum();right=mat[groups[j+1]].cumsum()
    for typ in order[:-1]:
     ax.plot([j+.24,j+1-.24],[left[typ],right[typ]],color='black',lw=.55)
     boundaries.append(dict(comparison=comp,dataset=ds,parent=parent,subtype_boundary=typ,left_group=groups[j],right_group=groups[j+1],left_cumulative=left[typ],right_cumulative=right[typ]))
   tick=[]
   for g in groups:
    q=a[a.group.eq(g)].iloc[0];tick.append(f'{g}\n{int(q.n_specimens)} specimens / {int(q.n_patients_recorded)} mapped IDs')
   title='All-cell major composition' if parent=='Atlas' else 'TME major composition' if parent=='TME' else 'Composition within '+parent
   ax.set(title=title+'\n'+('Common datasets, equal weights' if ds=='COMMON_DATASETS' else ds),ylabel='Ratio',ylim=(0,1),xticks=np.arange(len(groups)),xticklabels=tick,yticks=[0,.25,.5,.75,1],xlim=(-.45,max(len(groups)-1,1)+.45))
   if len(groups)==1:ax.text(.5,.96,'No within-dataset two-group contrast',transform=ax.transAxes,ha='center',fontsize=9)
   ax.legend(loc='center left',bbox_to_anchor=(1.02,.5),frameon=False,fontsize=9)
   fig.tight_layout();save(fig,ds,'stacked_'+comp+'_'+safe(parent))
pd.concat(rows).to_csv(T/'stacked_bar_summary.tsv',sep='\t',index=False)
pd.DataFrame(boundaries).to_csv(T/'stacked_bar_boundaries.tsv',sep='\t',index=False)
pd.DataFrame(skips).to_csv(T/'stacked_bar_unavailable.tsv',sep='\t',index=False)
# Tables contain every P/q; matrices show effect and q and explicitly identify unestimable entries.
r=pd.read_csv(T/'proportion_models_v2.tsv',sep='\t');r=r[r.analysis_variant.eq('primary_asin')]
cs=['Normal_Tumor','T1_T2','sample_MP10_group_allEpi','MP10_continuous','MP10_stage_complete_unadjusted','MP10_stage_adjusted']
for ds in ['ALL']+sorted(s.dataset.unique()):
 for lev,den in [('major','TME'),('subtype','parent')]:
  a=r[r.scope.eq(ds)&r.annotation_level.eq(lev)&r.denominator_type.eq(den)&r.comparison.isin(cs)]
  order=a[['major','subtype']].drop_duplicates().sort_values(['major','subtype']).subtype.tolist()
  eff=a.pivot(index='subtype',columns='comparison',values='coefficient').reindex(index=order,columns=cs);q=a.pivot(index='subtype',columns='comparison',values='q').reindex(index=order,columns=cs)
  vals=np.sign(eff.to_numpy())*-np.log10(np.maximum(q.to_numpy(),1e-12));fig,ax=plt.subplots(figsize=(14,max(4,len(order)*.43)))
  im=ax.imshow(np.ma.masked_invalid(vals),cmap='RdBu_r',vmin=-3,vmax=3,aspect='auto');ax.set_xticks(range(len(cs)),['Tumor - Normal','T2 - T1','high - low','MP10 continuous','Stage-complete\nunadjusted','Stage-complete\nadjusted'],rotation=25,ha='right');ax.set_yticks(range(len(order)),order)
  for i in range(len(order)):
   for j in range(len(cs)):
    v=eff.iloc[i,j];qq=q.iloc[i,j];lab='NE' if not np.isfinite(v) else f'{v:.2g}\nq={qq:.2g}'
    ax.text(j,i,lab,ha='center',va='center',fontsize=7,color='white' if np.isfinite(vals[i,j]) and abs(vals[i,j])>2 else 'black')
  ax.set_title(f'{ds}: {lev} / {den}\nCoefficients: asin(sqrt(raw proportion)); continuous effects per 1 original MP10 score unit',pad=14)
  fig.colorbar(im,ax=ax,shrink=.4,label='Direction × -log10(q), clipped at ±3');fig.tight_layout();save(fig,ds,'effects_'+lev+'_'+den)
print('Abundance figures completed',len(rows),'summary families')
