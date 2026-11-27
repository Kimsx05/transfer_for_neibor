from base_data import *
from scipy.spatial import cKDTree
from sklearn.decomposition import PCA
from sklearn.metrics import normalized_mutual_info_score
from scipy.stats import kruskal
import matplotlib.pyplot as plt
platform=sys.argv[1];samples=HD if platform=='hd' else XE;rng=np.random.default_rng(1729)
parts=[];recovery=[];crit=[];support=[];markers=[];rem=[]
for s in samples:
 d=paths(s);z=pd.read_parquet(d/'analysis_index.parquet');lat=np.load(d/'latent.npy');idx=np.flatnonzero(z.eligible_expression & np.isfinite(lat[:,0]));take=rng.choice(idx,min(4000,len(idx)),replace=False)
 q=z.iloc[take].copy();q['cell_id']=q.index
 for j in range(lat.shape[1]):q['latent_'+str(j)]=lat[take,j]
 parts.append(q)
 # The sorted per-sample graph uses original eligible-row order.
 e=z.index[z.niche_id.ne('unassigned_QC')]
 if not (d/'niche_graph_ids.tsv.gz').exists():pd.Series(e).to_csv(d/'niche_graph_ids.tsv.gz',header=False,index=False)
 g=sparse.load_npz(d/'niche_graph.npz');assert g.shape==(len(e),len(e)) and (g.diagonal()==0).all() and (g!=g.T).nnz==0
 for niche,unit in z.groupby('niche_id'):
  if niche=='unassigned_QC':continue
  support.append((s,niche,len(unit),len(unit)/len(e),bool(len(unit)>=50 and len(unit)>=.01*len(e))))
 for oldmajor,sub in z.groupby('old_major',dropna=False):
  for (new,level),n in sub.groupby(['major_identity','review_level']).size().items():recovery.append((s,oldmajor,new,level,n))
 if platform=='xenium':
  hi=z.old_major_support.ge(.8)&z.old_major_support.sub(z.old_major_runner).ge(.25)&~z.old_major.isin(['Unknown','Mixed'])
  for (o,n,r),num in z[hi&z.major_identity.ne(z.old_major)].groupby(['old_major','major_identity','conflict']).size().items():crit.append((s,o,n,r,num))
 else:
  rem.append(z.groupby(['source','old_major','other_major','major_identity','review_level']).size().rename('n').reset_index().assign(sample=s))
f=pd.concat(parts,ignore_index=True);lat=f.filter(regex='^latent_').values;pc=PCA(2,random_state=1729);proj=pc.fit_transform(lat)
f['PC1']=proj[:,0];f['PC2']=proj[:,1];f[['cell_id','sample','PC1','PC2','major_identity','niche_id']].to_parquet(O/platform/'latent_projection_sample.parquet')
tab=pd.DataFrame(support,columns=['sample','niche','n','fraction_eligible','robust_presence']);write(tab,platform+'_niche_recurrence.tsv')
write(pd.DataFrame(recovery,columns=['sample','old_major','new_major','new_level','n']),platform+'_major_transition.tsv')
if crit:write(pd.DataFrame(crit,columns=['sample','old_major','new_major','review_reason','n']),platform+'_old_confident_conflicts.tsv')
if rem:write(pd.concat(rem),platform+'_TACCO_comparison.tsv')
pal=json.loads((O/'shared'/('colors_'+platform+'.json')).read_text())
for col in ['sample','major_identity','n_counts','area_um2']:
 fig,ax=plt.subplots(figsize=(5.6,4.5))
 if col in ['n_counts','area_um2']:
  sc=ax.scatter(*proj.T,c=np.log1p(f[col]),s=1,cmap='viridis',alpha=.6,rasterized=True);fig.colorbar(sc,ax=ax,label='log(1 + '+col+')')
 else:
  for t,sub in f.groupby(col):ax.scatter(sub.PC1,sub.PC2,s=1,c=pal.get(t,None),label=t,rasterized=True); 
  ax.legend(markerscale=5,fontsize=6,frameon=False,loc='center left',bbox_to_anchor=(1,.5))
 ax.set_xlabel('scVI latent PC1');ax.set_ylabel('scVI latent PC2');ax.set_title(platform+' · '+col)
 for ext in ['png','pdf']:fig.savefig(O/'figures'/f'{platform}_latent_{col}.{ext}',bbox_inches='tight',dpi=155)
 plt.close(fig)
q=f[f.major_identity=='Epithelial'];write(q.groupby('sample').agg(eligible_balanced_draw=('sample','size'),PC1_q10=('PC1',lambda a:a.quantile(.1)),PC1_q90=('PC1',lambda a:a.quantile(.9)),PC2_q10=('PC2',lambda a:a.quantile(.1)),PC2_q90=('PC2',lambda a:a.quantile(.9))).reset_index(),platform+'_epithelial_dispersion.tsv')
m=f.niche_id.ne('unassigned_QC');qc=[]
for col in ['sample','major_identity','review_level','source','geometry_flag','n_counts','nFeature','area_um2']:
 v=f.loc[m,col]
 if col in ['n_counts','nFeature','area_um2']:v=pd.qcut(v,5,duplicates='drop')
 qc.append((col,normalized_mutual_info_score(f.loc[m,'niche_id'],v.astype(str))))
write(pd.DataFrame(qc,columns=['field','normalized_mutual_information_balanced_sample']),platform+'_niche_confounding.tsv')

