from common import *
from scipy.spatial import cKDTree
from scipy.ndimage import distance_transform_edt
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score,silhouette_score
from scipy.optimize import linear_sum_assignment
import joblib
from datetime import datetime

def match_min(a,b,K):
 tab=np.zeros((K,K),int);np.add.at(tab,(a,b),1);i,j=linear_sum_assignment(-tab);return float(min(tab[i,j]/np.maximum(1,tab.sum(1)[i]+tab.sum(0)[j]-tab[i,j])))
if __name__=='__main__':
 gate();sel=pd.read_csv(R/'tables/unit_selection.tsv',sep='\t');assert sel.loc[sel.selected,'unit'].iloc[0]=='cell','This implementation must be adapted if bin selected; never substitute counts for model weights silently'
 d=pd.read_parquet(R/'tables/Cervilla_cell_annotation.parquet');nodes=pd.read_parquet(A/'hd_epi_niche_dual_methods_v1/Cervilla/nodes.parquet');xy=d[['x_um','y_um']].to_numpy();nx=nodes[['x_um','y_um']].to_numpy();tree=cKDTree(nx);dist,nn=tree.query(xy,workers=1);d['component']=nodes.component.to_numpy()[nn];d['nearest_supplied_bin_um']=dist
 # Distance to observed tissue/supplied-grid boundary on the actual original array; no inferred expression mask.
 row=nodes.array_row.to_numpy(int);col=nodes.array_col.to_numpy(int);r0=row.min()-1;c0=col.min()-1;mask=np.zeros((row.max()-r0+2,col.max()-c0+2),bool);mask[row-r0,col-c0]=nodes.in_tissue.to_numpy()==1;edt=distance_transform_edt(mask)*8;d['distance_to_observed_boundary_um']=np.maximum(0,edt[row[nn]-r0,col[nn]-c0]-dist-4)
 types=sorted(set(d.type_supported)-{'Unknown'});typemap={t:i for i,t in enumerate(types)};ti=d.type_supported.map(typemap).fillna(-1).to_numpy(int);anchor=np.flatnonzero(d.epi_anchor);anc=d.iloc[anchor].copy().reset_index(drop=True);anc['global_index']=anchor
 raw_arrays={k:d[k+'_UMI'].to_numpy() for k in ['Epi','Endo','T']};tissue_array=nodes.in_tissue.to_numpy();modes=[];metrics=[];interpret=[];cfg={};Tlists=[]
 for radius in [50,200]:
  C=np.zeros((len(anc),len(types)),np.int32);nall=np.zeros(len(anc),int);cover=np.zeros(len(anc));observed=np.zeros(len(anc),int);tissue=np.zeros(len(anc),int)
  raw_epi=np.zeros(len(anc));raw_endo=np.zeros(len(anc));raw_t=np.zeros(len(anc));mp10neighbors=np.zeros(len(anc))
  for component in sorted(d.component.unique()):
   members=np.flatnonzero(d.component.to_numpy()==component);centers=np.flatnonzero(anc.component.to_numpy()==component);nd=np.flatnonzero(nodes.component.to_numpy()==component);kt=cKDTree(xy[members]);nt=cKDTree(nx[nd]);
   for start in range(0,len(centers),1000):
    z=centers[start:start+1000];near=kt.query_ball_point(anc.iloc[z][['x_um','y_um']].to_numpy(),radius+1e-7,workers=1);knownbins=nt.query_ball_point(anc.iloc[z][['x_um','y_um']].to_numpy(),radius+1e-7,workers=1)
    for ii,ns,nb in zip(z,near,knownbins):
     idx=members[ns];idx=idx[idx!=anchor[ii]];nall[ii]=len(idx);lab=ti[idx];C[ii]=np.bincount(lab[lab>=0],minlength=len(types));observed[ii]=len(nb);tissue[ii]=int(tissue_array[nd[nb]].sum());raw_epi[ii]=raw_arrays['Epi'][idx].sum();raw_endo[ii]=raw_arrays['Endo'][idx].sum();raw_t[ii]=raw_arrays['T'][idx].sum()
  nknown=C.sum(1);prop=np.divide(C,nknown[:,None],out=np.zeros_like(C,dtype=float),where=nknown[:,None]>0);cover=np.minimum(1,tissue*64/(np.pi*radius**2));valid=(nknown>=5)&(cover>=.25);out=anc[['cell_id','x_um','y_um','component','global_index','nCount_RNA','nFeature_RNA','epi_mixed_flag','mp10_epithelial_eligible','distance_to_observed_boundary_um']].copy();out['radius_um']=radius;out['n_neighbors']=nall;out['n_known']=nknown;out['known_fraction']=nknown/np.maximum(1,nall);out['supplied_positions']=observed;out['observed_tissue_positions']=tissue;out['observed_area_um2']=tissue*64;out['observed_area_fraction']=cover;out['neighborhood_eligible']=valid;out['neighborhood_reason']=np.where(nknown==0,'no_known_neighbors',np.where(nknown<5,'too_few_known_objects',np.where(cover<.25,'less_than_quarter_circle_observed','eligible')));out['raw_neighbor_Endo_UMI']=raw_endo;out['raw_neighbor_T_UMI']=raw_t;out['raw_neighbor_Epi_UMI']=raw_epi
  for k,t in enumerate(types):out['count_'+t]=C[:,k];out['P_'+t]=prop[:,k]
  # Exact nearest support distance within component, censored at observed boundary.
  for target in ['Endothelial','Treg']:
   near=np.full(len(out),np.nan)
   for comp in sorted(d.component.unique()):
    inds=np.flatnonzero((d.component==comp)&d.type_supported.eq(target));q=np.flatnonzero(out.component==comp)
    if len(inds):near[q]=cKDTree(xy[inds]).query(out.iloc[q][['x_um','y_um']],workers=1)[0]
   out['distance_'+target+'_um']=near;out['distance_'+target+'_censored']=np.isnan(near)|(near>out.distance_to_observed_boundary_um);out['distance_'+target+'_observed_um']=np.where(out['distance_'+target+'_censored'],np.nan,near)
  # Training centers balanced by pre-specified 400um blocks; all neighbors remain intact.
  rng=np.random.default_rng(123);block=(out.x_um//400).astype(int).astype(str)+'_'+(out.y_um//400).astype(int).astype(str);out['spatial_block']=block;T=[]
  for _,g in out[valid].groupby('spatial_block'):T.extend(rng.choice(g.index,min(150,len(g)),False))
  T=np.array(sorted(T),int);Z=np.sqrt(prop);Tlists.append(pd.DataFrame(dict(radius_um=radius,cell_id=out.iloc[T].cell_id)));hold=np.array([int(b.split('_')[0])*73856093+int(b.split('_')[1])*19349663 for b in block])%5==0;train_block=T[~hold[T]];small=rng.choice(T,min(2000,len(T)),False);models={};labels={}
  for K in [2,4,6,8]:
   km=KMeans(K,n_init=20,max_iter=300,random_state=123).fit(Z[T]);pred=km.predict(Z);models[K]=km;labels[K]=pred;ars=[];jacs=[]
   for seed in [456,789]:
    alt=KMeans(K,n_init=20,max_iter=300,random_state=seed).fit(Z[T]).predict(Z[T]);ars.append(adjusted_rand_score(pred[T],alt));jacs.append(match_min(pred[T],alt,K))
   bm=KMeans(K,n_init=20,max_iter=300,random_state=123).fit(Z[train_block]);bp=bm.predict(Z[T]);ba=adjusted_rand_score(pred[T],bp);bj=match_min(pred[T],bp,K);metric=dict(radius_um=radius,K=K,n_train=len(T),seed_ARI_min=min(ars),seed_min_cluster_Jaccard=min(jacs),block_ARI=ba,block_min_cluster_Jaccard=bj,silhouette=silhouette_score(Z[small],pred[small]),min_cluster_n=int(np.bincount(pred[valid],minlength=K).min()));metrics.append(metric);out[f'niche_K{K}_candidate']=np.where(valid,pred+1,0)
  eligible=[x for x in metrics if x['radius_um']==radius and x['K']>2 and x['seed_ARI_min']>=.8 and x['block_ARI']>=.6 and min(x['block_min_cluster_Jaccard'],x['seed_min_cluster_Jaccard'])>=.4];K=min([x['K'] for x in eligible],default=2);out['niche']=np.where(valid,labels[K]+1,0);out['niche_K2']=np.where(valid,labels[2]+1,0);out['niche_assignment']=np.where(valid,np.where(np.isin(np.arange(len(out)),T),'training','model_predict'),'unresolved');cfg[str(radius)]=dict(K=K,n_anchor=len(out),n_eligible=int(valid.sum()),n_training=len(T),types=types,quality_rule='>=5 known objects and >=25% observed circle area; no TME minimum',association_area_fraction=.5,feature='sqrt full known-type relative composition, including Epi',selected_at=datetime.now().isoformat(),selection_without_MP10=True)
  for k in sorted(out.niche.unique()):
   if k==0:continue
   ix=out.niche.eq(k);mean=prop[ix].mean(0);top=np.argsort(-mean)[:3];name=' / '.join(f'{types[j]} {mean[j]:.0%}'for j in top);interpret.append(dict(radius_um=radius,niche=k,n=int(ix.sum()),composition_description=name,mean_known_fraction=float(out.loc[ix,'known_fraction'].mean()),mean_observed_area_fraction=float(cover[ix].mean()),**{t:float(mean[j])for j,t in enumerate(types)}))
  out.to_parquet(R/'tables'/f'Cervilla_cell_neighborhood_r{radius}.parquet',index=False);joblib.dump(models[K],R/'models'/f'niche_r{radius}_K{K}.joblib');joblib.dump(models[2],R/'models'/f'niche_r{radius}_K2.joblib');print('NICHE',radius,cfg[str(radius)],flush=True)
 d.to_parquet(R/'tables/Cervilla_cell_formal_metadata.parquet',index=False);put(pd.DataFrame(metrics),'niche_candidate_metrics.tsv');put(pd.DataFrame(interpret),'niche_composition.tsv');put(pd.concat(Tlists),'niche_training_ids.tsv.gz');(R/'niche_frozen.json').write_text(json.dumps(cfg,indent=2));print('NICHE_FROZEN_NO_ASSOCIATION_INSPECTED',flush=True)
