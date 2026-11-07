from common import *
from scipy.cluster.hierarchy import linkage,cut_tree
from scipy.spatial.distance import pdist,squareform,cdist
from scipy.cluster.vq import kmeans2

REP=['niche_all_Epi_all','B1_R0_partial_r_full','B3_Epi_distance_partial_r_full','C1_CAF_full','cell_composition_whole_tissue']
def matrix_inputs():
 reg=read('endpoints.tsv.gz');u=read('unit_endpoints.tsv.gz');u=u[u.identity_level.eq('analysis_unit')&u.value.notna()].copy();mapping=read('canonical_programs.tsv');meta=read('unit_metadata.tsv').set_index('record_id');out={}
 def make(name,r,kind,filter_n=False):
  v=u[u.endpoint_id.isin(r.endpoint_id)].copy()
  if filter_n:v=v[v.n_min.ge(10)&v.n_other_min.ge(10)]
  v=v.merge(r[['endpoint_id','canonical_program','candidate_type']],on='endpoint_id',validate='m:1');v['feature']=v.canonical_program+' / '+v.candidate_type;mat=v.pivot(index='record_id',columns='feature',values='value');out[name]=(mat,kind)
 def version(r,ver):
  desired=set(mapping.full_program if ver=='full' else mapping.noMP10_program);return r[r.actual_program.isin(desired)]
 for mo in ['B1_R0','B1_R1','B1_R2','B3_Epi_distance','B3_tissue_distance']:
  for st in ['r','partial_r']:
   for ver in ['full','noMP10']:make(mo+'_'+st+'_'+ver,version(reg[reg.module.eq(mo)&reg.statistic.eq(st)],ver),'r')
 for mo in ['C1','C2']:
  r=reg[reg.module.eq(mo)&reg.score_method.eq('AddModuleScore_SCT_data')&reg.range_type.eq('common_support_0_1000')&reg.statistic.eq('r' if mo=='C1' else 'paired_delta_r')]
  for block in sorted(r.candidate_block.unique()):
   for ver in ['full','noMP10']:
    make(mo+'_'+block+'_'+ver,version(r[r.candidate_block.eq(block)&r.core],ver),'r' if mo=='C1' else 'identity',True)
    make(mo+'_'+block+'_expanded_'+ver,version(r[r.candidate_block.eq(block)],ver),'r' if mo=='C1' else 'identity',True)
 for st in ['r','partial_r']:make('C4_'+st,reg[reg.module.eq('C4')&reg.statistic.eq(st)],'r')
 for st in ['mean_difference_um','median_difference_um','matched_mean_difference_um','zero_fraction_difference']:make('C3_'+st,reg[reg.module.eq('C3')&reg.statistic.eq(st)],'identity')
 bg=read('unit_background.tsv.gz');bg=bg[bg.record_id.isin(meta.index)]
 names=['niche_all_Epi_all','niche_internal_Epi_all','niche_internal_Epi_P95','niche_front_only_all','niche_near_inside_all','niche_deep_all','UCell_internal_Epi_median','cell_composition_whole_tissue']
 names+=sorted(n for n in bg.input_family.unique() if n.startswith('UCell_P95'))
 for name in names:
  z=bg[bg.input_family.eq(name)];out[name]=(z.pivot(index='record_id',columns='feature',values='value'),'identity' if name.startswith('UCell_') else 'hellinger')
 z1=bg[bg.input_family.eq('niche_internal_Epi_P95')].pivot(index='record_id',columns='feature',values='value');z0=bg[bg.input_family.eq('niche_internal_Epi_nonP95')].pivot(index='record_id',columns='feature',values='value');out['niche_P95_minus_nonP95']=(z1-z0,'identity')
 return out,meta

def processed(mat,kind,datasets,meta,reference=None,alternate=False):
 ix=mat.index.intersection(meta.index[meta.dataset_id.isin(datasets)]);a=mat.loc[ix].copy();n_registered=len(a)
 if len(a)<3 or a.empty:return None
 cov=a.notna().mean();keep=cov.ge(.8)
 for ds in datasets:
  sub=a.loc[a.index.intersection(meta.index[meta.dataset_id.eq(ds)])]
  if not len(sub):return None
  keep&=sub.notna().mean().ge(.5)
 a=a.loc[:,keep];a=a.loc[a.notna().any(axis=1)]
 if len(a)<3 or a.shape[1]<2:return None
 missing=a.isna().mean(axis=1);x=a.fillna(a.median()).to_numpy(float)
 if kind=='r':x=transform(x,'r')
 if kind=='hellinger':x=np.sqrt(np.maximum(x,0))
 scale=x.std(0,ddof=1);vary=scale>1e-10;x=x[:,vary];features=a.columns[vary]
 if x.shape[1]<2:return None
 x-=x.mean(0)
 if kind!='hellinger' or alternate:
  sd=np.median(np.abs(x-np.median(x,axis=0)),axis=0)*1.4826 if alternate and kind!='hellinger' else x.std(0,ddof=1);sd=np.where(sd>1e-8,sd,x.std(0,ddof=1));x/=sd
 uu,ss,vt=np.linalg.svd(x,full_matrices=False);coords=uu[:,:2]*ss[:2];loads=vt[:2].T
 for j in range(2):
  sign=1
  if reference is not None:
   overlap=features.intersection(reference.index)
   if len(overlap):sign=np.sign(loads[features.get_indexer(overlap),j]@reference.loc[overlap,'PC'+str(j+1)].to_numpy()) or 1
  else:sign=np.sign(loads[np.argmax(np.abs(loads[:,j])),j]) or 1
  loads[:,j]*=sign;coords[:,j]*=sign
 return dict(X=x,index=a.index,features=features,coords=coords,loadings=loads,explained=ss[:2]**2/(ss**2).sum(),missing=missing.to_numpy(),n_registered=n_registered)

def sil(D,lab):
 vals=[]
 for i,k in enumerate(lab):
  same=np.flatnonzero(lab==k);same=same[same!=i]
  if not len(same):vals.append(0.);continue
  aa=D[i,same].mean();bb=min(D[i,lab==l].mean() for l in np.unique(lab) if l!=k);vals.append((bb-aa)/max(aa,bb) if max(aa,bb)>0 else 0.)
 return float(np.mean(vals))
def ari(a,b):
 ct=pd.crosstab(np.asarray(a),np.asarray(b)).to_numpy();choose=lambda x:x*(x-1)/2;t=choose(ct).sum();r=choose(ct.sum(1)).sum();c=choose(ct.sum(0)).sum();nn=choose(ct.sum());ex=r*c/nn if nn else 0;den=(r+c)/2-ex;return float((t-ex)/den) if den else 1.
def cluster(x):
 tree=linkage(x,method='ward');D=squareform(pdist(x));rows=[]
 for k in range(2,min(6,len(x)//3)+1):
  lab=cut_tree(tree,n_clusters=k).ravel();sz=np.bincount(lab);rows.append((k,sil(D,lab),lab,sz))
 return rows

def eta(y,ds):
 tot=((y-y.mean())**2).sum();return sum(len(y[ds==d])*(y[ds==d].mean()-y.mean())**2 for d in np.unique(ds))/tot if tot>0 else np.nan

def main():
 mats,meta=matrix_inputs();rows=[];coords=[];loads=[];stability=[];scopes=[];rng=np.random.default_rng(SEED)
 for name,(mat,kind) in mats.items():
  supported=sorted(meta.loc[mat.dropna(how='all').index.intersection(meta.index),'dataset_id'].unique());ref=processed(mat,kind,supported,meta)
  if ref is None:continue
  refl=pd.DataFrame(ref['loadings'],index=ref['features'],columns=['PC1','PC2']);scopes.append(dict(input_family=name,kind=kind,datasets='|'.join(supported),n_nonempty_combinations=2**len(supported)-1,n_units=len(ref['index']),n_features=len(ref['features'])))
  for comb in combos(supported):
   p=processed(mat,kind,comb,meta,refl)
   if p is None:
    rows.append(dict(input_family=name,combination=combo_id(comb),status='insufficient_comparable_features_or_units'));continue
   ds=meta.loc[p['index'],'dataset_id'].to_numpy();ntech=np.log1p(meta.loc[p['index'],'n_P95' if name.startswith(('C1','C2','C3','UCell_P95')) else 'n_Epi'].to_numpy());techr=stats.pearsonr(p['coords'][:,0],ntech).statistic if np.std(ntech)>0 else np.nan;missr=stats.pearsonr(p['coords'][:,0],p['missing']).statistic if np.std(p['missing'])>0 else np.nan
   cl=cluster(p['X']);base=dict(input_family=name,combination=combo_id(comb),n_units=len(p['index']),n_features=len(p['features']),n_no_retained_features=p['n_registered']-len(p['index']),mean_imputed_fraction=p['missing'].mean(),max_imputed_fraction=p['missing'].max(),PC1_variance=p['explained'][0],PC2_variance=p['explained'][1],PC1_dataset_eta2=eta(p['coords'][:,0],ds),PC1_log_coverage_r=techr,PC1_missing_fraction_r=missr,PC1_top_features='|'.join(str(p['features'][ii])+':'+format(p['loadings'][ii,0],'.4g') for ii in np.argsort(-np.abs(p['loadings'][:,0]))[:3]),PC2_top_features='|'.join(str(p['features'][ii])+':'+format(p['loadings'][ii,1],'.4g') for ii in np.argsort(-np.abs(p['loadings'][:,1]))[:3]),status='computed')
   if not cl:rows.append(base|dict(K=0,silhouette=np.nan,cluster_sizes='',n_rare_clusters=0,cluster_dataset_table='',status='PCA_only_n_too_small_for_K_ge2'))
   for k,sv,lab,sz in cl:rows.append(base|dict(K=k,silhouette=sv,cluster_sizes='|'.join(map(str,sz)),n_rare_clusters=int((sz<=2).sum()),cluster_dataset_table=';'.join(str(i)+':'+','.join(d+'='+str(sum((lab==i)&(ds==d))) for d in sorted(set(ds))) for i in range(k))))
  if name in REP:
   p=ref;choices=cluster(p['X']);k,sv,lab,sz=max(choices,key=lambda q:q[1]);z=pd.DataFrame(p['coords'],index=p['index'],columns=['PC1','PC2']);z['record_id']=z.index;z['input_family']=name;z['cluster']=lab;z['K']=k;z['missing_fraction']=p['missing'];coords.append(z.reset_index(drop=True));l=refl.reset_index(names='feature');l['input_family']=name;loads.append(l)
   for it in range(100):
    idx=np.sort(rng.choice(len(lab),max(k+1,int(.8*len(lab))),replace=False));sub=cut_tree(linkage(p['X'][idx],method='ward'),n_clusters=k).ravel();features=rng.choice(p['X'].shape[1],max(2,int(.8*p['X'].shape[1])),replace=False);flab=cut_tree(linkage(p['X'][:,features],method='ward'),n_clusters=k).ravel();stability.append(dict(input_family=name,check='unit_subset_80percent',value=ari(lab[idx],sub)));stability.append(dict(input_family=name,check='feature_subset_80percent',value=ari(lab,flab)))
   for ds in supported:
    ix=np.flatnonzero(meta.loc[p['index'],'dataset_id'].to_numpy()!=ds)
    if len(ix)>k:
     la=cut_tree(linkage(p['X'][ix],method='ward'),n_clusters=k).ravel();stability.append(dict(input_family=name,check='leave_dataset_'+ds,value=ari(lab[ix],la)))
   for seed in range(10):
    cent,la=kmeans2(p['X'],k,iter=100,minit='++',seed=SEED+seed);stability.append(dict(input_family=name,check='kmeans_same_K',value=ari(lab,la)))
   alt=processed(mat,kind,supported,meta,refl,True)
   if alt is not None:
    la=cut_tree(linkage(alt['X'],method='ward'),n_clusters=k).ravel();stability.append(dict(input_family=name,check='alternative_feature_scale',value=ari(lab,la)))
   for altname in [name.replace('_full','_noMP10'),name.replace('_partial_r','_r')]:
    if altname!=name and altname in mats:
     ma,ki=mats[altname];ap=processed(ma,ki,supported,meta,refl)
     if ap is not None:
      la=cut_tree(linkage(ap['X'],method='ward'),n_clusters=k).ravel();common=p['index'].intersection(ap['index']);stability.append(dict(input_family=name,check='alternate_input_'+altname,value=ari(lab[p['index'].get_indexer(common)],la[ap['index'].get_indexer(common)])))
  print('heterogeneity',name,len(supported),flush=True)
 save(rows,'heterogeneity_configurations.tsv.gz');save(scopes,'heterogeneity_input_scope.tsv');save(pd.concat(coords,ignore_index=True),'representative_coordinates.tsv');save(pd.concat(loads,ignore_index=True),'representative_loadings.tsv')
 st=pd.DataFrame(stability).groupby(['input_family','check']).value.agg(['count','mean','median','min','max']).reset_index();save(st,'stability_summary.tsv')
 # Block-balanced integration is fixed after separate representations: two PCs from each input, each block normalized to equal total variance.
 bases=['niche_all_Epi_all','B1_R0_partial_r_full','B3_Epi_distance_partial_r_full'];zs=[]
 for name in bases:
  ma,ki=mats[name];p=processed(ma,ki,sorted(meta.loc[ma.dropna(how='all').index,'dataset_id'].unique()),meta);zz=pd.DataFrame(p['coords'],index=p['index'],columns=[name+'_PC1',name+'_PC2']);zz/=np.sqrt((zz.var()).sum());zs.append(zz)
 joint=pd.concat(zs,axis=1,join='inner');p=processed(joint,'hellinger' if False else 'already_scaled',sorted(meta.loc[joint.index,'dataset_id'].unique()),meta)
 # preserve equal block variance: no additional per-column scaling
 X=joint.to_numpy();X-=X.mean(0);uu,ss,vt=np.linalg.svd(X,full_matrices=False);co=pd.DataFrame(uu[:,:2]*ss[:2],index=joint.index,columns=['PC1','PC2']);co['record_id']=co.index;save(co.reset_index(drop=True),'joint_balanced_coordinates.tsv');save(pd.DataFrame(vt[:2].T,index=joint.columns,columns=['PC1','PC2']).reset_index(names='feature'),'joint_balanced_loadings.tsv')
 print('heterogeneity complete',len(rows),flush=True)
if __name__=='__main__':main()
