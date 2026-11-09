from pathlib import Path
import sys,json,itertools
import numpy as np,pandas as pd
from scipy import ndimage
from sklearn.metrics import adjusted_rand_score,normalized_mutual_info_score
root=Path(__file__).resolve().parents[1];s=sys.argv[1];o=root/s;base=root.parent/'hd_bin8_comparison_v1'/s
read=lambda p:pd.read_csv(p,sep='\t')
d=read(base/'analysis_metadata.tsv.gz').set_index('bin_id',drop=False)
new=read(o/'type_metadata.tsv.gz').set_index('bin_id',drop=False)
for fn in ['spatial_metadata.tsv.gz','neighbor_statistics.tsv.gz']:
 z=read(o/fn).set_index('bin_id');assert z.index.equals(d.index);new=new.join(z)
for tag in ['l00','l02','l08']:
 z=read(o/f'{tag}_clusters.tsv.gz').set_index('bin_id');assert z.index.equals(d.index);new[f'banksy_{tag}_cluster']=z.cluster.astype(str)
assert new.index.equals(d.index)
for c in new.columns:
 if c not in d:d[c]=new[c]
for c in ['native_edge_distance_um','neighbor_max_um']:assert c in d
branches={'historical':'cluster_projected',**{t:f'banksy_{t}_cluster' for t in ['l00','l02','l08']}}
summaries=[];comps=[];ct=[];spmetrics=[];conf=[]
for tag,col in branches.items():
 d[col]=d[col].astype(str)
 for cl,z in d.groupby(col,sort=True):
  r=dict(branch=tag,cluster=cl,n=len(z),UMI=int(z.nCount_RNA.sum()),nCount_median=z.nCount_RNA.median(),nFeature_median=z.nFeature_RNA.median(),nFeature_q10=z.nFeature_RNA.quantile(.1),nFeature_q90=z.nFeature_RNA.quantile(.9),low_info_fraction=(z.nFeature_RNA<20).mean(),edge16_fraction=(z.native_edge_distance_um<=16).mean(),long_neighbor80_fraction=(z.neighbor_max_um>80).mean(),epi_fib_ge2_fraction=z.epi_fib_ge2.mean(),marker_conflict_fraction=z.marker_conflict.mean(),forward_conflict_n=(z.reference_conflict_group=='BLCA_Fibroblast_HPCA_Epithelial').sum())
  for c in [c for c in d if c.startswith('detected_')]:r[c+'_ge2_fraction']=(z[c]>=2).mean()
  r['Plasma_core_supported_ge2_fraction']=((z.detected_Plasma>=2)&(z.plasma_core_detected>0)).mean()
  for rt in ['blca','hpca']:
   labs=z[f'singler_{rt}_pruned_labels'].fillna('Unknown');co=z[f'coarse_{rt}'];r[f'{rt}_top']=labs.value_counts().index[0];r[f'{rt}_top_fraction']=labs.value_counts(normalize=True).iloc[0];r[f'{rt}_unknown_fraction']=(co=='Unknown').mean()
   for mode,k in [('major',f'singler_{rt}_pruned_labels'),('coarse',f'coarse_{rt}')]:
    for la,zz in z.groupby(z[k].fillna('Unknown')):comps.append(dict(branch=tag,cluster=cl,reference=rt,level=mode,label=la,n=len(zz),UMI=int(zz.nCount_RNA.sum()),within_cluster_fraction=len(zz)/len(z),within_cluster_UMI_fraction=zz.nCount_RNA.sum()/z.nCount_RNA.sum()))
  summaries.append(r)
 # Adjacent grid-edge agreement is descriptive expected smoothness, not accuracy.
 vals=d[col].astype('category').cat.codes.to_numpy()+1;g=np.zeros((int(d.array_row.max()+1),int(d.array_col.max()+1)),np.int32);g[d.array_row,d.array_col]=vals
 same=total=0
 for di,dj in [(1,0),(0,1)]:
  aa=g[di:,dj:];bb=g[:g.shape[0]-di or None,:g.shape[1]-dj or None];ok=(aa>0)&(bb>0);same+=((aa==bb)&ok).sum();total+=ok.sum()
 fragments=sum(ndimage.label(g==v,structure=np.ones((3,3)))[1] for v in np.unique(vals));spmetrics.append(dict(branch=tag,clusters=d[col].nunique(),adjacent_pairs=total,same_cluster_pairs=same,adjacent_agreement=same/total,connected_fragments=fragments))
 for gg,z in d.groupby('reference_conflict_group'):
  for cl,zz in z.groupby(col):conf.append(dict(branch=tag,conflict_group=gg,cluster=cl,n=len(zz),UMI=int(zz.nCount_RNA.sum()),within_group_fraction=len(zz)/len(z)))
for (ta,ca),(tb,cb) in itertools.combinations(branches.items(),2):
 z=d.groupby([ca,cb],observed=True).size().reset_index(name='n');z.columns=['cluster_A','cluster_B','n'];z['branch_A']=ta;z['branch_B']=tb;ct.append(z)
pd.DataFrame(summaries).to_csv(o/'cluster_summary.tsv',sep='\t',index=False);pd.DataFrame(comps).to_csv(o/'cluster_reference_composition.tsv',sep='\t',index=False);pd.DataFrame(spmetrics).to_csv(o/'spatial_structure_metrics.tsv',sep='\t',index=False);pd.concat(ct).to_csv(o/'cluster_crosstabs.tsv.gz',sep='\t',index=False);pd.DataFrame(conf).to_csv(o/'conflict_cluster_distribution.tsv',sep='\t',index=False)
pd.DataFrame([dict(branch_A=ta,branch_B=tb,ARI=adjusted_rand_score(d[ca],d[cb]),NMI=normalized_mutual_info_score(d[ca],d[cb])) for (ta,ca),(tb,cb) in itertools.combinations(branches.items(),2)]).to_csv(o/'cluster_agreement.tsv',sep='\t',index=False)
conf=[]
for (gg,band),z in d.groupby(['reference_conflict_group','banksy_information_band']):
 r=dict(group=gg,information_band=band,n=len(z),UMI=int(z.nCount_RNA.sum()),nCount_median=z.nCount_RNA.median(),nFeature_median=z.nFeature_RNA.median(),edge16_fraction=(z.native_edge_distance_um<=16).mean(),epi_fib_ge2_fraction=z.epi_fib_ge2.mean(),long_neighbor80_fraction=(z.neighbor_max_um>80).mean());r.update({c+'_ge2_fraction':(z[c]>=2).mean() for c in d if c.startswith('detected_')});conf.append(r)
pd.DataFrame(conf).to_csv(o/'conflict_information_summary.tsv',sep='\t',index=False)
# New fields only; spatial coordinates are retained to verify the attachment.
new['x_um']=d.x_um;new['y_um']=d.y_um
new.to_csv(o/'metadata.tsv.gz',sep='\t',index=False)
print(s,'SUMMARY_COMPLETE',len(d));print(pd.DataFrame(spmetrics).to_string(index=False))
