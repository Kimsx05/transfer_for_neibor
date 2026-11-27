from base_data import *
from sklearn.cluster import MiniBatchKMeans
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
platform=sys.argv[1];samples=XE if platform=='xenium' else HD
zlist=[];lat=[];markerlist=[];anchorlist=[];avail=[]
if platform=='hd':
 history=pd.read_parquet(A/'cell_annotation_v2_1/cell_annotations.parquet');history=history.set_index('cell_id',drop=False) if 'cell_id' in history else history
evaluation=pd.read_csv(O/'tables/reference_evaluation.tsv',sep='\t').set_index(['head','class'])
for s in samples:
 a=dataset(s);d=paths(s);ct=pd.read_parquet(d/'prediction.parquet');assert a.obs_names.equals(ct.index)
 z=a.obs[['x_um','y_um','nFeature']].copy();z['cell_id']=z.index;z['sample']=s;z['source']=a.obs.source if platform=='hd' else 'original_segmentation'
 z['n_counts']=a.obs.n_counts if platform=='hd' else a.obs.nGeneTranscripts;z['area_um2']=a.obs.area_um2 if platform=='hd' else a.obs.cell_area
 z['geometry_flag']=a.obs.geometry_flag if platform=='hd' else a.obs.area_outlier|a.obs.boundary_degenerate
 z['eligible_expression']=eligible(a,s)
 if platform=='hd':
  old=history.reindex(a.obs.old_id.values);z['old_id']=a.obs.old_id;z['old_major']=old.proposed_major.fillna('no_prior_label').values;z['old_label']=old.proposed_subtype.mask(old.proposed_subtype.eq('Unknown'),old.proposed_major).fillna('no_prior_label').values
  w=pd.read_parquet(d/'tacco_weights.parquet');finite=w.notna().all(1);values=w.fillna(-1).values;ix=np.argsort(-values,axis=1,kind='stable')
  z['other_major']=np.where(finite,w.columns[ix[:,0]],'no_fit');z['other_major_support']=np.where(finite,values[np.arange(len(z)),ix[:,0]],np.nan);z['tacco_second']=np.where(finite,w.columns[ix[:,1]],'no_fit')
  z['tacco_runner']=np.where(finite,values[np.arange(len(z)),ix[:,1]],np.nan);z['other_sub']=old.proposed_subtype.fillna('').values
 else:
  for n in ['major_candidate','major_support','major_second','major_runner','candidate','support','second','runner','flag','major','supported','marker_detected','anchor_detected','boundary_vertices','segmentation_method','nucleus_count']:z['old_'+n]=a.obs[n]
  z['old_label']=a.obs.supported;z['other_major']=a.obs.major_candidate;z['other_major_support']=a.obs.major_support;z['other_sub']=a.obs.candidate
 z=z.join(ct);zg=np.load(d/'latent.npy');lat.append(zg);zlist.append(z.reset_index(drop=True))
 def sig(gs):
  jj=a.var_names.get_indexer(gs);jj=jj[jj>=0];return a.X[:,jj].getnnz(1) if len(jj) else np.full(len(z),np.nan)
 markerlist.append(pd.DataFrame({k:sig(v) for k,v in M.items()}));anchorlist.append(pd.DataFrame({k:sig(str(v).split()) for k,v in tax.anchors.items()}))
 for k,gs in {**{k+' major':v for k,v in M.items()},**{k:str(v).split() for k,v in tax.anchors.items()}}.items():
  avail.append((s,k,' '.join(g for g in gs if g in a.var_names),' '.join(g for g in gs if g not in a.var_names)))
 del a
write(pd.DataFrame(avail,columns=['sample','identity','measured','not_measured']),platform+'_marker_availability.tsv')
z=pd.concat(zlist,ignore_index=True);mk=pd.concat(markerlist,ignore_index=True);ak=pd.concat(anchorlist,ignore_index=True);L=np.vstack(lat)
def pick(t,labels):
 out=np.full(len(t),np.nan)
 for k in t:
  i=(labels==k).fillna(False).to_numpy();out[i]=t[k].values[i]
 return out
ct=z.ct_major_first;other=z.other_major;im=pick(mk,ct);om=pick(mk,other);agree=ct.eq(other);good=z.eligible_expression;second=z.ct_major_margin>0
if platform=='hd':mix=(z.tacco_runner>=.25)&(z.other_major_support<.65)&(om>=1)&(pick(mk,z.tacco_second)>=1)
else:mix=(z.old_flag.eq('mixed')|((om>=2)&(im>=2)&~agree))
mixed=z.geometry_flag|mix
cand=ct.copy();choose_other=~agree&(om>=2)&(np.nan_to_num(im)==0)&(z.other_major_support>=.65);cand[choose_other]=other[choose_other]
q=pick(mk,cand);conflict=~agree&(im>=1)&(om>=1)&(z.other_major_support>=.60)&~choose_other
z['lineage_cluster']='unclustered';z['candidate_major']=cand;z['within_major_markers']=q;z['reference_agree']=agree
for k in sorted(cand.dropna().unique()):
 idx=np.flatnonzero((cand==k)&good);n=min(max(2,sum(tax.major==k)*2),max(1,len(idx)//100),14)
 if len(idx)<20:continue
 km=MiniBatchKMeans(n_clusters=n,random_state=1729,n_init=5,batch_size=1024);yy=km.fit_predict(L[idx]);z.loc[idx,'lineage_cluster']=[f'{k}.{v:02d}' for v in yy]
cl=z.lineage_cluster;cl_ok=pd.Series(q>=1).groupby(cl).transform('mean')>=.55
known=(good&~mixed&~conflict&((agree&((q>=1)|cl_ok))|((q>=2)&((z.ct_major_margin>=1)|choose_other))|((q>=1)&cl_ok&(z.ct_major_margin>=2)&(np.nan_to_num(om)==0))))
z['major_identity']=np.where(known,cand,np.where(mixed|conflict,'Mixed','Unknown'))
z['review_level']=np.where(known,'major_only',np.where(mixed|conflict,'mixed_or_suspect','major_unresolved'))
z['fine_identity']='';fine=z.ct_sub_first.fillna('');anc=pick(ak,fine)
for k in tax.index:
 par=tax.at[k,'major'];sel=known&(cand==par)&(fine==k)
 if not sel.any() or not ak[k].notna().any():continue
 refok=True
 if platform=='xenium':
  if (par,k) not in evaluation.index:continue
  row=evaluation.loc[(par,k)];refok=row.precision>=.55 and row.recall>=.3 and row.held_n>=5
 if not refok:continue
 cluster_evidence=(ak[k].fillna(0)>0).groupby(cl).transform('mean');base=(ak.loc[(cand==par)&good,k].fillna(0)>0).mean()
 agree_f=z.other_sub.eq(k)
 accept=sel&(z.ct_sub_margin>=.5)&(((anc>=1)&(agree_f|(cluster_evidence>=max(.08,base*.9))))|(agree_f&(cluster_evidence>=max(.25,base*1.15))&(q>=1)))
 z.loc[accept,'review_level']='fine_supported';z.loc[accept,'fine_identity']=k
z['display_label']=np.where(z.review_level.eq('fine_supported'),z.fine_identity,z.major_identity)
z['review_origin']='Codex_cluster_rules_v1';z['human_review']='pending';z['fine_anchor_count']=anc
z['conflict']=np.select([z.geometry_flag,mix,conflict,~good,~agree,z.ct_sub_first.ne(z.other_sub)],['geometry_outlier','mixed_component_or_old_mixed','cross_lineage_markers','low_expression','reference_disagreement','fine_disagreement'],default='none')
for s in samples:
 p=z[z['sample']==s].copy();p.set_index('cell_id').to_parquet(paths(s)/'annotation.parquet')
write(z.groupby(['sample','old_major','old_label','major_identity','display_label','review_level'],dropna=False).size().rename('n').reset_index(),platform+'_label_changes.tsv')
write(z.groupby(['sample','conflict','review_level']).size().rename('n').reset_index(),platform+'_conflicts.tsv')
tables=[]
for group,i in z.groupby('lineage_cluster').indices.items():
 sub=z.iloc[i];markers=mk.iloc[i].gt(0).mean().sort_values(ascending=False);labels=sub.display_label.value_counts(normalize=True)
 tables.append(dict(cluster=group,n=len(i),sample_n=sub['sample'].nunique(),first=labels.index[0],first_fraction=labels.iloc[0],reference_agree=sub.reference_agree.mean(),mixed_fraction=sub.review_level.eq('mixed_or_suspect').mean(),major_only_fraction=sub.review_level.eq('major_only').mean(),marker_families='; '.join(f'{k} {v:.2f}' for k,v in markers.head(4).items()),principal_conflict=sub.conflict.value_counts().index[0]))
write(pd.DataFrame(tables),platform+'_cluster_review.tsv')

