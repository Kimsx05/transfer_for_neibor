"""Expression-cluster and measured-marker review of model candidates; no physical neighbors."""
from common import *
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse,stats
from marker_rules import MAJOR,SUBTYPE,ANCHORS
np.random.seed(123);dictionary=pd.read_csv(ROOT/'metadata/reference_label_dictionary.tsv',sep='\t').set_index('reference_celltype_v2').major_label.to_dict()
markers=sorted(set(sum(MAJOR.values(),[])+sum(SUBTYPE.values(),[])+sum(ANCHORS.values(),[])));mi={g:i for i,g in enumerate(markers)}
for name in ['Cervilla','D1_Johnson']:
 out=ROOT/name;a=ad.read_h5ad(out/'cells.h5ad');o=a.obs;canon=a.var.reference_symbol.to_numpy();rows=[];cols=[]
 for j,g in enumerate(markers):
  for i in np.flatnonzero(canon==g):rows.append(i);cols.append(j)
 transform=sparse.csr_matrix((np.ones(len(rows),dtype='int32'),(rows,cols)),shape=(a.n_vars,len(markers)));mx=(a.X@transform).toarray();det=mx>0;tot=a.obs.n_counts.to_numpy();major_cluster_support=np.zeros(a.n_obs,bool);fine_cluster_support=np.zeros(a.n_obs,bool);cluster_reasons=np.repeat('',a.n_obs).astype(object);group_rows=[];marker_rows=[]
 def group_support(ix,genes,fraction=.03,mincells=3):
  jj=[mi[g] for g in genes if g in mi];n=len(ix);hits=det[ix][:,jj].sum(0);required=max(mincells,int(np.ceil(n*fraction)));return [g for g,h in zip(genes,hits) if h>=required],hits,required
 def program_fractions(ix,genes):
  jj=[mi[g] for g in genes if g in mi];n=det[ix][:,jj].sum(1);return float((n>=1).mean()),float((n>=2).mean())
 # Globally defined expression clusters: majority alone never promotes a label.
 for cluster,ix in o.groupby('expression_cluster',observed=True).indices.items():
  if cluster=='not_in_structure_core':continue
  n=len(ix);pred=o.predicted_major.to_numpy()[ix];freq=pd.Series(pred).value_counts(normalize=True)
  for major,purity in freq.items():
   found,hits,req=group_support(ix,MAJOR[major]);one,two=program_fractions(ix,MAJOR[major]);coherent=n>=30 and purity>=.70 and len(found)>=3 and one>=.30 and two>=.10
   if major=='B':coherent=coherent and any(g in found for g in ['CD79A','MS4A1','CD19'])
   # support only cells independently predicted into this lineage; not all members of cluster.
   if coherent:major_cluster_support[ix[pred==major]]=True
   group_rows.append(dict(level='global_major',cluster=str(cluster),candidate=major,n=n,candidate_fraction=float(purity),support=bool(coherent),marker_genes_supported=';'.join(found),min_marker_positive_cells=req,program_any_fraction=one,program_two_fraction=two))
 # Within-lineage clusters retain individual model candidates and require a substantial coherent program.
 for cluster,ix in o.groupby('lineage_cluster',observed=True).indices.items():
  if cluster in ['not_in_structure_core','not_clustered']:continue
  n=len(ix);major=cluster.split('::')[0];found,_,req=group_support(ix,MAJOR[major]);one,two=program_fractions(ix,MAJOR[major]);major_coherent=n>=20 and len(found)>=3 and one>=.30 and two>=.10
  if major=='B':major_coherent=major_coherent and any(g in found for g in ['CD79A','MS4A1','CD19'])
  if major_coherent:major_cluster_support[ix]=True
  pred=o.predicted_subtype.to_numpy()[ix];freq=pd.Series(pred).value_counts(normalize=True)
  group_rows.append(dict(level='within_lineage_major',cluster=cluster,candidate=major,n=n,candidate_fraction=1.,support=bool(major_coherent),marker_genes_supported=';'.join(found),min_marker_positive_cells=req,program_any_fraction=one,program_two_fraction=two))
  for subtype,purity in freq.items():
   if subtype not in SUBTYPE:continue
   found,hits,req=group_support(ix,SUBTYPE[subtype],fraction=.03,mincells=3);anchor=ANCHORS.get(subtype,[]);anchor_ok=not anchor or len(group_support(ix,anchor,fraction=.03,mincells=3)[0])>=1
   one,two=program_fractions(ix,SUBTYPE[subtype]);supported=major_coherent and n>=30 and purity>=.60 and len(found)>=2 and one>=.20 and two>=.05 and anchor_ok and dictionary.get(subtype)==major
   if subtype=='GC B':
    found,_,req=group_support(ix,SUBTYPE[subtype],fraction=.05,mincells=3);supported=major_coherent and n>=30 and purity>=.60 and len(found)>=2 and one>=.20 and two>=.05
   if supported:
    ii=ix[pred==subtype];fine_cluster_support[ii]=True;cluster_reasons[ii]='coherent_expression_cluster_and_observed_'+','.join(found)
   group_rows.append(dict(level='within_lineage_subtype',cluster=cluster,candidate=subtype,n=n,candidate_fraction=float(purity),support=bool(supported),marker_genes_supported=';'.join(found),min_marker_positive_cells=req,program_any_fraction=one,program_two_fraction=two))
  for gene,j in mi.items():marker_rows.append(dict(section=name,lineage_cluster=cluster,major=major,n=n,gene=gene,detected_n=int(det[ix,j].sum()),detection_fraction=float(det[ix,j].mean()),raw_count_sum=int(mx[ix,j].sum()),mean_log1p_10k=float(np.log1p(mx[ix,j]/np.maximum(tot[ix],1)*1e4).mean())))
 mp=o.predicted_major.to_numpy();fp=o.predicted_subtype.to_numpy();own_major=np.array([o['observed_markers__'+x.replace('/','_').replace(' ','_')].iloc[i] for i,x in enumerate(mp)])
 major_ok=o.major_calibrated.to_numpy()&((own_major>=2)|major_cluster_support)
 # Direct two-marker evidence is one route, not a universal requirement. B requires an anchor for the direct route.
 banchors=det[:,[mi[g] for g in ['CD79A','MS4A1','CD19']]].any(1);major_ok[(mp=='B')&~major_cluster_support]&=banchors[(mp=='B')&~major_cluster_support]
 direct_fine=np.zeros(a.n_obs,bool)
 for sub in SUBTYPE:
  ix=np.flatnonzero(fp==sub)
  if not len(ix):continue
  jj=[mi[g] for g in SUBTYPE[sub]];yes=det[ix][:,jj].sum(1)>=2
  if sub in ANCHORS:yes&=det[ix][:,[mi[g] for g in ANCHORS[sub]]].any(1)
  direct_fine[ix]=yes
 compatible=np.array([dictionary.get(f)==m for f,m in zip(fp,mp)]);fine_ok=major_ok&o.subtype_calibrated.to_numpy()&compatible&(direct_fine|fine_cluster_support)
 broad=np.isin(mp,['Epithelial','Plasma','Mast']);fine_ok[broad]=False;mixed=o.mixed_or_ambiguous.to_numpy()
 majorlabel=np.where(major_ok,mp,'Unknown').astype(object);final=np.array([m.replace('/','_')+'_unresolved' if ok else 'Unknown' for m,ok in zip(mp,major_ok)],dtype=object);status=np.where(major_ok,'subtype_unresolved','lineage_unresolved').astype(object);evidence=np.repeat('insufficient_reference_or_observed_lineage_evidence',a.n_obs).astype(object)
 evidence[major_ok&major_cluster_support]='reference_depth_calibrated_and_coherent_expression_cluster_major_markers';evidence[major_ok&~major_cluster_support]='reference_depth_calibrated_and_cell_observed_major_markers'
 final[major_ok&broad]=mp[major_ok&broad];status[major_ok&broad]='reference_broad_identity';final[fine_ok]=fp[fine_ok];status[fine_ok]='subtype_supported';evidence[fine_ok&fine_cluster_support]=cluster_reasons[fine_ok&fine_cluster_support];evidence[fine_ok&~fine_cluster_support]='depth_calibrated_lineage_subtype_with_cell_observed_marker_combination'
 final[mixed]='mixed_or_ambiguous';status[mixed]='mixed_or_ambiguous';evidence[mixed]=o.mixed_reason.to_numpy()[mixed]
 o['major_label']=majorlabel;o['final_label']=final;o['subtype_label']=final;o['annotation_status']=status;o['annotation_evidence']=evidence;o['major_cluster_support']=major_cluster_support;o['subtype_cluster_support']=fine_cluster_support;o['lineage_unresolved']=~major_ok&~mixed;o['subtype_unresolved']=major_ok&~fine_ok&~broad&~mixed;o['annotation_review_status']='computational_review_not_pathologist_validated'
 o['uncertainty_reason']=np.where(mixed,o.mixed_reason,np.where(~major_ok,'major_depth_calibration_or_coherent_expression_evidence_insufficient',np.where(broad,'reference_has_no_finer_identity',np.where(~fine_ok,'subtype_calibration_or_cluster_marker_evidence_insufficient',''))))
 a.uns['analysis_status']='cell_annotation_v2_review_complete'
 a.uns['annotation_rule_summary']='Candidate preserved; section-specific reference thinning calibration plus per-cell or coherent expression-cluster observed markers; no physical majority vote; broad epithelial malignancy unassigned; cluster evidence only helps independently consistent candidates.'
 a.write_h5ad(out/'cells.h5ad',compression='gzip');o.to_parquet(out/'cell_labels.parquet',compression='zstd');pd.DataFrame(group_rows).to_csv(out/'metadata/cluster_annotation_review.tsv',sep='\t',index=False);pd.DataFrame(marker_rows).to_csv(out/'metadata/cluster_marker_evidence.tsv',sep='\t',index=False)
 # Sparse observed-expression differential summaries for review. Only top25 per lineage cluster retained.
 x=a.X.astype('float32').tocsr();x=sparse.diags(1e4/np.maximum(tot,1)).dot(x).tocsr();x.data=np.log1p(x.data);de=[]
 for major,mirows in o[o.scvi_training_core].groupby('predicted_major',observed=True).groups.items():
  midx=o.index.get_indexer(mirows);xx=x[midx];n=len(midx)
  if n<100:continue
  total=np.asarray(xx.sum(0)).ravel();sq=xx.copy();sq.data**=2;total2=np.asarray(sq.sum(0)).ravel();del sq
  for cl in sorted(o.iloc[midx].lineage_cluster.unique()):
   sub=midx[o.iloc[midx].lineage_cluster.eq(cl).to_numpy()];n1=len(sub);n2=n-n1
   if min(n1,n2)<10:continue
   y=x[sub];s1=np.asarray(y.sum(0)).ravel();y2=y.copy();y2.data**=2;v1=np.maximum((np.asarray(y2.sum(0)).ravel()-s1*s1/n1)/(n1-1),0);s2=total-s1;v2=np.maximum((total2-np.asarray(y2.sum(0)).ravel()-s2*s2/n2)/(n2-1),0);delta=s1/n1-s2/n2;score=delta/np.sqrt(v1/n1+v2/n2+1e-8);idx=np.argsort(-score)[:25]
   for j in idx:de.append(dict(lineage=major,cluster=cl,n=n1,reference_n=n2,feature_id=a.var_names[j],gene_symbol=a.var.gene_symbol.iloc[j],mean_log_difference=float(delta[j]),descriptive_welch_score=float(score[j]),note='observed log-normalized expression; exploratory review ranking, not confirmatory testing'))
 pd.DataFrame(de).to_csv(out/'metadata/lineage_cluster_top_markers.tsv',sep='\t',index=False)
 summary=dict(section=name,retained=a.n_obs,low_information=int(o.low_information.sum()),status=o.annotation_status.value_counts().to_dict(),major=o.major_label.value_counts().to_dict(),supported_subtypes=o.loc[o.annotation_status.eq('subtype_supported'),'final_label'].astype(str).value_counts().to_dict(),low_information_status=o.loc[o.low_information,'annotation_status'].value_counts().to_dict(),scvi_core=int(o.scvi_training_core.sum()),umap_unmapped=int((~o.umap_mapped).sum()),cluster_supported_subtype_cells=int((fine_ok&fine_cluster_support&~mixed).sum()),cell_marker_supported_subtype_cells=int((fine_ok&~fine_cluster_support&~mixed).sum()),limitations='Computational suggested labels; no spatial voting; reference and cluster evidence are not independent ground truth; sparse populations may remain unresolved')
 jwrite(out/'metadata/annotation_summary.json',summary);print('ANNOTATION_REVIEW',summary,flush=True)
jwrite(ROOT/'metadata/annotation_rules.json',dict(seed=123,reference_version='cell2location_v2_frozen_20260916__scRNA32',major_markers=MAJOR,subtype_markers=SUBTYPE,anchors=ANCHORS,cluster_support='Major global >=30 cells and >=70% candidate with >=3 lineage genes; within-lineage >=20 cells with >=3 major genes. Both require >=30% cells with any lineage marker and >=10% with >=2. Fine >=30 cells, >=60% subtype candidates, >=2 subtype genes (each >=max(3,3%cluster)), appropriate anchor and >=20% cells with any subtype marker plus >=5% with>=2; GC uses5% and2GC genes plusB lineage anchor.',individual_assignment='Major calibrated and either own>=2 markers or coherent cluster; subtype calibrated and consistent and either own>=2 subtype markers+anchor or supported expression cluster. No automatic whole-cluster label transfer.',low_information='10-30 observed genes retained; shallow calibration applied; noncore not in graph/UMAP, latent/spatial kept',malignancy='not inferred',cell_markers_are_not_universal_veto=True))
