from base import *
from sklearn.metrics import adjusted_rand_score,normalized_mutual_info_score
from scipy.cluster.hierarchy import linkage,dendrogram
from scipy.spatial.distance import pdist,squareform
from scipy.optimize import linear_sum_assignment
import itertools,shutil,zipfile
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
C=OUT/'04_comparison';F=C/'figures';F.mkdir(exist_ok=True);m,f,ref=inputs();m=pd.read_parquet(OUT/'01_features/metadata.parquet');q=np.concatenate([f['q1'],f['q2']],1)
import subprocess
subprocess.run(['/usr/bin/python3.12',str(OUT/'code/reconcile_selection.py')],check=True)
configs=sorted((OUT/'03_main').iterdir());assert len(configs)==18
cfg={p.name:json.loads((p/'working_config.json').read_text()) for p in configs};labels={p.name:pd.read_parquet(p/'labels.parquet') for p in configs};mapping=pd.read_csv(OLD/'00_preflight/type_mapping.tsv',sep='\t');major=mapping.set_index('subtype').major.reindex(TYPES);majors=sorted(major.unique())

def palette(k):return plt.cm.hsv(np.arange(k)/max(k,1)*.94)
def save(fig,p,pdf=True):
 fig.savefig(str(p)+'.png',dpi=150,transparent=True)
 if pdf:fig.savefig(str(p)+'.pdf',transparent=True)
 plt.close(fig)
def discrete(ax,x,y,lab,title,k,size=2,legend=True):
 colors=palette(k)
 for c in range(k):
  ix=lab==c;ax.scatter(np.asarray(x)[ix],np.asarray(y)[ix],s=size,c=[colors[c]],label=str(c),rasterized=True)
 ix=lab<0
 if ix.any():ax.scatter(np.asarray(x)[ix],np.asarray(y)[ix],s=size,c='#bdbdbd',label='UNASSIGNED',rasterized=True)
 ax.set_title(title,fontsize=9)
 if legend:ax.legend(fontsize=5,markerscale=2,ncol=min(5,k+1),loc='upper right',framealpha=.65)

# Required master mapping includes the two target sections with no anchors.
over=pd.read_csv(ROOT/'spatial_ecology_top20_v2/05_07_shared_K12/v1/section_identity_overlay.tsv',sep='\t');cov=pd.read_csv(OUT/'00_preflight/batch_key_mapping.tsv',sep='\t');missing=over[~over.section_id.isin(cov.section_id)].copy()
for r in missing.itertuples():
 row=dict(dataset_id=r.dataset,section_id=r.section_id,sample_uid=r.dataset+'::'+r.section_id,n_common_anchors=0,n_reference_anchors=0,n_complete=0,original_sample_id=r.source_sample_id,specimen_uid_if_known='',patient_uid=r.patient_uid,analysis_unit_uid=r.analysis_unit_uid,analysis_unit_type=r.analysis_unit_type,source_patient_identity_known=r.source_patient_identity_known,sample_key_definition='dataset_id::section_id (independent section)');cov=pd.concat([cov,pd.DataFrame([row])],ignore_index=True)
assert len(cov)==63;tab(OUT/'00_preflight/batch_key_mapping.tsv',cov);tab(OUT/'00_preflight/coverage_by_dataset_sample.tsv',cov);shutil.copy2(OLD/'00_preflight/type_mapping.tsv',OUT/'01_features/type_mapping.tsv')

pairs=[];corr=[];complete=[];dictionary=[];supports=[];allprofiles=[];allsource=[];integs=[];samek=[];legacy=[];selections=[]
for a,b in itertools.combinations(labels,2):
 la,lb=labels[a],labels[b];common=(la.accepted_label>=0)&(lb.accepted_label>=0)
 for cohort,ix,col in [('reference_native',ref,'reference_native_label'),('all_candidate',np.arange(len(m)),'candidate_label'),('shared_accepted',np.flatnonzero(common),'accepted_label')]:
  pairs.append(dict(config1=a,config2=b,cohort=cohort,n=len(ix),k1=la[col].iloc[ix].nunique(),k2=lb[col].iloc[ix].nunique(),ARI=adjusted_rand_score(la[col].iloc[ix],lb[col].iloc[ix]),NMI=normalized_mutual_info_score(la[col].iloc[ix],lb[col].iloc[ix]),coverage1=float((la.accepted_label>=0).mean()),coverage2=float((lb.accepted_label>=0).mean())))
 cont=pd.crosstab(la.candidate_label,lb.candidate_label);aa,bb=linear_sum_assignment(-cont.to_numpy())
 for i,j in zip(aa,bb):corr.append(dict(config1=a,cluster1=int(cont.index[i]),config2=b,cluster2=int(cont.columns[j]),overlap_n=int(cont.iloc[i,j]),fraction_of_1=float(cont.iloc[i,j]/cont.iloc[i].sum()),fraction_of_2=float(cont.iloc[i,j]/cont.iloc[:,j].sum()),method='maximum-overlap matching; IDs remain config-specific'))
tab(C/'pairwise_ARI_NMI.tsv',pairs);tab(C/'niche_correspondence.tsv',corr)
for p in configs:
 name=p.name;cc=cfg[name];lab=labels[name];native=lab.candidate_label.to_numpy();k=cc['n_clusters'];pr=pd.read_csv(p/'profiles_R1_R2.tsv',sep='\t');sup=pd.read_csv(p/'cluster_support.tsv',sep='\t');sup.insert(0,'config_id',name);supports.append(sup);pr.insert(0,'config_id',name);allprofiles.append(pr);src=pd.read_csv(p/'profile_by_dataset_sample.tsv',sep='\t');src.insert(0,'config_id',name);allsource.append(src)
 cent=np.array([q[native==c].mean(0) for c in range(k)]);d=squareform(pdist(cent));np.fill_diagonal(d,np.inf)
 for c in range(k):
  vals=cent[c];glob=q.mean(0);enrich=np.log2((vals+1e-6)/(glob+1e-6));parts=[]
  for ri in range(2):
   off=ri*30;idx=[j for j in np.argsort(-enrich[off:off+30]) if vals[off+j]>=.03][:3];parts.append(', '.join(f'{TYPES[j]} {vals[off+j]:.1%} ({enrich[off+j]:+.2f} log2)' for j in idx))
  s=sup[sup.cluster==c].iloc[0];dictionary.append(dict(config_id=name,cluster=c,n=int(s.n),R1_relative_enrichment=parts[0],R2_relative_enrichment=parts[1],shared_across_datasets=bool(s.n_dataset>=3 and s.max_dataset_fraction<.8),n_dataset=int(s.n_dataset),n_sample=int(s.n_sample),n_source_known_patients=int(s.n_source_known_patients),n_analysis_units=int(s.n_analysis_units),dataset_dominated=bool(s.dataset_dominated),sample_dominated=bool(s.sample_dominated),small=bool(s.small),redundant_profile=bool(d[c].min()<.03),nearest_profile_distance=float(d[c].min()),interpretation='relative environmental enrichment; not pure cell type or validated anatomical interface'))
 for masklabel,mask in [('all',np.ones(len(m),bool)),('complete_rings',m.complete_rings.to_numpy())]:
  for c in range(k):
   ix=(native==c)&mask;pp=q[ix].mean(0) if ix.any() else np.full(60,np.nan);complete.append(dict(config_id=name,subset=masklabel,cluster=c,n=int(ix.sum()),denominator=int(mask.sum()),fraction=float(ix.sum()/mask.sum()),accepted_fraction=float((lab.accepted_label.to_numpy()[mask]>=0).mean()),profile_L2_from_all=float(np.linalg.norm(pp-cent[c]))))
 # Major proportions are sums of frozen subtype proportions, not marginal major q05.
 mp=[]
 for c in range(k):
  for ri in range(2):
   for maj in majors:mp.append(dict(cluster=c,ring=['R1','R2'][ri],major=maj,proportion=float(cent[c,ri*30:(ri+1)*30][major.to_numpy()==maj].sum())))
 tab(p/'profiles_major.tsv',mp)
 fig,ax=plt.subplots(2,1,figsize=(14,max(6,k*.5)))
 for ri,a in enumerate(ax):
  im=a.imshow(cent[:,ri*30:(ri+1)*30],aspect='auto',cmap='viridis');a.set_xticks(range(30),TYPES,rotation=70,ha='right',fontsize=6);a.set_yticks(range(k),[str(i) for i in range(k)]);a.set_title(name+' original '+['R1','R2'][ri]+' subtype composition');fig.colorbar(im,ax=a,shrink=.6)
 fig.tight_layout();save(fig,F/(name+'__profiles'))
 fig,axs=plt.subplots(2,2,figsize=(14,9));counts=np.bincount(native,minlength=k);axs[0,0].bar(range(k),counts,color=palette(k));axs[0,0].set_title('Candidate sizes');axs[0,1].bar(range(k),[float((lab.accepted_label[native==c]>=0).mean()) for c in range(k)],color=palette(k));axs[0,1].set_ylim(0,1);axs[0,1].set_title('Accepted / candidate per niche')
 for ax,key in zip(axs[1],['dataset_id','sample_uid']):
  t=pd.crosstab(native,m[key],normalize='index');bottom=np.zeros(k)
  for j,col in enumerate(t.columns):ax.bar(range(k),t[col],bottom=bottom,color=palette(len(t.columns))[j],label=col);bottom+=t[col].to_numpy()
  ax.set_title(key+' contribution');ax.legend(fontsize=4,ncol=2 if key=='sample_uid' else 1,loc='upper left',bbox_to_anchor=(1,1))
 fig.suptitle(name);fig.tight_layout();save(fig,F/(name+'__size_source_coverage'))
 fig,ax=plt.subplots(figsize=(8,5));dendrogram(linkage(cent,method='average'),labels=[str(c) for c in range(k)],ax=ax);ax.set_title(name+' original q1/q2 centroid hierarchy');save(fig,F/(name+'__hierarchy'))
 selections.append(f"{name}: K={k}, parameter={cc['k_or_resolution']}, embedding silhouette={cc['silhouette_embedding']:.4f}, median conditional stability ARI={cc['stability_ARI_median']:.4f}, accepted={cc['accepted_fraction']:.2%}, source-dominated niches (dataset/sample)={cc['dataset_dominance_summary']}/{cc['sample_dominance_summary']}. Rule: {cc['selection_rule']}. Negative silhouette or source domination limits recommendation; no identity inferred from cluster IDs.")

tab(C/'niche_dictionary.tsv',dictionary);tab(C/'dataset_sample_support.tsv',pd.concat(supports));tab(C/'profiles_R1_R2_all.tsv',pd.concat(allprofiles));tab(C/'profile_by_dataset_sample_all.tsv',pd.concat(allsource));tab(C/'complete_ring_sensitivity.tsv',complete);(C/'parameter_selection_rationale.md').write_text('\n\n'.join(selections)+'\n')
reps={p.resolve().parent for p in configs}|{p for p in (OUT/'02_R12_pilot').iterdir() if p.is_dir() and (p/'COMPLETE.json').exists()}
for od in sorted(reps):
 rid=od.name;integ=pd.read_csv(od/'integration_metrics.tsv',sep='\t');integ.insert(0,'representation_id',rid);integs.append(integ)
 if (od/'sameK_diagnostics.tsv').exists():samek.append(pd.read_csv(od/'sameK_diagnostics.tsv',sep='\t'))
 import fcntl
 ul=open(od/'umap.lock','w');fcntl.flock(ul,fcntl.LOCK_EX)
 if not (od/'umap_reference.npy').exists():
  import umap
  log('UMAP '+rid);z=np.load(od/'embedding.npy',mmap_mode='r');u=umap.UMAP(n_neighbors=30,min_dist=.3,random_state=123,n_jobs=1).fit_transform(z[ref]);np.save(od/'umap_reference.npy',u)
 fcntl.flock(ul,fcntl.LOCK_UN);ul.close()
 u=np.load(od/'umap_reference.npy');fig,axs=plt.subplots(2,2,figsize=(16,12))
 for ax,key in zip(axs[0],['dataset_id','sample_uid']):
  cats=pd.Categorical(m.iloc[ref][key]);discrete(ax,u[:,0],u[:,1],cats.codes,rid+' '+key,len(cats.categories),size=1,legend=False);handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor=palette(len(cats.categories))[i],markersize=3,label=str(c)) for i,c in enumerate(cats.categories)];ax.legend(handles=handles,fontsize=4,ncol=3 if key=='sample_uid' else 1,loc='upper right')
 for ax,alg in zip(axs[1],['KM','LEIDEN']):
  lab=np.load(od/alg/'reference_native_labels.npy');discrete(ax,u[:,0],u[:,1],lab,rid+' '+alg,len(set(lab)),1)
 fig.suptitle('Reference-only UMAP; shared axes within representation; axes NOT matched across representations');fig.tight_layout();save(fig,F/(rid+'__UMAP'))
tab(C/'integration_metrics.tsv',pd.concat(integs));tab(C/'sameK_diagnostics.tsv',pd.concat(samek))
for ring in ['R1','R2','R12']:
 models,ls,rr=joblib.load(OLD/'04_clustering'/ring/'candidates.joblib');new=np.load(OUT/'03_main'/f'NONE_{ring}_KM'/'reference_native_labels.npy')
 for k in [2,6]:
  old=ls[f'KM_{k}'];none_rep=(OUT/'03_main'/f'NONE_{ring}_KM').resolve().parent;_,nl,_=joblib.load(none_rep/'candidates.joblib');same=nl[f'KM_{k}'];legacy.append(dict(ring=ring,legacy_K=k,comparison='legacy full-Hellinger vs NONE PCA same-K',ARI=adjusted_rand_score(old,same),NMI=normalized_mutual_info_score(old,same)));legacy.append(dict(ring=ring,legacy_K=k,comparison='legacy vs NONE chosen working K',working_K=len(set(new)),ARI=adjusted_rand_score(old,new),NMI=normalized_mutual_info_score(old,new)))
tab(C/'legacy_vs_NONE_PCA.tsv',legacy)
# Every target section, including zero-anchor sections; physical image coordinates preserve orientation.
master=pd.read_parquet(OLD/'03_features/anchor_master.parquet',columns=['dataset_id','section_id','array_row','array_col','spot_uid']);coords=pd.read_csv(OLD/'00_preflight/full_image_tissue_coordinates.csv.gz',usecols=['section','barcode','imagerow','imagecol']);mm=m.merge(coords.rename(columns={'section':'section_id'}),on=['section_id','barcode'],how='left',validate='one_to_one');assert len(mm)==len(m)
for r in cov.itertuples():
 s=r.section_id;ds=r.dataset_id;bd=C/'by_dataset'/ds;bd.mkdir(parents=True,exist_ok=True);ix=np.flatnonzero(m.section_id==s);bg=coords[coords.section==s]
 for ring in ['R1','R2','R12']:
  dest=bd/(s+'__'+ring+'__six_modes')
  if Path(str(dest)+'.png').exists():continue
  fig,axs=plt.subplots(2,3,figsize=(16,10))
  for i,alg in enumerate(['KM','LEIDEN']):
   for j,mode in enumerate(['NONE','DATASET','SAMPLE']):
    name=f'{mode}_{ring}_{alg}';ax=axs[i,j];ax.scatter(bg.imagecol,bg.imagerow,s=1,c='#eeeeee',rasterized=True);discrete(ax,mm.iloc[ix].imagecol,mm.iloc[ix].imagerow,labels[name].accepted_label.to_numpy()[ix],f'{name} | accepted {int((labels[name].accepted_label.to_numpy()[ix]>=0).sum())}/{len(ix)}',cfg[name]['n_clusters'],size=3);ax.set_aspect('equal');ax.invert_yaxis();ax.set_xlabel('image column');ax.set_ylabel('image row')
  fig.suptitle(f'{ds} / {s}: accepted labels; pale gray = physical tissue not eligible; gray = UNASSIGNED');fig.tight_layout();save(fig,dest,pdf=False)
 # Separate candidate and accepted coverage summary on fixed same cohort for all six R12 labels.
 fig,axs=plt.subplots(2,3,figsize=(16,10))
 for ax,name in zip(axs.ravel(),[f'{mode}_R12_{alg}' for alg in ['KM','LEIDEN'] for mode in ['NONE','DATASET','SAMPLE']]):
  ax.scatter(bg.imagecol,bg.imagerow,s=1,c='#eeeeee',rasterized=True);discrete(ax,mm.iloc[ix].imagecol,mm.iloc[ix].imagerow,labels[name].candidate_label.to_numpy()[ix],name+' candidates',cfg[name]['n_clusters'],3);ax.set_aspect('equal');ax.invert_yaxis()
 fig.suptitle(s+' candidate labels (before confidence rejection)');fig.tight_layout();save(fig,bd/(s+'__R12__candidates'),pdf=False)
# Trade-off rather than mixing-only ranking.
i=pd.concat(integs);fig,axs=plt.subplots(1,2,figsize=(12,5))
for key,ax in zip(['dataset_id','sample_uid'],axs):
 for r in i[i.source_key==key].itertuples():ax.scatter(1-r.same_source_neighbor_fraction,r.raw_Hellinger_neighbor_mean);ax.annotate(r.representation_id,(1-r.same_source_neighbor_fraction,r.raw_Hellinger_neighbor_mean),fontsize=5)
 ax.set_xlabel('cross-source neighbor fraction');ax.set_ylabel('original R12 Hellinger distance of neighbors');ax.set_title(key+' mixing vs composition separation')
fig.tight_layout();save(fig,F/'mixing_vs_composition')
log('Comparisons and figures complete')
import subprocess
subprocess.run(['/usr/bin/python3.12',str(OUT/'code/additional_diagnostics.py')],check=True)
subprocess.run(['/usr/bin/python3.12',str(OUT/'code/plot_running.py')],check=True)
subprocess.run(['/usr/bin/python3.12',str(OUT/'code/enrichment_figures.py')],check=True)

# Final scope and semantic checks use only niche outputs.
checks=[]
for name,lab in labels.items():
 assert lab.spot_uid.tolist()==m.spot_uid.tolist() and len(lab)==129632
 assert (lab.reference_native_label.to_numpy()[ref]==np.load(OUT/'03_main'/name/'reference_native_labels.npy')).all()
 assert (lab.candidate_label.to_numpy()[ref]==lab.reference_native_label.to_numpy()[ref]).all()
 accepted=~lab.low_support & ~lab.distance_outlier;assert ((lab.accepted_label>=0)==accepted).all();assert (lab.loc[accepted,'accepted_label']==lab.loc[accepted,'candidate_label']).all()
 assert np.isclose(cfg[name]['accepted_fraction'],accepted.mean());assert cfg[name]['n_reference']==20000
 assert len(pd.read_csv(OUT/'03_main'/name/'stability.tsv',sep='\t'))==5
 checks.append(dict(config_id=name,n=len(lab),native_reference_preserved=True,accepted_semantics_valid=True,niche_columns_only=not any('mp10' in c.lower() for c in lab.columns),status='PASS'))
tab(C/'validation.tsv',checks)
summary=pd.DataFrame(list(cfg.values()));summary['recommendation_status']=np.where((summary.stability_ARI_median>=.8)&(summary.silhouette_embedding>0)&(summary.silhouette_raw_Hellinger>0)&(summary.convergence_status!='iteration_limit'),'niche_working_candidate','exploratory_with_limitations');tab(C/'eighteen_scheme_comparison.tsv',summary)
R=OUT/'reports';R.mkdir(exist_ok=True)
methods='''# Niche methods and scope

The current run stops at niche. No MP10 score columns, gene programs, source-validation result values or downstream unblinding scripts were used. Metadata are read using an explicit whitelist. Existing Seurat objects were opened solely to extract image slots, coordinates and scale factors; molecular assays and metadata scores were not accessed by the image extraction code. Region review does not block niche completion.

Frozen accepted input: cell2location v2_30type subtype q05, 30 factors in frozen order. Center gate is posterior mean Epithelial / all30 >=0.25, explicitly approved by the user. Exact physical hex ring1 and ring2, minimum 4/8 valid observations, full rings 6/12; missing-composition physical tissue nodes retained in topology. Fixed 129632 common anchors, 108345 full-ring anchors and 20000 reference anchors. Target 63 sections, 61 with common anchors/reference. sample_uid = dataset_id::section_id, never patient ID. Patient identity and analysis-unit evidence remain distinct.

Reused v2 Hellinger features and fitted PCA; source features were saved float32 in v2, converted to float64 before original PCA model transform. No further standardization or whitening. Each ring has its own minimal 95% PCA. NONE completely skips Harmony. R Harmony 2.0.5 matrix API fits all common anchors jointly using one batch covariate; this is not independent holdout representation learning. WEAK theta=1 lambda=5; STRONG theta=2 lambda=1; sigma=.1, internal nclust=100, ncores=1, seed123, max_iter20, early_stop=TRUE, max.iter.cluster20, epsilon.cluster1e-3, epsilon.harmony1e-4. Only true iteration-limit branches receive one finite max50 rerun. Version 2.0.5 corrected coordinates use getZcorr(), not the historical Z_corr field. Objective traces and manifests accompany each representation. Internal nclust is not final niche K.

R12 five-representation pilot precedes global strength selection separately for DATASET and SAMPLE. Strong correction requires an absolute same-source-neighbor reduction >.02, <=5% raw neighbor Hellinger distance increase, <=.02 loss of both algorithms' stability and coverage, and <=5% increase in raw composition dispersion; otherwise weak is retained. These conservative operational criteria are project choices, not literature-derived universal cutoffs. Identical chosen strengths extend to R1/R2 and both algorithms. All pilot branches retained.

K-means works on the respective PCA/Harmony coordinates, n_init50, max_iter300, Lloyd, seed123; finite1000 continuation only if needed. Working candidates K6/8/10/12; K2/4 are coarse diagnostics. Conservative K6 rule and per-configuration measured evidence are in parameter_selection_rationale.md. Leiden uses exact Euclidean30NN including self, binary intersection/(60-intersection) shared-neighbor Jaccard, edges >1/15, undirected/no selfloops; RBConfigurationVertexPartition, weighted, n_iterations10. Grid .005/.01/.02/.05/.1; finite extensions only when all candidates fall below6 or above12. Candidates with positive embedding and original-Hellinger silhouettes are preferred whenever available before considering a moderate K range. A working K outside6–12 is reported when low-resolution moderate solutions fail this separation check. Initial NONE low-resolution maps remain archived as diagnostics. Source-dominated results remain exploratory. Graph components/isolates and degree summaries are recorded. Neither desired biological naming nor UMAP determines parameters.

Finalist stability seeds124–128 reuses exactly the same v2 80% reference UID subsets, refitting clustering only; PCA/Harmony remain fixed. ARI, matched overlap, K and original-profile displacement are reported. This is conditional stability, not patient validation. R12 additionally removes 20% of section-specific 8x8 array blocks. Sparse block mixing keeps the first fixed reference anchor in each8x8block. These diagnostics reduce local spatial redundancy but do not create independent biological replicates.

Reference native labels are preserved. Candidate = nearest centroid for K-means or 15NN vote for Leiden, with reference native override. Common transfer = 15NN for both. Diagnostics exclude each reference's self-neighbor. Ties use highest count, smallest mean distance among label voters, then numeric ID. Support for K-means is the vote share of its centroid candidate, not another class's majority vote. Accepted requires support>=.6 and nearest-reference distance<=reference leave-self-out Q99; otherwise -1 / UNASSIGNED. Full candidate, reference-native and shared-accepted comparisons use explicit matched denominators.

Interpretations always use original q1/q2 point compositions; major profiles are sums of subtype proportions, not marginal major q05. Enrichment is log2((niche mean+1e-6)/(global common-anchor mean+1e-6)). Dataset>80%, sample>50%, size<.5% and centroid distance<.03 are descriptive flags. No sample is removed for dominance. Niche IDs are config-specific; matching tables provide explicit correspondence. Composition strata for conditional mixing are frozen dominant mean(q1,q2) subtype, not pathological labels. Entropies for dataset and sample must not be compared directly across category counts. Random pairing is descriptive only.

UMAP uses the fixed reference only, n_neighbors30,min_dist.3,seed123; common coordinates within a representation, different axes across representations. All section maps preserve original image-row/column orientation, equal aspect and physical-tissue background; rejected anchors remain gray. Embedded H&E is the actual existing image at stored resolution. Automatic boundary is the perimeter of ROI member spots, not a claimed lesion boundary. No final pathology or interface is inferred.

References: [Harmony original method](https://www.nature.com/articles/s41592-019-0619-0), [official implementation](https://github.com/immunogenomics/harmony/blob/master/R/ui.R), [BANKSY Harmony vignette](https://prabhakarlab.github.io/Banksy/articles/batch-correction.html). The BANKSY vignette supports integration of spatial embeddings, not this project's gate/ring/strength parameters. [Integration benchmark](https://www.nature.com/articles/s41592-021-01336-8) and [ENKTL paper](https://www.nature.com/articles/s41375-026-03118-7) were supplied in the request but full pages were not retrievable in this session; no unverified specific clustering settings are attributed to them.
'''
(R/'METHODS.md').write_text(methods)
lines=['# Niche v3 实际完成情况',f'18/18 主配置已完成：NONE、DATASET、SAMPLE × R1、R2、R12 × KM、Leiden。另保留全部 R12 五表示 pilot。共同 anchors=129632，完整双环=108345，reference=20000；63 张目标切片中 61 张有 anchors。','## 主要结果','|配置|群数|稳定性 ARI 中位数|标签接受率|状态|','|---|---:|---:|---:|---|']
for r in summary.itertuples():lines.append(f'|{r.config_id}|{r.n_clusters}|{r.stability_ARI_median:.3f}|{r.accepted_fraction:.1%}|{r.recommendation_status}|')
lines+=['## 三种整合模式的收益与代价',(OUT/'02_R12_pilot/pilot_selection.md').read_text(),'|R12 表示|同 dataset 邻居|同 sample 邻居|原始 Hellinger 邻居距离|','|---|---:|---:|---:|']
for rid,g in i[i.representation_id.str.endswith('R12')].groupby('representation_id'):
 gg=g.set_index('source_key');lines.append(f"|{rid}|{gg.loc['dataset_id','same_source_neighbor_fraction']:.1%}|{gg.loc['sample_uid','same_source_neighbor_fraction']:.1%}|{gg.raw_Hellinger_neighbor_mean.iloc[0]:.4f}|")
lines+=['来源混合不能单独证明批次效应被正确消除；patient/clinical 差异也可能是真实组成差异。请结合原始组成距离、来源支持、稳定性和拒绝率评估，暂不宣称唯一最佳方案。','## 工作粒度和组成解释','v2 的 K2 是粗粒度结构，原 Leiden 38–45 群偏细。本轮 K2/4 仅作粗诊断，工作 K-means 从 K6 起评估；Leiden 从更低分辨率评估。负轮廓系数或强来源依赖的配置被保留为探索结果，不命名为已验证的病理 niche。详细逐配置选择依据见比较目录。','[niche_dictionary.tsv](../04_comparison/niche_dictionary.tsv) 给出每个 niche 的 R1/R2 相对富集、来源共享、真实患者证据与来源主导标记。SPP1 TAM、CAF、T/NK 等只按原始组成描述，不推断纯细胞群或机制。','## 结果入口','- [18方案比较](../04_comparison/eighteen_scheme_comparison.tsv)\n- [图形目录](../04_comparison/figures/)\n- [逐数据集空间图](../04_comparison/by_dataset/)\n- [ROI 操作说明](../05_user_review/WHAT_USER_NEEDS_TO_DO.md)\n- [轻量回传包](review.zip)','## 你现在需要做什么','先查看 ROI 操作说明中的五个重点区域，填写 ROI_review_to_fill.tsv 的 reviewed_label、review_comment、reviewer、review_date。允许 uncertain 和 mixed_needs_split；混合区域可返回 polygon/spot成员表，不强迫整个区域单一标签。人工判断不影响本轮 niche 完成。','## 停止范围与续跑','本轮止于 niche。MP10 分析、解盲、来源验证全部 DEFERRED_BY_USER，留待后续授权。未运行旧 v2 resume.sh。','```bash\nbash "/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_harmony_v3/code/resume_niche_only.sh"\n```','续跑命令复用已完成 checkpoint，仅生成/补全本轮 niche 产物，不触发 MP10。']
# Concrete, descriptive R12 interpretation; no biological ground-truth labels.
explain=['## 建议优先复核的工作结果','K-means 可先以 DATASET_WEAK 的 R12 K6 作为较保守的整合工作视图，同时保留 NONE K6 作为未校正基线；SAMPLE_WEAK 的 R12 Leiden K9 可作为细分环境的补充视图。它们是不同用途的工作候选，不宣称唯一最佳；具体稳定性、覆盖率和组成代价见上表。','NONE R2/R12 Leiden 的低分辨率中等群数解在原始组成上分离不足，因此工作结果采用有限网格中的25/27群，R1保留有正轮廓证据的6群。它们仍有来源集中限制；v2的38–45群不被直接当作生物学真值。','|SAMPLE_WEAK R12 Leiden ID|占共同anchors|R1 相对富集|R2 相对富集|支持数据集|最大单dataset占比|','|---|---:|---|---|---:|---:|']
nd=pd.DataFrame(dictionary)
for r in nd[nd.config_id=='SAMPLE_R12_LEIDEN'].itertuples():
 sp=pd.concat(supports);ss=sp[(sp.config_id==r.config_id)&(sp.cluster==r.cluster)].iloc[0];explain.append(f'|{r.cluster}|{r.n/len(m):.1%}|{r.R1_relative_enrichment}|{r.R2_relative_enrichment}|{r.n_dataset}|{ss.max_dataset_fraction:.1%}|')
explain+=['这些是原始 q05 双环组成的相对富集描述。多数群仍有大量 Epithelial 成分；不能把8% SPP1 TAM 或8% Matrix_CAF称为纯髓系/纯CAF niche。R1/R2总体相近，部分群如中性粒细胞富集群在R2更高，尚不能据此声称组织界面或作用方向。','来源主导群详见 niche_dictionary.tsv；明确共享需结合多来源支持而非要求每群包含全部7个dataset。已知患者数与匿名/假设分析单位分别记录。']
lines.extend(explain)
(R/'RUN_SUMMARY.md').write_text('\n\n'.join(lines)+'\n')
(R/'SCOPE_AND_ACCEPTANCE.md').write_text('# Scope and acceptance\n\n18 main configurations validated. 5 R12 pilot representations retained. Scope: niche only. Region review independent. MP10 analysis, unblinding and source validation DEFERRED_BY_USER. See validation.tsv for label semantics, reference and denominator checks. No unique best integration declared without pathology. NONE Leiden low-resolution choices were revised after negative silhouette and source-concentration evidence; archived diagnostics retained and fixed-representation stability recomputed for the selected finite-grid alternatives. Initial Harmony export compatibility error was corrected using installed 2.0.5 getZcorr(); affected runs repeated with identical parameters. All final results derive from successful corrected runs.\n')
review=['# Review index','[运行总结](RUN_SUMMARY.md)','[方法](METHODS.md)','[18方案比较](../04_comparison/eighteen_scheme_comparison.tsv)','[用户 ROI 操作说明](../05_user_review/WHAT_USER_NEEDS_TO_DO.md)','## Figures']+[f'[{p.stem}](../04_comparison/figures/{p.name})' for p in sorted(F.glob('*.png'))]
review+=['## Section maps']+[f'[{p.parent.name}/{p.stem}](../04_comparison/by_dataset/{p.parent.name}/{p.name})' for p in sorted((C/'by_dataset').glob('*/*.png'))];(R/'REVIEW_INDEX.md').write_text('\n\n'.join(review)+'\n')
manifest=dict(status='NICHE_COMPLETE',stop_after='niche',mp10_authorized_this_run=False,mp10_values_loaded=False,mp10_analysis_executed=False,source_validation_status='DEFERRED_BY_USER',region_review_required_for_niche_completion=False,region_final_frozen=False,n_main=18,n_R12_pilot_representations=5,selected_harmony_parameters=json.loads((OUT/'02_R12_pilot/selected_harmony_parameters.json').read_text()),configs=summary.to_dict('records'),input_hashes=pd.read_csv(OUT/'00_preflight/input_hashes.tsv',sep='\t').to_dict('records'),code_hashes={str(p):sha(p) for p in (OUT/'code').glob('*') if p.is_file()},label_hashes={name:sha(OUT/'03_main'/name/'labels.parquet') for name in labels},PCA_hashes={p.name:sha(p) for p in (OUT/'01_features').glob('*PCA.joblib')},reference_hash=sha(OUT/'01_features/reference_anchor_ids.tsv'),batch_mapping_hash=sha(OUT/'00_preflight/batch_key_mapping.tsv'),limitations=['conditional fixed-representation stability, not patient validation','Harmony fits all anchors jointly, not independent holdout representation','stored v2 features float32, projected float64','source mixing may erase real sample-level biology','embedded H&E resolution; final pathology pending','Leiden source-specific/negative-silhouette branches exploratory'],deviations=['Harmony 2.0.5 export API adapted to getZcorr(); affected initial exports rerun','ENKTL and integration-benchmark full pages unavailable; no unverified parameters adopted','NONE Leiden low-resolution working choice revised using positive composition separation; initial maps retained as diagnostics, global strength choice unaffected'])
js(OUT/'frozen_niche_manifest.json',manifest)
# Lightweight review package: no raw counts, labels, corrected embeddings or fitted model objects.
priority=pd.read_csv(OUT/'05_user_review/ROI_priority_review.tsv',sep='\t').head(5);paths=[]
for sub in ['reports','00_preflight','04_comparison','05_user_review']:
 paths.extend(p for p in (OUT/sub).glob('*') if p.is_file() and p.suffix in ['.md','.tsv','.json'])
paths+=list(F.glob('*.png'))+[OUT/'frozen_niche_manifest.json']
paths += [Path(p) for p in priority.evidence_figure if Path(p).is_file() and str(p).startswith(str(OUT))]
paths += [p for rep in (OUT/'02_R12_pilot').iterdir() if rep.is_dir() for p in [rep/'objective.tsv',rep/'harmony_manifest.json'] if p.exists()]
paths += [OUT/'02_R12_pilot'/p for p in ['pilot_comparison.tsv','pilot_selection.md','selected_harmony_parameters.json']]
with zipfile.ZipFile(R/'review.zip','w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(set(paths)):z.write(p,str(p.relative_to(OUT)))
for ds in sorted((C/'by_dataset').iterdir()):
 with zipfile.ZipFile(R/(ds.name+'__spatial_figures.zip'),'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for p in ds.glob('*.png'):z.write(p,str(p.relative_to(OUT)))
with zipfile.ZipFile(R/'ROI_evidence_all.zip','w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in (OUT/'05_user_review/evidence').glob('*.png'):z.write(p,str(p.relative_to(OUT)))
log('NICHE_COMPLETE; scope stops before MP10')
