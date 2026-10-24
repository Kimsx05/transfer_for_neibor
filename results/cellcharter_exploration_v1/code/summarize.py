from common import *
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex,LogNorm
plt.rcParams.update({'figure.facecolor':'none','savefig.transparent':True,'font.size':8})
o=pd.read_parquet(ROOT/'cellcharter_assignments.parquet');original,z,xy=load();assert o.cell_id.equals(original.cell_id);K=json.loads((ROOT/'metadata/K_selection.json').read_text())['K'];palette={f'niche{k}':to_hex(plt.get_cmap('tab20')(k)) for k in range(K)};js(ROOT/'metadata/niche_palette.json',palette);versions=['PRIMARY','EXPANDED_LOWINFO','NO_CURRENT_CONFLICT'];summ=[];comp=[];fine=[];stability=[]
for version in versions:
 for (sec,niche),v in o[o[version+'_included']].groupby(['section_id',version+'_niche'],observed=True):
  total=int((o.section_id.eq(sec)&o[version+'_included']).sum());known=~v.interpretation_major.isin(['Unknown','Disputed']);row=dict(version=version,section_id=sec,niche=niche,n=len(v),section_fraction=len(v)/total,interpretable_n=int(known.sum()),interpretable_coverage=float(known.mean()),median_genes=float(v.n_genes.median()),median_counts=float(v.n_counts.median()),Unknown_fraction=float(v.interpretation_major.eq('Unknown').mean()),Disputed_fraction=float(v.interpretation_major.eq('Disputed').mean()))
  for c in ['low_information','review_strong_competition','working_secondary_signal','mixed_or_ambiguous','morphology_identity_uncertain','artifact_review_flag']:row[c+'_fraction']=float(v[c].mean())
  row['isolated_n']=int(v[version+'_context'].eq('isolated_self_expression_only').sum());summ.append(row)
  for lab,n in v.interpretation_major.value_counts().items():comp.append(dict(version=version,section_id=sec,niche=niche,major=lab,n=n,all_denominator=len(v),fraction_all=n/len(v),interpretable_denominator=int(known.sum()),fraction_interpretable=n/known.sum() if lab not in ['Unknown','Disputed'] and known.sum() else np.nan))
  for lab,n in v.interpretation_display.value_counts().items():fine.append(dict(version=version,section_id=sec,niche=niche,display=lab,n=n,fraction_all=n/len(v),meaning='Inherited subtype if available, otherwise parent; Unknown/Disputed explicit'))
for version in versions[1:]:
 common=o.PRIMARY_included&o[version+'_included']
 for sec in ['Cervilla','D1_Johnson']:
  s=o[common&o.section_id.eq(sec)]
  for niche in sorted(s.PRIMARY_niche.unique()):
   a=s.PRIMARY_niche.eq(niche);b=s[version+'_niche'].eq(niche);stability.append(dict(version=version,section_id=sec,niche=niche,primary_common_n=int(a.sum()),sensitivity_common_n=int(b.sum()),retained_same_mapped_n=int((a&b).sum()),retention=float((a&b).sum()/a.sum()),Jaccard=float((a&b).sum()/(a|b).sum())))
summary=pd.DataFrame(summ);composition=pd.DataFrame(comp);summary.to_csv(ROOT/'metadata/niche_summary.tsv',sep='\t',index=False);composition.to_csv(ROOT/'metadata/major_composition.tsv',sep='\t',index=False);pd.DataFrame(fine).to_csv(ROOT/'metadata/subtype_parent_composition.tsv',sep='\t',index=False);pd.DataFrame(stability).to_csv(ROOT/'metadata/per_niche_sensitivity.tsv',sep='\t',index=False)
# Descriptive full-feature measured markers. Only K x genes summaries are dense.
markers=[];panels=json.loads((BASE/'joint_umap_major_review/metadata/review_panels.json').read_text())['panels'];known=set(g for p in panels.values() for gs in p.values() for g in gs)
for sec in ['Cervilla','D1_Johnson']:
 a=ad.read_h5ad(BASE/'cell_annotation_v2'/sec/'cells.h5ad');idx=pd.Index(o.cell_id).get_indexer(a.obs_names);sub=o.iloc[idx];symbols=np.where(a.var.reference_symbol.astype(str).to_numpy()!='',a.var.reference_symbol.astype(str),a.var.gene_symbol.astype(str));genes=pd.Index(pd.unique(symbols));mapx=sparse.csr_matrix((np.ones(a.n_vars),(np.arange(a.n_vars),genes.get_indexer(symbols))),shape=(a.n_vars,len(genes)));x=(a.X@mapx).tocsr();sel=np.flatnonzero(sub.PRIMARY_included);x=x[sel];labels=sub.iloc[sel].PRIMARY_niche.to_numpy();norm=sparse.diags(1e4/np.maximum(np.asarray(x.sum(1)).ravel(),1)).dot(x).tocsr();norm.data=np.log1p(norm.data);total=np.asarray(norm.sum(0)).ravel();detall=np.asarray((x>0).sum(0)).ravel()
 for niche in sorted(pd.unique(labels)):
  ii=np.flatnonzero(labels==niche);n=len(ii);mean=np.asarray(norm[ii].mean(0)).ravel();det=np.asarray((x[ii]>0).mean(0)).ravel();other=(total-mean*n)/max(len(sel)-n,1);od=(detall-det*n)/max(len(sel)-n,1);diff=(mean-other)/np.log(2);eligible=(det>=.1)&(diff>0)&~genes.str.match('^(MT-|RPL|RPS)');top=np.flatnonzero(eligible)[np.argsort(-diff[eligible])[:15]];wanted=np.union1d(top,np.flatnonzero(genes.isin(known)))
  for j in wanted:markers.append(dict(section_id=sec,niche=niche,n=n,gene=genes[j],detection=float(det[j]),mean_log1p10k=float(mean[j]),difference_log2=float(diff[j]),other_detection=float(od[j]),representative=bool(j in top)))
 del a,x,norm
m=pd.DataFrame(markers);m.to_parquet(ROOT/'metadata/measured_marker_evidence.parquet',compression='zstd');m[m.representative].to_csv(ROOT/'metadata/representative_markers.tsv',sep='\t',index=False)
hops=pd.read_parquet(ROOT/'metadata/hop_sample.parquet');hs=hops.groupby(['version','section_id','hop']).agg(sampled_cells=('cell_id','size'),median_neighbors=('neighbor_count','median'),p95_neighbors=('neighbor_count',lambda x:x.quantile(.95)),median_max_um=('max_distance_um','median'),p95_max_um=('max_distance_um',lambda x:x.quantile(.95))).reset_index();hs.to_csv(ROOT/'metadata/hop_scale_summary.tsv',sep='\t',index=False)
def save(fig,name):fig.savefig(ROOT/'figures'/name,dpi=180,bbox_inches='tight');plt.close(fig)
def spatial(ax,ix,col,title):
 for lab in sorted(o.iloc[ix][col].dropna().unique()):
  jj=ix[o.iloc[ix][col].eq(lab).to_numpy()];ax.scatter(xy[jj,0],xy[jj,1],s=.35,c=palette[lab],linewidths=0,label=f'{lab}: {len(jj):,}',rasterized=True)
 ax.set_aspect('equal');ax.invert_yaxis();ax.set_title(title);ax.set_xlabel('µm');ax.set_ylabel('µm');ax.legend(markerscale=5,fontsize=6,bbox_to_anchor=(1,1),loc='upper left',frameon=False)
for version in versions:
 fig,axs=plt.subplots(1,2,figsize=(17,7))
 for ax,sec in zip(axs,['Cervilla','D1_Johnson']):
  ix=np.flatnonzero(o[version+'_included']&o.section_id.eq(sec));spatial(ax,ix,version+'_niche',f'{version} | {sec} n={len(ix):,}')
 fig.subplots_adjust(wspace=.6);save(fig,version+'_two_sections.png')
for sec in ['Cervilla','D1_Johnson']:
 ix=np.flatnonzero(o.PRIMARY_included&o.section_id.eq(sec));fig,ax=plt.subplots(figsize=(9,7));spatial(ax,ix,'PRIMARY_niche',sec+' PRIMARY');save(fig,sec+'_PRIMARY.png');fig,axs=plt.subplots(2,2,figsize=(12,11))
 for ax,col in zip(axs.ravel(),['interpretation_major','n_genes','n_counts','review_strong_competition']):
  if col=='interpretation_major':vals=np.select([o.iloc[ix][col].eq('Disputed'),o.iloc[ix][col].eq('Unknown')],[2,1],default=0);im=ax.scatter(*xy[ix].T,c=vals,s=.3,vmin=0,vmax=2,cmap='coolwarm',linewidths=0);fig.colorbar(im,ax=ax,ticks=[0,1,2],label='0 Known / 1 Unknown / 2 Disputed')
  elif col=='review_strong_competition':im=ax.scatter(*xy[ix].T,c=o.iloc[ix][col].astype(int),s=.3,cmap='coolwarm',vmin=0,vmax=1,linewidths=0);fig.colorbar(im,ax=ax,ticks=[0,1])
  else:im=ax.scatter(*xy[ix].T,c=o.iloc[ix][col],s=.3,cmap='viridis',norm=LogNorm(),linewidths=0);fig.colorbar(im,ax=ax,label='log color')
  ax.set_title(sec+' '+col);ax.set_aspect('equal');ax.invert_yaxis();ax.set_xlabel('µm');ax.set_ylabel('µm')
 save(fig,sec+'_quality_diagnostics.png')
 # Complete major composition, including unknown/disputed in the denominator.
 c=composition[(composition.version=='PRIMARY')&(composition.section_id==sec)];p=c.pivot(index='niche',columns='major',values='fraction_all').fillna(0);ax=p.plot.bar(stacked=True,figsize=(11,5),colormap='tab20');ax.set_ylabel('Fraction of all niche cells');ax.legend(bbox_to_anchor=(1,1));save(ax.figure,sec+'_major_composition.png')
 mm=m[m.section_id.eq(sec)];gs=list(dict.fromkeys(mm[mm.representative].sort_values('difference_log2',ascending=False).groupby('niche').head(4).gene));gs=[g for g in gs if mm[mm.gene.eq(g)].niche.nunique()==mm.niche.nunique()][:48];v=mm[mm.gene.isin(gs)].pivot(index='niche',columns='gene',values='mean_log1p10k').reindex(columns=gs).fillna(0);fig,ax=plt.subplots(figsize=(max(10,.3*len(gs)),max(3,.4*K)));im=ax.imshow(v,aspect='auto',cmap='viridis');ax.set_xticks(range(len(gs)),gs,rotation=90);ax.set_yticks(range(len(v)),v.index);ax.set_title(sec+' measured marker mean log1p10k');fig.colorbar(im,ax=ax);save(fig,sec+'_marker_heatmap.png')
curve=pd.read_csv(ROOT/'metadata/autok_stability.tsv',sep='\t');fig,ax=plt.subplots(figsize=(8,4))
for stage,v in curve.groupby('stage'):ax.errorbar(v.K,v.stability,yerr=v.sd,marker='o',label=stage)
ax.axvline(K,c='red',ls='--',label=f'selected K={K}');ax.set_xlabel('K');ax.set_ylabel('Official stability (Fowlkes–Mallows)');ax.legend();save(fig,'autok_stability.png')
print('SUMMARIES_FIGURES_COMPLETE',flush=True)
