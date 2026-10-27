from common import *
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import re
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
s=read('samples_v2.tsv');p=read('proportions_v2.tsv.gz').merge(s[['sample_uid','sample_MP10_group_allEpi','MP10_sample_allEpi','n_parent_B','T_category','tissue_status']],on='sample_uid');r=read('proportion_models_v2.tsv');colors=dict(read('cell_colors.tsv')[['cell_type','color']].values);colors.update(Epi_MP10_low='#4575b4',Epi_MP10_high='#d73027');pd.DataFrame({'cell_type':list(colors),'color':list(colors.values())}).to_csv(T/'cell_colors.tsv',sep='\t',index=False)
def safe(x):return re.sub(r'[^a-zA-Z0-9_-]+','_',str(x))
def figsave(fig,ds,name):
 dest=O/'figures'/('overall' if ds=='ALL' else 'by_dataset/'+ds);dest.mkdir(parents=True,exist_ok=True)
 for ext in ['png','pdf']:fig.savefig(dest/(name+'.'+ext),dpi=300,transparent=True,bbox_inches='tight')
 plt.close(fig)
# Old versus new eligibility, exact affected specimen list in source table.
v=read('old_vs_new_samples.tsv');fig,axs=plt.subplots(1,2,figsize=(12,4));ds=sorted(s.dataset.unique());xx=np.arange(len(ds));old=v.groupby('dataset').old_MP10_eligible.sum().reindex(ds);new=v.groupby('dataset').MP10_eligible.sum().reindex(ds);axs[0].bar(xx-.18,old,.36,label='v1 malignant-like');axs[0].bar(xx+.18,new,.36,label='v2 all epithelial');axs[0].set(xticks=xx,xticklabels=ds,ylabel='Eligible tumor specimens',title='Eligibility changes');axs[0].tick_params(axis='x',rotation=25);axs[0].legend()
for d,z in v[v.old_MP10_eligible&v.MP10_eligible].groupby('dataset'):axs[1].scatter(z.old_MP10_sample,z.MP10_sample_allEpi,label=d,marker='o',s=np.where(z.group_flipped,80,25),facecolors='none')
axs[1].set(xlabel='v1 mean malignant-like MP10',ylabel='v2 mean all-epi MP10',title='Larger open points: sample-group flips');axs[1].legend(fontsize=7);fig.tight_layout();figsave(fig,'ALL','old_vs_new_eligibility')
# Every original proportion point retained, including low B and zeros; MP10 only eligible new frozen grouping.
focus=[('B','major','TME'),('Plasma','major','TME')]+[(x,'subtype','parent') for x in ['GC B','Memory B','Naive B','IFN-stimulated B']]
for ds in ['ALL']+sorted(s.dataset.unique()):
 z=p if ds=='ALL' else p[p.dataset.eq(ds)];fig,axs=plt.subplots(2,6,figsize=(22,8))
 for j,(ct,lev,den) in enumerate(focus):
  a=z[z.subtype.eq(ct)&z.level.eq(lev)&z.denominator_type.eq(den)&z.sample_MP10_group_allEpi.notna()];rng=np.random.default_rng(123)
  for d,aa in a.groupby('dataset'):
   color=plt.get_cmap('tab10')(sorted(s.dataset.unique()).index(d));low=aa.n_parent_B.lt(30) if den=='parent' else np.zeros(len(aa),bool);pos=aa.sample_MP10_group_allEpi.map({'low':0,'high':1}).to_numpy()+rng.uniform(-.1,.1,len(aa))
   for jj,(_,row) in enumerate(aa.iterrows()):axs[0,j].scatter(pos[jj],row.proportion,facecolors='none' if np.asarray(low)[jj] else color,edgecolors=[color],s=30);axs[1,j].scatter(row.MP10_sample_allEpi,row.proportion,facecolors='none' if np.asarray(low)[jj] else color,edgecolors=[color],s=30)
  axs[0,j].set(title=f'{ct} / {den}',xticks=[0,1],xticklabels=['sample low','sample high'],ylabel='Raw proportion',ylim=(-.03,1.03) if den=='parent' else None);axs[1,j].set(xlabel='All-epi specimen MP10',ylabel='Raw proportion')
 fig.suptitle(ds+': specimen points; hollow B-subtype points = total B <30; no zero-cell samples removed');fig.tight_layout();figsave(fig,ds,'B_Plasma_MP10_points')
 # B composition with actual B counts, all specimens and raw zeros where parent exists.
 a=z[z.level.eq('subtype')&z.denominator_type.eq('parent')&z.major.eq('B')];ids=s.loc[s.dataset.eq(ds) if ds!='ALL' else np.ones(len(s),bool)].sort_values(['dataset','sample_uid']);mat=a.pivot(index='sample_uid',columns='subtype',values='proportion').reindex(ids.sample_uid);fig,(ax,bx)=plt.subplots(2,1,figsize=(max(10,len(ids)*.24),6),gridspec_kw={'height_ratios':[3,1]},sharex=True);bottom=np.zeros(len(ids))
 for ct in sorted(mat.columns):val=mat[ct].fillna(0).to_numpy();ax.bar(np.arange(len(ids)),val,bottom=bottom,color=colors[ct],label=ct,width=.85);bottom+=val
 for i,n in enumerate(ids.n_parent_B):
  if n==0:ax.text(i,.4,'NA',rotation=90,ha='center',fontsize=6)
 bx.scatter(np.arange(len(ids)),ids.n_parent_B,c='black',s=12);bx.axhline(30,color='red',lw=.6);bx.set_yscale('symlog',linthresh=1);bx.set(ylabel='Total B',xticks=np.arange(len(ids)),xticklabels=ids['orig.ident']);bx.tick_params(axis='x',rotation=90,labelsize=6);ax.set(ylabel='Raw B composition',ylim=(0,1),title=ds+'; no B denominator = NA');ax.legend(bbox_to_anchor=(1,1));fig.tight_layout();figsave(fig,ds,'B_composition_with_coverage')
 # B parent coverage/transform sensitivities and stage adjustments, all coefficients and CI.
 rr=r[r.scope.eq(ds)&r.major.eq('B')&r.denominator_type.eq('parent')&r.analysis_variant.isin(['primary_asin','parent30_asin','parent30_logit'])];cts=['GC B','Memory B','Naive B','IFN-stimulated B'];cs=['sample_MP10_group_allEpi','MP10_continuous','MP10_stage_complete_unadjusted','MP10_stage_adjusted'];fig,axs=plt.subplots(4,4,figsize=(16,11))
 for i,ct in enumerate(cts):
  for j,co in enumerate(cs):
   aa=rr[rr.subtype.eq(ct)&rr.comparison.eq(co)].set_index('analysis_variant');ax=axs[i,j]
   for k,var in enumerate(['primary_asin','parent30_asin','parent30_logit']):
    if var not in aa.index:continue
    b=aa.loc[var];val=b.coefficient
    if np.isfinite(val):ax.errorbar(val,k,xerr=[[val-b.CI_low],[b.CI_high-val]],fmt='o',capsize=2)
    else:ax.text(.02,k,'NE',transform=ax.get_yaxis_transform())
   ax.axvline(0,c='gray',lw=.5);ax.set(yticks=range(3),yticklabels=['primary asin','B>=30 asin','B>=30 logit'] if j==0 else [],ylim=(-.6,2.6),title=ct+'\n'+co.replace('MP10_',''))
 fig.suptitle(ds+': effect and 95% CI; logit and asin have DIFFERENT scales; continuous per 1 MP10');fig.tight_layout();figsave(fig,ds,'B_sensitivity_effects')
# Limited B marker check: specimen-equal mean and detection; missing vs zero distinguished.
b=read('B_marker_sample_support.tsv');markers=read('B_marker_gene_coverage.tsv').gene.tolist();a=b.groupby(['dataset','subtype','gene']).agg(mean=('mean','mean'),detection=('detection','mean'),n_specimens=('sample_uid','nunique')).reset_index();save(a,'B_marker_dataset_support.tsv');groups=a[['dataset','subtype']].drop_duplicates().sort_values(['dataset','subtype']);fig,ax=plt.subplots(figsize=(15,9));mx=a['mean'].max()
for i,row in enumerate(groups.itertuples(index=False)):
 aa=a[a.dataset.eq(row.dataset)&a.subtype.eq(row.subtype)].set_index('gene').reindex(markers)
 ax.scatter(np.arange(len(markers)),np.full(len(markers),i),s=aa.detection*130,c=aa['mean'],cmap='viridis',norm=Normalize(0,mx),edgecolors='gray',linewidths=.2)
ax.set(yticks=range(len(groups)),yticklabels=[f'{x.dataset}: {x.subtype} (cells={int(b.loc[b.dataset.eq(x.dataset)&b.subtype.eq(x.subtype), "n_cells"].sum()/len(markers))}; specimens={b.loc[b.dataset.eq(x.dataset)&b.subtype.eq(x.subtype), "sample_uid"].nunique()})' for x in groups.itertuples()],xticks=range(len(markers)),xticklabels=markers,title='Frozen B marker panel; specimen-equal expression; size = detection; no relabeling');ax.tick_params(axis='x',rotation=60);fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,mx),cmap='viridis'),ax=ax,label='Mean RNA:data');fig.tight_layout();figsave(fig,'ALL','B_marker_comparability')
print('Focus figures done')
