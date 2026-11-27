from lib import *
import argparse,matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
p=argparse.ArgumentParser();p.add_argument('--version',default='primary');arg=p.parse_args()
out=O if arg.version=='primary' else O/'versions'/arg.version
figdir=out/'figures';figdir.mkdir(exist_ok=True,parents=True);idx=[]
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'font.family':'DejaVu Sans','savefig.bbox':'tight'})
palette=['#27789a','#c87941','#925c91','#aaa8a1'];bn=['All E','ENN + ENS'];groups=G+['uncertain'];core=read(O/'classification/feature_selection.tsv').query("decision=='include'").feature.tolist()
source_colors=dict(zip(sorted(M.dataset_id.unique()),plt.cm.tab10.colors))
def dsname(s):return s.split('__')[0]
def get(tab,loc=out):return read(loc/tab)
def plot(w=5.6,h=4):return plt.subplots(figsize=(w,h),layout='constrained')
def export(fig,key,title,src,note):
    for ext in ['pdf','png']:fig.savefig(figdir/(key+'.'+ext),dpi=300)
    plt.close(fig);idx.append(dict(figure=key,title=title,plot_data=src,caption=note))
def legend(ax):ax.legend(frameon=False,fontsize=7)
def forest(a,xlabel,title,key,src):
    fig,ax=plot(6.5,max(3,.36*len(a)+1));a=a.reset_index(drop=True)
    for i,r in a.iterrows():
        if np.isfinite(r.effect):
            ax.scatter(r.effect,i,color='#286c85',s=22)
            if np.isfinite(r.low):ax.plot([r.low,r.high],[i,i],color='#286c85',lw=1)
    ax.set_yticks(np.arange(len(a)),a.label);ax.invert_yaxis();ax.axvline(0,c='#777',lw=.7);ax.set_xlabel(xlabel);ax.set_title(title)
    export(fig,key,title,src,'Points are equal-source means. Whiskers are 95% intervals across sources when available; see support columns.')
meta=get('tables/all_section_levels.tsv').set_index('section_id')
tum=meta[~meta.source_context.eq('adjacent_normal_source')]
pc=get('classification/pca_coordinates.tsv',O).set_index('section_id').loc[tum.index]
x=get('classification/immune_input.tsv.gz',O).set_index('section_id')
for color in ['source']+core:
    fig,ax=plot()
    if color=='source':
        for ds,a in tum.groupby('dataset_id'):ax.scatter(pc.loc[a.index,'PC1'],pc.loc[a.index,'PC2'],s=22,color=source_colors[ds],label=dsname(ds),alpha=.8)
        legend(ax)
    else:
        art=ax.scatter(pc.PC1,pc.PC2,c=x.loc[pc.index,color],cmap='viridis',s=25);fig.colorbar(art,ax=ax,label=color)
    ax.set(xlabel='PC1 (39.7%)',ylabel='PC2 (21.1%)',title='Immune architecture · '+color)
    export(fig,'PCA_'+color,'Immune PCA: '+color,'classification/pca_coordinates.tsv; classification/immune_input.tsv.gz','Tumor-source sections; 4 PCs entered model selection. Colors use original metric units.')
k=get('classification/model_selection.tsv',O);fig,ax=plot();ax.plot(k.K,k.BIC,'o-',color='#317b84');ax.set(xticks=k.K,xlabel='Number of components',ylabel='BIC',title='Mixture model selection')
export(fig,'GMM_BIC','GMM model selection','classification/model_selection.tsv','The models share one PC space and covariance regularization.')
kde=get('plot_data/KDE.tsv.gz',O)
for metric in core:
    fig,ax=plot();xs=x.loc[tum.index,metric].dropna()
    if metric in json.loads((O/'classification/natural_model.json').read_text())['log1p']:xs=np.log1p(xs)
    ax.hist(xs,bins='sturges',density=True,color='#d6deda',edgecolor='white')
    for f,a in kde[kde.metric==metric].groupby('bandwidth_factor'):ax.plot(a.x,a.density,label=str(f)+' × Scott',lw=1.1)
    legend(ax);ax.set(xlabel=metric+(' (log1p)' if metric in ['whole_T_density','CLQ_E_T_R1','window500_E_T_fraction_IQR'] else ''),ylabel='Density')
    export(fig,'density_'+metric,metric,'plot_data/KDE.tsv.gz; classification/immune_input.tsv.gz','Each tumor section contributes once. Histograms and sensitivity bandwidths use the same transformed scale.')
c=x.loc[tum.index,core].corr(method='spearman');fig,ax=plot(7,6);im=ax.imshow(c,vmin=-1,vmax=1,cmap='RdBu_r');ax.set_xticks(range(len(core)),core,rotation=65,ha='right',fontsize=7);ax.set_yticks(range(len(core)),core,fontsize=7);fig.colorbar(im,ax=ax,label='Spearman correlation');export(fig,'core_correlations','Core feature correlations','classification/correlations.tsv.gz','Tumor sections only. Complete matrix for every immune candidate is supplied.')
cur=get('plot_data/group_curves.tsv.gz');effects=get('tables/group_effects.tsv.gz');paired=get('tables/paired_summary.tsv');tests=get('tables/main_six_tests.tsv');lev=get('tables/level_composition_summary.tsv.gz')
for b,title in zip(B,bn):
    for m in ['score','centered','P95']:
        fig,(ax,count)=plt.subplots(2,1,figsize=(6.5,4.8),height_ratios=[3,1],sharex=True,layout='constrained')
        for g,col in zip(G,palette):
            a=cur[(cur.boundary==b)&(cur.group==g)&(cur.metric==m)&(cur.n_units>0)].sort_values('center')
            sm=cur[(cur.boundary==b)&(cur.group==g)&(cur.metric==(m if m=='P95' else m+'_smooth'))&(cur.n_units>0)].sort_values('center')
            if not len(a):ax.plot([],[],color=col,label=g+' (0)');continue
            ax.scatter(a.center,a.effect,s=10,color=col,alpha=.7,label=g)
            count.scatter(a.center,a.n_units,s=9,color=col)
            for _,v in sm.groupby((sm.center.diff().ne(100)|(sm.center*sm.center.shift()<0)).cumsum()):
                ax.plot(v.center,v.effect,color=col,lw=1.1);ax.fill_between(v.center,v.low,v.high,color=col,alpha=.12)
        for aa in [ax,count]:aa.axvline(0,color='#333',lw=.8)
        ax.set(ylabel={'score':'MP10 score','centered':'Centered MP10 score','P95':'MP10-high proportion'}[m],title=title);legend(ax);count.set(xlabel='Signed distance (μm)',ylabel='Units');count.set_ylim(bottom=0)
        export(fig,b+'_'+m,title+' | '+m,'plot_data/group_curves.tsv.gz','Observed 100 μm bins with ≥5 points per section. Lines use neighboring bins on the same side for continuous scores. Shading is a 95% interval over contributing source means. Positive-side E points occur only in EIN/EIS.')
    for ep in ['score','adjusted','P95_pp']:
        a=effects.query("boundary==@b and variant=='native' and width==200 and endpoint==@ep and scope=='overall'").copy();a['label']=a.group+' (n='+a.n_units.astype(str)+')'
        forest(a,'B1 − B2 (percentage points)' if ep=='P95_pp' else 'B1 − B2 (MP10 score)',title+' · '+ep,b+'_effects_'+ep,'tables/group_effects.tsv.gz')
    a=effects.query("boundary==@b and endpoint=='score' and group=='all' and scope=='overall'").copy();a['label']=a.variant.str.replace('_',' ')+' · '+a.width.astype(str)+' μm'
    forest(a,'B1 − B2 (MP10 score)',title+' · sensitivity',b+'_sensitivity','tables/group_effects.tsv.gz')
    sec=get('tables/section_effects.tsv.gz',O).query("boundary==@b and variant=='native' and width==200 and eligible").merge(tum.reset_index(),on='section_id').sort_values(['dataset_id','score'])
    fig,ax=plot(7.5,8.5);i=np.arange(len(sec));ax.scatter(sec.score,i,c=sec.dataset_id.map(source_colors).tolist(),s=10);ax.set_yticks(i,sec.section_id,fontsize=5);ax.axvline(0,c='#777',lw=.7);ax.set_xlabel('B1 − B2 (MP10 score)');ax.set_title(title+' · evaluable tumor sections');export(fig,b+'_sections',title+' | per section','tables/section_effects.tsv.gz','Each observed section contributes one descriptive estimate.')
a=paired.query("endpoint=='paired_difference' and group=='all' and scope=='overall'").copy();a['label']=a.domain+' (n='+a.n_units.astype(str)+')';forest(a,'D (ENN + ENS) − D (all E)','Paired boundaries','paired_comparison','tables/paired_summary.tsv')
a=tests.query("g1=='infiltrated-like' and g2=='non-infiltrated-like'").copy();a['label']=a.boundary;forest(a,'Difference in D (MP10 score)','Infiltrated-like minus non-infiltrated-like','direct_comparison','tables/main_six_tests.tsv')
hl=get('tables/highlow_boundary.tsv').query("scope=='overall'").copy();hl['label']=hl.boundary.str.replace('E_without_I_enriched','ENN + ENS')+' · '+hl.exposure_rule.str.replace('exclude_B1_.*','exclude B1',regex=True);forest(hl,'D (high) − D (low)','Sample MP10 and boundary effect','highlow_D','tables/highlow_boundary.tsv')
for col in ['mean_MP10','median_MP10']:
    fig,ax=plot();rng=np.random.default_rng(SEED)
    for i,g in enumerate(groups):
        a=tum[tum.provisional_architecture_group.eq(g)]
        ax.scatter(i+rng.uniform(-.14,.14,len(a)),a[col],s=17,c=[source_colors[d] for d in a.dataset_id],alpha=.8)
    ax.set_xticks(range(len(groups)),['Infiltrated','Non-infiltrated','Mixed','Uncertain'],rotation=15);ax.set(ylabel='MP10 score',title='Section '+col.split('_')[0]+' in all E')
    export(fig,'overall_'+col,'Section '+col,'tables/all_section_levels.tsv','Tumor sections. Inference balances units and sources in level_composition_tests.')
for kind in ['region','R12','within_all_E','score_in']:
    a=lev[(lev.scope=='overall')&lev.measure.str.startswith(kind+'_')&lev.group.isin(groups)&lev.effect.notna()].copy()
    if kind=='region':a=a[~a.measure.eq('region_Unknown')]
    piv=a.pivot(index='group',columns='measure',values='effect');piv.columns=[str(x)[len(kind)+1:] for x in piv]
    fig,ax=plot(6.4,3.8)
    if kind=='score_in':piv.T.plot.bar(ax=ax,color=[palette[groups.index(x)] for x in piv.index],rot=0);ax.set_ylabel('MP10 score')
    else:piv.plot.barh(stacked=True,ax=ax,cmap='Set3');ax.set_xlabel('Region proportion' if kind!='R12' else 'Assigned niche proportion')
    ax.set_title(kind.replace('_',' '));ax.legend(frameon=False,fontsize=6,ncol=3,loc='upper left',bbox_to_anchor=(1.0,1.))
    export(fig,'composition_'+kind,kind,'tables/level_composition_summary.tsv.gz','Equal-weighted source means among units in each group. Classified tissue, all E or assigned R12 domain as indicated.')
pr=get('tables/highlow_architecture_proportions.tsv');pr=pr[pr.den>0].groupby(['highlow','architecture']).fraction.mean().unstack().reindex(columns=groups,fill_value=0)
fig,ax=plot();pr.plot.barh(stacked=True,ax=ax,color=palette);ax.set_xlabel('Architecture proportion');ax.legend(frameon=False,fontsize=7);export(fig,'highlow_proportions','Architecture by MP10 exposure','tables/highlow_architecture_proportions.tsv','One unit vote divided across its labeled sections, then equal source weights; ties retained.')
# Frozen representative, same selection as v1 (coverage and tissue structure).
sec='GSE319536_BC01';d=load('spots/'+sec+'.parquet');seg=load('boundary/'+sec+'.parquet')
for b,title in zip(B,bn):
    fig,ax=plot(6,5.3);ax.scatter(d.x_um,d.y_um,c='#e5e5e1',s=4,lw=0);a=d[d.new_epi_domain];art=ax.scatter(a.x_um,a.y_um,c=a.MP10,cmap='magma',s=5,lw=0);v=seg[seg.boundary==b];ax.add_collection(LineCollection(v[['x1','y1','x2','y2']].to_numpy().reshape(-1,2,2),colors=np.where(v.contact_type=='to_EIN_EIS','#59b29c','#5093ba'),linewidths=.55));ax.set_aspect('equal');ax.invert_yaxis();fig.colorbar(art,ax=ax,label='MP10 score');ax.set(title=sec+' · '+title,xlabel='x (μm)',ylabel='y (μm)');export(fig,'map_'+b,title+' and MP10','spots/GSE319536_BC01.parquet; boundary/GSE319536_BC01.parquet','Grey marks observed tissue outside the scoring domain. Blue edges contact non-E, green edges contact EIN/EIS.')
save(idx,out/'FIGURE_INDEX.tsv')

