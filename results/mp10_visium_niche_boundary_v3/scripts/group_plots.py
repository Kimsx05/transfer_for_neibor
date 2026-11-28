from common3 import *
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from down import locked
theme();cp={'C1':'#31708b','C2':'#db8247','C3':'#845898'};C=OUT/'classification';T=OUT/'tables';F=OUT/'figures/groups'
lab=locked();labs=lab[lab.tumor].set_index('section_id');rows=[]
def export(fig,key,tab,note):
    for ext in ['pdf','png']:fig.savefig(F/(key+'.'+ext),dpi=300,bbox_inches='tight')
    plt.close(fig);rows.append(dict(figure=key,pdf='figures/groups/'+key+'.pdf',png='figures/groups/'+key+'.png',data=tab,caption=note))
def pieces(frame,x):
    frame=frame.sort_values(x)
    return frame.groupby(frame[x].diff().ne(100).cumsum())
def one():return plt.subplots(figsize=(6.2,4.3),layout='constrained')
def tags(b,dom):return ('All E' if b==B[0] else 'ENN + ENS')+' · '+{'E_inner':'E interior','nonE_outer':'non-E exterior','EIN_EIS_outer':'EIN/EIS exterior'}[dom]
def xaxis(ax):ax.set_xlabel('Signed distance (μm)');ax.invert_xaxis();ax.axvline(0,c='#555',ls=':',lw=.7)
bins=read(O/'tables/section_bins.tsv.gz').merge(labs[['group']],left_on='section_id',right_index=True)
curve=read(T/'descriptive_group_curves.tsv.gz');diff=read(T/'model_differences.tsv.gz')
dsnames=sorted(labs.dataset_id.unique());dscol=dict(zip(dsnames,plt.cm.tab10.colors))
for (b,dom),t in bins.groupby(['boundary','domain']):
    g=curve[(curve.boundary==b)&(curve.domain==dom)];tit=tags(b,dom)
    for ep,field in [('raw','median'),('centered','centered')]:
        yl='MP10 score' if ep=='raw' else 'Centered MP10 score';fig,ax=one()
        for (sec,grp),a in t[t.n>=5].groupby(['section_id','group']):
            for _,v in pieces(a,'center'):ax.plot(v.center,v[field],c=cp[grp],alpha=.17,lw=.65,marker='.',ms=1.6)
        for grp,a in g[g.scale==ep].groupby('group'):
            for _,v in pieces(a,'center'):ax.plot(v.center,v['mean'],color=cp[grp],lw=1.8,marker='.',ms=3)
        ax.set(ylabel=yl,title=tit);xaxis(ax);ax.legend(handles=[Line2D([0],[0],color=c,label=k) for k,c in cp.items()],frameon=False)
        export(fig,b+'_'+dom+'_'+ep+'_sections', 'tables/section_bins.tsv.gz; tables/descriptive_group_curves.tsv.gz','Thin = tumor section median including unknown unit, thick = equal sections within unit, equal units within source and equal available sources. Domains plotted separately.')
        key=labs.sort_values(['group','dataset_id','stability']);cols=np.arange(t.center.min(),t.center.max()+1,100)[::-1]
        arr=t[t.n>=5].pivot(index='section_id',columns='center',values=field).reindex(index=key.index,columns=cols)
        fig=plt.figure(figsize=(8,8.4),layout='constrained');axg=fig.add_gridspec(1,4,width_ratios=[.23,.23,.23,8]);axes=[fig.add_subplot(axg[0,i]) for i in range(4)]
        for label,vals,colors,aa in [('C',key.group.map({k:i for i,k in enumerate(G)}),list(cp.values()),axes[0]),('Src',key.dataset_id.map({k:i for i,k in enumerate(dsnames)}),list(dscol.values()),axes[1])]:
            aa.imshow(np.array(vals)[:,None],cmap=ListedColormap(colors),aspect='auto',vmin=0,vmax=len(colors)-1);aa.set_title(label,fontsize=7);aa.axis('off')
        axes[2].imshow(key.stability.to_numpy()[:,None],cmap='Greys',vmin=0,vmax=1,aspect='auto');axes[2].set_title('Freq.',fontsize=7);axes[2].axis('off')
        ax=axes[3];ma=arr.to_numpy();kwargs=dict(cmap='viridis') if ep=='raw' else dict(cmap='RdBu_r',vmin=-np.nanmax(abs(ma)),vmax=np.nanmax(abs(ma)))
        im=ax.imshow(ma,aspect='auto',interpolation='nearest',**kwargs);ax.set_yticks(range(len(key)),key.index,fontsize=4.6);ax.yaxis.tick_right();ticks=np.unique(np.linspace(0,len(cols)-1,min(8,len(cols))).astype(int));ax.set_xticks(ticks,cols[ticks].astype(int));ax.set_xlabel('Signed distance (μm)');ax.set_title(tit);fig.colorbar(im,ax=ax,location='bottom',shrink=.5,label=yl);fig.legend(handles=[Line2D([0],[0],c=dscol[x],lw=4,label=x.split('__')[0]) for x in dsnames],loc='outside lower center',ncol=4,frameon=False,fontsize=7)
        export(fig,b+'_'+dom+'_'+ep+'_heatmap','tables/section_bins.tsv.gz; classification/section_classification.tsv','All 60 tumor sections; white cells lack an eligible bin. Annotation strips indicate PAM cluster, source and assignment frequency.')
        fig,ax=one();n=0
        for (g1,g2),d in diff[(diff.boundary==b)&(diff.domain==dom)&(diff.variant=='primary')&(diff.scale==ep)&diff.contrast.str.count(':').eq(1)].groupby(['g1','g2']):
            d=d.sort_values('center');col=list(cp.values())[n];ax.plot(d.center,d['diff'],color=col,label=g1+' − '+g2);ax.fill_between(d.center,d.lo,d.hi,color=col,alpha=.13);n+=1
        if n:
            ax.axhline(0,color='#555',lw=.7);ax.set(ylabel=yl+' difference',title=tit);xaxis(ax);ax.legend(frameon=False)
            export(fig,b+'_'+dom+'_'+ep+'_differences','tables/model_differences.tsv.gz','Pair-specific common sources and contiguous 70% distance support. Model contrasts with pointwise 95% whole-unit bootstrap intervals; no simultaneous band.')
        else:plt.close(fig)
pc=read(C/'PCA_scores.tsv').set_index('section_id').join(labs[['group','dataset_id','medoid']]).dropna();var=read(C/'PCA_variance.tsv')
for by in ['group','dataset_id']:
    fig,ax=one()
    for name,v in pc.groupby(by):ax.scatter(v.PC1,v.PC2,c=[(cp if by=='group' else dscol)[name]],s=26,label=name.split('__')[0])
    vv=pc[pc.medoid.eq(True)];ax.scatter(vv.PC1,vv.PC2,marker='o',s=90,facecolors='none',edgecolors='#222',lw=1);ax.set(xlabel=f'PC1 ({100*var.variance.iloc[0]:.1f}%)',ylabel=f'PC2 ({100*var.variance.iloc[1]:.1f}%)');ax.legend(frameon=False,fontsize=8)
    export(fig,'PCA_'+by,'classification/PCA_scores.tsv; classification/PCA_variance.tsv','PCA of family-weighted X_pam. Rings mark medoids. Tumor discovery sections.')
f=read(C/'bootstrap_runs.tsv.gz');fig,ax=one();sl=f[(f.rule=='primary')&(f.status=='ok')];kk=list(sl.groupby(['method','K']));ax.boxplot([g.ARI for _,g in kk],tick_labels=[m+' '+str(k) for (m,k),g in kk]);ax.set_ylabel('Adjusted Rand index');export(fig,'bootstrap_ARI','classification/bootstrap_runs.tsv.gz','100 unit-level refits; both variable selection and scaling re-estimated.')
gm=read(C/'GMM_models.tsv');fig,ax=one()
for cov,v in gm[gm.reg==.001].groupby('cov'):ax.plot(v.K,v.BIC,'o-',label=cov)
ax.set(xlabel='Number of components',ylabel='BIC',xticks=[1,2,3,4]);ax.legend(frameon=False);export(fig,'GMM_selection','classification/GMM_models.tsv','Identical unweighted standardized input and 20 starts for each fitted model.')
levels=read(T/'section_levels.tsv');fig,ax=one();rng=np.random.default_rng(52)
for j,grp in enumerate(G):
    v=levels[levels.tumor&levels.group.eq(grp)];ax.scatter(j+rng.uniform(-.12,.12,len(v)),v.mean_MP10,s=20,c=[dscol[n] for n in v.dataset_id])
ax.set(xticks=range(3),xticklabels=G,ylabel='Mean MP10 score');export(fig,'overall_MP10','tables/section_levels.tsv','Each tumor section shown once. Source-balanced unit-level contrasts are in auxiliary_tests.')
save(rows,OUT/'audit/group_figure_index.tsv')

