from lib import *
import argparse,matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('--version',default='primary');arg=p.parse_args();out=O if arg.version=='primary' else O/'versions'/arg.version
figdir=out/'figures';index=read(out/'FIGURE_INDEX.tsv');rows=[]
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
def emit(f,key,title,src,note):
    for ext in ['pdf','png']:f.savefig(figdir/(key+'.'+ext),dpi=300,bbox_inches='tight')
    plt.close(f);rows.append(dict(figure=key,title=title,plot_data=src,caption=note))
cl=read(out/'classification/section_classification.tsv');cl=cl[cl['mode']=='tumor_discovery'].sort_values('F_entry')
art=dict(zip(G+['uncertain'],['#27789a','#c87941','#925c91','#888888']))
f,a=plt.subplots(figsize=(5.4,4.3))
for g,v in cl.groupby('provisional_architecture_group'):
    a.scatter(v.F_entry,v.F_non,label=g,c=art[g],s=24,alpha=.75)
cfg=json.loads((out/'classification/classification_lock.json').read_text())['config']
a.axvline(cfg['mixed_area_share'],c='#777',ls=':',lw=.75);a.axhline(cfg['mixed_area_share'],c='#777',ls=':',lw=.75);a.set(xlabel='Infiltrated-like share of eligible E area',ylabel='Non-infiltrated-like share of eligible E area');a.legend(frameon=False,fontsize=7)
emit(f,'supported_areas','Local pattern coexistence','classification/section_classification.tsv','Only spatially supported pattern areas count. Dotted lines indicate minimum mixed share, with additional coverage and window-count requirements.')
st=read(O/'classification/cluster_stability.tsv.gz');f,a=plt.subplots(figsize=(5.4,4));a.boxplot([v.ARI for _,v in st.groupby('K')],tick_labels=[str(n) for n in sorted(st.K.unique())],showfliers=True);a.set(xlabel='Number of components',ylabel='Resampled adjusted Rand index')
emit(f,'GMM_stability','Unit-resampled model stability','classification/cluster_stability.tsv.gz','80 unit-block resamples per K in the same original PC space. K=1 agreement is automatic.')
miss=read(O/'classification/missingness.tsv.gz').set_index('section_id');f,a=plt.subplots(figsize=(8,4));im=a.imshow(miss.astype(int).to_numpy(),cmap='Greys',vmin=0,vmax=1,aspect='auto',interpolation='nearest');a.set_xlabel('Candidate metric (ordered as in source table)');a.set_ylabel('Tumor section');a.set_title('Structural missingness');emit(f,'candidate_missingness','Missing candidate values','classification/missingness.tsv.gz','Black = missing, white = observed; exact row and column labels are in the source table.')
pts=read(out/'plot_data/association_points.tsv.gz');ds=sorted(M.dataset_id.unique());cc=dict(zip(ds,plt.cm.tab10.colors))
for metric in ['penetration_T_200','whole_T_density']:
    for end in ['D_all_E','D_E_without_I_enriched','mean_MP10']:
        v=pts[(pts.metric==metric)&(pts.endpoint==end)];f,a=plt.subplots(figsize=(6.2,4))
        for name,s in v.groupby('dataset_id'):
            a.scatter(s.x,s.y,s=13,alpha=.75,label=name.split('__')[0],color=cc[name])
            if len(s)>=3 and s.x.nunique()>1:
                xx=np.array([s.x.min(),s.x.max()]);a.plot(xx,np.polyval(np.polyfit(s.x,s.y,1),xx),lw=1,c=cc[name])
        a.set(xlabel='T-cell penetration (log2 density ratio)' if metric=='penetration_T_200' else 'T-cell density (estimated cells/mm²)',ylabel='B1 − B2 (MP10 score)' if end.startswith('D_') else 'Mean MP10 score',title=end.replace('D_E_without_I_enriched','ENN + ENS').replace('D_all_E','All E').replace('mean_MP10','Overall level'))
        a.legend(frameon=False,fontsize=7);emit(f,'association_'+metric+'_'+end,metric+' · '+end,'plot_data/association_points.tsv.gz','Unit-level points; source-specific linear display limited to observed x. Inference is by Spearman association.')
allrows=pd.DataFrame(rows);save(pd.concat([index[~index.figure.isin(allrows.figure)],allrows]),out/'FIGURE_INDEX.tsv')

