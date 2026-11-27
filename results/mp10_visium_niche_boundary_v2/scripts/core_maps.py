from lib import *
import argparse,matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('--version',default='primary');arg=p.parse_args();out=O if arg.version=='primary' else O/'versions'/arg.version
sec='GSE319536_BC01';d=load(f'spots/{sec}.parquet');d['tissue_component_id']=d.block
valid=d.expression_QC.to_numpy()&np.isfinite(d.T_f);adj=tissue_graph(d,valid);deg=np.asarray(adj.sum(1)).ravel();ok=deg>0;W=sparse.diags(1/np.maximum(deg,1))@adj
t=d.T_f.to_numpy();e=d.E_f.to_numpy()
def standardized(a):return (a-a[ok].mean())/a[ok].std()
zT=standardized(t);zE=standardized(e);lag=W@t
w=load(f'local/{sec}.parquet').query("level=='window500'").set_index('local_id')
a=d[['section_id','barcode','x_um','y_um']].copy()
a['T_density']=d.T_ab/AREA
a['P_T_window']=d.window500.map(w.P_T.where((w.E_n>=5)&(w.per_n>=4)))
a['neighbor_T_relative']=np.where(ok,lag/t[ok].mean(),np.nan)
a['local_Moran_T']=np.where(ok,zT*(W@zT),np.nan)
a['E_T_window_fraction']=d.window500.map(w.E_T_fraction.where(w.E_n>=5))
a['Lee_local_product']=np.where(ok,(W@zE)*(W@zT),np.nan)
a['mixed_E_window_share']=d.window500.map(d.region_main.isin(['EIS','EIN']).groupby(d.window500).sum()/d.all_E.groupby(d.window500).sum().replace(0,np.nan))
a['I_fraction']=d.I_f
save(a,out/'plot_data/core_maps_BC01.tsv.gz');idx=read(out/'FIGURE_INDEX.tsv');rows=[]
labels=['T density (estimated cells/mm²)','Local T penetration (log2 density ratio)','Neighbor T / global T fraction','Local Moran I (total T)','T fraction in E (500 μm window)','Local contribution to Lee L (E, T)','Immune-enriched fraction of E area','All immune fraction']
for metric,label in zip(a.columns[4:],labels):
    fig,ax=plt.subplots(figsize=(5.3,4.6),layout='constrained');ax.scatter(a.x_um,a.y_um,c='#dedede',s=4,lw=0);q=a[metric].notna();im=ax.scatter(a.x_um[q],a.y_um[q],c=a.loc[q,metric],s=4,lw=0,cmap='viridis');fig.colorbar(im,ax=ax,label=label);ax.invert_yaxis();ax.set_aspect('equal');ax.set(xlabel='x (μm)',ylabel='y (μm)',title=sec)
    key='spatial_'+metric
    for ex in ['pdf','png']:fig.savefig(out/'figures'/(key+'.'+ex),dpi=300,bbox_inches='tight')
    plt.close(fig);rows.append(dict(figure=key,title=label,plot_data='plot_data/core_maps_BC01.tsv.gz',caption='Spatial input underlying the section-level metric; E-weighted summaries, window dispersion and global correlations are computed separately. Grey denotes missing local support.'))
j=pd.DataFrame(rows);save(pd.concat([idx[~idx.figure.isin(j.figure)],j]),out/'FIGURE_INDEX.tsv')

