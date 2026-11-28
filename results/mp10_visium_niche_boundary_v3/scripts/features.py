from common3 import *
from scipy.sparse.csgraph import connected_components
import shapely
spec=[('500',500,0,0),('1000',1000,0,0),('x250',500,250,0),('y250',500,0,250),('xy250',500,250,250)]
def quant(x,n,p):
    if not len(x):return np.nan
    i=np.argsort(x,kind='stable');return x[i][np.searchsorted(np.cumsum(n[i]),p*n.sum(),side='left')]
def section(s):
    cols=['barcode','x_um','y_um','array_row','array_col','block','all_E','region_main','expression_QC','total_abundance','E_f','T_f','I_f','T_ab','I_ab','UMI','nGene','dist_all_E','all_E_component']
    d=spot(s,cols);valid=(d.expression_QC&np.isfinite(d.total_abundance)&(d.total_abundance>0)).to_numpy();E=valid&d.all_E.to_numpy()
    adj=graph(d,valid);deg=np.asarray(adj.sum(1)).ravel();has=deg>0;W=sparse.diags(1/np.maximum(deg,1))@adj
    cross=sparse.triu(graph(d,E),k=1).tocoo();ei,ej=cross.row,cross.col
    near=(d.dist_all_E.gt(0)&d.dist_all_E.le(200)&valid).to_numpy()
    # Associate each exterior point with the closest observed contact within its exterior component.
    a=sparse.triu(graph(d,d.region_main.ne('Unknown').to_numpy()),k=1).tocoo()
    mask=d.all_E.to_numpy();same=mask[a.row]==mask[a.col]
    _,comp=connected_components(sparse.coo_matrix((np.ones(sum(same)),(a.row[same],a.col[same])),shape=adj.shape),directed=False)
    owner=np.full(len(d),-1);edges=pd.read_parquet(V/'boundary'/f'{s}.parquet').query("boundary=='all_E'");ids=pd.Index(d.barcode)
    for cc,sg in edges.groupby('outside_component'):
        ids2=np.flatnonzero(near&(comp==cc))
        if len(ids2):
            tree=shapely.STRtree(shapely.linestrings(sg[['x1','y1','x2','y2']].to_numpy().reshape(-1,2,2)))
            jj=tree.nearest(shapely.points(d.loc[ids2,['x_um','y_um']].to_numpy()));owner[ids2]=ids.get_indexer(sg.inside_barcode.iloc[jj])
    assert np.all(owner[near]>=0)
    locs=[];metrics=[];alllinks=[]
    for v,width,dx,dy in spec:
        cell=d.block.astype(str)+'|'+np.floor((d.x_um-dx)/width).astype(int).astype(str)+'|'+np.floor((d.y_um-dy)/width).astype(int).astype(str)
        links=pd.DataFrame({'a':cell.iloc[ei].to_numpy(),'b':cell.iloc[ej].to_numpy()});links=links[links.a!=links.b]
        ll=np.sort(links[['a','b']].to_numpy(),axis=1);links=pd.DataFrame(ll,columns=['a','b']).value_counts().reset_index(name='E_edges')
        alllinks.append(links.assign(section_id=s,grid=v))
        z=pd.DataFrame({'window':cell,'E_n':E.astype(int),'obs_n':valid.astype(int),'E_total':np.where(E,d.total_abundance,0),'lag_n':(has&E).astype(int),'x_E':np.where(E,d.x_um,0),'y_E':np.where(E,d.y_um,0),'E_fraction_sum':np.where(E,d.E_f,0)})
        for line in ['T','I']:
            z['ab_'+line]=np.where(E,d[line+'_ab'],0);z['lag_'+line]=np.where(has&E,W@d[line+'_f'],0)
        a=z.groupby('window').sum();a=a[a.E_n>0].copy();outer=d.loc[near,['T_ab','I_ab','total_abundance']].assign(window=cell.iloc[owner[near]].to_numpy());og=outer.groupby('window')
        a['outer_n']=og.size().reindex(a.index,fill_value=0);a['outer_total']=og.total_abundance.sum().reindex(a.index,fill_value=0)
        a['E_area']=a.E_n*area;a['x']=a.x_E/a.E_n;a['y']=a.y_E/a.E_n;a['eligible']=a.E_n>=5
        sx=a.index[a.eligible];ix=pd.Index(sx);ed=links[links.a.isin(ix)&links.b.isin(ix)]
        rr=ix.get_indexer(ed.a);cc=ix.get_indexer(ed.b);wa=sparse.coo_matrix((np.ones(2*len(rr)),(np.r_[rr,cc],np.r_[cc,rr])),shape=(len(ix),len(ix))).tocsr();nd=np.asarray(wa.sum(1)).ravel();w=sparse.diags(1/np.maximum(nd,1))@wa
        for line in ['T','I']:
            total=d[line+'_ab'][valid].sum();inner=d[line+'_ab'][E].sum();frac=d[line+'_f'].to_numpy();e=d.E_f.to_numpy();lag=W@frac
            a['frac_'+line]=a['ab_'+line]/a.E_total;a['density_'+line]=a['ab_'+line]/a.E_area;a['neighbor_'+line]=a['lag_'+line]/a.lag_n.replace(0,np.nan)
            a['outer_ab_'+line]=og[line+'_ab'].sum().reindex(a.index,fill_value=0);od=a['outer_ab_'+line]/(a.outer_n.replace(0,np.nan)*area)
            a['penetration_'+line]=np.log2((1+a['density_'+line])/(1+od)).where(a.outer_n>=4)
            f=a.loc[sx,'frac_'+line].to_numpy();n=a.loc[sx,'E_n'].to_numpy()
            mor=np.nan;non=nd>0
            if non.sum()>1 and f[non].std()>0:
                zz=(f-f[non].mean())/f[non].std();mor=np.mean(zz[non]*(w@zz)[non])
            qs=[quant(f,n,j) for j in [.1,.25,.5,.75,.9]]
            metrics.append(dict(section_id=s,grid=v,axis=line,A=qs[2],B=divide(divide(inner,total),divide(E.sum(),valid.sum())),C=divide(divide((e[has]*lag[has]).sum(),e[has].sum()),frac[has].mean()),D=qs[3]-qs[1],E=mor,F=qs[0],G=qs[4],P=np.log2((1+divide(inner,E.sum()*area))/(1+divide(d.loc[near,line+'_ab'].sum(),near.sum()*area))),whole_ab=total,E_ab=inner,total_area=valid.sum()*area,E_area=E.sum()*area,valid_E_area=n.sum()*area,win_n=len(sx),nonisolated_n=non.sum(),islands=d.loc[E,'all_E_component'].nunique(),median_island=d.loc[E,'all_E_component'].value_counts().median(),UMI=d.loc[valid,'UMI'].median(),nGene=d.loc[valid,'nGene'].median(),mixed_E_share=d.region_main.isin(['EIN','EIS']).sum()/E.sum(),outer_coverage=divide(a.loc[a.eligible&a.outer_n.ge(4),'E_n'].sum(),n.sum()),usable_area_fraction=divide(n.sum(),E.sum()),whole_fraction=divide(total,d.loc[valid,'total_abundance'].sum()),density_median=quant(a.loc[sx,'density_'+line].to_numpy(),n,.5)))
        locs.append(a.reset_index().assign(section_id=s,grid=v))
    return metrics,locs,alllinks
if __name__=='__main__':
    results=[];windows=[];links=[]
    for s in M.index:
        m,l,e=section(s);results+=m;windows+=l;links+=e
    z=pd.DataFrame(results);save(z,O/'cache/immune_features.tsv');pd.concat(windows).to_parquet(O/'cache/immune_windows.parquet',index=False);save(pd.concat(links),O/'cache/immune_links.tsv.gz')
    prev=read(V/'classification/immune_input.tsv.gz').set_index('section_id')
    chk=[]
    for axis in ['T','I']:
        now=z[(z.grid=='500')&(z.axis==axis)].set_index('section_id')
        for k,old in [('B','area_corrected_entry_'+axis),('C','CLQ_E_'+axis+'_R1'),('P','penetration_'+axis+'_200')]:
            delta=(prev[old]-now[k]).abs();chk.append(dict(axis=axis,feature=k,max_abs_difference=delta.max(),comparable=delta.notna().sum()))
    save(chk,O/'audit/v2_metric_concordance.tsv')

