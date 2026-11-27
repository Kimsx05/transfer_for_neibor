from lib import *
import matplotlib.pyplot as plt
links=[]
for sec in M.index:
    d=pd.read_parquet(O/'spots'/f'{sec}.parquet',columns=['barcode','x_um','y_um','array_row','array_col','block','all_E','T_f','T_ab','window500','window1000'])
    d['tissue_component_id']=d.block
    g=sparse.triu(tissue_graph(d,d.all_E.to_numpy()),k=1).tocoo()
    for w in [500,1000]:
        v=d['window'+str(w)].to_numpy();a=v[g.row];b=v[g.col];pairs=sorted(set(tuple(sorted((i,j))) for i,j in zip(a,b) if i!=j))
        links.extend(dict(section_id=sec,level='window'+str(w),a=i,b=j) for i,j in pairs)
    if sec in ['GSE319536_BC01','GSE171351_A1','GSE238145_B13']:
        fig,ax=plt.subplots(1,2,figsize=(7,3))
        for a,k in zip(ax,['all_E','T_f']):a.scatter(d.x_um,d.y_um,c=d[k],s=3,lw=0,cmap='viridis');a.set_aspect('equal');a.set_title(k);a.axis('off')
        fig.tight_layout();fig.savefig(O/'classification'/f'{sec}_blind.jpg',dpi=110);plt.close(fig)
save(links,'classification/window_links.tsv.gz')
