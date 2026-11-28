from common3 import *
from scipy.sparse.csgraph import connected_components
rows=[]
for s in M.index:
    d=spot(s,['region_main','expression_QC','MP10','array_row','array_col','block','all_E',B[1],'dist_'+B[0],'dist_'+B[1]])
    valid=d.region_main.notna()&d.region_main.ne('Unknown');adj=graph(d,valid.to_numpy())
    for b in B:
        m=d[b].to_numpy();edges=adj.tocoo();ok=m[edges.row]==m[edges.col]
        sub=sparse.csr_matrix((edges.data[ok],(edges.row[ok],edges.col[ok])),shape=adj.shape)
        _,cc=connected_components(sub,directed=False)
        di=d['dist_'+b];v=valid&d.expression_QC&np.isfinite(d.MP10)&di.notna()
        z=pd.DataFrame(dict(block=d.block[v],side_component=cc[v],domain=np.where(d.loc[v,b],'E_inner',np.where(d.loc[v,'all_E'],'EIN_EIS_outer','nonE_outer')),center=np.where(m[v],-1,1)*(np.maximum(0,np.ceil(abs(di[v])/100-1e-10)-1)*100+50)))
        n=z.value_counts().rename('n').reset_index();n['boundary']=b;n['section_id']=s;rows.append(n)
out=pd.concat(rows,ignore_index=True)
cols=['section_id','boundary','domain','center']
ref=read(O/'tables/section_bins.tsv.gz').set_index(cols).n.sort_index()
np.testing.assert_array_equal(out.groupby(cols).n.sum().sort_index().to_numpy(),ref.to_numpy())
save(out,O/'tables/bin_side_component_support.tsv.gz')

