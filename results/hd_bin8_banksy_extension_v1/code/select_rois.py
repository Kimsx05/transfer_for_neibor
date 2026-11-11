"""Select 256um squares from input geometry/information, before predictions or BANKSY."""
from pathlib import Path
import pandas as pd,numpy as np,sys,json
r=Path(__file__).resolve().parents[1];s=sys.argv[1];o=r/s
if (o/'fixed_rois.tsv').exists():sys.exit()
d=pd.read_csv(o/'analysis_metadata.tsv.gz',sep='\t');grid=pd.DataFrame({'ix':np.floor(d.x_um/256).astype(int),'iy':np.floor(d.y_um/256).astype(int),'nf':d.nFeature_RNA,'umi':d.nCount_RNA})
a=grid.groupby(['ix','iy']).agg(n=('nf','size'),median_nf=('nf','median'),UMI=('umi','sum')).reset_index();a=a[a.n>=128].copy();rows=[];chosen=[]
for name,qt in [('low_information',.1),('typical_information',.5),('high_information',.9)]:
 target=a.median_nf.quantile(qt);a['dist']=(a.median_nf-target).abs();c=a.sort_values(['dist','n','ix','iy'],ascending=[True,False,True,True]);q=c[~c.index.isin(chosen)].iloc[0];chosen.append(q.name);rows.append(dict(roi=name,x=int(q.ix)*256,y=int(q.iy)*256,X=(int(q.ix)+1)*256,Y=(int(q.iy)+1)*256,selection=f'256um grid window closest to section window median nFeature quantile {qt}; >=128 QC bins; before reference/BANKSY',n_preselection=int(q.n),median_nf_preselection=q.median_nf))
# Fourth window covers a sparse boundary of the supplied/QC support; not claimed anatomical edge.
c=a[~a.index.isin(chosen)].sort_values(['n','ix','iy']);q=c.iloc[0];rows.append(dict(roi='sparse_input_edge',x=int(q.ix)*256,y=int(q.iy)*256,X=(int(q.ix)+1)*256,Y=(int(q.iy)+1)*256,selection='Lowest occupied 256um window with >=128 QC bins; input/expression support edge, not established tissue edge',n_preselection=int(q.n),median_nf_preselection=q.median_nf))
pd.DataFrame(rows).to_csv(o/'fixed_rois.tsv',sep='\t',index=False);print(s,'ROI_FROZEN',rows)
