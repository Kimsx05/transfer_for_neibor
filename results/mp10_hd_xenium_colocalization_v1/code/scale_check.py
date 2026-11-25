from base import *
from geometry import mask_for,pixel,edges
from spatial_analysis import enrich
s='Cervilla';d=O/'hd'/s;mask,aff,b=mask_for(s)
ann=pd.read_parquet(A/'hd_epi_niche_dual_methods_v1'/s/'nodes.parquet',columns=['bin_id','type_blca'])
b=b.merge(ann,on='bin_id',validate='one_to_one');b['label']=b.type_blca.fillna('Unknown')
xy=b[['x_um','y_um']].to_numpy();center=xy[((xy-np.median(xy,axis=0))**2).sum(1).argmin()]
use=(np.abs(xy-center)<=512).all(1);small=b[use].copy();x=xy[use];rr=pixel(x,aff);cc=small.component.to_numpy('i4')
adj=edges(x,x,rr,rr,mask,cc,cc,50,True);enrich(s,d,adj,small.label.to_numpy(),50,'scale8_major')
z=pd.read_csv(d/'units.tsv.gz',sep='\t');members=pd.read_parquet(d/'members.parquet').merge(b[['bin_id','label']],on='bin_id',validate='one_to_one')
tab=pd.crosstab(members.region_id,members.label);win=tab.idxmax(1).mask(tab.eq(tab.max(1),axis=0).sum(1)>1,'Mixed').reindex(z.region_id)
use=(np.abs(z[['x_um','y_um']].to_numpy()-center)<=512).all(1);g=sparse.load_npz(d/'graph_50.npz')[use][:,use]
enrich(s,d,g,win.to_numpy()[use],50,'scale16_major')
put(pd.DataFrame([{'sample':s,'x_center_um':center[0],'y_center_um':center[1],'width_um':1024,'n_8':len(x),'n_16':sum(use),'rule':'prespecified section, window centered on observed bin closest to marginal median coordinates; BLCA SingleR major votes at both scales, ties Mixed'}]),d/'scale_window.tsv')

