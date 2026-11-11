"""Input-support edge distance: explicitly treat outside the observed grid as unsupported.
Does not change any spatial neighbor, count, QC, model or previous result file.
"""
from pathlib import Path
import pandas as pd,numpy as np
from scipy import ndimage
r=Path(__file__).resolve().parents[1];old=r.parent/'hd_bin8_banksy_v1';base=r.parent/'hd_bin8_comparison_v1';rows=[]
for s in ['Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2']:
 new=s not in ['Cervilla','D1_Johnson'];o=r/s;b=o if new else base/s
 d=pd.read_csv(b/'bin_metadata.tsv.gz',sep='\t');f='in_tissue_original' if new else 'in_tissue';raw=d[d[f]==1] if s in ['Cervilla','D1_Johnson','UC_SCC'] else d[d.supplied_bin]
 grid=np.zeros(tuple((d[['array_row','array_col']].max()+1).astype(int)),bool);grid[raw.array_row,raw.array_col]=True
 dist=ndimage.distance_transform_edt(np.pad(grid,1,constant_values=False))[1:-1,1:-1]*8
 sp=pd.read_csv((o if new else old/s)/'spatial_metadata.tsv.gz',sep='\t');q=d.set_index('bin_id').loc[sp.bin_id];v=dist[q.array_row,q.array_col];prev=sp.native_edge_distance_um.to_numpy();change=np.abs(v-prev)>1e-6
 rows.append(dict(section=s,n=len(v),distance_changed_bins=int(change.sum()),old_edge16_bins=int((prev<=16).sum()),corrected_edge16_bins=int((v<=16).sum()),edge16_membership_changed_bins=int(((prev<=16)!=(v<=16)).sum()),changed_bins_UMI=int(q.nCount_RNA.to_numpy()[change].sum()),neighbor_graph_changed=False))
 if new:
  sp['native_edge_distance_um']=v;sp.to_csv(o/'spatial_metadata.tsv.gz',sep='\t',index=False)
 else:
  supp=pd.read_csv(o/'supplement_metadata.tsv.gz',sep='\t');assert supp.bin_id.equals(sp.bin_id);supp['ext_input_edge_distance_um']=v;supp.to_csv(o/'supplement_metadata.tsv.gz',sep='\t',index=False)
pd.DataFrame(rows).to_csv(r/'tables/input_edge_distance_correction.tsv',sep='\t',index=False);print(pd.DataFrame(rows).to_string(index=False))
