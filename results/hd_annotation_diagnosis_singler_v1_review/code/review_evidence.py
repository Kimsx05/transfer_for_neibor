from pathlib import Path
import pandas as pd,numpy as np,json
root=Path(__file__).resolve().parents[1];old=root.parent/'hd_major_annotation_v1'
for s in ['Cervilla','D1_Johnson']:
 o=root/s;d=pd.read_csv(o/'diagnosis.tsv.gz',sep='\t',low_memory=False).set_index('cell_id');e=pd.read_csv(old/s/'tables/marker_evidence.tsv.gz',sep='\t').set_index('cell_id');assert d.index.equals(e.index)
 u=d.final_major.eq('Unknown');b=d.singler_blca_pruned_labels;h=d.singler_hpca_pruned_labels
 groups={'BLCA_Fib_HPCA_Epi':b.eq('Fibroblast')&h.eq('Epithelial'),'Both_Epi':b.eq('Epithelial')&h.eq('Epithelial'),'Both_Fib':b.eq('Fibroblast')&h.eq('Fibroblast'),'BLCA_Epi_HPCA_Fib':b.eq('Epithelial')&h.eq('Fibroblast'),'BLCA_Mural':b.eq('Mural'),'BLCA_Plasma':b.eq('Plasma'),'BLCA_B':b.eq('B'),'BLCA_T_NK':b.eq('T/NK')}
 rows=[];cl=[];ro=[];roi=pd.read_csv(old/s/'spatial_assets/rois.tsv',sep='\t');tr=json.loads((old/s/'spatial_assets/transform.json').read_text());crop=tr['fullres_crop_xyxy'];sc=tr['local_to_model_scale_xy'];mx=(d.nucleus_x_fullres-crop[0]+.5)*sc[0]-.5;my=(d.nucleus_y_fullres-crop[1]+.5)*sc[1]-.5
 for group,mask in groups.items():
  hit=u&mask;g=d[hit];eg=e[hit]
  if not len(g):continue
  both=(eg.detected_Epithelial>=2)&(eg.detected_Fibroblast>=2)
  rows.append(dict(group=group,n=len(g),UMI=int(g.nCount_RNA.sum()),Epi_ge2=int((eg.detected_Epithelial>=2).sum()),Fib_ge2=int((eg.detected_Fibroblast>=2).sum()),Epi_Fib_both_ge2=int(both.sum()),Epi_Fib_both_fraction=both.mean(),original_conflict=int(g.diag_fail_conflict.sum()),median_nFeature=g.nFeature_RNA.median()))
  cc=g.groupby('cluster_projected').agg(n=('object_id','size'),UMI=('nCount_RNA','sum')).reset_index();cc['group']=group;cl.append(cc)
  for _,r in roi.iterrows():
   loc=(mx>=r.x)&(mx<r.X)&(my>=r.y)&(my<r.Y);ro.append(dict(group=group,roi=r.roi,n=int((hit&loc).sum()),ROI_all_objects=int(loc.sum()),ROI_original_Unknown=int((loc&u).sum()),UMI=int(d.loc[hit&loc,'nCount_RNA'].sum())))
 pd.DataFrame(rows).to_csv(o/'focused_cell_evidence.tsv',sep='\t',index=False);pd.concat(cl).to_csv(o/'focused_groups_by_cluster.tsv',sep='\t',index=False);pd.DataFrame(ro).to_csv(o/'focused_groups_by_ROI.tsv',sep='\t',index=False)
 rs=[]
 for _,bnd in roi.iterrows():
  loc=(mx>=bnd.x)&(mx<bnd.X)&(my>=bnd.y)&(my<bnd.Y);zz=d[loc]
  rs.append(dict(roi=bnd.roi,n=len(zz),median_nFeature=zz.nFeature_RNA.median(),median_nCount=zz.nCount_RNA.median(),BLCA_top3=json.dumps(zz.singler_blca_pruned_labels.value_counts().head(3).to_dict()),HPCA_top3=json.dumps(zz.singler_hpca_pruned_labels.value_counts().head(3).to_dict())))
 pd.DataFrame(rs).to_csv(o/'roi_information_candidates.tsv',sep='\t',index=False)
 print('FOCUSED_CELL_EVIDENCE',s)
# Reuse original filtered cluster marker tables; absence is not zero expression.
rows=[];availability=[]
for s in ['Cervilla','D1_Johnson']:
 g=pd.read_csv(root/s/'focused_groups_by_cluster.tsv',sep='\t');cs=g[g.group.eq('BLCA_Fib_HPCA_Epi')].sort_values('n',ascending=False).head(5).cluster_projected
 m=pd.read_csv(old/s/'tables/cluster_markers.tsv.gz',sep='\t')
 for k in cs:
  z=m[m.cluster.eq(k)].sort_values('log2FC',ascending=False).head(12).copy();availability.append(dict(section=s,cluster=k,available_filtered_markers=len(z),meaning='original pct>=0.1, pct_difference>=0.05, log2FC>0, top50; empty does not imply absent expression'));z.insert(0,'section',s);rows.append(z)
pd.concat(rows).to_csv(root/'tables/focused_original_cluster_markers.tsv',sep='\t',index=False);pd.DataFrame(availability).to_csv(root/'tables/focused_cluster_marker_availability.tsv',sep='\t',index=False)
