from pathlib import Path
import numpy as np,pandas as pd,json
from scipy.spatial import cKDTree
from scipy import ndimage
R=Path(__file__).resolve().parents[1];A=R.parent;E=A/'hd_bin8_banksy_extension_v1'
sections=['Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2']
types={'blca':['Epithelial','T/NK','Myeloid','Fibroblast','Endothelial','Mast','B','Plasma','Mural'],'hpca':['Epithelial','T/NK','Myeloid','Fibroblast','Endothelial','B']}
(R/'type_mapping.json').write_text(json.dumps(types,indent=2))
for s in sections:
 o=R/s;o.mkdir(exist_ok=True)
 if (o/'r200.parquet').exists():continue
 base=A/'hd_bin8_comparison_v1' if s in sections[:2] else E
 d=pd.read_csv(base/s/'bin_metadata.tsv.gz',sep='\t');d['section']=s
 q=pd.read_csv(E/s/'comparison_metadata.tsv.gz',sep='\t');m=pd.read_csv(base/s/'analysis_metadata.tsv.gz',sep='\t');assert q.bin_id.is_unique and m.bin_id.is_unique
 d=d.merge(q.drop(columns=['nCount_RNA','nFeature_RNA']),on='bin_id',how='left',validate='one_to_one')
 d['analysis_qc']=d.bin_id.isin(m.bin_id);d['supplied']=True
 d['in_tissue_original']=d['in_tissue'] if 'in_tissue' in d else d['in_tissue_original']
 d['epi_marker_detected_n']=(d.ext_Keratin_extended_detected+d.ext_Urothelial_controls_detected).fillna(0).astype(int)
 d['epi_marker_supported']=d.epi_marker_detected_n>=2
 for ref in types:
  d['type_'+ref]=d['singler_'+ref+'_pruned_labels'].where(d['singler_'+ref+'_pruned_labels'].isin(types[ref]))
 d['epi_candidate_any']=(d.type_blca=='Epithelial')|(d.type_hpca=='Epithelial')
 d['epi_anchor']=d.analysis_qc&d.epi_candidate_any&d.epi_marker_supported
 d['epi_mixed_flag']=d.epi_marker_supported&(d.ext_Collagen_stroma_detected>=2)
 d['reference_conflict']=d.coarse_blca.ne(d.coarse_hpca)&d.type_blca.notna()&d.type_hpca.notna()
 shape=tuple((d[['array_row','array_col']].max()+1).astype(int));grid=np.zeros(shape,bool)
 support=d if s.startswith('MPBC') else d[d.in_tissue_original==1]
 grid[support.array_row,support.array_col]=True;lab,n=ndimage.label(grid,np.ones((3,3)))
 d['component']=lab[d.array_row,d.array_col]
 # Coordinates already verified by upstream; all original supplied positions retained.
 d.to_parquet(o/'nodes.parquet',index=False)
 anchor=d[d.epi_anchor].copy();anchor.to_parquet(o/'anchor_manifest.parquet',index=False)
 for radius in [50,200]:
  if (o/f'r{radius}.parquet').exists():continue
  out=anchor[['bin_id','section','x_um','y_um','component','nFeature_RNA','nCount_RNA','epi_mixed_flag']].copy().set_index('bin_id')
  offsets=np.arange(-int(radius//8),int(radius//8)+1);theory=int(sum((i*8)**2+(j*8)**2<=radius**2+1e-7 for i in offsets for j in offsets))-1
  out['theoretical_grid_neighbors']=theory
  for comp,ac in anchor.groupby('component'):
   assert comp>0
   full=d[d.component==comp];pool=full[full.analysis_qc];xy=ac[['x_um','y_um']].to_numpy();ids=ac.bin_id
   def count(z):
    return cKDTree(z[['x_um','y_um']].to_numpy()).query_ball_point(xy,radius+1e-7,return_length=True,workers=4) if len(z) else np.zeros(len(ac),int)
   ns=count(full)-1;nq=count(pool)-1
   out.loc[ids,'n_supplied_component']=ns;out.loc[ids,'n_analysis_neighbors']=nq
   out.loc[ids,'low_information_fraction']=np.divide(count(pool[pool.nFeature_RNA<20])-(ac.nFeature_RNA<20).to_numpy(),nq,out=np.full(len(ac),np.nan),where=nq>0)
   out.loc[ids,'reference_conflict_fraction']=np.divide(count(pool[pool.reference_conflict])-ac.reference_conflict.to_numpy(),nq,out=np.full(len(ac),np.nan),where=nq>0)
   for ref,cats in types.items():
    cs=[]
    for typ in cats:
     v=count(pool[pool['type_'+ref]==typ])-(ac['type_'+ref]==typ).to_numpy();out.loc[ids,ref+'_C_'+typ]=v;cs.append(v)
    total=np.sum(cs,axis=0);out.loc[ids,ref+'_n_known']=total
    out.loc[ids,ref+'_known_fraction']=np.divide(total,nq,out=np.full(len(ac),np.nan),where=nq>0)
  out['input_coverage_component']=out.n_supplied_component/theory
  for ref,cats in types.items():
   for typ in cats:out[ref+'_P_'+typ]=out[ref+'_C_'+typ]/out[ref+'_n_known'].replace(0,np.nan)
  out.reset_index().to_parquet(o/f'r{radius}.parquet',index=False)
  print(s,radius,len(out),'anchors',out.filter(regex='n_known').describe().to_dict(),flush=True)
print('PREPARE_COMPLETE',flush=True)
