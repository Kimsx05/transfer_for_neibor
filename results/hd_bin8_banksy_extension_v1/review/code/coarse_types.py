from pathlib import Path
import pandas as pd,numpy as np,json,sys
from scipy.signal import convolve2d
from PIL import Image
root=Path(__file__).resolve().parents[1];base=root;hd=root.parents[1];s=sys.argv[1];o=root/s
mp={'Epithelial':'Epithelial','Fibroblast':'Mesenchymal','Mural':'Mesenchymal','Endothelial':'Endothelial',**{k:'Immune' for k in ['T/NK','Myeloid','B','Plasma','Mast']}};pd.DataFrame(list(mp.items()),columns=['original','coarse']).to_csv(root/'tables/coarse_mapping.tsv',sep='\t',index=False)
d=pd.read_csv(base/s/'analysis_metadata.tsv.gz',sep='\t',low_memory=False);z=pd.DataFrame(dict(bin_id=d.bin_id));z['banksy_information_band']=pd.cut(d.nFeature_RNA,[-np.inf,19,99,499,np.inf],labels=['<20','20-99','100-499','>=500']);summary=[];composition=[];cross=[]
for ref in ['blca','hpca']:
 for stage,suffix in [('raw','labels'),('pruned','pruned_labels')]:
  col=f'singler_{ref}_{suffix}';new=f'coarse_{ref}'+('_raw' if stage=='raw' else '');unknown=d[col].isna();assert d.loc[~unknown,col].isin(mp).all();z[new]=d[col].map(mp).fillna('Unknown');c=d.assign(label=z[new]).groupby('label',dropna=False).agg(n=('bin_id','size'),UMI=('nCount_RNA','sum')).reset_index();c['reference']=ref;c['stage']=stage;c['fraction']=c.n/len(d);c['UMI_fraction']=c.UMI/d.nCount_RNA.sum();composition.append(c)
b=d.singler_blca_pruned_labels;h=d.singler_hpca_pruned_labels;both=b.notna()&h.notna();covered=both&b.isin(['Epithelial','Fibroblast','T/NK','Myeloid','B','Endothelial']);fine=b.fillna('Unknown')!=h.fillna('Unknown');coarse=z.coarse_blca!=z.coarse_hpca
for name,mask in [('all_QC_bins',np.ones(len(d),bool)),('both_candidates',both),('both_candidates_original_HPCA_coverage',covered)]:
 for level,dis in [('original_major',fine),('four_coarse',coarse)]:
  hit=mask&dis;summary.append(dict(denominator=name,level=level,n=int(mask.sum()),UMI=int(d.loc[mask,'nCount_RNA'].sum()),disagree_n=int(hit.sum()),disagree_fraction=float(hit.sum()/mask.sum()),disagree_UMI=int(d.loc[hit,'nCount_RNA'].sum()),disagree_UMI_fraction=float(d.loc[hit,'nCount_RNA'].sum()/d.loc[mask,'nCount_RNA'].sum())))
 zz=d.loc[mask,['bin_id','nCount_RNA']].assign(BLCA=z.coarse_blca[mask],HPCA=z.coarse_hpca[mask]);cc=zz.groupby(['BLCA','HPCA']).agg(n=('bin_id','size'),UMI=('nCount_RNA','sum')).reset_index();cc['denominator']=name;cc['fraction']=cc.n/mask.sum();cc['UMI_fraction']=cc.UMI/zz.nCount_RNA.sum();cross.append(cc)
fwd=(b=='Fibroblast')&(h=='Epithelial');rev=(b=='Epithelial')&(h=='Fibroblast');z['reference_conflict_group']=np.select([fwd,rev,(b=='Epithelial')&(h=='Epithelial'),(z.coarse_blca=='Mesenchymal')&(z.coarse_hpca=='Mesenchymal'),(z.coarse_blca=='Immune')&(z.coarse_hpca=='Immune'),~both],['BLCA_Fibroblast_HPCA_Epithelial','BLCA_Epithelial_HPCA_Fibroblast','Agree_Epithelial','Agree_Mesenchymal','Agree_Immune','Not_both_evaluable'],default='Other');z.to_csv(o/'type_metadata.tsv.gz',sep='\t',index=False);pd.DataFrame(summary).to_csv(o/'coarsening_comparison.tsv',sep='\t',index=False);pd.concat(composition).to_csv(o/'coarse_composition.tsv',sep='\t',index=False);pd.concat(cross).to_csv(o/'coarse_crosstabs.tsv',sep='\t',index=False)

print(s,'TYPES_COMPLETE',len(d))
