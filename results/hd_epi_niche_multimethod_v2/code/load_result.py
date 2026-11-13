from pathlib import Path
import pandas as pd

def load_niche_metadata(root,section):
 root=Path(root);a=pd.read_parquet(root/'tables/frozen_anchors.parquet');a=a[a.section==section].set_index('bin_id');out=pd.DataFrame({'niche_v2_anchor':True,'niche_v2_x_um':a.x_um,'niche_v2_y_um':a.y_um},index=a.index)
 index=pd.read_csv(root/'tables/assignment_index.tsv',sep='\t')
 for _,v in index.iterrows():
  d=pd.read_parquet(root/v.path);d=d[d.section==section].set_index('bin_id');assert d.index.equals(a.index)
  cols=['candidate','label','source','reason','confidence','distance','evaluable'];out=out.join(d[cols].add_prefix('niche_v2_'+v.tag+'_'))
 return out

def attach_niches(obj,root,section):
 d=load_niche_metadata(root,section);assert obj.obs_names.is_unique and d.index.isin(obj.obs_names).all();assert not set(d.columns)&set(obj.obs.columns),'No existing field is overwritten'
 obj.obs=obj.obs.join(d);obj.obs['niche_v2_anchor']=obj.obs.niche_v2_anchor.fillna(False).astype(bool);return obj
if __name__=='__main__':
 import sys
 load_niche_metadata(sys.argv[1],sys.argv[2]).rename_axis('bin_id').reset_index().to_csv(sys.argv[3],sep='\t',index=False)
