from pathlib import Path
import pandas as pd

def load_niche_metadata(root,section):
 root=Path(root);m=pd.read_parquet(root/section/'anchor_manifest.parquet');out=m[['bin_id','epi_anchor','epi_candidate_any','epi_marker_supported','epi_marker_detected_n','epi_mixed_flag']].set_index('bin_id')
 for ref in ['blca','hpca']:
  for r in [50,200]:
   p=root/f'{ref}_r{r}'/'labels.parquet'
   if not p.exists():continue
   d=pd.read_parquet(p);d=d[d.section==section].set_index('bin_id');assert d.index.is_unique and d.index.equals(out.index)
   out=out.join(d.drop(columns='section').add_prefix(f'niche_{ref}_r{r}_'))
 return out

def attach_niches(obj,root,section):
 """AnnData: original X/layers untouched; non-anchor labels remain missing."""
 m=load_niche_metadata(root,section);assert obj.obs_names.is_unique and m.index.isin(obj.obs_names).all()
 assert not set(m.columns)&set(obj.obs.columns),'Will not overwrite existing metadata'
 obj.obs=obj.obs.join(m);obj.obs['epi_anchor']=obj.obs.epi_anchor.fillna(False)
 return obj
if __name__=='__main__':
 import sys
 load_niche_metadata(sys.argv[1],sys.argv[2]).rename_axis('bin_id').reset_index().to_csv(sys.argv[3],sep='\t',index=False)
