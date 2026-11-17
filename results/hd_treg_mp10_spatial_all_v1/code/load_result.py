from pathlib import Path
import pandas as pd,json
ROOT=Path(__file__).resolve().parents[1]
def load_section(section,root=ROOT,include_members=False,include_pairs=False):
 root=Path(root);p=root/section
 if section not in ['Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2']:raise ValueError(section)
 out={k:pd.read_parquet(p/f'{k}_regions.parquet')for k in ['T32','T16','Epi64']}
 out['inputs']=json.loads((p/'prepared.json').read_text())
 if include_members:out['members']={k:pd.read_parquet(p/f'{k}_members.parquet')for k in ['T32','T16','Epi64']}
 if include_pairs:out['pairs']=pd.read_parquet(p/'region_pairs.parquet')
 for k in ['T32','T16','Epi64']:assert out[k].region_id.is_unique
 return out

def bin_membership(section,root=ROOT):
 p=Path(root)/section;b=pd.read_parquet(p/'bin_evidence.parquet',columns=['bin_id','x_um','y_um','in_tissue_original','eligible_pool','epi_support'])
 for k in ['T32','T16','Epi64']:
  m=pd.read_parquet(p/f'{k}_members.parquet').rename(columns={'region_id':'trmp_'+k+'_region_id'});assert m.bin_id.is_unique;b=b.merge(m,on='bin_id',how='left',validate='one_to_one')
 return b

def attach_membership(adata,section,root=ROOT):
 d=bin_membership(section,root).set_index('bin_id');ids=pd.Index(adata.obs_names);assert ids.is_unique and ids.isin(d.index).all()
 for c in d.filter(like='trmp_').columns:
  if c in adata.obs:raise ValueError('Existing field would be overwritten: '+c)
  adata.obs[c]=d.loc[ids,c].to_numpy()
 return adata

if __name__=='__main__':
 import argparse
 a=argparse.ArgumentParser();a.add_argument('section');a.add_argument('kind',choices=['T32','T16','Epi64','membership']);a.add_argument('output');z=a.parse_args();d=bin_membership(z.section)if z.kind=='membership'else load_section(z.section)[z.kind];d.to_csv(z.output,sep='\t',index=False,na_rep='NA')
