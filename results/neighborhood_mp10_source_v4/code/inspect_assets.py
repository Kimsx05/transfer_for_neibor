from base import *
import h5py,anndata as ad
p=ROOT/'cell2location_full_run/v2_30type/02_sections/CanDis_2026_Sample1'
a=ad.read_h5ad(p/'spatial_posterior.h5ad',backed='r');print('adata',a.shape,'layers',list(a.layers),'uns',list(a.uns),'obs',list(a.obs.columns),'var',list(a.var.columns));print('unsmod',[(k,type(v).__name__) for k,v in a.uns.items()]);
for k in ['mod','model','cell2location']:
 if k in a.uns:
  print(k,list(a.uns[k]));
  for kk,v in a.uns[k].items():print(kk,list(v) if isinstance(v,dict) else getattr(v,'shape',str(v)[:300]))
a.file.close()
import torch
v=torch.load(p/'model/model.pt',map_location='cpu',weights_only=False);print('model_keys',list(v));d=v.get('attr_dict',{});print('attrs',list(d));
for k in ['samples','samples_','factor_names_','cell_state_df_','adata_manager']:
 if k in d:
  q=d[k];print(k,list(q) if isinstance(q,dict) else getattr(q,'shape',str(q)[:400]));
  if isinstance(q,dict):print('subkeys',[(kk,list(vv) if isinstance(vv,dict) else getattr(vv,'shape',type(vv).__name__)) for kk,vv in q.items()])
print('statekeys',list(v.get('model_state_dict',{}))[:12]);print('registry',d.get('registry_',{}))
