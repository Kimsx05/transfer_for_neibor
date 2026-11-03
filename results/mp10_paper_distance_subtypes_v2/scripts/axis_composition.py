from plots import *
rows=[]
for sec in M.section_id:
 d=pd.read_parquet(O/'spots'/f'{sec}.parquet')
 cols=[c for c in d if c.endswith('_fraction') and not c.startswith('fraction__')]+[k+'_within_lineage' for k in list(TARGETS)[:7]]
 for axis in ['local_interface']+[k+'_P95' for k in list(TARGETS)[:7]]:
  x=d.local_interface_um if axis=='local_interface' else d[axis+'_distance']
  b=pd.cut(x,[-np.inf,-500,-200,0,100,200,500,1000,np.inf])
  for band,g in d.groupby(b,observed=True):
   rows.append(dict(section_id=sec,axis=axis,bin=str(band),n=len(g),**g[cols].mean().to_dict()))
a=pd.DataFrame(rows).merge(M[['section_id','dataset_id','source_context']],on='section_id')
write(a,O/'tables/axis_composition.tsv.gz')
z=a[a.source_context.str.startswith('tumor')&a.axis.eq('local_interface')]
cols=[c for c in ['Epithelial_fraction','Fibroblast_fraction','Myeloid_fraction','T_NK_fraction','Endothelial_fraction','Mural_fraction','B_fraction','Plasma_fraction'] if c in z]
m=z.groupby(['dataset_id','bin'])[cols].mean().groupby('bin').mean()
heat(m,'Lineage composition / signed local interface / equal datasets','interface_composition',1)
