from base import *
d=O/'hd/Cervilla';keys=['from','to'];cols=['n_from','n_to','observed','expected','z_score','observed_expected']
v=[]
for n in [8,16]:
 t=pd.read_csv(d/f'scale{n}_major_50_matrix.tsv.gz',sep='\t')
 t=t[t['from']<t.to][keys+cols].set_index(keys).add_suffix('_'+str(n));v.append(t)
z=v[0].join(v[1],how='outer');known=~z.index.to_frame(index=False).isin(['Unknown','Mixed']).any(axis=1).to_numpy()
valid=known&(z[['n_from_8','n_to_8','n_from_16','n_to_16']].min(axis=1,skipna=False)>=10)&z[['z_score_8','z_score_16']].notna().all(axis=1)
z['direction']=np.where(valid,np.where(np.sign(z.z_score_8)==np.sign(z.z_score_16),'same','different'),'insufficient')
put(z.reset_index(),d/'scale_direction_50.tsv')
print(z.direction.value_counts().to_dict())

