from common import *
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
versions=['EXPANDED_LOWINFO','NO_CURRENT_CONFLICT']
sections=['Cervilla','D1_Johnson']
o=pd.read_parquet(ROOT/'cellcharter_assignments.parquet')
p=ROOT/'metadata/per_niche_sensitivity.tsv'
t=pd.read_csv(p,sep='\t')
K=json.loads((ROOT/'metadata/K_selection.json').read_text())['K']
labels=[f'niche{i}' for i in range(K)]
fig,axes=plt.subplots(2,2,figsize=(16,14),constrained_layout=True)
for row,v in enumerate(versions):
 for col,sec in enumerate(sections):
  s=o[o.section_id.eq(sec)&o.PRIMARY_included&o[v+'_included']]
  counts=pd.crosstab(s.PRIMARY_niche,s[v+'_niche']).reindex(index=labels,columns=labels,fill_value=0)
  fractions=counts.div(counts.sum(axis=1).replace(0,np.nan),axis=0)
  for niche in counts.index[counts.sum(axis=1).gt(0)]:
   dest=fractions.loc[niche].idxmax()
   mask=t.version.eq(v)&t.section_id.eq(sec)&t.niche.eq(niche)
   assert mask.sum()==1
   t.loc[mask,'largest_destination']=dest
   t.loc[mask,'largest_destination_fraction']=fractions.loc[niche,dest]
  counts.to_csv(ROOT/f'metadata/{v}_{sec}_mapped_crosstab.tsv',sep='\t')
  ax=axes[row,col]
  im=ax.imshow(fractions.to_numpy(),vmin=0,vmax=1,cmap='Blues')
  ax.set_xticks(range(K),labels,rotation=90)
  ax.set_yticks(range(K),[f'{n} (n={counts.loc[n].sum():,})' for n in labels])
  ax.set_xlabel('Sensitivity: joint-matched display label')
  ax.set_ylabel('PRIMARY: shared cells only')
  ax.set_title(f'{sec} | {v}')
  for i in range(K):
   for j in range(K):
    value=fractions.iloc[i,j]
    if pd.notna(value) and value>=.1:
     ax.text(j,i,f'{value:.0%}',ha='center',va='center',fontsize=8,color='white' if value>.55 else 'black')
fig.colorbar(im,ax=axes.ravel().tolist(),label='Fraction of each PRIMARY niche',shrink=.6)
fig.savefig(ROOT/'figures/sensitivity_transitions.png',dpi=160,transparent=True)
plt.close(fig)
t.to_csv(p,sep='\t',index=False)
print('SENSITIVITY_REVIEW_COMPLETE',len(t))
