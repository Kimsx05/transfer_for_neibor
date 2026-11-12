from pathlib import Path
import numpy as np,pandas as pd,json
R=Path(__file__).resolve().parents[1];sections=['Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2'];types=json.loads((R/'type_mapping.json').read_text());rows=[]
for ref,cats in types.items():
 for radius in [50,200]:
  o=R/f'{ref}_r{radius}';z=np.load(o/'B_training.npz');p=z['P'];y=z['labels'];T=z['T'];q=pd.concat([pd.read_parquet(R/s/f'r{radius}.parquet') for s in sections],ignore_index=True);P=q[[ref+'_P_'+t for t in cats]].to_numpy();lab=pd.read_parquet(o/'labels.parquet');native=lab.B_assignment_source=='native_train'
  # All exact-profile ties are used, never an arbitrary subset of 15 identical neighbors.
  u,inv,cnt=np.unique(p,axis=0,return_inverse=True,return_counts=True);tab=pd.crosstab(inv,y);distinct=(tab>0).sum(1);votes=tab.max(1)/tab.sum(1);pred=tab.idxmax(1);lookup={tuple(v):i for i,v in enumerate(u)};ix=np.array([lookup.get(tuple(v),-1) for v in P]);exact=ix>=0;nontrain=exact&~native.to_numpy()&lab.evaluable.to_numpy()
  lab['B_exact_profile_n_labels']=np.where(exact,distinct.to_numpy()[np.maximum(ix,0)],0);lab['B_exact_profile_training_n']=np.where(exact,cnt[np.maximum(ix,0)],0)
  lab.loc[nontrain,'B_candidate']=pred.to_numpy()[ix[nontrain]];lab.loc[nontrain,'B_support']=votes.to_numpy()[ix[nontrain]];lab.loc[nontrain,'B_distance']=0
  ok=nontrain&(lab.B_support>.5);bad=nontrain&~ok;lab.loc[ok,'B_niche']=lab.loc[ok,'B_candidate'];lab.loc[ok,'B_assignment_source']='external_projection';lab.loc[ok,'B_reason']='exact_profile_all_ties_majority';lab.loc[bad,'B_niche']='unresolved';lab.loc[bad,'B_assignment_source']='unresolved';lab.loc[bad,'B_reason']='identical_composition_split_native_labels'
  assert np.array_equal(lab.loc[native,'B_niche'].to_numpy(),y)
  lab['squidpy_niche_train']=lab.B_niche.where(native)
  lab['squidpy_niche_projected']=lab.B_candidate.where(~native&lab.evaluable)
  lab.to_parquet(o/'labels.parquet',index=False)
  cover=[]
  for method in ['A','B']:
   for sec,g in lab.groupby('section'):
    for state,n in g[method+'_assignment_source'].value_counts().items():cover.append(dict(section=sec,method=method,state=state,n=n,anchors=len(g),evaluable=int(g.evaluable.sum())))
  pd.DataFrame(cover).to_csv(o/'coverage.tsv',sep='\t',index=False)
  rows.append(dict(reference=ref,radius=radius,training_profiles=len(u),profiles_split_across_native_labels=int((distinct>1).sum()),training_bins_in_split_profiles=int(cnt[distinct.to_numpy()>1].sum()),pure_epi_train=int((p[:,0]==1).sum()),pure_epi_native_labels=len(set(y[p[:,0]==1])),external_exact_ambiguous_unresolved=int(bad.sum()),native_labels_preserved=True))
pd.DataFrame(rows).to_csv(R/'tables/identical_composition_diagnostic.tsv',sep='\t',index=False)
