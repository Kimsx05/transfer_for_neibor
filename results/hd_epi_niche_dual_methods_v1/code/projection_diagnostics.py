from pathlib import Path
import numpy as np,pandas as pd,json,joblib
from sklearn.neighbors import NearestNeighbors
R=Path(__file__).resolve().parents[1];sections=['Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2'];types=json.loads((R/'type_mapping.json').read_text());rows=[]
for ref,cats in types.items():
 for radius in [50,200]:
  o=R/f'{ref}_r{radius}';z=np.load(o/'B_training.npz');P=z['P'];T=z['T'];y=z['labels'];d=pd.concat([pd.read_parquet(R/s/f'r{radius}.parquet') for s in sections],ignore_index=True).iloc[T]
  keep=((d.x_um//400+d.y_um//400).astype(int)%5!=0).to_numpy();idx=np.flatnonzero(~keep);nn=NearestNeighbors(n_neighbors=15,n_jobs=4).fit(P[keep]);dd,jj=nn.kneighbors(P[idx]);ty=y[keep];pred=[];vote=[]
  lookup={}
  for profile,label in zip(P[keep],ty):
   key=tuple(profile);lookup.setdefault(key,{});lookup[key][label]=lookup[key].get(label,0)+1
  for profile,ds,js in zip(P[idx],dd,jj):
   exact=lookup.get(tuple(profile))
   if exact:
    u=np.array(sorted(exact));c=np.array([exact[k] for k in u]);n=c.sum()
   else:
    u,c=np.unique(ty[js],return_counts=True);n=len(js)
   pred.append(sorted(zip(-c,u))[0][1]);vote.append(c.max()/n)
  rules={'A_max_distance':float(np.quantile(joblib.load(o/'A_model.joblib').transform(np.sqrt(P)).min(1),.995)*1.5+1e-8),'B_max_distance':float(np.quantile(dd[:,0],.995)*1.5+1e-8),'B_vote_strictly_above':.5,'projection_representation':{'A':'sqrt(P)','B':'P'},'tie_rule':'lexical among maximum vote labels; all exact-profile training ties when present; otherwise k15','native_labels_never_replaced':True}
  (o/'projection_rules.json').write_text(json.dumps(rules,indent=2));rows.append(dict(reference=ref,radius=radius,heldout_blocks_n=len(idx),B_label_transfer_agreement=float(np.mean(np.asarray(pred)==y[idx])),B_vote_median=float(np.median(vote)),B_vote_le_half=float(np.mean(np.asarray(vote)<=.5)),note='Existing native label transfer consistency, not independent biological accuracy'))
pd.DataFrame(rows).to_csv(R/'tables/projection_diagnostics.tsv',sep='\t',index=False)
