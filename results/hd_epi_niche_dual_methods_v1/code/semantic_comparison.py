from pathlib import Path
import pandas as pd,numpy as np
R=Path(__file__).resolve().parents[1];c=pd.read_csv(R/'tables/niche_composition.tsv',sep='\t',dtype={'niche':str});rows=[];labels={}
for ref in ['blca','hpca']:
 for radius in [50,200]:
  d=pd.read_parquet(R/f'{ref}_r{radius}/labels.parquet');z=c[(c.reference==ref)&(c.radius==radius)&(c.method=='A')&(c.niche!='unresolved')];high=z.loc[z.Epithelial.idxmax(),'niche'];d['semantic']=np.where(d.A_niche=='unresolved','unresolved',np.where(d.A_niche==high,'Higher_Epi_candidate_fraction','Epi_mesenchymal_composite'));labels[ref,radius]=d
  for s,g in d.groupby('section'):
   for label,n in g.semantic.value_counts().items():rows.append(dict(section=s,reference=ref,radius=radius,interpretation=label,n=n,fraction=n/len(g)))
pd.DataFrame(rows).to_csv(R/'tables/A_semantic_counts.tsv',sep='\t',index=False);pairs=[]
for kind in ['scale','reference']:
 for val in (['blca','hpca'] if kind=='scale' else [50,200]):
  a=labels[val,50] if kind=='scale' else labels['blca',val];b=labels[val,200] if kind=='scale' else labels['hpca',val]
  for s in a.section.unique():
   ok=(a.section==s)&a.evaluable&b.evaluable&(a.semantic!='unresolved')&(b.semantic!='unresolved');pairs.append(dict(comparison=kind,reference_or_radius=val,section=s,n=int(ok.sum()),semantic_agreement=float((a.loc[ok,'semantic']==b.loc[ok,'semantic']).mean())))
pd.DataFrame(pairs).to_csv(R/'tables/A_semantic_paired_agreement.tsv',sep='\t',index=False)
