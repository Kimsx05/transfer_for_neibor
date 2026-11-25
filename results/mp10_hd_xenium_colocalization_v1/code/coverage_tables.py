from base import *
def coverage(s):
 hd=s in HD;d=O/('hd' if hd else 'xenium')/s;z=pd.read_parquet(d/'labels.parquet')
 tax=pd.read_csv(O/'shared/taxonomy.tsv',sep='\t').drop_duplicates('label')
 focus=['T/NK','Endothelial','Fibroblast','Myeloid'];maj=sorted(tax.major.unique());fine=sorted(tax.loc[tax.major.isin(focus),'label'])
 resolved=[a for a in maj if a not in focus]+fine
 vocab={'major':maj+['Unknown','Mixed'],'subtype':resolved+[a+' (unresolved)' for a in focus]+['Unknown','Mixed'],'chol':[a for a in resolved if a!='Epithelial']+['Chol-high Epi','other Epi']}
 out=[];raw=[]
 for r in ([50] if hd else [50,200]):
  for k in (['major','subtype'] if hd else ['major','subtype','chol']):
   p=d/f'{k}_{r}_matrix.tsv.gz'
   if not p.exists():continue
   t=pd.read_csv(p,sep='\t');nodes=t[['from','n_from']].drop_duplicates().set_index('from').n_from
   labels=sorted(set(vocab[k])|set(nodes.index));n=nodes.reindex(labels,fill_value=0);n.index.name='label'
   for lab in labels:out.append((s,k,r,lab,n[lab],n[lab]>0))
   ix=pd.MultiIndex.from_product([labels,labels],names=['from','to']);v=t.set_index(['from','to']).reindex(ix).reset_index();v['n_from']=v['from'].map(n);v['n_to']=v.to.map(n);v['estimability']=np.select([(v.n_from==0)|(v.n_to==0),v.null_sd.eq(0)],['absent class','constant permutation statistic'],default='estimated')
   v['matrix']=k;v['radius_um']=r;raw.append(v)
 put(pd.DataFrame(out,columns=['sample','matrix','radius_um','label','n_nodes','represented']),d/'matrix_nodes.tsv')
 pd.concat(raw,ignore_index=True).to_parquet(d/'matrix_full_taxonomy.parquet',index=False)
 # Report support against the original underlying marker panels including zero-observation classes.
 rows=[]
 for t in tax.itertuples():
  if t.major not in focus:continue
  hit=z.supported==t.label;candidate=z.candidate==t.label
  rows.append((t.label,t.major,int(candidate.sum()),int(hit.sum()),int(z.major.eq(t.major).sum()),int((z.major.eq(t.major)&z.flag.eq('major_only')).sum()),z.loc[hit,'support'].median(),z.loc[hit,'marker_detected'].median()))
 put(pd.DataFrame(rows,columns=['subtype','major','n_candidate','n_supported','n_major','n_parent_only','median_supported_weight','median_detected_markers']),d/'subtype_support.tsv')
if __name__=='__main__':
 for s in sys.argv[1:]:coverage(s)

