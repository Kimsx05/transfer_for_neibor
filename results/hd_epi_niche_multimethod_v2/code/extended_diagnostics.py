from common import *
import joblib
from sklearn.metrics import adjusted_rand_score as ari,adjusted_mutual_info_score as ami

def run():
 a=anchors();T=train_index(a);index=pd.read_csv(R/'tables/assignment_index.tsv',sep='\t');duplicates=[];sizes=[];ys={};coarse={};cross=[]
 for _,row in index.iterrows():
  o=R/Path(row.path).parent;z=pd.read_parquet(R/row.path);j=json.loads((o/'complete.json').read_text());tag=row.tag;ys[tag]=z.label.to_numpy();y=ys[tag][T];r=int(j['radius'])
  if j['method']=='CellCharter':
   x=np.concatenate([np.load(R/'models'/f'{s}_expression_neighbors_r{r}.npz')['X'] for s in SECTIONS]);x=x[T,:10] if j['representation']=='scVI' else x[T,10:]
  elif (o/'training.npz').exists():x=np.load(o/'training.npz')['X']
  elif j['representation']=='WEAK':x=np.load(R/'models'/f'harmony_{j["reference"]}_r{r}.npz')['X'][T]
  else:
   _,P,_=composition(j['reference'],r);pre=joblib.load(o/'preprocessor.joblib');x=pre['PCA'].transform(np.sqrt(P[T]))[:,:pre['rank']]
  u,inv,n=np.unique(x,axis=0,return_inverse=True,return_counts=True);lo=np.full(len(u),y.max()+1);hi=np.full(len(u),-1);np.minimum.at(lo,inv,y);np.maximum.at(hi,inv,y);split=lo!=hi;duplicates.append(dict(tag=tag,actual_input_rows=len(x),unique_full_inputs=len(u),rows_with_duplicate_input=int(n[n>1].sum()),identical_full_inputs_split=int(split.sum()),rows_in_split_inputs=int(n[split].sum()),input_dimension=x.shape[1]))
  for (s,lab),q in z.groupby(['section','label']):sizes.append(dict(tag=tag,section=s,label=int(lab),n=len(q),fraction_section_anchors=len(q)/int((a.section==s).sum()),native_train=int((q.source=='native_train').sum()),assigned=int((q.label>=0).sum())))
  if (o/'K2.npz').exists():
   f=np.load(o/'K2.npz');v=np.full(len(a),-1);v[f['S']]=f['labels'];coarse[tag]=v
 pd.DataFrame(duplicates).to_csv(R/'tables/identical_input_diagnostics.tsv',sep='\t',index=False);pd.DataFrame(sizes).to_csv(R/'tables/niche_sizes_by_section.tsv',sep='\t',index=False)
 comparisons=pd.read_csv(R/'tables/method_comparisons.tsv',sep='\t')
 for (left,right),q in comparisons.groupby(['left','right']):
  x,y=ys[left],ys[right]
  for scope,ix in [('common_native_T',T),('common_assigned',np.flatnonzero((x>=0)&(y>=0)))]:
   ct=pd.crosstab(x[ix],y[ix]);ct=ct.stack();ct=ct[ct>0]
   for (ll,rr),n in ct.items():cross.append(dict(left=left,right=right,scope=scope,left_label=int(ll),right_label=int(rr),n=int(n),denominator=len(ix)))
 pd.DataFrame(cross).to_csv(R/'tables/partition_cross_tables.tsv',sep='\t',index=False)
 k2=[]
 for tag,y in coarse.items():k2.append(dict(left=tag+'_K2',right=tag+'_working',scope='common_native_T',n=len(T),ARI=ari(y[T],ys[tag][T]),AMI=ami(y[T],ys[tag][T])))
 for ref in TYPES:
  for r in [50,200]:
   aa=f'hellinger_kmeans_{ref}_r{r}_NONE';bb=f'CellCharter_both_r{r}_scVI'
   if aa in coarse and bb in coarse:k2.append(dict(left=aa+'_K2',right=bb+'_K2',scope='common_native_T',n=len(T),ARI=ari(coarse[aa][T],coarse[bb][T]),AMI=ami(coarse[aa][T],coarse[bb][T])))
 pd.DataFrame(k2).to_csv(R/'tables/common_K2_comparisons.tsv',sep='\t',index=False)
 # The common four-type layer is explanatory only; missing fine types are not zero evidence.
 mapping={'Epithelial':'Epithelial','Fibroblast':'Mesenchymal','Mural':'Mesenchymal','Endothelial':'Endothelial','T/NK':'Immune','Myeloid':'Immune','B':'Immune','Plasma':'Immune','Mast':'Immune'}
 cm=pd.read_csv(R/'tables/niche_composition.tsv',sep='\t');four=[]
 for _,v in cm.iterrows():
  row={k:v[k] for k in ['tag','reference','radius','label','n']}
  for major in ['Epithelial','Mesenchymal','Endothelial','Immune']:row[major]=sum(v[t] for t in TYPES[v.reference] if mapping[t]==major)
  four.append(row)
 pd.DataFrame(four).to_csv(R/'tables/niche_common_four_composition.tsv',sep='\t',index=False)
 fits=[]
 for variant in ['scVI','scVI_without_Top20']:
  o=R/'models/expression'/variant;j=json.loads((o/'complete.json').read_text());h=pd.read_csv(o/'reference_validation_elbo.tsv',sep='\t').iloc[:,1];j['reference_validation_ELBO_final']=float(h.iloc[-1]);j['reference_validation_ELBO_last10_relative_change']=float((h.iloc[-1]-h.iloc[max(0,len(h)-10)])/abs(h.iloc[max(0,len(h)-10)]));j['reference_epochs']=len(h);fits.append(j)
 pd.DataFrame(fits).to_csv(R/'tables/expression_fit_diagnostics.tsv',sep='\t',index=False)
 print('EXTENDED_DIAGNOSTICS_COMPLETE',flush=True)
if __name__=='__main__':run()
