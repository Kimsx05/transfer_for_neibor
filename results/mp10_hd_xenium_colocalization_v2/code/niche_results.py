from common import *
from spatial import groups,shuffle_effect,boot
from scipy.spatial.distance import jensenshannon
from sklearn.metrics import adjusted_rand_score
def collect(s):
 v,d=folders(s);p=OUT/'niche'/s;ct=pd.read_csv(p/'niche_input.tsv.gz',sep='\t').set_index('id').ct;z=labels(s);idt=pd.read_parquet(d/'identity.parquet');res=16 if s in HD else 'cell';ms=['MP10'] if s in HD else STATES
 e=pd.read_parquet(d/f'centers_{res}.parquet');n=pd.read_parquet(d/f'neighborhood_{res}_50.parquet').set_index('id').loc[e.index];depth=z['nCount' if s in HD else 'transcript_counts']
 imm=idt.grade.eq('supported');agg=[];comp=[];metrics=[];idx_niche={}
 for name,fp in [('MENDER',p/'MENDER.parquet'),('BuildNicheAssay',p/'BuildNicheAssay.tsv.gz'),('CellCharter all',p/'CellCharter_all.parquet'),('CellCharter masked',p/'CellCharter_masked.parquet')]:
  if not fp.exists():continue
  nc=(pd.read_csv(fp,sep='\t').set_index('id') if name=='BuildNicheAssay' else pd.read_parquet(fp)).niche.reindex(z.index);assert nc.notna().all()
  idx_niche[name]=nc
  for k in sorted(set(nc)):
   sl=nc==k;ec=sl.loc[e.index].values;qc=n.valid.to_numpy();yi=n.supported.values/e.area_50.values*1e6;yt=n.supported.values/np.where(n.reliable_T>0,n.reliable_T,np.nan)
   for lab,ct0 in ct[sl].value_counts().items():comp.append((s,name,k,lab,ct0,ct0/sl.sum()))
   for m in ms:
    h=e[m+'_high'].to_numpy(bool);valid=ec&qc
    for end,a in [('density',yi),('Treg/T',yt)]:
     hh=valid&h&np.isfinite(a);ll=valid&~h&np.isfinite(a)
     agg.append(dict(sample=s,method=name,niche=k,state=m,endpoint=end,units=sl.sum(),epi=ec.sum(),high=(ec&h).sum(),high_fraction=h[ec].mean() if ec.any() else np.nan,mean_score=e.loc[ec,m].mean(),valid_high=hh.sum(),valid_low=ll.sum(),mean_high=np.mean(a[hh]) if hh.any() else np.nan,mean_low=np.mean(a[ll]) if ll.any() else np.nan,targets=imm[sl].sum(),reliable_T=idt.reliable_T[sl].sum(),mean_depth=depth[sl].mean(),unknown_mixed=ct[sl].isin(['Unknown','Mixed']).mean(),mean_coverage=n.coverage[ec].mean()))
  if name=='MENDER':
   for k in sorted(set(nc)):
    take=nc.loc[e.index].eq(k).values
    if not take.any():continue
    y=e.loc[take];hh=y[[m+'_high' for m in ms]].to_numpy(bool);vv=np.column_stack([np.where(n.valid,yi,np.nan)[take],np.where(n.valid,yt,np.nan)[take]])
    pr=shuffle_effect(hh,np.ascontiguousarray(vv),groups(y,400),4999,CFG['seed']);bb=groups(y,400,False)
    for j,m in enumerate(ms):
     for c,end in enumerate(['density','Treg/T']):
      nu=pr[1:,j,c];ob=pr[0,j,c];mn=np.mean(nu);sd=np.std(nu);nH=(hh[:,j]&np.isfinite(vv[:,c])).sum();nL=(~hh[:,j]&np.isfinite(vv[:,c])).sum();ok=(sd>1e-10 and np.isfinite(ob) and nH>=5 and nL>=5);lo,hi=boot(vv[:,c],hh[:,j],bb) if ok else [np.nan]*2
      metrics.append(dict(sample=s,platform='HD' if s in HD else 'Xenium',state_definition=m,unit_resolution=res,identity_definition='supported',radius_or_contact_um=50,endpoint=end,niche=k,effect=ob,CI_low=lo,CI_high=hi,null_mean=mn,null_sd=sd,P=(1+(abs(nu-mn)>=abs(ob-mn)-1e-12).sum())/5000 if ok else np.nan,family='within MENDER niche',targets=int(imm.sum()),valid_high=nH,valid_low=nL,reason='estimable' if ok else 'sparse or degenerate',block_um=400,version=CFG['version'],human_review_status=CFG['human_review_status']))
 put(pd.DataFrame(metrics),p/'within_niche_effects.tsv');put(pd.DataFrame(agg),p/'niche_overlay.tsv')
 put(pd.DataFrame(comp,columns=['sample','method','niche','identity','units','fraction']),p/'composition.tsv')
 put(pd.DataFrame([(a,b,adjusted_rand_score(idx_niche[a],idx_niche[b])) for a,b in itertools.combinations(idx_niche,2)],columns=['a','b','ARI']),p/'method_agreement.tsv')
def match():
 comp=pd.concat([pd.read_csv(OUT/'niche'/s/'composition.tsv',sep='\t') for s in HD+XE]);tab=comp[comp.method=='MENDER'].pivot_table(index=['sample','niche'],columns='identity',values='fraction',fill_value=0);out=[]
 for key,x in tab.iterrows():
  for s in (HD if key[0] in HD else XE):
   if s==key[0]:continue
   ds=tab.loc[s].apply(lambda r:jensenshannon(r,x,base=2),axis=1);ix=ds.idxmin();out.append((*key,s,ix,ds.min(),x.idxmax(),tab.loc[(s,ix)].idxmax()))
 put(pd.DataFrame(out,columns=['sample','niche','comparison_sample','closest_niche','JS_distance','dominant','comparison_dominant']),OUT/'reports/niche_similarity.tsv')
if __name__=='__main__':
 for s in sys.argv[1:] or HD+XE:collect(s)
 if not sys.argv[1:]:match()
