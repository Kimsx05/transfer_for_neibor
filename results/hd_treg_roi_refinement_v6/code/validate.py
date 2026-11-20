from classifier import *
m,x,sy=reference();split=pd.read_parquet(R/'reference_split.parquet');assert np.array_equal(m.cell_id,split.cell_id)
assert split.groupby('group_id').role.nunique().max()==1
test=split.role.eq('evaluation').to_numpy();x=x[test];m=m[test].reset_index(drop=True);mod=joblib.load(R/'model.joblib');assert np.array_equal(sy,mod[0])
rng=np.random.default_rng(123);results=[];summary=[]
for s in sections:
 d=pd.read_parquet(R/s/'cells_geometry.parquet');depth=d.loc[(d.direct_width>0)&(d.nCount>0),'nCount'].to_numpy()
 z=thin(x,rng.choice(depth,size=len(m)),rng);p=predict(mod,z);old=classify(z,sy)
 def width(gs):return np.asarray((z[:,np.isin(sy,gs)]>0).sum(1)).ravel()
 t=width(T);di=width(direct);k=width(core[1:]);f=width(['FOXP3']);lib=np.asarray(z.sum(1)).ravel()
 common=(lib>=25)&(di>0)&(p.p_T>=.5);a=common&(f>0);b=common&(f==0)&(k>0)&p.p_Treg_given_T.ge(.5)
 mixed=(t>=2)&((width(epi)>=2)|(width(fib)>=2)|(width(mye)>=2))
 prev=old.pred.eq('Treg')&(t>=2)&(k+f>=1)&((old.pTreg>=.95)|((f>0)&(k+f>=2)&old.pTreg.ge(.5)))&~mixed
 for name,v in [('v5',prev),('A',a),('B',b),('main',a|b),('strict',a&(k>0)),('extended',(lib>=25)&(di>0)&(p.p_T>=.3)&((f>0)|((k>0)&(p.p_Treg_given_T>=.3))))]:
  tmp=pd.DataFrame({'observed':v,'label':m.label,'donor_or_sample':m.group_id,'UMI':lib});tmp['depth']=pd.cut(lib,[-1,24,99,299,np.inf]).astype(str)
  for keys in [['label'],['label','depth']]:
   zz=tmp.groupby(keys)['observed'].agg(n='size',positive='sum').reset_index();zz['section']=s;zz['rule']=name;results.append(zz)
  y=m.label.eq('Treg');tp=(y&v).sum();fp=(~y&v).sum();fn=(y&~v).sum();summary.append(dict(section=s,rule=name,TP=tp,FP=fp,FN=fn,negatives=(~y).sum(),precision=tp/(tp+fp) if tp+fp else np.nan,recall=tp/y.sum(),FPR=fp/(~y).sum(),groups=m.group_id.nunique(),N=len(y),T_candidate_depth_median=np.median(depth)))
tab(pd.concat(results),'validation_by_competitor.tsv');tab(pd.DataFrame(summary),'reference_validation.tsv')
