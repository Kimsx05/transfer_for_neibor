from common5 import *
import importlib.util
from scipy.stats import spearmanr
from scipy.signal import fftconvolve
import warnings
spec=importlib.util.spec_from_file_location('historical',V/'code/common.py');hist=importlib.util.module_from_spec(spec);spec.loader.exec_module(hist)
def score(x,sy):
 assert x.shape[1]==18085 and all(np.sum(sy==g)==1 for g in hist.TOP)and list(hist.TOP)==top
 selected=np.array([np.flatnonzero(sy==g)[0]for g in top]);ranks=hist.selected_ranks(x.indptr,x.indices,x.data,x.shape[1],selected)
 return hist.mp10(pd.DataFrame(ranks,columns=top))
def check(s):
 # Rebuild eight old Epi64 count rows directly from frozen native8 members.
 old=pd.read_parquet(V/s/'Epi64_regions.parquet',columns=['region_id','MP10','measurement_eligible']);old=old[old.measurement_eligible].iloc[::max(1,old.measurement_eligible.sum()//8)].head(8).reset_index(drop=True)
 m=pd.read_parquet(V/s/'Epi64_members.parquet');m=m[m.region_id.isin(old.region_id)]
 with h5py.File(js(V/s/'prepared.json')['counts_path'])as f:
  g=f['matrix'];sy=g['features/name'][:].astype(str);ix=pd.Index(g['barcodes'][:].astype(str)).get_indexer(m.bin_id);assert ix.min()>=0;ptr=g['indptr'][:];data=np.zeros((len(old),len(sy)),np.int64);target=pd.Index(old.region_id).get_indexer(m.region_id)
  for i,l in zip(ix,target):
   a,b=ptr[i:i+2];np.add.at(data[l],g['indices'][a:b],g['data'][a:b])
 sc=score(sparse.csr_matrix(data),sy);old['recomputed']=sc;assert np.allclose(sc,old.MP10,atol=1e-7);return old
def run(s):
 out=R/s;d=pd.read_parquet(out/'cells_pre_MP10.parquet');X=sparse.load_npz(out/'counts.npz');sy=rd(out/'genes.tsv').symbol.to_numpy();old=check(s);put(old,out/'historical_score_check.tsv')
 excluded=(d.roi.isin([2,3,4]) if s=='Cervilla'else pd.Series(False,index=d.index));d['image_review_exclusion']=excluded
 d['epi_eligible']=d.epithelial&~excluded&d.in_D
 d['score_ok']=d.epi_eligible&d.nFeature.ge(50);ok=d.score_ok.to_numpy();v=score(X[ok],sy);d['MP10_cell_score']=np.nan;d.loc[ok,'MP10_cell_score']=v
 highrows=[]
 for col,q in [('H05',.95),('H10',.9)]:
  cut=np.quantile(v,q,method='linear');d[col]=pd.Series(pd.NA,index=d.index,dtype='boolean');d.loc[ok,col]=v>cut;highrows.append(dict(section=s,definition=col,threshold=cut,high=int((v>cut).sum()),eligible=len(v),fraction=np.mean(v>cut),ties_at_cut=(v==cut).sum(),zeros=(v==0).sum(),score_depth_spearman=spearmanr(v,d.nCount[ok]).statistic,score_feature_spearman=spearmanr(v,d.nFeature[ok]).statistic))
 d['group']=np.select([d.status.eq('supported_Treg'),d.status.eq('supported_other_T'),d.FOXP3.gt(0)&d.core_width.ge(2)&d.T_width.ge(2)&~d.mixed&~d.geometry_review&~d.image_confirmed_invalid&d.nuclear_bins.gt(0)],['Treg','other_T','RNA_candidate'],default='other')
 d.loc[excluded|~d.in_D,'group']='other';d.to_parquet(out/'cells_scored.parquet',index=False);put(pd.DataFrame(highrows),out/'thresholds.tsv')
 cent=d[d.group.ne('other')].copy();ep=d[d.epi_eligible].copy();b=pd.read_parquet(V/s/'bin_evidence.parquet',columns=['x_um','y_um','component','eligible_pool']);rows=[]
 for c,cs in cent.groupby('component'):
  e=ep[ep.component==c].reset_index(drop=True);b0=b[(b.component==c)&b.eligible_pool]
  te=cKDTree(e[['x','y']])if len(e)else None;tc=cKDTree(b0[['x_um','y_um']])if len(b0)else None
  for r in [50,200]:
   xy=cs[['x','y']].to_numpy();nb=te.query_ball_point(xy,r)if te else [[] for _ in xy];coverage=tc.query_ball_point(xy,r,return_length=True)*64/(np.pi*r*r)if tc else np.zeros(len(xy))
   for z,ix,frac in zip(cs.itertuples(),nb,coverage):
    ev=e.iloc[ix];sc=ev[ev.score_ok];rows.append(dict(section=s,roi=z.roi,id=z.id,x=z.x,y=z.y,group=z.group,radius=r,nCount=z.nCount,component=c,coverage=min(1,frac),epi=len(ev),scored=len(sc),scoreable_fraction=len(sc)/len(ev)if len(ev)else np.nan,H05_n=sc.H05.sum(),H10_n=sc.H10.sum(),H05_fraction=sc.H05.mean()if len(sc)else np.nan,H10_fraction=sc.H10.mean()if len(sc)else np.nan,mean_MP10=sc.MP10_cell_score.mean(),eligible=(c>=0 and len(sc)>0 and frac>=.5),epi_ids='|'.join(map(str,sc.id))))
 nb=pd.DataFrame(rows);nb.to_parquet(out/'neighborhoods.parquet',index=False);effects=[]
 for r in [50,200]:
  for group in ['Treg','RNA_candidate']:
   a=nb[nb.radius.eq(r)&nb.eligible].copy();a['block']=(a.x//400).astype(int).astype(str)+'_'+(a.y//400).astype(int).astype(str)
   a['opp']=np.digitize(a.scored,[5,20]);a['Tdepth']=(a.nCount>=100).astype(int);keys=['component','block','opp','Tdepth']
   for metric in ['H05_fraction','H10_fraction','mean_MP10']:
    cases=a[a.group==group];control=a[a.group=='other_T'];bars=control.groupby(keys)[metric].agg(['mean','size']).reset_index();p=cases.merge(bars,on=keys)
    delta=(p[metric]-p['mean']).astype(float);blockdelta=delta.groupby(p.block).mean();lo=hi=np.nan
    if len(blockdelta)>=5:
     rng=np.random.default_rng(123);sims=rng.choice(blockdelta.to_numpy(),size=(2000,len(blockdelta)),replace=True).mean(axis=1);lo,hi=np.quantile(sims,[.025,.975])
    effects.append(dict(section=s,radius=r,group=group,metric=metric,all_cases=int((cent.group==group).sum()),all_controls=int((cent.group=='other_T').sum()),eligible_cases=len(cases),eligible_controls=len(control),matched_centers=len(p),matched_blocks=len(blockdelta),unadjusted=cases[metric].mean()-control[metric].mean(),matched_block_mean_difference=blockdelta.mean(),lo=lo,hi=hi,reason='no supported case'if not (cent.group==group).any()else 'no eligible neighborhood'if not len(cases)else 'no comparable controls'if not len(p)else 'fewer than five matched blocks'if len(blockdelta)<5 else 'descriptive spatial block interval'))
 put(pd.DataFrame(effects),out/'effects.tsv')
if __name__=='__main__':
 for s in sections:run(s)
