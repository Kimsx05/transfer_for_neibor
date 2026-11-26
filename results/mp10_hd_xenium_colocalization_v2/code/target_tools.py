from common import *
from niche import data
import squidpy as sq
def run(s):
 hd=s in HD;v,d,z,i,ct,xy,a,p=data(s);res=16 if hd else 'cell';e=pd.read_parquet(d/f'centers_{res}.parquet');ee=e.index;keys=['MP10'] if hd else STATES
 g=sparse.load_npz(v/'graph_50.npz');a.obsp['spatial_connectivities']=g;mid=xy[np.argmin(((xy-np.median(xy,axis=0))**2).sum(1))];sub=(abs(xy-mid)<=320).all(1);out=[]
 gr=sparse.load_npz(d/f'epi_{res}_200.npz');ii,jj=gr.nonzero();dd=np.linalg.norm(e[['x_um','y_um']].values[ii]-xy[jj],axis=1);tf=i.grade.eq('supported').values;allowed=ee.isin(z.index);curves=[]
 for key in keys:
  h=e[key+'_high'].to_numpy(bool)
  for r in [25,50,100,150,200]:
   for b in [True,False]:
    take=(dd<=r)&(h[ii]==b);den=take.sum();n=(take&tf[jj]).sum();curves.append((s,key,r,'high' if b else 'low',den,n,(n/den)/tf.mean() if den and tf.any() else np.nan))
  if hd:continue
  lab=z.major.where(z.major.isin(['Unknown','Mixed']),'other known').copy();lab[i.grade.eq('supported')]='Treg';lab.loc[ee]=np.where(h,'high','low');a.obs['labels']=pd.Categorical(lab)
  rr=sq.gr.nhood_enrichment(a,cluster_key='labels',n_perms=1000,n_jobs=1,rng=CFG['seed'],copy=True,handle_nan='keep',show_progress_bar=False);cats=a.obs.labels.cat.categories
  for u in ['high','low']:
   a1=cats.get_loc(u);b1=cats.get_loc('Treg');n0=lab.eq(u).sum();n1=lab.eq('Treg').sum();exp=g.nnz*(n0/len(z))*(n1/(len(z)-1));out.append((s,key,u,rr.counts[a1,b1],exp,rr.zscore[a1,b1],rr.counts[a1,b1]/exp))
  ss=a[sub].copy();out0=sq.gr.co_occurrence(ss,cluster_key='labels',interval=np.array([0,25,50,100,150,200.]),copy=True);np.savez_compressed(d/f'squidpy_target_native_{key}.npz',types=ss.obs.labels.cat.categories.to_numpy(dtype=str),values=out0[0],radius=out0[1])
 put(pd.DataFrame(curves,columns=['sample','state','r','group','edges','Treg_edges','ratio_to_unit_prevalence']),d/'fixed_target_curve.tsv')
 if out:put(pd.DataFrame(out,columns=['sample','state','group','observed','expected_global','Z','OE']),d/'squidpy_state.tsv')
if __name__=='__main__':
 for s in sys.argv[1:] or HD+XE:run(s)
