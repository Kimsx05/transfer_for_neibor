from common import *
from numba import njit
from scipy.spatial import cKDTree
@njit
def simulate(x,y,mark,buf,rad,seed):
 n=len(x);L=mark.shape[1];reps=999;result=np.zeros((reps,L));mx=int(np.ceil(np.ptp(x)/rad))+3;my=int(np.ceil(np.ptp(y)/rad))+3
 for p in range(reps):
  tile=((x-x.min()+p*rad/reps)/rad).astype(np.int64)*my+((y-y.min()+p*rad/reps)/rad).astype(np.int64)
  # O(N) stable counting sort within the offset square tiles.
  offs=np.concatenate((np.zeros(1,np.int64),np.cumsum(np.bincount(tile,minlength=mx*my))));pos=offs[:-1].copy();order=np.zeros(n,np.int64)
  for i in range(n):
   j=tile[i];order[pos[j]]=i;pos[j]+=1
  for j in range(len(pos)):
   a=offs[j];b=offs[j+1];m=b-a
   if m==0:continue
   np.random.seed(seed+p);perm=np.random.permutation(m)
   for i in range(m):
    if buf[order[a+i]]:
     src=order[a+perm[i]]
     for k in range(L):result[p,k]+=mark[src,k]
 return result
def run(s):
 v,d=folders(s);t=pd.read_csv(d/'crawdad_input.tsv.gz',sep='\t');xy=t[['x_um','y_um']].values;tree=cKDTree(xy[t.Treg]);near=tree.query(xy)[0]<=50;keys=['MP10'] if s in HD else STATES;mark=(t[keys].values==1).astype('i4');out=[]
 # Fixed component focal centers. HD can carry high and Treg component simultaneously.
 # Exclude all focal IDs consistently from the observed/expected neighborhood.
 buf=near&~t.Treg.to_numpy();obs=mark[buf].sum(0);nb=buf.sum()
 for r in [100,200,400,800]:
  null=simulate(xy[:,0],xy[:,1],mark,buf,r,1);pool=(obs+null)/(2*nb);den=np.sqrt(pool*(1-pool)*2/nb);z=np.divide((obs-null)/nb,den,out=np.full_like(null,np.nan),where=den>0)
  for k,key in enumerate(keys):
   nu=null[:,k];mn=nu.mean();sd=nu.std();valid=sd>1e-9 and t.Treg.sum()>1
   out.append((s,key,r,50,999,int(t.Treg.sum()),nb,obs[k],mn,np.nanmean(z[:,k]),(1+(abs(nu-mn)>=abs(obs[k]-mn)-1e-10).sum())/1000 if valid else np.nan))
 put(pd.DataFrame(out,columns=['sample','state','tile_um','neighDist','perms','focal','unique_buffer_nodes','observed_high','expected_high','mean_proportion_Z','P']),d/'CRAWDAD_tile_sensitivity.tsv')
if __name__=='__main__':
 for s in sys.argv[1:] or HD+XE:run(s)
