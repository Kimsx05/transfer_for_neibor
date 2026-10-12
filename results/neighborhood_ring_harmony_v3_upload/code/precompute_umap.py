from base import *
import umap,fcntl
_,_,ref=inputs();done=set()
while len(done)<11:
 paths=list((OUT/'02_R12_pilot').glob('*/embedding.npy'))+list((OUT/'01_features').glob('*/embedding.npy'))
 for p in paths:
  od=p.parent
  if od.name in done:continue
  lock=open(od/'umap.lock','w');fcntl.flock(lock,fcntl.LOCK_EX)
  try:
   if not (od/'umap_reference.npy').exists():
    log('UMAP '+od.name);z=np.load(p,mmap_mode='r');u=umap.UMAP(n_neighbors=30,min_dist=.3,random_state=123,n_jobs=1).fit_transform(z[ref]);np.save(od/'umap_reference.npy',u)
   done.add(od.name)
  finally:fcntl.flock(lock,fcntl.LOCK_UN);lock.close()
 if len(done)<11:time.sleep(5)
log('All 11 representation UMAPs cached')
