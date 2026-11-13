from common import *
import fit_composition as fc
from sklearn.decomposition import PCA

def main():
 gate();ref='blca';r=50;a=anchors();T=train_index(a);d,P,ok=composition(ref,r);c=[]
 for s in SECTIONS:
  z=np.load(R/'models'/f'{s}_rings.npz');ix=[list(z['columns']).index(ref+'::'+t) for t in TYPES[ref]];c.append(z['counts'][:,:2,ix])
 c=np.concatenate(c).astype(float);known=c.sum(2,keepdims=True);v=np.divide(c,known,out=np.zeros_like(c),where=known>0).reshape(len(a),-1);target=float(np.median(v[T].sum(1)));v=np.log1p(np.divide(v*target,v.sum(1,keepdims=True),out=np.zeros_like(v),where=v.sum(1,keepdims=True)>0));pca=PCA(n_components=v.shape[1]-1,svd_solver='full').fit(v[T]);x=pca.transform(v)
 original=json.loads((R/'models/MENDER_adapted_blca_r50_ring_counts/complete.json').read_text());fc.CFG['resolutions']=[original['selected_parameter']]
 fc.run(x,ok,'MENDER_adapted',ref,r,'per_ring_normalized',{'PCA':pca,'target_total':target,'normalization':'known proportions per ring; empty ring zero with missing support retained; normalize_total then log1p','feature_order':'ring-major known types','include_self':False,'fixed_resolution':original['selected_parameter']})
 p=R/'models/MENDER_adapted_blca_r50_per_ring_normalized/complete.json';j=json.loads(p.read_text());j['branch_role']='conditional_quality_sensitivity';j['trigger']='Main pure-Epi clusters14/5/3/22 had average input coverage1.00/.77/.74/.67; investigate ring support rather than invent biology';p.write_text(json.dumps(j,indent=2))
if __name__=='__main__':
 import fcntl
 with (R/'logs/mender_ring_sensitivity.lock').open('w') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX)
  if not (R/'models/MENDER_adapted_blca_r50_per_ring_normalized/complete.json').exists():main()
