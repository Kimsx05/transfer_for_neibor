from refit import *
import umap
cid=sys.argv[1];p=OUT/'review'/cid;a=np.load(p/'reference_embedding.npz');z=umap.UMAP(n_neighbors=30,min_dist=.3,random_state=123,n_jobs=1).fit_transform(a['reference_coordinates']);npz(p/'umap.npz',xy=z,reference_index=a['reference_index'],labels=a['reference_labels'])
