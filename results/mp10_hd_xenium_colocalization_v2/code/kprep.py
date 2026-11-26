from common import *
from geometry import pixel
for s in sys.argv[1:] or HD+XE:
 v,d=folders(s);wi=np.load(v/'window.npz');mask=wi['mask'];mask.astype(np.uint8).ravel(order='F').tofile(d/'k_mask.bin');pd.Series(mask.shape).to_csv(d/'k_shape.txt',index=False,header=False)
 a=pd.read_csv(d/'crawdad_input.tsv.gz',sep='\t');rc=pixel(a[['x_um','y_um']].values,wi['affine']);a['x']=rc[:,1]*8;a['y']=rc[:,0]*8
 put(a[['id','x','y','Treg']+[k for k in ['MP10']+STATES if k in a]],d/'k_points.tsv.gz')
