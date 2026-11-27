from base_data import *
from scipy import ndimage as ndi
from skimage.transform import resize
from stardist.models import StarDist2D
import matplotlib.pyplot as plt
import h5py
from matplotlib.collections import LineCollection
from skimage.segmentation import find_boundaries

model=StarDist2D(None,name='fluo',basedir=str(O/'shared'))
model.keras_model.predict=lambda x,**kw: [y.numpy() for y in model.keras_model(x,training=False)]
def segment(gy,primary,covered,bounds,offset):
 r,c,h,w=bounds; pad=30
 y=max(0,r-pad);x=max(0,c-pad); sub=np.s_[y:r+h+pad,x:c+w+pad]
 raw=gy[sub]; p=primary[sub]; ok=covered[sub]
 out=np.zeros((min(h,gy.shape[0]-r),min(w,gy.shape[1]-c)),np.int32)
 if not np.any(ok&(p==0)&(raw>0)):return out,offset
 # Sample fluorescence-like expression at 1 um, with a fixed section-wide intensity scale.
 im=ndi.gaussian_filter(resize(raw,(raw.shape[0]*2,raw.shape[1]*2),order=1,preserve_range=True,anti_aliasing=False).astype('f4'),1.5)
 seg,_=model.predict_instances(im,prob_thresh=.10,nms_thresh=.3,verbose=False,n_tiles=(1,1))
 lab=seg[1::2,1::2]; bad=np.unique(lab[(p>0)|~ok]);lab[np.isin(lab,bad[bad>0])]=0
 core=lab[r-y:r-y+out.shape[0],c-x:c-x+out.shape[1]];labels=np.unique(core);labels=labels[labels>0]
 if len(labels):
  # Discard nuclei crossing the core tile seam rather than allow fragmented/duplicated objects.
  total=np.bincount(lab.ravel()); cnt=np.bincount(core.ravel(),minlength=len(total)); take=labels[(cnt[labels]==total[labels])&(cnt[labels]>=3)&(cnt[labels]<=70)]
  mapping=np.zeros(len(total),np.int32);mapping[take]=np.arange(len(take))+offset+1;out=mapping[core];offset+=len(take)
 return out,offset

def run(s,mode):
 d=paths(s);d.mkdir(exist_ok=True);old=H/s
 b=pd.read_parquet(old/'inputs/bin_to_cell.parquet');r=b.array_row.values;c=b.array_col.values
 shape=(int(r.max())+1,int(c.max())+1);cov=np.zeros(shape,bool);cov[r,c]=b.image_covered
 prim=np.zeros(shape,'i4');prim[r,c]=b.cell_label;g=np.zeros(shape,'f4');g[r,c]=b.raw_counts
 cap=np.quantile(b.raw_counts[b.raw_counts>0],.995);g=np.minimum(g/cap,1)
 pos=pd.read_parquet(src(s).parent/'spatial/tissue_positions.parquet',filters=[('in_tissue','==',1)],columns=['array_row','array_col','pxl_col_in_fullres','pxl_row_in_fullres']).iloc[::157]
 tx=np.linalg.lstsq(np.c_[pos.array_row,pos.array_col,np.ones(len(pos))],pos[['pxl_col_in_fullres','pxl_row_in_fullres']],rcond=None)[0]
 t=json.loads((old/'metadata/transform.json').read_text());mpp=t['native_mpp']; trans=np.array(t['local_to_model_scale_xy']);origin=t['fullres_crop_xyxy'][:2]
 def grid(xy):return (xy-tx[2])@np.linalg.inv(tx[:2])
 if mode=='pilot':
  boxes=[]
  for name,box in t['representative_model_bboxes'].items():
   xy=(np.mean(np.reshape(box,(2,2)),axis=0)+.5)/trans-.5+origin
   boxes.append((name,grid(xy)))
  oldlabels=pd.read_parquet(old/'cell_labels.parquet')
  for name,sel in [('immune',oldlabels.predicted_major.eq('T/NK')),('ambiguous',oldlabels.low_information)]:
   q=oldlabels.loc[sel,['x_um','y_um']];tiles=np.floor(q/128).astype(int);k=tiles.value_counts().index[0];boxes.append((name,grid((np.array(k)+.5)*128/mpp)))
  rows=[];im=np.load(old/'inputs/model_image.npy',mmap_mode='r')
  for name,rc in boxes:
   y,x=np.maximum(0,np.rint(rc).astype(int)-42);sl=np.s_[y:y+84,x:x+84];p=prim[sl];sec,end=segment(g,prim,cov,(y,x,84,84),int(prim.max()))
   R,C=np.mgrid[y:y+84,x:x+84];pts=np.c_[R.ravel(),C.ravel(),np.ones(R.size)]@tx;mp=((pts-origin+.5)*trans-.5).reshape(84,84,2)
   # View source image on its native model-image frame; overlay bin assignment locations in that frame.
   x0,y0=np.floor(mp.min((0,1))).astype(int);x1,y1=np.ceil(mp.max((0,1))).astype(int)+1
   for route,lab in [('old',p),('HE',p),('joint',np.where(p,p,sec))]:
    fig,ax=plt.subplots(figsize=(4,4));ax.imshow(im[max(0,y0):y1,max(0,x0):x1],extent=(max(0,x0),x1,y1,max(0,y0)))
    for arr,col in [(p,'#24c6c5'),(sec if route=='joint' else np.zeros_like(p),'#ff8c00')]:
     v=find_boundaries(arr,mode='inner')&(arr>0);ax.scatter(mp[:,:,0][v],mp[:,:,1][v],s=2,color=col,alpha=.8)
    ax.set_title(s+' '+name+' '+route);ax.axis('off');fig.tight_layout()
    for ext in ['png','pdf']:fig.savefig(O/'figures'/f'{s}_{name}_{route}.{ext}',dpi=130)
    plt.close(fig)
    rows.append(dict(sample=s,roi=name,row=y,col=x,core_bins=84,buffer_bins=30,route=route,objects=len(np.unique(lab[lab>0])),secondary=len(np.unique(sec[sec>0])) if route=='joint' else 0,unassigned_umi=int(b.raw_counts[(r>=y)&(r<y+84)&(c>=x)&(c<x+84)].sum()-np.rint(0)),assigned_bins=int((lab>0).sum())))
    rows[-1]['unassigned_umi']=int(b.loc[(r>=y)&(r<y+84)&(c>=x)&(c<x+84)&(lab[np.clip(r-y,0,83),np.clip(c-x,0,83)]==0),'raw_counts'].sum())
  write(pd.DataFrame(rows),s+'_ROI.tsv');return
 assert (O/'shared/reconstruction_freeze.md').exists()
 sec=np.zeros(shape,'i4');offset=int(prim.max());pmax=offset
 for y in range(0,shape[0],224):
  for x in range(0,shape[1],224):
   block,offset=segment(g,prim,cov,(y,x,224,224),offset);sec[y:y+block.shape[0],x:x+block.shape[1]]=block
  print(s,y,offset-pmax,flush=True)
 final=np.maximum(prim,sec)[r,c];assert np.array_equal(final[b.cell_label>0],b.cell_label[b.cell_label>0])
 xx,ids,genes=tenx(src(s));assert np.array_equal(b.barcode,ids) and np.array_equal(xx.sum(1).A1,b.raw_counts)
 keep=final>0;u,inv=np.unique(final[keep],return_inverse=True);agg=sparse.csr_matrix((np.ones(keep.sum(),'i4'),(inv,np.flatnonzero(keep))),shape=(len(u),len(b)));raw=(agg@xx).tocsr();raw.eliminate_zeros()
 assert np.array_equal(raw.sum(0),xx[keep].sum(0))
 n=agg.sum(1).A1;xy=np.c_[(agg@r.astype(float))/n,(agg@c.astype(float))/n,np.ones(len(n))]@tx*mpp
 obs=pd.DataFrame(dict(source=np.where(u<=pmax,'HE_primary','GEX_secondary'),label=u,n_bins=n,area_um2=n*4,n_counts=raw.sum(1).A1,nFeature=raw.getnnz(1),x_um=xy[:,0],y_um=xy[:,1]),index=[f'{s}__retry__{i}' for i in u])
 oldobs=pd.read_parquet(old/'cell_labels.parquet');obs['old_id']=[f'{s}__v2__cell_{i:07d}' if i<=pmax else '' for i in u]
 obs['geometry_flag']=(n>110)|((u<=pmax)&~obs.old_id.isin(oldobs.index)&(obs.nFeature>=10))
 b['final_label']=final;b[['barcode','source_bin_index','final_label','cell_label','array_row','array_col','image_covered']].to_parquet(d/'members.parquet',index=False)
 a=ad.AnnData(raw,obs=obs,var=pd.DataFrame({'feature_id':genes.id.values,'symbol':genes.name.values},index=genes.name.values));a.var_names_make_unique();a.write_h5ad(d/'cells.h5ad',compression='lzf')
 write(obs.groupby('source').agg(objects=('label','size'),umi=('n_counts','sum'),median_UMI=('n_counts','median'),median_genes=('nFeature','median'),median_bins=('n_bins','median')).reset_index(),s+'_reconstruction.tsv')
 (d/'conservation.json').write_text(json.dumps(dict(source_umi=int(xx.sum()),assigned_umi=int(raw.sum()),unassigned_umi=int(xx[~keep].sum()),source_bins=len(ids),assigned_bins=int(keep.sum()),per_gene_check=True,cap=cap,grid_to_fullres=tx.tolist())))
for s in HD:run(s,sys.argv[1])

