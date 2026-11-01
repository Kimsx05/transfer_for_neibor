"""Coarse author-figure reference, registered ONLY using paired H&E images.
No barcode truth, independent validation or spot-level accuracy is inferred.
"""
import os,sys,json
from pathlib import Path
R=Path(__file__).resolve().parents[3];O=R/'results/mp10_pathology_regions_v3_3';PRE=R/'results/mp10_pathology_regions_v3_2';OLD=R/'results/mp10_pathology_regions_v3_1'
sys.path.insert(0,str(R/'results/mp10_pathology_regions_v3/scripts'))
from run_v3 import np,pd,M,Image,plt,Patch
import cv2
from scipy.ndimage import distance_transform_edt,label,uniform_filter
cv2.setNumThreads(1);cv2.setRNGSeed(20261030)
LAYOUT=[['BC01','BC11','BC02'],['BC05','BC06','BC07'],['BC08','BC09','BC10'],['BC12','BC03','BC04'],['BC15','BC13','BC16'],['BC17','BC18','BC14'],['BC19','BC20','BC22'],['BC21']]
POS={s:(r,c) for r,row in enumerate(LAYOUT) for c,s in enumerate(row)}
# Pixel centers of published legend swatches, sampled from the existing S1 raster.
NAMES=['Tumor Cell','Tumor Stroma','Immune Cell','Necrosis','TLS/TIL/LA','Vessel','Urothelium','Normal Stroma','Muscle','Others']
SWATCH=[(475,1516),(475,1538),(475,1559),(625,1516),(625,1538),(625,1559),(749,1516),(749,1538),(749,1559),(865,1559)]
def palette(im):
 # Robust median of 3x3 center swatch pixels.
 return np.array([np.median(im[y-1:y+2,x-1:x+2].reshape(-1,3),axis=0) for x,y in SWATCH])
def fit(src,tgt):
 # Extract descriptors at the paper's effective resolution (no artificial extra detail).
 sc=320/max(src.shape[:2]);srcs=cv2.resize(src,None,fx=sc,fy=sc,interpolation=cv2.INTER_AREA);tt=cv2.resize(tgt,None,fx=2,fy=2,interpolation=cv2.INTER_CUBIC)
 sift=cv2.SIFT_create(nfeatures=5000,contrastThreshold=.008,edgeThreshold=20)
 k1,a=sift.detectAndCompute(cv2.cvtColor(srcs,cv2.COLOR_RGB2GRAY),None);k2,b=sift.detectAndCompute(cv2.cvtColor(tt,cv2.COLOR_RGB2GRAY),None)
 if a is None or b is None:return None,{'status':'insufficient_HE_features'}
 pairs=cv2.BFMatcher().knnMatch(a,b,k=2);good=[x for x,y in pairs if x.distance<.75*y.distance]
 if len(good)<6:return None,{'status':'insufficient_HE_matches','n_matches':len(good)}
 p=np.array([k1[x.queryIdx].pt for x in good],np.float32)/sc;q=np.array([k2[x.trainIdx].pt for x in good],np.float32)/2
 A,mask=cv2.estimateAffinePartial2D(p,q,method=cv2.RANSAC,ransacReprojThreshold=2.,maxIters=10000,confidence=.999,refineIters=50)
 if A is None:return None,{'status':'HE_registration_failed','n_matches':len(good)}
 keep=mask.ravel().astype(bool);err=np.linalg.norm(np.c_[p,np.ones(len(p))]@A.T-q,axis=1);pts=q[keep];coverage=cv2.contourArea(cv2.convexHull(pts))/np.prod(tgt.shape[:2]) if len(pts)>2 else 0
 warped=cv2.warpAffine(src,A,(tgt.shape[1],tgt.shape[0]),borderValue=(255,255,255));g1=cv2.cvtColor(warped,cv2.COLOR_RGB2GRAY);g2=cv2.cvtColor(tgt,cv2.COLOR_RGB2GRAY);maskim=(g1<245)|(g2<245);corr=np.corrcoef(g1[maskim],g2[maskim])[0,1]
 reliable=keep.sum()>=8 and keep.mean()>=.35 and np.quantile(err[keep],.95)<=2.5 and coverage>=.12 and corr>=.5
 info={'status':'reliable_coarse_HE_registration' if reliable else 'registration_unreliable_no_label_transfer','n_matches':len(good),'n_inliers':int(keep.sum()),'inlier_fraction':float(keep.mean()),'inlier_error_p95_author_px':float(np.quantile(err[keep],.95)),'inlier_spatial_coverage':float(coverage),'HE_gray_correlation':float(corr),'source_to_author_HE_affine':A.tolist(),'rotation_degrees':float(np.degrees(np.arctan2(A[1,0],A[0,0]))),'scale_author_px_per_source_px':float(np.sqrt(abs(np.linalg.det(A[:,:2])))),'reflection':'not detected by similarity model; unreliable fits not transferred','registration_basis':'H&E features only; spatial labels not used'}
 return (A,p,q,keep,warped),info

def run(sec,allim,colors):
 sample=sec.split('_')[-1];r,c=POS[sample];p=O/'registration'/sec;p.mkdir(exist_ok=True)
 # Independent exact panel extents: map and H&E have their own crop origins.
 ax=[90,438,785][c];hx=[247,597,943][c];y=57+202*r;h=155;w=155
 annot=allim[y:y+h,ax:ax+w];he=allim[y:y+h,hx:hx+w];src=np.array(Image.open(M.loc[sec,'image_path']).convert('RGB'))
 result,info=fit(src,he);info.update(section_id=sec,author_sample=sample,source='SpatialBC S1; figshare file 69527807',annotation_crop=[ax,y,w,h],HE_crop=[hx,y,w,h],reference_form='author-figure-derived COARSE regional reference, not original barcode truth')
 d=pd.read_parquet(OLD/'by_dataset'/M.loc[sec,'dataset_id']/sec/'spot_regions_v3_1.parquet');ref=np.full(len(d),'unavailable',object);classlabel=np.full(len(d),'unavailable',object);margin=np.full(len(d),np.nan);mapped=np.full((len(d),2),np.nan)
 fig,axs=plt.subplots(1,4,figsize=(16,5));axs[0].imshow(annot);axs[0].set_title('Published S1 compartments');axs[1].imshow(he);axs[1].set_title('Published paired H&E')
 if result is not None:
  A,p1,p2,keep,warped=result;axs[2].imshow(he);axs[2].imshow(warped,alpha=.5);axs[2].scatter(p2[keep,0],p2[keep,1],s=8,facecolors='none',edgecolors='lime');axs[2].set_title('H&E registration / inliers')
  # Matching map and H&E printed extents may differ by <1px; normalize crop grids.
  mapped=np.c_[d[['x','y']].to_numpy(),np.ones(len(d))]@A.T
  # Identify only saturated flat author palette colors; original background H&E / gray Others excluded.
  dist=np.linalg.norm(annot[:,:,None,:].astype(float)-colors[None,None,:,:],axis=3);cl=dist.argmin(2);exact=dist.min(2)<28
  # Macro categories combine Tumor Cell+Tumor Stroma. Reference Others/necrosis/immune/muscle are not a normal/tumour binary.
  macro=np.zeros((h,w),np.uint8);macro[exact&np.isin(cl,[0,1])]=1;macro[exact&np.isin(cl,[6,7])]=2
  # Use a coarse neighbourhood instead of eroding every antialiasing/white pinhole.
  # The spot centre must itself have an explicit author tumour/normal colour.
  # Others never inherit a label; opposite macro colours in the uncertainty radius veto transfer.
  masks=[];buffer=max(4.,float(info['inlier_error_p95_author_px'])*2+2.);radius=int(np.ceil(buffer));win=2*radius+1
  fractions=[uniform_filter((macro==code).astype(float),size=win,mode='constant') for code in [1,2]]
  for code in [1,2]:
   mask=(macro==code);lab,n=label(mask);sz=np.bincount(lab.ravel());large=sz[lab]>=60
   consensus=(fractions[code-1]>=.70)&(fractions[2-code]<=.02)
   # The entire registration-uncertainty radius must stay within the coarse consensus region.
   safe=(distance_transform_edt(consensus)>=buffer)&mask&large
   masks.append(safe)
  if info['status']=='reliable_coarse_HE_registration':
   for i,(x,y1) in enumerate(mapped):
    xx,yy=int(round(x)),int(round(y1))
    if 0<=xx<w and 0<=yy<h:
     classlabel[i]=NAMES[cl[yy,xx]] if exact[yy,xx] else 'boundary_mixed_or_background'
     if classlabel[i]=='Others':ref[i]='Others_not_assigned'
     else:ref[i]='reference_edge_or_unassigned'
     for code,safe in enumerate(masks,1):
      if safe[yy,xx]:ref[i]='author_coarse_inside' if code==1 else 'author_coarse_outside';margin[i]=fractions[code-1][yy,xx]
  info['exclusion_buffer_author_px']=buffer;info['minimum_same_class_patch_author_px']=60
  axs[3].imshow(annot,alpha=.5);col={'author_coarse_inside':'#b61535','author_coarse_outside':'#0066c0'}
  for k,v in col.items():sel=ref==k;axs[3].scatter(mapped[sel,0],mapped[sel,1],s=5,c=v,label=k)
  axs[3].legend(fontsize=6,loc='lower center',bbox_to_anchor=(.5,-.12));axs[3].set_title('Buffered coarse reference only')
 else:
  axs[2].text(.05,.5,info['status'],transform=axs[2].transAxes);axs[3].imshow(annot)
 for a in axs:a.axis('off')
 fig.suptitle(sec+' | '+info['status']+' | no barcode-level accuracy claim');fig.tight_layout();fig.savefig(p/'registration_check.png',dpi=300);fig.savefig(p/'registration_check.pdf');plt.close(fig)
 pd.DataFrame({'section_id':sec,'barcode':d.barcode,'author_coarse_reference':ref,'author_panel_class_at_center':classlabel,'author_reference_neighborhood_fraction':margin,'author_x':mapped[:,0],'author_y':mapped[:,1]}).to_csv(p/'coarse_reference.tsv.gz',sep='\t',index=False)
 info['n_coarse_inside']=int((ref=='author_coarse_inside').sum());info['n_coarse_outside']=int((ref=='author_coarse_outside').sum());(p/'registration.json').write_text(json.dumps(info,indent=2));print(sec,info,flush=True);return info
if __name__=='__main__':
 im=np.array(Image.open(PRE/'sources/image1.png').convert('RGB'));colors=palette(im);print('PALETTE',dict(zip(NAMES,colors.tolist())),flush=True)
 rows=[run(sec,im,colors) for sec in M.index if sec.startswith('GSE319536_')];pd.DataFrame(rows).drop(columns=['source_to_author_HE_affine'],errors='ignore').to_csv(O/'tables/author_registration_summary.tsv',sep='\t',index=False)
