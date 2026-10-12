from base import *
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.spatial import ConvexHull
from PIL import Image
od=OUT/'05_user_review';ev=od/'evidence';roi=pd.read_csv(OLD/'01_regions/ROI_priority_review.tsv',sep='\t');roi=roi.sort_values('n_common_anchors',ascending=False).reset_index(drop=True)
master=pd.read_parquet(OLD/'03_features/anchor_master.parquet',columns=['spot_uid','section_id','barcode','array_row','array_col','epithelial_fraction','copykat_prediction_original','cnv_score_original','roi_id','common_anchor'])
images=[];evidence=[]
for s,rr in roi.groupby('section_id',sort=False):
 d=ev/s;ip=d/'HE_embedded.png';g=master[master.section_id==s].copy();has=ip.exists();reg='MISSING';scale=np.nan
 if has:
  im=np.asarray(Image.open(ip));co=pd.read_csv(d/'image_coordinates.csv');info=json.loads((d/'registration.json').read_text());sf=info['scale_factors'];g=g.merge(co[['barcode','imagerow','imagecol']],on='barcode',validate='one_to_one');
  # Seurat Visium image slots contain low resolution image. Verify coordinates against both known scales, never infer rotation.
  candidates=[(k,float(v)) for k,v in sf.items() if k in ['lowres','hires'] and isinstance(v,(int,float))]
  checks=[]
  for key,v in candidates:
   xx=g.imagecol*v;yy=g.imagerow*v;inside=((xx>=0)&(xx<im.shape[1])&(yy>=0)&(yy<im.shape[0])).mean();extent=max(xx.max()/im.shape[1],yy.max()/im.shape[0]);checks.append((key,v,inside,extent))
  if checks:
   checks.sort(key=lambda a:(a[2]<.99,abs(a[3]-.9)));key,scale,inside,extent=checks[0];reg='IN_BOUNDS_REQUIRES_VISUAL_CONFIRMATION' if inside>=.99 else 'COORDINATE_BOUNDS_WARNING';g['px']=g.imagecol*scale;g['py']=g.imagerow*scale
  else:has=False;reg='MISSING_SCALE'
 images.append(dict(section_id=s,HE=str(ip) if has else '',availability='embedded image' if has else 'missing usable HE',registration_status=reg,scale_key=key if has else '',scale_used=scale,coordinate_file=str(d/'image_coordinates.csv') if has else '',required_missing='original full-resolution H&E and pathology annotation' if has else 'H&E image + scalefactors + tissue_positions linked by barcode',source='ST/all_ST_combined_list_SeuratObject_updated.qs image slots',pathology_status='no finalized pathology supplied'))
 for r in rr.itertuples():
  mask=g.roi_id==r.roi_id;members=g[mask];dest=ev/(r.roi_id+'.png');q='请结合 H&E 判断：肿瘤上皮、非肿瘤上皮、其他组织，还是仍不确定？如 ROI 混合，请拆分并返回 spot 成员或 polygon；不要以 diploid 判正常。'
  if has and len(members):
   fig,axs=plt.subplots(2,3,figsize=(16,10));xs=members.px.to_numpy();ys=members.py.to_numpy();pad=20;box=(max(0,xs.min()-pad),min(im.shape[1],xs.max()+pad),max(0,ys.min()-pad),min(im.shape[0],ys.max()+pad))
   for ax in axs.ravel():ax.imshow(im);ax.set_aspect('equal');ax.set_xlim(box[0],box[1]);ax.set_ylim(box[3],box[2]);ax.set_xticks([]);ax.set_yticks([])
   axs[0,0].set_title('Embedded H&E | same-coordinate ROI crop')
   axs[0,1].scatter(g.px,g.py,s=1,c='cyan',alpha=.35);axs[0,1].scatter(xs,ys,s=2,c='red',alpha=.5);axs[0,1].set_title('Registered spots: red = exact ROI members')
   # Boundary is exact member perimeter on the physical array, not a convex lesion claim.
   lookup={(int(a),int(b)) for a,b in zip(members.array_row,members.array_col)};boundary=[any((int(a)+da,int(b)+db) not in lookup for da,db in [(0,-2),(0,2),(-1,-1),(-1,1),(1,-1),(1,1)]) for a,b in zip(members.array_row,members.array_col)]
   axs[0,2].scatter(xs[boundary],ys[boundary],s=4,facecolors='none',edgecolors='yellow',linewidths=.5);axs[0,2].set_title('ROI member perimeter | not pathology border')
   sc=axs[1,0].scatter(g.px,g.py,c=g.epithelial_fraction,s=4,cmap='viridis',vmin=0,vmax=1);fig.colorbar(sc,ax=axs[1,0],shrink=.65);axs[1,0].set_title('Center c2l mean Epithelial / all30')
   calls=g.copykat_prediction_original.fillna('missing').astype(str);cat=np.where(calls.eq('aneuploid'),'aneuploid',np.where(calls.eq('diploid'),'diploid',np.where(calls.eq('missing'),'missing','unknown/low.conf')))
   for c,color in [('aneuploid','#cc3311'),('diploid','#0077bb'),('unknown/low.conf','#eeaa33'),('missing','#bbbbbb')]:
    mm=cat==c;axs[1,1].scatter(g.loc[mm,'px'],g.loc[mm,'py'],s=4,c=color,label=c)
   axs[1,1].legend(fontsize=7);axs[1,1].set_title('Existing CopyKAT calls')
   sc=axs[1,2].scatter(g.px,g.py,c=g.cnv_score_original,s=4,cmap='magma');fig.colorbar(sc,ax=axs[1,2],shrink=.65);axs[1,2].set_title('Existing MSQ score | author scale unverified')
   fig.suptitle(f'{r.roi_id}\nanchors={r.n_common_anchors}; {reg}; pathology pending',fontsize=11);fig.tight_layout();fig.savefig(dest,dpi=150,transparent=True);plt.close(fig)
  else:dest=Path(r.evidence_figure)
  evidence.append(dict(roi_id=r.roi_id,dataset_id=r.dataset_id,section_id=s,n_common_anchors=r.n_common_anchors,n_spots=r.n_spots,parent_roi_id=r.roi_id,roi_version='cnv_candidate_v1',evidence_figure=str(dest),raw_HE=str(ip) if has else '',registration_status=reg,mean_epi=r.mean_epi,aneuploid_n=r.aneuploid_n,diploid_n=r.diploid_n,question=q,impact='later region stratification only; niche completion and clustering unchanged',reviewed_label='',review_comment='',reviewer='',review_date=''))
result=pd.DataFrame(evidence).sort_values('n_common_anchors',ascending=False);tab(od/'ROI_priority_review.tsv',result);tab(od/'ROI_review_to_fill.tsv',result);tab(od/'image_availability.tsv',images)
text=['# 你现在需要做什么','本轮 niche 的计算不等待 ROI 判断。请先复核下列 5 个高优先级区域，再按需要处理其余 128 项。133 项是计算连通块，不等同于独立病灶；原 8257 个连通块无需逐一人工复核。','打开 [待填表](ROI_review_to_fill.tsv)，保留 roi_id，填写 reviewed_label、review_comment、reviewer、review_date。允许标签 tumor_epithelium、non_tumor_epithelium、other、mixed_needs_split、uncertain。不要填写不存在的病理依据。','若区域混合，不要求给整个 ROI 一个标签：在 review_comment 写明需拆分，另返回 polygon 坐标（所附图像像素坐标）或 spot_uid 成员表，保留 parent_roi_id；后续据此生成稳定 child UID 和新 ROI version。']
for r in result.head(5).itertuples():text.append(f'## {r.roi_id}\n共同 anchors：{r.n_common_anchors}；mean Epi={r.mean_epi:.3f}；aneuploid={r.aneuploid_n}，diploid={r.diploid_n}。\n\n[证据图]({r.evidence_figure}) · [原嵌入 H&E]({r.raw_HE})\n\n明确问题：{r.question}\n配准状态：{r.registration_status}；请确认 spot 与组织边缘吻合。')
text+=['## 缺少什么','[逐切片图像清单](image_availability.tsv) 列明可用的嵌入 H&E、坐标和比例因子。已有图像属于 Seurat 保存的分辨率，不能代替原始全分辨率 H&E；尚缺病理标注。若无法判断，请针对对应 section_id 提供原始 H&E、配套 scalefactors_json.json / tissue_positions 文件及病理标注来源，不需要重做 CopyKAT。','## 如何返回','将填写后的 ROI_review_to_fill.tsv 保存到此目录或上传，必要时附 ROI 拆分 polygon/spot_uid 表。空白表示尚未复核，uncertain 表示已看证据仍不能判定。表内不预填人工姓名或日期。返回后只进入后续区域/界面确认；不会自动触发 MP10、解盲或来源验证。']
(od/'WHAT_USER_NEEDS_TO_DO.md').write_text('\n\n'.join(text)+'\n');log('ROI evidence complete '+str(len(result)))
