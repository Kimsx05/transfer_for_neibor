from common import *
from PIL import Image,ImageDraw
from geometry import pixel
Image.MAX_IMAGE_PIXELS=None
def run():
 notes=pd.read_csv(OUT/'review/visual_notes.tsv',sep='\t');rows=[]
 for s in HD+XE:
  v,d=folders(s);p=paths(s);nt=notes[notes['sample']==s].set_index('id');nt['path']=['review/'+s+'/cards/'+id0+'.png' for id0 in nt.index]
  for fn in [d/'identity.parquet',p/'candidates.parquet']:
   e=pd.read_parquet(fn);assert nt.index.isin(e.index).all();e['review_window_note']=e.index.map(nt.window_observation);e.loc[nt.index,'review_mode']='table_and_image';e.loc[nt.index,'evidence_path']=nt.path;e.to_parquet(fn)
  if s not in HD:continue
  hu=pd.read_csv(p/'human_review.tsv',sep='\t');z=labels(s)
  if s!='UC_SCC':
   folder=imagedir(s)/'recovered';cfg=json.loads((folder/'recovered_image_coordinates.json').read_text());fn=Path(cfg['image_path']);xx,yy=cfg['source_fullres_crop_xyxy'][:2];sc=1/cfg['microns_per_pixel']
  else:
   fn=imagedir(s)/'tissue_hires_image.png';cfg=sf(s);xx=yy=0;sc=cfg['tissue_hires_scalef']/cfg['microns_per_pixel']
  img=Image.open(fn);dst=p/'histology';dst.mkdir(exist_ok=True)
  aff=np.load(v/'window.npz')['affine'];inv=np.linalg.inv(aff[:2]);off=np.array([[-1,-1],[-1,1],[1,1],[1,-1],[-1,-1]])@inv
  for id0 in hu.id:
   row=z.loc[id0];cx,cy=row.x_um,row.y_um;box=((cx-55)*sc-xx,(cy-55)*sc-yy,(cx+55)*sc-xx,(cy+55)*sc-yy);im=img.crop(tuple(round(j) for j in box)).convert('RGB').resize((440,440))
   # Draw grid footprint from global origin, not an invented cell boundary.
   center=(row[['array_row_grid','array_col_grid']].values.astype(float)*2+.5-aff[2])@inv;poly=(center+off-[cx-55,cy-55])*4
   draw=ImageDraw.Draw(im);draw.line([tuple(q) for q in poly],fill='#008fb5',width=2);draw.line((10,426,90,426),fill='black',width=3);im.save(dst/(id0+'.png'));im.save(dst/(id0+'.pdf'),resolution=200)
   rows.append((s,id0,str(fn),sc,xx,yy,*box))
  idx=p/'index.html';src=idx.read_text();pos=''.join('<li><a href="histology/'+q+'.png">'+q+'</a></li>' for q in hu.id)
  if 'Full-resolution crops' not in src:idx.write_text(src+'<h3>Full-resolution crops (UC_SCC uses deposited hires)</h3><ol>'+pos+'</ol>')
 put(pd.DataFrame(rows,columns=['sample','id','image','px_per_um','origin_x','origin_y','left_px','top_px','right_px','bottom_px']),OUT/'audit/registered_images.tsv')
if __name__=='__main__':run()
