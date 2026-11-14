import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import pathlib,pandas as pd,numpy as np,fitz
from PIL import Image,ImageDraw
B=pathlib.Path('results/mp10_hypothesis_followup_v1');fm=pd.read_csv(B/'tables/figure_index.tsv',sep='\t');checks=[];thumbs=[]
for row in fm.itertuples():
 for filename in row.files.split(';'):
  p=B/'figures'/filename;d=dict(figure=row.figure,file=str(p),direct_data_exists=(B/row.data).exists(),passed=p.exists())
  if p.suffix=='.pdf':
   doc=fitz.open(p);page=doc[0];spans=[s for bl in page.get_text('dict')['blocks'] if 'lines' in bl for ln in bl['lines'] for s in ln['spans']];outside=[s['text'] for s in spans if s['bbox'][0]<-1 or s['bbox'][1]<-1 or s['bbox'][2]>page.rect.width+1 or s['bbox'][3]>page.rect.height+1];d.update(pages=len(doc),raster_images=len(page.get_images()),outside_text=' | '.join(outside),fonts=';'.join(sorted({s['font'] for s in spans})));d['passed'] &= not outside and not d['raster_images'] and bool(spans)
   pix=page.get_pixmap(matrix=fitz.Matrix(.65,.65),alpha=False)
   if not any(n in p.stem for n in ['panel_only','legend_only']):
    im=Image.frombytes('RGB',(pix.width,pix.height),pix.samples);im.thumbnail((510,370));can=Image.new('RGB',(540,405),'white');can.paste(im,((540-im.width)//2,28));ImageDraw.Draw(can).text((8,8),row.figure,fill='black');thumbs.append(can)
   doc.close()
  else:
   with Image.open(p) as im:
    dpi=im.info.get('dpi',(0,0));d.update(width=im.width,height=im.height,dpi=dpi[0]);d['passed'] &= 299<=dpi[0]<=301 and im.width>=int(row.width_inches*299)
  d['passed'] &= d['direct_data_exists'];checks.append(d)
for i in range(0,len(thumbs),6):
 sh=Image.new('RGB',(1620,810),'#eeeeee')
 for j,im in enumerate(thumbs[i:i+6]):sh.paste(im,((j%3)*540,(j//3)*405))
 sh.save(B/'qc'/f'contact_sheet_{i//6+1:02d}.png')
x=pd.DataFrame(checks);x.to_csv(B/'qc/render_audit.tsv',sep='\t',index=False);print('Render',x.passed.sum(),'/',len(x), 'figures',len(fm));print(x[~x.passed].to_string(index=False));assert x.passed.all()
