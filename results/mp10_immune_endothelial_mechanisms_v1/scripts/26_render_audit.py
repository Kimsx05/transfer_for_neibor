"""Audit all delivered figure parts and render contact sheets. No biological data are read."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import pathlib,sys,json
import pandas as pd
from PIL import Image,ImageOps,ImageDraw
import fitz
B=pathlib.Path('results/mp10_immune_endothelial_mechanisms_v1');Q=B/'qc';Q.mkdir(exist_ok=True)
files=[B/'config/figure_manifest.tsv',B/'config/figure_manifest_extended.tsv',B/'config/figure_manifest_integration.tsv'];fm=pd.concat([pd.read_csv(x,sep='\t') for x in files if x.exists()],ignore_index=True).drop_duplicates(['figure_id','part'],keep='last');fm.to_csv(B/'config/figure_manifest_all.tsv',sep='\t',index=False)
checks=[];thumbs=[]
for row in fm.itertuples():
 p=pathlib.Path(row.file);entry=dict(figure_id=row.figure_id,part=row.part,file=str(p),exists=p.exists(),direct_data_exists=pathlib.Path(row.data_file).exists())
 if p.suffix=='.pdf':
  doc=fitz.open(p);page=doc[0];rect=page.rect;spans=[s for b in page.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans']];outside=[s['text'] for s in spans if s['bbox'][0]<-1 or s['bbox'][1]<-1 or s['bbox'][2]>rect.width+1 or s['bbox'][3]>rect.height+1];entry.update(pages=len(doc),raster_images=len(page.get_images()),outside_text=' | '.join(outside),fonts=';'.join(sorted({s['font'] for s in spans})),render_ok=True)
  pix=page.get_pixmap(matrix=fitz.Matrix(.5,.5),alpha=False)
  if not any(x in p.stem for x in ['legend_only','panel_only']):
   im=Image.frombytes('RGB',[pix.width,pix.height],pix.samples);im.thumbnail((400,290));canvas=Image.new('RGB',(430,330),'white');canvas.paste(im,((430-im.width)//2,25));ImageDraw.Draw(canvas).text((8,5),row.figure_id,fill='black');thumbs.append(canvas)
  doc.close()
 elif p.suffix=='.png':
  with Image.open(p) as im:entry.update(width_px=im.width,height_px=im.height,dpi_x=im.info.get('dpi',(0,0))[0],dpi_y=im.info.get('dpi',(0,0))[1],render_ok=True)
 checks.append(entry)
qa=pd.DataFrame(checks);qa.to_csv(Q/'render_audit.tsv',sep='\t',index=False)
for i in range(0,len(thumbs),9):
 sheet=Image.new('RGB',(1290,990),'#eeeeee')
 for j,t in enumerate(thumbs[i:i+9]):sheet.paste(t,((j%3)*430,(j//3)*330))
 sheet.save(Q/f'contact_sheet_{i//9+1:02d}.png')
lines=['# Figure index','', 'All figures have direct table data. PDF is vector; PNG is 300 dpi. Where a legend exists, full/panel/legend PDFs are retained. See `qc/render_audit.tsv` for automated rendering checks.','']
for fig,z in fm.groupby('figure_id',sort=False):
 lines+=['## '+fig,'',z.note.iloc[0] if pd.notna(z.note.iloc[0]) else '', '',f'Direct data: `{pathlib.Path(z.data_file.iloc[0]).relative_to(B.resolve()) if pathlib.Path(z.data_file.iloc[0]).is_absolute() else z.data_file.iloc[0]}`','']
 for r in z.itertuples():lines.append(f'- {r.part}: `figures/{pathlib.Path(r.file).name}`')
 lines.append('')
(B/'FIGURE_INDEX.md').write_text('\n'.join(lines))
assert qa.exists.all() and qa.direct_data_exists.all()
assert not qa.get('outside_text',pd.Series(dtype=str)).fillna('').str.len().gt(0).any(),'Text outside PDF page; inspect render audit'
assert qa.loc[qa.file.str.endswith('.png'),'dpi_x'].between(299,301).all()
assert qa.loc[qa.file.str.endswith('.pdf'),'raster_images'].eq(0).all(),'PDF contains raster plot elements'
print(f'Rendered and audited {fm.figure_id.nunique()} figures / {len(fm)} files; {len(thumbs)} full-figure thumbnails')
