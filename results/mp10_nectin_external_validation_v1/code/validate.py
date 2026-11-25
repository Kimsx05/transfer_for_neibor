import csv,gzip,math,json,statistics,os,sys
from pathlib import Path
import fitz
from PIL import Image,ImageDraw
D=Path(__file__).resolve().parents[1];T=D/'tables';checks=[]
def read(p):
 with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return list(csv.DictReader(f,delimiter='\t'))
def ck(label,b):
 checks.append((label,bool(b)))
 if not b:raise ValueError(label)
for a,n in [('TCGA',407),('GSE32894',308),('GSE31684',93),('GSE13507',165)]:
 d=read(T/'bulk'/a/'scores.tsv');z=read(T/'bulk'/a/'core.tsv');ck(a+' unique specimen and patient',len(d)==n==len(set(r['patient_sample_id'] for r in d))==len(set(r['GSM'] for r in d)))
 ck(a+' 9 core 20 family 3 auxiliary',len(z)==9 and len(read(T/'bulk'/a/'family.tsv'))==20 and len(read(T/'bulk'/a/'Treg.tsv'))==3)
 for s in z:
  x=[float(v[s['x']]) for v in d];y=[float(v[s['y']]) for v in d];rr=statistics.correlation(x,y);ck(s['analysis_id'],abs(rr-float(s['r']))<1e-9 and int(s['n'])==n and float(s['lo'])<=rr<=float(s['hi']))
for a,line,n in [('GSE236742','HCV29',9),('GSE143803','KYSE150',4),('GSE143803','TE5',4),('GSE267018','HeLa',18)]:
 p=T/'perturbation'/a/line;m=read(p/'sample_mapping.tsv');ck(line+' GSM and groups',n==len(m)==len(set(v['GSM'] for v in m)) and len(set(v['matrix_column'] for v in m))==n)
 f=read(p/'family.tsv');ck(line+' family retained',all(len([w for w in f if w['contrast']==c])==5 for c in set(w['contrast'] for w in f)))
mods=read(T/'competition/models.tsv');ck('18 singleton HC3 NE',all(v['status'].startswith('NE') and v['p']=='NA' and v['n']=='18' for v in mods if v['exposure']=='H05_NECTIN2' and v['population']=='primary'))
idx=read(D/'FIGURE_INDEX.tsv');tiles=[]
for r in idx:
 key=r['figure'];dat=read(D/r['data']);ck(key+' data',bool(dat))
 with fitz.open(D/r['pdf']) as pdf:
  page=pdf[0];ck(key+' single vector page, fonts and labels',len(pdf)==1 and len(page.get_text().strip())>5 and len(page.get_fonts())>0)
  pix=page.get_pixmap(matrix=fitz.Matrix(.6,.6),alpha=False);ck(key+' renders',len(set(pix.samples[::97]))>2)
  if any(key.endswith(x) for x in ['core_overview','core_heatmap','paired_effects','HCV29/family','HeLa/FBS/family','forest_SREBF1_NECTIN2']):
   img=Image.frombytes('RGB',(pix.width,pix.height),pix.samples);img.thumbnail((390,235));tiles.append((key,img))
 with Image.open(D/r['png']) as im:ck(key+' 300 dpi',im.width>1000 and 298<=im.info['dpi'][0]<=302)
for a in ['GSE32894','GSE31684','GSE13507']:
 ds=[r for r in idx if r['figure'].startswith('bulk/'+a+'/')];ck(a+' independent presentation',len(ds)==13 and (D/'reports'/(a+'.md')).is_file())
ck('nine external REML',len([r for r in read(T/'meta_external.tsv') if r['family_id']=='core' and r['k']=='3'])==9)
sheet=Image.new('RGB',(800,255*((len(tiles)+1)//2)),'white');draw=ImageDraw.Draw(sheet)
for j,(key,im) in enumerate(tiles):
 x=(j%2)*400;y=(j//2)*255;sheet.paste(im,(x,y+18));draw.text((x+2,y+1),key,fill='black')
sheet.save(D/'qc/figure_review.png')
with (D/'qc/checks.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['check','passed']);w.writerows(checks)
summary=dict(status='complete',core_external_tests=27,external_specimens=566,figures=len(idx),checks=len(checks),virtual=json.loads((D/'qc/virtual_interface.json').read_text())['status'],optional_ChIP='not_started: standalone processed peaks unavailable',H05_HC3='NE: unit leverage')
(D/'completion_manifest.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2));print(summary)

