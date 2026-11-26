from common import *
from graph import polygons
from geometry import pixel
from PIL import Image,ImageDraw,ImageFont
import html
from scipy.spatial import cKDTree
Image.MAX_IMAGE_PIXELS=None
JS="""let num=document.getElementById('idx'),sel=document.getElementById('items'),cv=document.querySelector('canvas'),ctx=cv.getContext('2d');let tile={},imgs=new Image(),loaded={};
for(let i=0;i<index.length;i++){let op=new Option(index[i][0]+' / '+index[i][3],i);sel.add(op)}
if(cfg.img)imgs.src=cfg.img;
async function show(){const r=index[Number(sel.value)],cx=r[1],cy=r[2];num.value=sel.selectedIndex+1;let files=[];for(let i=Math.floor((cx-54)/256);i<=Math.floor((cx+54)/256);i++)for(let j=Math.floor((cy-54)/256);j<=Math.floor((cy+54)/256);j++){let k=i+'_'+j;if(!exists.includes(k))continue;if(!loaded[k])loaded[k]=new Promise(ok=>{let e=document.createElement('script');e.src='tiles/'+k+'.js';e.onload=ok;document.body.appendChild(e)});files.push([k,loaded[k]])}await Promise.all(files.map(a=>a[1]));ctx.clearRect(0,0,cv.width,cv.height);let rows=files.flatMap(a=>tile[a[0]]).filter(p=>Math.abs(p[0]-cx)<55&&Math.abs(p[1]-cy)<55);
for(let k=0;k<8;k++){let xx=(k%4)*200,yy=Math.floor(k/4)*216;ctx.save();ctx.translate(xx,yy+20);ctx.beginPath();ctx.rect(0,0,200,192);ctx.clip();ctx.fillStyle='#f5f4f0';ctx.fillRect(0,0,200,192);if(!k&&imgs.complete&&imgs.naturalWidth)ctx.drawImage(imgs,(cx-50)*cfg.sc,(cy-48)*cfg.sc,100*cfg.sc,96*cfg.sc,0,0,200,192);
for(let p of rows){let poly=p[2],v=p[3+k-1];ctx.beginPath();poly.forEach((q,j)=>{let x=100+2*(q[0]-cx),y=96+2*(q[1]-cy);if(j)ctx.lineTo(x,y);else ctx.moveTo(x,y)});ctx.closePath();ctx.strokeStyle=k?'#bdbdbd':'#596879';ctx.lineWidth=.55;if(k){ctx.fillStyle=v==null?'#bfb3ca':v>0?'rgba(137,17,61,'+(0.35+0.65*Math.min(1,Math.log1p(v)/Math.log1p(k===7?800:5)))+')':'#fafafa';ctx.fill()}ctx.stroke()}if(cfg.hd){ctx.strokeStyle='#0c5371';ctx.lineWidth=1.2;for(let p of rows){let box=p[10];if(!box)continue;ctx.beginPath();box.forEach((q,j)=>{let x=100+2*(q[0]-cx),y=96+2*(q[1]-cy);if(j)ctx.lineTo(x,y);else ctx.moveTo(x,y)});ctx.closePath();ctx.stroke()}}ctx.strokeStyle='#0375a0';ctx.lineWidth=1.6;ctx.strokeRect(96,92,8,8);ctx.fillStyle='black';ctx.fillRect(9,182,40,2);ctx.restore();ctx.font='14px sans-serif';ctx.fillStyle='black';ctx.fillText(cfg.names[k],xx+4,yy+15)}}
sel.onchange=show;num.onchange=()=>{sel.selectedIndex=Math.min(index.length-1,Math.max(0,num.value-1));show()};imgs.onload=show;show();"""
def run(s):
 hd=s in HD;q=paths(s);v,d=folders(s);z=labels(s);idt=pd.read_parquet(d/'identity.parquet');hu=pd.read_csv(q/'human_review.tsv',sep='\t')
 ids=idt.index[idt.candidate_union|idt.index.isin(hu.id)];idx=idt.loc[ids,['x_um','y_um','grade']].reset_index();idx['window_index']=np.arange(1,len(idx)+1);put(idx,q/'window_index.tsv')
 keys=CORE+['pan-T','competition','UMI'];names=['H&E + grid' if hd else 'Segmentation',*CORE,'pan-T','NK / cytotoxic','UMI']; img=None;sc=1
 if hd:
  fn=imagedir(s)/'tissue_hires_image.png';im=Image.open(fn);cfg=sf(s);sc=cfg['tissue_hires_scalef']/cfg['microns_per_pixel'];img=im.convert('RGB')
  p=json.loads((v/'conservation.json').read_text());x,gids,feat=tenx(p['input']);gs=feat.name.values
  b=pd.read_parquet(ROOT/'highres_VisiumHD/analysis/hd_epi_niche_dual_methods_v1'/s/'nodes.parquet',columns=['bin_id','x_um','y_um','array_row','array_col']).set_index('bin_id').loc[gids]
  aff=np.load(v/'window.npz')['affine'];inv=np.linalg.inv(aff[:2]);offs=np.array([[-.5,-.5],[-.5,.5],[.5,.5],[.5,-.5]])
  rc=b[['array_row','array_col']].values;xy=b[['x_um','y_um']].values;ps=xy[:,None,:]+(offs@inv)[None,:,:];center=((rc//2)*2+.5-aff[2])@inv;boxes=center[:,None,:]+2*(offs@inv)[None,:,:]
 else:
  x,gs,gids=expression(s);xy=z[['x_um','y_um']].values;ps=[np.asarray(poly.exterior.coords) for poly in polygons(s,z.index)];boxes=None
 a=raw(x,gs,gids,CORE+PAN+['NKG7','GNLY','CD8A']);a['pan-T']=a[PAN].sum(axis=1,min_count=1);a['competition']=a[['NKG7','GNLY','CD8A']].sum(axis=1,min_count=1);a['UMI']=x.sum(1).A1;vals=a[keys].values
 xy=np.round(xy,2);tree=cKDTree(xy);td=q/'tiles';td.mkdir(exist_ok=True);wanted=set()
 for cx,cy in idx[['x_um','y_um']].values:
  for i in range(int((cx-55)//256),int((cx+55)//256)+1):
   for j in range(int((cy-55)//256),int((cy+55)//256)+1):wanted.add((i,j))
 valid=[];xyq=np.floor(xy/256).astype(int)
 for ij,rr in pd.DataFrame(xyq).groupby([0,1]).indices.items():
  if ij not in wanted:continue
  rows=[]
  for i in rr:
   vs=[None if not np.isfinite(n) else int(n) for n in vals[i]]
   row=[*xy[i],np.round(ps[i],1).tolist(),*vs]
   if hd:row.append(np.round(boxes[i],1).tolist())
   rows.append(row)
  k='_'.join(map(str,ij));valid.append(k);(td/(k+'.js')).write_text('tile['+json.dumps(k)+']='+json.dumps(rows,separators=(',',':')))
 config={'hd':hd,'img':os.path.relpath(fn,q) if hd else None,'sc':sc,'names':names}
 (q/'atlas.html').write_text('<meta charset="utf-8"><title>Blinded review '+s+'</title><h2>'+s+' · provisional</h2><label>Window <input id="idx" type="number" min="1" value="1"></label> <select id="items"></select><br><canvas width="800" height="440"></canvas><p>Center marked in blue. Scale bar 20 µm. x right, y down. White = zero, purple-gray = unmeasured. Pan-T and competing signals are sums in the displayed original unit. UMI has a separate scale.</p><script>const index='+idx[['id','x_um','y_um','grade']].to_json(orient='values')+',cfg='+json.dumps(config)+',exists='+json.dumps(valid)+';'+JS+'</script>')
 ff=q/'cards';ff.mkdir(exist_ok=True)
 for id0 in hu.id:
  cx,cy=z.loc[id0,['x_um','y_um']];near=tree.query_ball_point([cx,cy],77);out=Image.new('RGB',(900,496),'white');draw=ImageDraw.Draw(out)
  for k,title in enumerate(names):
   panel=Image.new('RGB',(216,208),'#fafafa');pn=ImageDraw.Draw(panel)
   if k==0 and hd:panel=img.crop(((cx-54)*sc,(cy-52)*sc,(cx+54)*sc,(cy+52)*sc)).resize(panel.size);pn=ImageDraw.Draw(panel)
   for i in near:
    poly=[(108+2*(x0-cx),104+2*(y0-cy)) for x0,y0 in ps[i]];vv=vals[i,k-1] if k else 0
    fill=None if not k else '#c2afc7' if np.isnan(vv) else tuple(int(250*(1-t)+c*t) for c in [135,15,55]) if (t:=min(1,np.log1p(vv)/np.log1p(800 if k==7 else 5))) else '#ffffff'
    pn.polygon(poly,fill=fill,outline='#a3a9ae')
   if hd:
    for i in near:pn.line([(108+2*(x0-cx),104+2*(y0-cy)) for x0,y0 in np.vstack([boxes[i],boxes[i][0]])],fill='#5b727f',width=1)
   pn.rectangle((104,100,112,108),outline='#0083af',width=2);pn.line((7,201,47,201),fill='black',width=2)
   xx=k%4*225;yy=26+k//4*230;out.paste(panel,(xx,yy+17));draw.text((xx+4,yy),title,fill='black')
  draw.text((4,4),s+'   '+str(id0)+'   |   Provisional, blinded   |   scale 20 um, x right / y down',fill='black');out.save(ff/(id0+'.png'));out.save(ff/(id0+'.pdf'),resolution=145)
 (q/'index.html').write_text('<meta charset="utf-8"><h2>'+s+' human review</h2><a href="atlas.html">All candidate windows</a> · <a href="human_review.tsv">Review form</a><ol>'+''.join('<li><a href="cards/'+i+'.png">'+html.escape(i)+'</a> '+st+'</li>' for i,st in zip(hu.id,hu.sampling_stratum))+'</ol>')
 print(s,len(idx),len(hu),flush=True)
if __name__=='__main__':
 for s in sys.argv[1:] or HD+XE:run(s)
