from PIL import Image, ImageDraw
from pathlib import Path
import pandas as pd
r=Path(__file__).resolve().parents[1]
out=r/'figures'/'review_sheets';out.mkdir(exist_ok=True)
for s,start in [('Cervilla',44),('D1_Johnson',0)]:
    files=sorted((r/'figures/details').glob(s+'*.jpg'))[start:]
    for j in range(0,len(files),4):
        page=Image.new('RGB',(1920,1950),'white')
        for k,f in enumerate(files[j:j+4]):
            src=Image.open(f)
            rows=(src.height-24)//67
            for i in range(rows*20):
                col=i%20;row=i//20
                patch=src.crop((54*col+2,24+67*row+2,54*col+50,24+67*row+50))
                page.paste(patch,((k%2*20+col)*48,30+(k//2*20+row)*48))
        ImageDraw.Draw(page).text((5,5),' | '.join(f.stem for f in files[j:j+4])+'    Native 0.5 um/pixel   |   24 um per field',fill='black')
        page.save(out/f'{s}_{j//4+1:02d}.jpg',quality=90)
