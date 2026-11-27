from base_data import *
files=sorted((O/'figures').glob('*.png'));rows=[]
for f in files:
 p=f.with_suffix('.pdf');assert p.exists() and min(f.stat().st_size,p.stat().st_size)>1000
 rows.append((f.stem,str(f),str(p)))
write(pd.DataFrame(rows,columns=['figure','png','pdf']),'figures.tsv')
priority=[]
for s in HD:
 a=pd.read_csv(O/'tables'/(s+'_GEX_examples.tsv'),sep='\t')
 for v in a.itertuples():
  priority.append((s,v.original_cell_id,'GEX '+v.selection,'有无独立核；与腔缘或相邻深染核的关系',str(O/'figures'/(s+'_GEX_'+v.selection+'.png'))))
for s,name in [('A2R2','conflict'),('A1R4','unresolved'),('A1R6','unresolved')]:
 a=pd.read_csv(O/'tables/xenium_review_ROI.tsv',sep='\t');r=a[(a['sample']==s)&(a.rationale==name)].iloc[0]
 priority.append((s,f'{r.x0_um},{r.y0_um} µm + 140 µm square',name,'分割边界与竞争谱系证据能否支持大类',str(O/'figures'/(s+'_ROI_'+name+'_boundaries.png'))))
write(pd.DataFrame(priority,columns=['sample','id_or_window','selection','review_question','image']),'priority_review.tsv')

