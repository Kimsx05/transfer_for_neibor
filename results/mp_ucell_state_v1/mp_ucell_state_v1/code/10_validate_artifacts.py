from pathlib import Path
import struct,csv,hashlib,sys
R=Path('/home/data/t070721/codex_workspace/Bladder Metabolism/mp_ucell_state_v1')
rows=[]
for p in sorted(R.rglob('*.png')):
 with p.open('rb') as f:
  sig=f.read(8);width=height=mode=None;dpi_x=dpi_y=None;end=False
  while True:
   h=f.read(8)
   if len(h)<8:break
   n,t=struct.unpack('>I4s',h)
   if t in [b'IHDR',b'pHYs']:
    d=f.read(n)
    if t==b'IHDR':width,height,bit,mode,_,_,_=struct.unpack('>IIBBBBB',d)
    else:
     x,y,unit=struct.unpack('>IIB',d)
     if unit==1:dpi_x=x*.0254;dpi_y=y*.0254
   else:f.seek(n,1)
   f.read(4)
   if t==b'IEND':end=True;break
 pdf=p.with_suffix('.pdf');pdfok=False
 if pdf.exists():
  with pdf.open('rb') as f:
   start=f.read(5);f.seek(max(0,pdf.stat().st_size-2048));tail=f.read();pdfok=start==b'%PDF-' and b'%%EOF' in tail
 ok=sig==b'\x89PNG\r\n\x1a\n' and end and mode in [4,6] and dpi_x is not None and abs(dpi_x-300)<.1 and pdfok
 rows.append(dict(path=str(p),width=width,height=height,png_color_type=mode,dpi_x=dpi_x,dpi_y=dpi_y,has_alpha=mode in [4,6],png_complete=end,pdf_complete=pdfok,pass_check=ok))
with (R/'00_input_audit/figure_file_validation.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
# Every required main view for all five datasets and ALL.
required=['cluster','MP10','all14_raw_UCell','dataset','sample','source','annotation','original_cluster','nFeature_counts','nCount_counts','percent.mt','MP10_distribution','cluster14MP_raw_heatmap','cluster14MP_GLOBALz_heatmap','sample_cluster_proportions','cluster_composition_dataset','cluster_composition_source']
instances={'FULL_Z':'03_full_z','FULL_RAW':'04_full_raw','BAL_Z_seed42':'05_balanced/BAL_Z_seed42','BAL_RAW_seed42':'05_balanced/BAL_RAW_seed42'}
datasets=['CNP0000460','GSE326225','HRA003620','InhouseData','PRJNA662018']
missing=[]
for nm,d in instances.items():
 for ds in ['ALL']+datasets:
  base=R/d/'figures' if ds=='ALL' else R/'08_figures_by_dataset'/ds/nm
  for v in required:
   for ext in ['png','pdf']:
    p=base/f'{nm}__{ds}__res0.4__{v}.{ext}'
    if not p.exists():missing.append(str(p))
for seed in [43,44]:
 for scale in ['Z','RAW']:
  nm=f'BAL_{scale}_seed{seed}'
  for v in ['cluster','MP10']:
   for ext in ['png','pdf']:
    p=R/'05_balanced'/nm/'figures'/f'{nm}__ALL__res0.4__{v}.{ext}'
    if not p.exists():missing.append(str(p))
(R/'00_input_audit/missing_required_figures.txt').write_text('\n'.join(missing))
print('PNG/PDF pairs:',len(rows),'invalid:',sum(not x['pass_check'] for x in rows),'missing required:',len(missing))
if missing or not all(x['pass_check'] for x in rows):sys.exit(1)
(R/'logs/ARTIFACT_VALIDATION_COMPLETE').write_text('complete\n')
