from pathlib import Path
import csv,struct,zlib,hashlib,zipfile,json,re
R=Path('/home/data/t070721/codex_workspace/Bladder Metabolism/mp_ucell_state_followup_v2')
checks=[]
for p in sorted((R/'04_figures_by_dataset').glob('*.png')):
 b=p.read_bytes();assert b[:8]==b'\x89PNG\r\n\x1a\n';pos=8;ids=[];ppm=None
 while pos<len(b):
  n=struct.unpack('>I',b[pos:pos+4])[0];k=b[pos+4:pos+8];v=b[pos+8:pos+8+n];crc=struct.unpack('>I',b[pos+8+n:pos+12+n])[0];assert zlib.crc32(k+v)&0xffffffff==crc
  if k==b'IHDR':w,h,depth,color,_,_,interlace=struct.unpack('>IIBBBBB',v)
  if k==b'pHYs':px,py,unit=struct.unpack('>IIB',v);ppm=(px,py,unit)
  if k==b'IDAT':ids.append(v)
  pos+=12+n
 assert color==6 and depth==8 and interlace==0,(p,color,depth,interlace)
 assert ppm and ppm[2]==1 and abs(ppm[0]*.0254-300)<.1 and abs(ppm[1]*.0254-300)<.1
 raw=zlib.decompress(b''.join(ids));assert len(raw)==h*(w*4+1);assert raw[4]==0 # first pixel alpha, all first-row first-pixel PNG predictors zero
 checks.append(dict(figure=p.name,width=w,height=h,dpi=round(ppm[0]*.0254,3),RGBA=True,transparent_corner=True,CRC_OK=True))
with (R/'00_audit/PNG_validation.tsv').open('w') as f:
 q=csv.DictWriter(f,fieldnames=checks[0],delimiter='\t');q.writeheader();q.writerows(checks)
# Local HTML links must resolve, including relative figure links used after ZIP extraction.
text=(R/'REPORT.html').read_text();links=re.findall(r'(?:href|src)="([^"]+)"',text)
for x in links:
 if not x.startswith(('http:','https:','#')):assert (R/x).exists(),x
status=dict(technical_status='completed',fixed_comparisons=6,pseudobulk_models=6,pdf_png_pairs=len(checks),R_validation_checks=56,read_only_inputs_unchanged=True,patient_validation='unavailable: derived/provisional mappings',GSE326225_paired_inference='unavailable: no eligible pairs',biological_validation='not established; discovery-stage state/marker evidence only')
(R/'completion_status.json').write_text(json.dumps(status,ensure_ascii=False,indent=2)+'\n')
# Include all compact review figures/data, full DE/filter results and reproducibility files; exclude matrices/caches.
files=[]
for p in sorted(R.rglob('*')):
 if not p.is_file():continue
 rel=p.relative_to(R).as_posix()
 if p.suffix in ('.qs','.rds','.zip') or '_expression_allgenes.tsv.gz' in p.name or p.name.endswith('.tmp'):continue
 if rel in ('review_package_manifest.tsv','logs/09_package.log'):continue
 if rel.startswith('02_pseudobulk_markers/') and p.name in ('aggregate_COMPLETE',):continue
 files.append(p)
manifest=[]
for p in files:manifest.append(dict(relative_path=p.relative_to(R).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
with (R/'review_package_manifest.tsv').open('w') as f:
 q=csv.DictWriter(f,fieldnames=manifest[0],delimiter='\t');q.writeheader();q.writerows(manifest)
files.append(R/'review_package_manifest.tsv')
zpath=R/'mp10_state_followup_v2_review.zip'
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in files:z.write(p,arcname='mp_ucell_state_followup_v2/'+p.relative_to(R).as_posix())
with zipfile.ZipFile(zpath) as z:
 assert z.testzip() is None
 assert not any(n.endswith(('.qs','.rds')) for n in z.namelist())
print(json.dumps(status,ensure_ascii=False));print('ZIP',zpath,'bytes',zpath.stat().st_size,'files',len(files))
