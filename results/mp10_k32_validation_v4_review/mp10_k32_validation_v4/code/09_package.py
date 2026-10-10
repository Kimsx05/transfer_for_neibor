from pathlib import Path
import csv,struct,zlib,hashlib,zipfile,json,re,datetime,html
R=Path('/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_k32_validation_v4')
exec((R/'code/00_check_inputs.py').read_text())
checks=[]
for p in sorted(R.rglob('*.png')):
 b=p.read_bytes();assert b[:8]==b'\x89PNG\r\n\x1a\n';pos=8;chunks=[];ppm=None
 while pos<len(b):
  n=struct.unpack('>I',b[pos:pos+4])[0];k=b[pos+4:pos+8];v=b[pos+8:pos+8+n];crc=struct.unpack('>I',b[pos+8+n:pos+12+n])[0];assert zlib.crc32(k+v)&0xffffffff==crc
  if k==b'IHDR':w,h,depth,color,_,_,interlace=struct.unpack('>IIBBBBB',v)
  if k==b'pHYs':ppm=struct.unpack('>IIB',v)
  if k==b'IDAT':chunks.append(v)
  pos+=n+12
 assert color==6 and depth==8 and interlace==0,(str(p),color,depth,interlace)
 assert ppm and ppm[2]==1 and abs(ppm[0]*.0254-300)<.1
 raw=zlib.decompress(b''.join(chunks));assert len(raw)==h*(w*4+1);assert raw[4]==0,(str(p),'opaque corner')
 checks.append(dict(figure=str(p.relative_to(R)),width=w,height=h,dpi=round(ppm[0]*.0254,3),RGBA=True,transparent_corner=True,CRC_OK=True))
def wt(p,rows):
 with (R/p).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
wt('00_audit/PNG_validation.tsv',checks)
for x in re.findall(r'(?:href|src)="([^"]+)"',(R/'REPORT.html').read_text()):
 x=html.unescape(x)
 if not x.startswith(('http:','https:','#')):assert (R/x).exists(),x
accept=list(csv.DictReader((R/'FINAL_ACCEPTANCE.tsv').open(),delimiter='\t'));assert all(x['passed']=='TRUE' for x in accept)
# Source identity has not changed; historical original aggregate/object inputs remain outside the archive.
status=dict(status='COMPLETED',completed_at=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),candidate='K32',partition='FULL_Z_r1.00',cluster='C18',n_full=123699,n_K32=3282,models_reused=3,models_new=3,models_not_testable=0,figure_pairs=len(checks),acceptance_checks=len(accept),acceptance_assertions=sum(int(x['n_assertions']) for x in accept),patient_sensitivity='NOT_EVALUABLE: no independently verified mapping',scope='K32 single-cell Supplement review only; no re-clustering/spatial/CNV/prognosis',inputs='read-only sources unchanged')
(R/'COMPLETED.json').write_text(json.dumps(status,ensure_ascii=False,indent=2)+'\n')
gc=[]
for line in (R/'logs/02_C14_aggregate.log').read_text().splitlines():
 if line.startswith(('Ncells','Vcells')):
  try:gc.append((line.split()[0],float(line.split()[-1])))
  except ValueError:pass
heap=sum(max([v for k,v in gc if k==cell] or [0]) for cell in ['Ncells','Vcells']);size=sum(p.stat().st_size for p in R.rglob('*') if p.is_file() and p.suffix!='.zip')
wt('00_audit/resource_usage.tsv',[dict(metric='output_bytes_excluding_zip',value=size,meaning='Includes new C14 aggregates and new sample-expression caches, excluded from review package'),dict(metric='observed_R_allocated_heap_peak_MiB',value=heap,meaning='Aggregation gc maximum Ncells+Vcells, not OS RSS'),dict(metric='original_RNA_loads',value=1,meaning='Single sequential load to fill missing C14 sparse aggregates'),dict(metric='maximum_data_threads',value=2,meaning='Heavy models sequential, BLAS and fgsea one thread')])
files=[p for p in sorted(R.rglob('*')) if p.is_file() and p.suffix not in ['.qs','.zip'] and '__pycache__' not in str(p) and p.name not in ['review_package_manifest.tsv','09_package.log'] and not p.name.endswith('.tmp')]
manifest=[dict(relative_path=p.relative_to(R).as_posix(),bytes=p.stat().st_size,sha256=hashlib.file_digest(p.open('rb'),'sha256').hexdigest()) for p in files];wt('review_package_manifest.tsv',manifest);files.append(R/'review_package_manifest.tsv');zp=R/'mp10_k32_validation_v4_review.zip'
with zipfile.ZipFile(zp,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in files:z.write(p,arcname='mp10_k32_validation_v4/'+p.relative_to(R).as_posix())
with zipfile.ZipFile(zp) as z:
 assert z.testzip() is None
 assert not any(n.endswith('.qs') for n in z.namelist())
 for m in ['full/rest','full/C17','full/C7','full/C14','exclude_P07T/rest','exclude_P07T/C14']:
  assert 'mp10_k32_validation_v4/02_rna_DEG/'+m+'/DE_annotated.tsv.gz' in z.namelist()
  assert 'mp10_k32_validation_v4/02_rna_DEG/'+m+'/function/enrichment_annotated.tsv.gz' in z.namelist()
print(json.dumps(status,ensure_ascii=False));print('ZIP',zp,'bytes',zp.stat().st_size,'MiB',round(zp.stat().st_size/1024**2,2),'files',len(files));print('R allocated heap MiB',heap,'output bytes',size)
