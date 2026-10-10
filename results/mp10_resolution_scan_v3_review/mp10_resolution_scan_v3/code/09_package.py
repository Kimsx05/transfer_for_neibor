from pathlib import Path
import csv,struct,zlib,hashlib,zipfile,json,re,datetime
R=Path('/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_resolution_scan_v3')
checks=[]
for p in sorted((R/'07_figures').glob('*.png')):
 b=p.read_bytes();assert b[:8]==b'\x89PNG\r\n\x1a\n';pos=8;chunks=[];ppm=None
 while pos<len(b):
  n=struct.unpack('>I',b[pos:pos+4])[0];k=b[pos+4:pos+8];v=b[pos+8:pos+8+n];crc=struct.unpack('>I',b[pos+8+n:pos+12+n])[0];assert zlib.crc32(k+v)&0xffffffff==crc
  if k==b'IHDR':w,h,depth,color,_,_,interlace=struct.unpack('>IIBBBBB',v)
  if k==b'pHYs':ppm=struct.unpack('>IIB',v)
  if k==b'IDAT':chunks.append(v)
  pos+=n+12
 assert color==6 and depth==8 and interlace==0,(str(p),color,depth,interlace)
 assert ppm and ppm[2]==1 and abs(ppm[0]*.0254-300)<.1
 raw=zlib.decompress(b''.join(chunks));assert len(raw)==h*(w*4+1);assert raw[4]==0,(str(p),'opaque first pixel')
 checks.append(dict(figure=p.name,width=w,height=h,dpi=round(ppm[0]*.0254,3),RGBA=True,transparent_corner=True,CRC_OK=True))
with (R/'00_audit/PNG_validation.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=checks[0],delimiter='\t');w.writeheader();w.writerows(checks)
report=(R/'REPORT.html').read_text()
for x in re.findall(r'(?:href|src)="([^"]+)"',report):
 if not x.startswith(('http:','https:','#')):assert (R/x).exists(),x
with (R/'FINAL_ACCEPTANCE.tsv').open() as f:accept=list(csv.DictReader(f,delimiter='\t'))
assert all(x['status']=='PASS' for x in accept)
status=dict(status='COMPLETED',completed_at=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),partitions=24,candidates=38,comparisons=108,figure_pairs=len(checks),acceptance_checks=len(accept),scope='Single-cell review complete; no downstream spatial/CNV/prognosis',read_only_inputs='unchanged')
(R/'COMPLETED.json').write_text(json.dumps(status,ensure_ascii=False,indent=2)+'\n')
# Actual R allocated-heap maximum is available in aggregation gc output; distinguish from RSS.
gc=[]
for line in (R/'logs/03_aggregate.log').read_text().splitlines():
 if line.startswith(('Ncells','Vcells')):
  fields=line.split()
  try:gc.append((fields[0],float(fields[-1])))
  except ValueError:pass
heap=sum(max([v for k,v in gc if k==cell] or [0]) for cell in ['Ncells','Vcells'])
size=sum(p.stat().st_size for p in R.rglob('*') if p.is_file() and p.suffix!='.zip')
resources=[dict(metric='output_bytes_before_review_zip',value=size,meaning='All v3 outputs including pseudobulk/sample-expression caches'),dict(metric='observed_R_allocated_heap_peak_MiB',value=heap,meaning='Aggregation gc maximum Ncells plus Vcells, not OS RSS'),dict(metric='estimated_RNA_process_peak_GiB',value='30-40',meaning='Pre-run estimate from same original Seurat in v2'),dict(metric='BLAS_threads',value=1,meaning='Sequential heavy models; no replicated original RNA workers')]
with (R/'00_audit/resource_usage.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=resources[0],delimiter='\t');w.writeheader();w.writerows(resources)
# Compact report package keeps every full DE/enrichment result, selected figures and all small review tables.
keepfig=set()
with (R/'07_figures/figure_index.tsv').open() as f:
 for row in csv.DictReader(f,delimiter='\t'):
  n=row['figure']
  if n.startswith('ALL__') or any(t in n for t in ['marker_DotPlot','sample_precision_recall','candidate_source_composition','candidate_PR','RNA__candidate_highlights']):keepfig.add(n)
files=[]
for p in sorted(R.rglob('*')):
 if not p.is_file():continue
 rel=p.relative_to(R).as_posix()
 if p.suffix in ('.qs','.zip') or '__pycache__' in rel or p.name.endswith('.tmp'):continue
 if rel in ('review_package_manifest.tsv','logs/09_package.log'):continue
 if rel.startswith('07_figures/') and p.suffix in ('.pdf','.png') and p.stem not in keepfig:continue
 if rel.startswith('07_figures/') and p.name.endswith('_plotdata.tsv.gz') and p.name[:-16] not in keepfig:continue
 # Complete rank/filter tables remain useful evidence; full per-gene descriptive arrays and large caches are excluded.
 if p.name=='descriptive_allgenes.tsv.gz':continue
 files.append(p)
manifest=[dict(relative_path=p.relative_to(R).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files]
with (R/'review_package_manifest.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=manifest[0],delimiter='\t');w.writeheader();w.writerows(manifest)
files.append(R/'review_package_manifest.tsv');zp=R/'mp10_resolution_scan_v3_review.zip'
with zipfile.ZipFile(zp,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in files:z.write(p,arcname='mp10_resolution_scan_v3/'+p.relative_to(R).as_posix())
with zipfile.ZipFile(zp) as z:
 assert z.testzip() is None
 assert not any(n.endswith('.qs') for n in z.namelist())
print(json.dumps(status,ensure_ascii=False));print('Review ZIP',zp,'bytes',zp.stat().st_size,'files',len(files));print('Observed allocated heap MiB',heap,'output bytes',size)
