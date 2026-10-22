"""One sequential local gzip scan: recover only images/metadata, verify gzip CRC to EOF."""
from pathlib import Path
import tarfile,gzip,json,shutil,hashlib,time
ROOT=Path(__file__).resolve().parents[1]
archive=ROOT.parent/'highres_ST/highres_spatial_audit_v3_GSE300692/00_input_audit/extracted_source/GSM9066380_D1_Johnson.tar.gz'
rows=[]; start=time.time()
with gzip.open(archive,'rb') as gz:
 with tarfile.open(fileobj=gz,mode='r|') as tf:
  for m in tf:
   rows.append({'name':m.name,'bytes':m.size,'type':m.type.decode(),'linkname':m.linkname})
   n=m.name
   wanted=('/outs/spatial/' in n and '/.' not in n) or n.endswith(('web_summary.html','metrics_summary.csv','probe_set.csv'))
   # Keep only bounded, explicitly relevant regular members; no bulk expression duplicates.
   if not wanted or not m.isfile() or m.size>100_000_000:continue
   dst=ROOT/'data/GSE300692'/n
   if dst.exists():continue
   assert dst.resolve().is_relative_to((ROOT/'data/GSE300692').resolve())
   dst.parent.mkdir(parents=True,exist_ok=True); tmp=dst.with_name(dst.name+'.part')
   with tf.extractfile(m) as src,tmp.open('wb') as out:shutil.copyfileobj(src,out,8<<20)
   assert tmp.stat().st_size==m.size;tmp.replace(dst);print('RECOVERED',n,m.size,flush=True)
 # tar EOF is earlier than gzip EOF; consume remainder to check gzip trailer.
 while gz.read(8<<20):pass
(ROOT/'metadata/thor_archive_members.json').write_text(json.dumps(rows,indent=2))
(ROOT/'logs/thor_integrity.json').write_text(json.dumps({'archive':str(archive),'bytes':archive.stat().st_size,'gzip_crc_to_eof':'PASS','members':len(rows),'elapsed_s':time.time()-start},indent=2))
print('GZIP_CRC_PASS',len(rows),time.time()-start,flush=True)
