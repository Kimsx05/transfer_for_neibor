"""Validated HTTP ranges and resumable minimal input recovery; no whole-record download."""
from pathlib import Path
import requests,json,io,zipfile,hashlib,time,re,sys,struct,zlib,binascii
ROOT=Path(__file__).resolve().parents[1]
S=requests.Session();S.trust_env=False;S.proxies={'https':'http://127.0.0.1:17890','http':'http://127.0.0.1:17890'}
def log(obj):
 with (ROOT/'logs/transfers.jsonl').open('a') as f:f.write(json.dumps(obj)+'\n')
 print(obj,flush=True)
def get_range(url,start,end,total):
 for attempt in range(3):
  try:
   r=S.get(url,headers={'Range':f'bytes={start}-{end}','Accept-Encoding':'identity'},stream=True,timeout=(15,90))
   expected=f'bytes {start}-{end}/{total}'
   if r.status_code!=206 or r.headers.get('Content-Range')!=expected:
    status=dict(url=url,status=r.status_code,content_range=r.headers.get('Content-Range'),expected=expected)
    r.close();raise RuntimeError('RANGE_NOT_HONORED '+str(status))
   b=r.content
   if len(b)!=end-start+1:raise IOError('range length mismatch')
   log(dict(url=url,start=start,end=end,bytes=len(b),content_range=expected,status=206));return b
  except Exception as e:
   log(dict(url=url,error=str(e),attempt=attempt))
   if attempt==2:raise
   time.sleep(2)
class RemoteZip(io.RawIOBase):
 def __init__(self,url,size):self.url=url;self.size=size;self.pos=0;self.cache=b'';self.begin=0
 def seekable(self):return True
 def readable(self):return True
 def tell(self):return self.pos
 def seek(self,offset,whence=0):
  self.pos=offset if whence==0 else self.pos+offset if whence==1 else self.size+offset
  return self.pos
 def read(self,n=-1):
  if n<0:n=self.size-self.pos
  n=min(n,self.size-self.pos)
  if n<=0:return b''
  if not(self.begin<=self.pos and self.pos+n<=self.begin+len(self.cache)):
   self.begin=self.pos;end=min(self.size-1,self.pos+max(n,1<<20)-1)
   self.cache=get_range(self.url,self.pos,end,self.size)
  b=self.cache[self.pos-self.begin:self.pos-self.begin+n];self.pos+=len(b);return b

def cervilla():
 meta=json.loads((ROOT/'metadata/zenodo_record.txt').read_text())
 entry=next(x for x in meta['files'] if x['key']=='visiumhd_Bladder.zip')
 url=entry['links']['self'];size=entry['size']
 # The first request is an explicit tiny range probe; reject any full body fallback.
 get_range(url,size-22,size-1,size)
 z=zipfile.ZipFile(RemoteZip(url,size));members=[dict(name=i.filename,bytes=i.file_size,compressed_bytes=i.compress_size,crc32=f'{i.CRC:08x}') for i in z.infolist()]
 (ROOT/'metadata/cervilla_zip_members.json').write_text(json.dumps(members,indent=2))
 print('ZIP_MEMBERS',len(members),'saved to metadata/cervilla_zip_members.json',flush=True)
 for i in z.infolist():
  n=i.filename
  wanted=('square_002um/' in n and n.endswith(('filtered_feature_bc_matrix.h5','tissue_positions.parquet','scalefactors_json.json','tissue_hires_image.png'))) or (i.file_size<20_000_000 and n.endswith(('metrics_summary.csv','web_summary.html','probe_set.csv','.mro','.json')) and 'analysis/' not in n)
  if not wanted or i.is_dir():continue
  dst=ROOT/'data/Cervilla'/n
  assert dst.resolve().is_relative_to((ROOT/'data/Cervilla').resolve())
  if dst.exists() and dst.stat().st_size==i.file_size:
   crc=0
   with dst.open('rb') as f:
    while b:=f.read(8<<20):crc=binascii.crc32(b,crc)
   if crc==i.CRC:continue
   raise IOError('Existing member has incorrect CRC: '+str(dst))
  dst.parent.mkdir(parents=True,exist_ok=True);tmp=dst.with_name(dst.name+'.part');compressed=dst.with_name(dst.name+'.deflate.part')
  # Central directory gives compressed size, while the local header gives payload offset.
  # Store only the active compressed member to permit byte-exact resume after interruption.
  header=get_range(url,i.header_offset,i.header_offset+29,size)
  if header[:4]!=b'PK\x03\x04':raise IOError('Invalid ZIP local header')
  name_len,extra_len=struct.unpack_from('<HH',header,26);data_start=i.header_offset+30+name_len+extra_len
  offset=compressed.stat().st_size if compressed.exists() else 0
  if offset>i.compress_size:raise IOError('Oversized partial compressed member')
  with compressed.open('ab') as out:
   while offset<i.compress_size:
    end=min(i.compress_size,offset+(8<<20));b=get_range(url,data_start+offset,data_start+end-1,size);out.write(b);out.flush();offset=end
  if i.compress_type not in (zipfile.ZIP_DEFLATED,zipfile.ZIP_STORED):raise ValueError('Unsupported ZIP method')
  decoder=zlib.decompressobj(-15) if i.compress_type==zipfile.ZIP_DEFLATED else None
  crc=0;h=hashlib.sha256();written=0
  with compressed.open('rb') as src,tmp.open('wb') as out:
   while b:=src.read(8<<20):
    raw=decoder.decompress(b) if decoder else b;written+=len(raw)
    if written>i.file_size:raise IOError('Uncompressed member exceeds declared size')
    out.write(raw);crc=binascii.crc32(raw,crc);h.update(raw)
   if decoder:
    raw=decoder.flush();written+=len(raw);out.write(raw);crc=binascii.crc32(raw,crc);h.update(raw)
    if not decoder.eof:raise IOError('Incomplete deflate stream')
  if written!=i.file_size or crc!=i.CRC:raise IOError('Member size/CRC mismatch')
  tmp.replace(dst);compressed.unlink()
  log(dict(member=n,local=str(dst),bytes=i.file_size,sha256=h.hexdigest(),zip_member_crc='PASS',whole_zip_md5_verified=False,removed_own_temporary=str(compressed)))

def geo():
 url='https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9040nnn/GSM9040999/suppl/GSM9040999_square_008um_filtered_feature_bc_matrix.h5'
 r=S.head(url,timeout=40);r.raise_for_status();size=int(r.headers['Content-Length'])
 (ROOT/'metadata/geo_h5_headers.json').write_text(json.dumps(dict(r.headers),indent=2))
 dst=ROOT/'data/GSE299573/GSM9040999/square_008um/filtered_feature_bc_matrix.h5';tmp=dst.with_name(dst.name+'.part')
 if dst.exists() and dst.stat().st_size==size:return
 dst.parent.mkdir(parents=True,exist_ok=True)
 start=tmp.stat().st_size if tmp.exists() else 0
 # Request a bounded chunk and verify Content-Range on every resumed chunk.
 with tmp.open('ab') as out:
  while start<size:
   end=min(size-1,start+(8<<20)-1);b=get_range(url,start,end,size);out.write(b);out.flush();start=end+1
 assert tmp.stat().st_size==size;tmp.replace(dst);log(dict(local=str(dst),bytes=size,download='complete'))
if __name__=='__main__':
 for name,func in [('Cervilla',cervilla),('GSE299573',geo)]:
  try:func()
  except Exception as e:log(dict(dataset=name,error=repr(e)))
