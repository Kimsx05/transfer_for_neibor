"""Incremental small-metadata checks only; reuse previous search and archive evidence."""
from common import *
import requests,io,json,re,h5py,zipfile,struct,concurrent.futures,csv
P=HD/'input_recovery_followup/public_metadata';dest=ROOT/'metadata/additional_inputs';dest.mkdir(exist_ok=True)
fixed=['HMGCS1','MSMO1','INSIG1','FDFT1','SQLE','FDPS','IDI1','ACSS2','ELOVL5','LDLR','MVD','SCD','HMGCR','ACAT2','SLC20A1','ACSL3','ACSL1','ADIPOR2','MXD1','NPC1']
class HTTPSFile(io.RawIOBase):
 def __init__(self,url,total=None,base=0,length=None):
  self.url=url;self.s=requests.Session();self.s.trust_env=False;self.base=base;self.pos=0;self.begin=0;self.cache=b'';self.transferred=0
  if total is None:
   x=self.s.get(url,headers={'Range':'bytes=0-0'},timeout=(10,30));x.raise_for_status();assert x.status_code==206;total=int(x.headers['Content-Range'].split('/')[-1])
  self.total=total;self.size=length or total
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.pos
 def seek(self,n,whence=0):self.pos=n if whence==0 else self.pos+n if whence==1 else self.size+n;return self.pos
 def readinto(self,b):d=self.read(len(b));b[:len(d)]=d;return len(d)
 def read(self,n=-1):
  n=min(self.size-self.pos,n if n>=0 else self.size-self.pos)
  if n<=0:return b''
  if not self.begin<=self.pos<=self.pos+n<=self.begin+len(self.cache):
   a=self.base+self.pos;b=min(self.base+self.size-1,a+max(n,65536)-1)
   for attempt in range(3):
    try:
     proxy=None if attempt==0 else {'https':'http://10.10.10.215:7897'};x=self.s.get(self.url,headers={'Range':f'bytes={a}-{b}','Accept-Encoding':'identity'},proxies=proxy,timeout=(10,35));assert x.status_code==206 and x.headers.get('Content-Range')==f'bytes {a}-{b}/{self.total}';assert len(x.content)==b-a+1;break
    except Exception:
     if attempt==2:raise
   self.begin=self.pos;self.cache=x.content;self.transferred+=len(self.cache)
  d=self.cache[self.pos-self.begin:self.pos-self.begin+n];self.pos+=len(d);return d

def panel(acc):
 try:
  t=(P/(acc+'.txt')).read_text();url=re.search(r'!Sample_supplementary_file = (\S+cell_feature_matrix.h5)',t).group(1).replace('ftp://','https://');f=HTTPSFile(url)
  with h5py.File(f,'r') as h:
   g=h['matrix/features'];names=g['name'].asstr()[:];types=g['feature_type'].asstr()[:];genes=set(names[types=='Gene Expression']);out=dict(sample=acc,url=url,remote_bytes=f.total,bytes_read=f.transferred,gene_targets=len(genes),top20_present=[z for z in fixed if z in genes],top20_missing=[z for z in fixed if z not in genes],counts_downloaded=False)
  jwrite(dest/(acc+'_features.json'),out);print('PANEL',acc,len(genes),len(out['top20_present']),f.transferred,flush=True);return out
 except Exception as e:return dict(sample=acc,error=repr(e))
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:panels=list(ex.map(panel,['GSM9625559','GSM9625560','GSM9625561','GSM9625562']))
jwrite(dest/'prime5k_panel_checks.json',panels)
try:
 m=json.loads((P/'zenodo_17644360.txt').read_text());zmeta=m['files'][0];f=HTTPSFile(zmeta['links']['self'],zmeta['size']);outer=zipfile.ZipFile(f);result=[]
 for name in ['LCCC_553_Xenium/cells.zarr.zip','LCCC_553_Xenium/transcripts.zarr.zip']:
  info=outer.getinfo(name);assert info.compress_type==zipfile.ZIP_STORED;f.seek(info.header_offset);header=f.read(30);n,e=struct.unpack_from('<HH',header,26);offset=info.header_offset+30+n+e;inner=HTTPSFile(f.url,f.total,base=offset,length=info.file_size)
  with zipfile.ZipFile(inner) as z:
   members=[dict(name=i.filename,bytes=i.file_size,compressed_bytes=i.compress_size) for i in z.infolist()];jwrite(dest/(Path(name).name+'_members.json'),members);small={}
   for i in z.infolist():
    if i.file_size<80000 and i.filename.endswith(('.zattrs','.zarray','.zgroup')):
     if len(small)>=40:break
     small[i.filename]=json.loads(z.read(i))
   result.append(dict(container=name,outer_member_compression='stored; byte offsets verified from local header',inner_members=len(members),metadata=small,selected_bytes_read=inner.transferred))
 jwrite(dest/'UNC_nested_zarr_inspection.json',result);print('UNC_NESTED_METADATA_DONE',[(z['container'],z['inner_members']) for z in result],flush=True)
except Exception as e:jwrite(dest/'UNC_nested_zarr_error.json',dict(error=repr(e)));print('UNC',repr(e),flush=True)
