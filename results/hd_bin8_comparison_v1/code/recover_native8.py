from pathlib import Path
import sys,json,zipfile,struct,zlib,binascii,subprocess,re,time
root=Path(__file__).resolve().parents[1];hd=root.parents[1];sys.path.insert(0,str(hd/'input_recovery_followup'));import https_utils;https_utils.ROOT=root/'logs';from https_utils import RemoteZip,get_range
url='https://zenodo.org/records/17999961/files/visiumhd_Bladder.zip';size=3362702746;base='visiumhd_Bladder/binned_outputs/square_008um/'
for rel in ['filtered_feature_bc_matrix.h5','spatial/tissue_positions.parquet']:
 dst=root/'Cervilla/native8'/rel;dst.parent.mkdir(parents=True,exist_ok=True)
 if dst.exists():continue
 with zipfile.ZipFile(RemoteZip(url,size)) as z:i=z.getinfo(base+rel)
 h=get_range(url,i.header_offset,i.header_offset+29,size);nl,el=struct.unpack_from('<HH',h,26);start=i.header_offset+30+nl+el;tmp=dst.with_suffix('.compressed.part');offset=tmp.stat().st_size if tmp.exists() else 0
 attempts=0
 with tmp.open('ab') as out:
  while offset<i.compress_size:
   end=min(i.compress_size,offset+4*1024*1024);chunk=dst.with_suffix('.chunk');header=dst.with_suffix('.headers')
   r=subprocess.run(['curl','-sS','-L','--max-time','40','--max-filesize','5000000','--range',f'{start+offset}-{start+end-1}','https://zenodo.org/api/records/17999961/files/visiumhd_Bladder.zip/content','-o',str(chunk),'-D',str(header)])
   h=header.read_text() if header.exists() else '';matches=re.findall(r'content-range: bytes (\d+)-(\d+)/(\d+)',h,re.I);got=chunk.stat().st_size if chunk.exists() else 0
   if got and matches and tuple(map(int,matches[-1]))==(start+offset,start+end-1,size):
    assert got<=end-offset;out.write(chunk.read_bytes());out.flush();offset+=got;attempts=0;print(rel,offset,i.compress_size,flush=True)
   else:
    attempts+=1
    if attempts>=12:raise RuntimeError('Repeated HTTP range failure')
    time.sleep(3)
   chunk.unlink(missing_ok=True);header.unlink(missing_ok=True)
 dec=zlib.decompressobj(-15);crc=0;written=0
 with tmp.open('rb') as src,dst.with_suffix('.part').open('wb') as out:
  while b:=src.read(4*1024*1024):
   v=dec.decompress(b);out.write(v);crc=binascii.crc32(v,crc);written+=len(v)
  v=dec.flush();out.write(v);crc=binascii.crc32(v,crc);written+=len(v)
 assert dec.eof and written==i.file_size and crc==i.CRC;dst.with_suffix('.part').replace(dst);tmp.unlink()
 print('NATIVE_MEMBER_OK',rel,written,flush=True)
