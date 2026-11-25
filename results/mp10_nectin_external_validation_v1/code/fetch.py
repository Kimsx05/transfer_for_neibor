from pathlib import Path
import subprocess,gzip,csv,os,hashlib,sys,urllib.parse,concurrent.futures
D=Path(__file__).resolve().parents[1]; I=D/'inputs'
accs=['GSE32894','GSE31684','GSE13507','GSE236742','GSE143803','GSE267018']
def get(url,name):
 p=I/name
 if p.exists() and p.stat().st_size:return p
 tmp=p.with_suffix(p.suffix+f'.{os.getpid()}.part')
 for base in [url,url.replace('ftp.ncbi.nlm.nih.gov','ftp.ncbi.nlm.nih.gov')]:
  r=subprocess.run(['/usr/bin/curl','-fLsS','--retry','1','--retry-all-errors','--connect-timeout','9','--max-time','1200','-o',str(tmp),'-w','%{http_code}',base],capture_output=True)
  if r.returncode==0:
   h=tmp.open('rb').read(200)
   if name.endswith('.gz') and h[:2]!=b'\x1f\x8b':continue
   if b'<html' in h.lower() and not name.endswith('.html'):continue
   tmp.replace(p);break
 with (D/'logs/download.tsv').open('a') as f:csv.writer(f,delimiter='\t').writerow([name,url,r.stdout.decode(),p.stat().st_size if p.exists() else 0,hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else '',r.stderr.decode()[-160:]])
 return p
def series(a):
 u=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{a[:-3]}nnn/{a}';p=get(u+'/matrix/'+a+'_series_matrix.txt.gz',a+'_series_matrix.txt.gz')
 return p
if __name__=='__main__':
 if not sys.argv[1:]:
  for a in accs:series(a)
 else:
  for u in sys.argv[1:]:get(u,u.split('/')[-1])
