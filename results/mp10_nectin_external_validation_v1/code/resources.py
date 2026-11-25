from fetch import *
import re
def run():
 for a in ['GSE236742','GSE143803','GSE267018']:
  u=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{a[:-3]}nnn/{a}'
  for sub in ['suppl','matrix']:
   ix=D/'logs'/f'{a}_{sub}.html'
   if not ix.exists():subprocess.run(['/usr/bin/curl','-fsSL','--retry','2',u+'/'+sub+'/', '-o',str(ix)],check=True)
   for name in re.findall('href="([^"]+)"',ix.read_text()):
    if name.startswith(a) and ('GPL24676_' in name or a=='GSE143803' and 'series_matrix' in name or any(s in name for s in ['Rawcounts','readCounts','hela_fpkm'])):get(u+'/'+sub+'/'+name,name)
 for a in accs[:3]:
  with gzip.open(I/(a+'_series_matrix.txt.gz'),'rt') as f:
   for line in f:
    if line.startswith('!Series_supplementary_file'):
     u=line.strip().split('\t')[1].strip('"').replace('ftp://','https://');n=u.split('/')[-1]
     if 'non-normalized' in n or 'preprocess_metadata' in n or 'clinical' in n:get(u,n)
    if line.startswith('!Sample_'):break
 for a in ['GPL6947','GPL6102']:
  get(f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={a}&targ=self&form=text&view=full',a+'.soft')
if __name__=='__main__':run()
