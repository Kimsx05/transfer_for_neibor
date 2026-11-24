from pathlib import Path
import lzma
import subprocess as sp, os, tempfile, json, gzip, zipfile, hashlib, sys, fcntl, csv, time
O=Path(__file__).resolve().parents[1]
C=Path.cwd()/'cache/mp10_public_resources_v1'
items=[]
def add(d,n,u): items.append((d,n,u))
u='https://d1hp0ufb0xxisr.cloudfront.net/CTD2/Broad/CTRPv2.1_2016_pub_NatChemBiol_12_109.zip'
add('CTRPv2.1',u.split('/')[-1],u)
add('DepMap_24Q4_v1','manifest.json','https://api.figshare.com/v2/articles/27993248/versions/1')
for n,i in [('Model.csv',51065297),('README.txt',51065795),('CRISPRGeneEffect.csv',51064667),('OmicsExpressionProteinCodingGenesTPMLogp1.csv',51065489)]:
 add('DepMap_24Q4_v1',n,'https://ndownloader.figshare.com/files/'+str(i))
for a,ns in {'GSE182638':['raw_readcounts_allsamples.txt.gz'],'GSE126696':['scores_read_count_matrix_786O_ML210_time_course.txt.gz','CP0041_reference_20160112.csv.gz'],'GSE175862':['Bernareggi_et_al.__CRISPR_Screen_Hits_and_scores.csv.gz'],'GSE175863':['CHMP2A_RNAseq_Abundance_and_Differential_Expression.xlsx'],'GSE149933':['mTKO_library.txt.gz','mTKO_readCounts_4May20.txt.gz','mVal_library.txt.gz','mVal_readCounts_4May20.txt.gz']}.items():
 for n in ns:
  name=a+'_'+n
  add(a,name,'https://ftp.ncbi.nlm.nih.gov/geo/series/'+a[:-3]+'nnn/'+a+'/suppl/'+name)
 name=a+'_series_matrix.txt.gz'
 add(a,name,'https://ftp.ncbi.nlm.nih.gov/geo/series/'+a[:-3]+'nnn/'+a+'/matrix/'+name)
add('Brunello_2016','broadgpp-brunello-library-contents.txt','https://media.addgene.org/cms/filer_public/8b/4c/8b4c89d9-eac1-44b2-bb2f-8fea95672705/broadgpp-brunello-library-contents.txt')
add('HMGCS1_2026','41467_2026_74022_MOESM4_ESM.xlsx','https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-026-74022-y/MediaObjects/41467_2026_74022_MOESM4_ESM.xlsx')
add('Lawson2020','tree.json','https://api.github.com/repos/NMikolajewicz/Lawson2020/git/trees/master?recursive=1')
def validate(p,n,d):
 if n.endswith(('.zip','.xlsx')):
  with zipfile.ZipFile(p) as z: assert z.testzip() is None
 elif n.endswith('.gz'):
  with gzip.open(p) as g:
   while g.read(1048576): pass
 elif n.endswith('.json'): json.loads(p.read_bytes())
 else:
  with p.open('rb') as f: assert not f.read(256).lstrip().startswith(b'<')
 if d=='DepMap_24Q4_v1' and n!='manifest.json':
  m=json.loads((p.parent/'manifest.json').read_bytes()); meta=next(v for v in m['files'] if v['name']==n)
  if p.suffix=='.xz':
   h=hashlib.md5();size=0
   with lzma.open(p) as f:
    for block in iter(lambda:f.read(4<<20),b''):h.update(block);size+=len(block)
   assert h.hexdigest()==meta['computed_md5'] and size==meta['size']
  else:assert p.stat().st_size==meta['size'] and digest(p,'md5')==meta['computed_md5']
def digest(p,a):
 h=hashlib.new(a)
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4<<20),b''): h.update(b)
 return h.hexdigest()
if __name__=='__main__':
 rows=[]
 for d,n,u in items:
  dest=C/d/n; dest.parent.mkdir(parents=True,exist_ok=True); tmp=None
  try:
   if not dest.exists() and (dest.parent/(n+'.xz')).exists():dest=dest.parent/(n+'.xz')
   if dest.exists(): validate(dest,n,d); state='reused'
   else:
    with tempfile.NamedTemporaryFile(delete=False,dir=dest.parent,prefix='.'+n+'.') as f: tmp=Path(f.name)
    r=sp.run(['curl','-fLsS','--retry','1','--connect-timeout','12','--max-time','600','-o',str(tmp),u],capture_output=True)
    if r.returncode: raise RuntimeError(r.stderr.decode()[-250:])
    validate(tmp,n,d); os.replace(tmp,dest);state='downloaded'
   rows.append([d,n,u,state,dest.stat().st_size,digest(dest,'sha256')])
  except Exception as e: rows.append([d,n,u,'failed',0,repr(e)])
  finally:
   if tmp and tmp.exists(): tmp.unlink()
  print(n,rows[-1][3],str(rows[-1][-1])[:150],flush=True)
  with (O/'tables/resources.tsv').open('w') as f:
   w=csv.writer(f,delimiter='\t'); w.writerow(['resource','file','url','state','bytes','sha256_or_error']);w.writerows(rows)

