"""Fetch GEO processed data only. Cache a verified file by atomic rename under a lock."""
from pathlib import Path
import sys,gzip,tarfile,hashlib,subprocess,tempfile,csv,fcntl,os
out=Path(__file__).resolve().parents[1]
cache=Path('cache/mp10_public_resources_v1')
geo='https://ftp.ncbi.nlm.nih.gov/geo'
resources=[]
for acc in ['GSE149573','GPL23038','GSE165257','GPL16570','GSE139325','GSE243021']:
    series=acc.startswith('GSE')
    url=f"{geo}/{'series' if series else 'platforms'}/{acc[:-3]}nnn/{acc}"
    if series:
        resources.append((acc,url+'/matrix/'+acc+'_series_matrix.txt.gz','series'))
    else:
        resources.append((acc,url+'/soft/'+acc+'_self.soft.gz','platform'))
    if acc=='GSE139325':resources.append((acc,url+'/suppl/GSE139325_RAW.tar','tar'))
    if acc=='GSE243021':resources.append((acc,url+'/suppl/GSE243021_Spl.VAT.WT.KO.RNAseq.rpkm.detected.csv.gz','csv'))
resources.append(('MGI','https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt','orthology'))
def verify(p,fmt):
    if fmt=='tar':
        with tarfile.open(p) as tar:
            members=tar.getnames()
            assert all(any(x in n for n in members) for x in ['GSM4138181','GSM4138186'])
            for m in tar.getmembers():
                if m.isfile():
                    with gzip.open(tar.extractfile(m),'rt') as f: f.read()
        return ','.join(members)
    if fmt=='orthology':
        h=p.open().readline().strip()
        assert 'Taxon' in h and 'Symbol' in h
        return h
    with gzip.open(p,'rt',errors='replace') as f:d=f.read()
    assert '<html' not in d[:300].lower()
    if fmt=='csv':
        assert 'VAT' in d.splitlines()[0]
        return d.splitlines()[0]
    key='!platform_table_begin' if fmt=='platform' else '!Sample_geo_accession'
    assert key in d
    return d.split(key,1)[1].strip().splitlines()[0]
pins={r['file']:r['sha256'] for r in csv.DictReader((out/'code/resource_hashes.tsv').open(),delimiter='\t')}
manifest=out/'resource_manifest.tsv'
fields=['accession','url','http','format','bytes','sha256','columns_or_members','path','status']
def record(row):
    new=not manifest.exists()
    with manifest.open('a') as f:
        w=csv.writer(f,delimiter='\t'); 
        if new:w.writerow(fields)
        w.writerow(row)
for acc,url,fmt in resources:
    if len(sys.argv)>1 and acc not in sys.argv[1:]:continue
    name=url.split('/')[-1]; folder=cache/acc/'snapshot_20261124';folder.mkdir(parents=True,exist_ok=True)
    p=folder/name
    with (folder/(name+'.lock')).open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        urls=[url]
        if fmt=='platform':urls=[f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&acc={acc}&form=text&view=full']
        if acc=='GPL16570':urls.extend(['https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16570&targ=self&form=text&view=full',f'{geo}/platforms/GPL16nnn/GPL16570/annot/GPL16570.annot.gz'])
        if fmt in ['csv','tar']:urls.append(f'https://www.ncbi.nlm.nih.gov/geo/download/?acc={acc}&format=file&file={name}')
        for u in urls:
            status='cached';http='cached';temp=None
            try:
                if p.exists(): columns=verify(p,fmt)
                else:
                    fd,n=tempfile.mkstemp(dir=folder,suffix='.part');os.close(fd);temp=Path(n)
                    r=subprocess.run(['/usr/bin/curl','-sSL','--retry','2','--retry-all-errors','--retry-delay','2','--connect-timeout','12','--max-time','150','-w','%{http_code}','-o',str(temp),u],capture_output=True,text=True)
                    http=r.stdout
                    if r.returncode or http!='200':raise ValueError(f'curl {r.returncode} HTTP {http}')
                    if fmt=='platform' and temp.open('rb').read(12).startswith(b'^PLATFORM'):
                        raw=temp.read_bytes();temp.write_bytes(gzip.compress(raw,mtime=0))
                    columns=verify(temp,fmt)
                    assert hashlib.sha256(temp.read_bytes()).hexdigest()==pins[name], 'snapshot checksum mismatch'
                    os.replace(temp,p);status='downloaded'
                assert hashlib.sha256(p.read_bytes()).hexdigest()==pins[name], 'cached checksum mismatch'
                record([acc,u,http,fmt,p.stat().st_size,hashlib.sha256(p.read_bytes()).hexdigest(),columns,p,status])
                print(acc,fmt,status,flush=True);break
            except Exception as e:
                record([acc,u,http,fmt,temp.stat().st_size if temp and temp.exists() else '', '', '',p,'FAILED '+repr(e)])
                print(acc,fmt,repr(e),flush=True)
            finally:
                if temp and temp.exists():temp.unlink()
