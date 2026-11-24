"""Finite official downloads, structural validation, atomic publish and SHA256."""
import csv,fcntl,gzip,hashlib,os,pathlib,sys,tarfile,tempfile,urllib.parse,zipfile
import requests
from xml.etree import ElementTree as ET
P=pathlib.Path(__file__).resolve().parents[1]
C=P.parents[1]/'cache/mp10_public_resources_v1'
log=P/'logs/resource_attempts.tsv'
def record(row):
    with log.open('a') as f:
        w=csv.writer(f,delimiter='\t')
        if f.tell()==0:w.writerow(['resource','file','url','status','bytes','sha256','detail'])
        w.writerow(row)
def fetch(acc,name,urls,kind):
    folder=C/acc; folder.mkdir(parents=True,exist_ok=True); dest=folder/name
    def valid(path):
        if kind in ('matrix','softgz'):
            with gzip.open(path,'rb') as f:
                head=f.readline(); assert len(head)>5
                if kind=='matrix':assert b'\t' in head
                while f.read(1024*512):pass
        elif kind=='tar':
            with tarfile.open(path) as f:assert len(f.getnames())>=20
        elif kind=='soft':assert b'^SAMPLE' in path.read_bytes()
        elif kind=="pdf":assert path.read_bytes().startswith(b"%PDF-") and b"%%EOF" in path.read_bytes()[-2048:]
        elif kind=='xlsx':
            with zipfile.ZipFile(path) as z:assert 'xl/workbook.xml' in z.namelist() and z.testzip() is None
        elif kind=='xml':assert ET.parse(path).getroot().tag in ('article','pmc-articleset')
        elif kind=='zip':
            with zipfile.ZipFile(path) as z:assert z.testzip() is None
    def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
    with (folder/('.'+name+'.lock')).open('a') as lk:
        fcntl.flock(lk,fcntl.LOCK_EX)
        if dest.exists():
            try:valid(dest); record([acc,name,'local','reused',dest.stat().st_size,digest(dest),'validated']);return dest
            except Exception as e:record([acc,name,'local','invalid_existing',0,'',repr(e)]);return
        for url in urls:
            tmp=None
            try:
                with requests.get(url,stream=True,timeout=(7,20),headers={'User-Agent':'Mozilla/5.0 (research GEO data access)'}) as r:
                    r.raise_for_status()
                    with tempfile.NamedTemporaryFile(dir=folder,delete=False) as f:
                        tmp=pathlib.Path(f.name)
                        for chunk in r.iter_content(524288):f.write(chunk)
                valid(tmp); sha=digest(tmp); os.replace(tmp,dest)
                record([acc,name,url,'downloaded',dest.stat().st_size,sha,'validated']);return dest
            except Exception as e:record([acc,name,url,'failed',0,'',repr(e)])
            finally:
                if tmp and tmp.exists():tmp.unlink()
if __name__=='__main__':
    files={'GSE174785':'GSE174785_RAW.tar','GSE201466':'GSE201466_Fig1RNAseq.txt.gz','GSE207787':'GSE207787_NormalizedCounts_PDH_SREBF2.txt.gz','GSE207919':'GSE207919_siKLF6_Normalizedcounts.txt.gz','GSE121770':'GSE121770_htseq.all.counts.txt.gz'}
    for acc,name in files.items():
        url='https://ftp.ncbi.nlm.nih.gov/geo/series/'+acc[:-3]+'nnn/'+acc
        fetch(acc,name,[url+'/suppl/'+name,'https://www.ncbi.nlm.nih.gov/geo/download/?acc='+acc+'&format=file&file='+urllib.parse.quote(name,safe='').replace('_','%5F'),'https://www.ncbi.nlm.nih.gov/geo/download/?acc='+acc+'&format=file&file='+name], 'tar' if name.endswith('tar') else 'matrix')
        fetch(acc,acc+'_family.soft.gz',[url+'/soft/'+acc+'_family.soft.gz'],'softgz')
        if not (C/acc/(acc+'_family.soft.gz')).exists():
            fetch(acc,acc+'.soft',['https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+acc+'&targ=all&form=text&view=full'],'soft')
    for pm in ['PMC7997658','PMC8432967','PMC10407908','PMC10515159','PMC6521885']:
        fetch(pm,'article.xml',['https://www.ebi.ac.uk/europepmc/webservices/rest/'+pm+'/fullTextXML'],'xml')
    fetch('eLife.63326','elife-63326-fig1-data1.xlsx',['https://cdn.elifesciences.org/articles/63326/elife-63326-fig1-data1-v2.xlsx','https://cdn.elifesciences.org/articles/63326/elife-63326-fig1-data1-v1.xlsx'],'xlsx')

