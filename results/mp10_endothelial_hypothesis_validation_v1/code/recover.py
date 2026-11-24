from fetch import *
import re
def ranged(acc,name,url,mode):
    folder=C/acc;folder.mkdir(parents=True,exist_ok=True);dest=folder/name
    with (folder/('.'+name+'.lock')).open('a') as lk:
        fcntl.flock(lk,fcntl.LOCK_EX)
        if dest.exists():return
        data=bytearray()
        try:
            while True:
                start=len(data);end=start+199999
                for attempt in range(3):
                    try:
                        r=requests.get(url+f'?download=1&range={start}&v={attempt}',headers={'Range':f'bytes={start}-{end}'},timeout=(6,12));r.raise_for_status();break
                    except requests.RequestException:
                        if attempt==2:raise
                r.raise_for_status()
                if r.status_code==200 and not start:data.extend(r.content);break
                a,b,n=map(int,re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)',r.headers['Content-Range']).groups())
                assert a==start and b-a+1==len(r.content)
                data.extend(r.content)
                if len(data)==n:break
            if mode=='gz':assert b'\t' in gzip.decompress(data)[:7000]
            elif mode=='pdf':assert data.startswith(b'%PDF-') and b'%%EOF' in data[-2048:]
            else:
                import io
                with zipfile.ZipFile(io.BytesIO(data)) as z:assert z.testzip() is None
            with tempfile.NamedTemporaryFile(dir=folder,delete=False) as f:f.write(data)
            os.replace(f.name,dest);record([acc,name,url,'downloaded',len(data),hashlib.sha256(data).hexdigest(),'range; validated'])
        except Exception as e:record([acc,name,url,'failed',len(data),'',repr(e)])
if __name__=='__main__':
    for a,n in [('GSE207787','NormalizedCounts_PDH_SREBF2.txt.gz'),('GSE207919','siKLF6_Normalizedcounts.txt.gz'),('GSE121770','htseq.all.counts.txt.gz')]:
        name=a+'_'+n;ranged(a,name,'https://ftp.ncbi.nlm.nih.gov/geo/series/'+a[:-3]+'nnn/'+a+'/suppl/'+name,'gz')
    ranged('eLife.63326','elife-63326-fig1-data1.xlsx','https://cdn.elifesciences.org/articles/63326/elife-63326-fig1-data1-v2.xlsx','xlsx')
    fetch('PMC8432967','mmc2.xlsx',['https://pmc.ncbi.nlm.nih.gov/articles/instance/8432967/bin/mmc2.xlsx','https://pmc-oa-opendata.s3.us-east-1.amazonaws.com/PMC8432967.1/mmc2.xlsx','https://www.ebi.ac.uk/biostudies/files/S-EPMC8432967/mmc2.xlsx'],'xlsx')
    
