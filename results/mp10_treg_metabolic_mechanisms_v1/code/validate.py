from pathlib import Path
import csv,gzip,hashlib,os,sys
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
import numpy as np
from scipy.stats import t as student
from PIL import Image
import fitz
out=Path(__file__).resolve().parents[1];tab=out/'tables';checks=[]
def ck(label,ok):
    if not ok:raise AssertionError(label)
    checks.append((label,'pass'))
def rows(p):
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return list(csv.DictReader(f,delimiter='\t'))
def bh(ps):
    p=np.array(ps); order=np.argsort(p);s=p[order]*len(p)/np.arange(1,len(p)+1)
    q=p.copy();q[order]=np.minimum(1,np.minimum.accumulate(s[::-1])[::-1]);return q
meta=rows(tab/'coverage.tsv')
ck('40 specimens 39 recorded IDs 13204 cells',len(meta)==40 and len({z['patient_uid_recorded'] for z in meta})==39 and sum(int(z['n_cells']) for z in meta)==13204)
joins=rows(tab/'sample_mapping.tsv')
taken=[r for r in joins if r['decision']=='included']
ck('35 unique GSM, explicit pools',len(taken)==35 and len({r['GSM'] for r in taken})==35)
for k,g in {'SCAP':(4,3),'FASN':(3,3),'TCR':(5,5),'CD36':(3,3),'SREBF2':(3,3)}.items():
    ck(k+' RNA group sizes',tuple(sum(r['contrast']==k and r['arm']==a for r in taken) for a in ['control','case'])==g)
proj=rows(tab/'specimen_scores.tsv.gz');models=rows(tab/'model_status.tsv')
klist=['SCAP','FASN','TCR','CD36','SREBF2']
for mask in ['Full','M20','Downstream']:
    ps=[];qs=[]
    for k in klist:
        d=sorted([r for r in proj if r['program']==k and r['mask']==mask and r['main_patient']=='TRUE' and int(r['n_cells'])>=50],key=lambda r:r['sample_uid'])
        ck(mask+' '+k+' membership',len(d)==len({r['patient_uid_recorded'] for r in d})==32)
        ds=sorted({r['dataset'] for r in d})
        x=np.array([[1,*[int(r['dataset']==v) for v in ds[1:]],float(r['Q50_z'])] for r in d])
        y=np.array([float(r['score']) for r in d])
        cf=np.linalg.lstsq(x,y,rcond=None)[0];df=len(d)-x.shape[1]
        se=np.sqrt(np.linalg.inv(x.T@x)[-1,-1]*sum((y-x@cf)**2)/df)
        quant=student.ppf(.975,df);p=2*student.sf(abs(cf[-1]/se),df)
        obs=next(r for r in models if r['setting']=='M0046' and r['adjustment']=='total' and r['mask']==mask and r['program']==k)
        ck(mask+' '+k+' independent beta/CI/p',np.allclose([cf[-1],cf[-1]-quant*se,cf[-1]+quant*se,p],[float(obs[t]) for t in ['beta','lo','hi','p']],rtol=1e-8,atol=1e-11))
        ps.append(p);qs.append(float(obs['q']))
    ck(mask+' independent BH',np.allclose([*bh(ps[:3]),*bh(ps[3:])],qs))
routes=rows(tab/'recovered_sender_receiver.tsv')
ck('18 main routes HC3 NE; 17 sensitivity fitted',all((r['reason']=='HC3_undefined_unit_leverage' and r['n_specimens']=='18') if r['setting']=='main' else (r['reason']=='estimated' and r['n_specimens']=='17') for r in routes))
pairs=[r for r in rows(tab/'gene_mapping.tsv.gz') if r['status']=='one_to_one']
ck('MGI both directions unique',len(pairs)==len({r['mouse'] for r in pairs})==len({r['human'] for r in pairs}))
pins={r['file']:r['sha256'] for r in rows(out/'code/resource_hashes.tsv')}
visited=set()
for r in rows(out/'resource_manifest.tsv'):
    if r['status'] not in ('cached','downloaded'):continue
    p=Path(r['path'])
    if p in visited:continue
    ck(p.name+' SHA256',hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256']==pins[p.name]);visited.add(p)
ck('all 9 pinned resources verified',len(visited)==len(pins)==9)
for f in rows(out/'FIGURE_INDEX.tsv'):
    p=out/f['PDF'];png=out/f['PNG'];t=out/f['data']
    with fitz.open(p) as doc:ck(f['figure']+' PDF',len(doc)==1 and doc[0].rect.width>0 and bool(doc[0].get_text().strip()))
    with Image.open(png) as im:ck(f['figure']+' PNG 300dpi',im.width>900 and min(im.info.get('dpi',(0,0)))>299)
    ck(f['figure']+' plot rows',len(rows(t))>0)
for name in ['METHODS.md','RESULTS.md','registration.md']:
    ck(name+' populated',(out/name).stat().st_size>1000)
with (out/'logs/validation.tsv').open('w') as f:
    wr=csv.writer(f,delimiter='\t');wr.writerow(['check','status']);wr.writerows(checks)
print('PASS',len(checks),'checks',sys.version.split()[0],np.__version__)
