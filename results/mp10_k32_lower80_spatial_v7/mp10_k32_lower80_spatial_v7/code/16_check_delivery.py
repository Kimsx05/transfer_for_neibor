"""Read-only delivery checks; deliberately outside package self-hash cycle."""
from common import *
import re,zipfile
from PIL import Image

review=R/'07_review';checks=[]
def ck(name,ok,detail=''):
 checks.append({'check':name,'passed':bool(ok),'detail':detail})
for base in [R,review]:
 text=(base/'REPORT.html').read_text()
 for src in re.findall(r'<img[^>]+src="([^"]+)"',text):ck(str(base.name)+' image '+src,(base/src).is_file())
for p in (review/'06_figures').rglob('*.png'):
 with Image.open(p) as im:
  ck('PNG alpha '+str(p.relative_to(review)),'A' in im.getbands())
  dpi=im.info.get('dpi',(0,0));ck('PNG 300dpi '+str(p.relative_to(review)),all(abs(x-300)<1 for x in dpi))
 stem=p.with_suffix('');ck('figure companion files '+str(stem.relative_to(review)),all(Path(str(stem)+s).exists() for s in ['.pdf','_caption.txt','_plot_data.tsv.gz']))
manifest=pd.read_csv(review/'PACKAGE_MANIFEST.tsv',sep='\t')
for row in manifest.itertuples():
 p=review/row.path;ck('package SHA '+row.path,p.is_file() and p.stat().st_size==row.bytes and sha(p)==row.sha256)
with zipfile.ZipFile(review/'review.zip') as z:
 ck('ZIP CRC',z.testzip() is None)
 members=z.namelist();ck('ZIP unique paths',len(members)==len(set(members)))
 ck('ZIP no large raw arrays or models',not any(Path(x).suffix in ['.h5ad','.pt','.ckpt','.npy','.npz'] for x in members))
 ck('ZIP all manifest members',set(manifest.path).issubset(members))
 for row in manifest.itertuples():ck('ZIP byte equality '+row.path,hashlib.sha256(z.read(row.path)).hexdigest()==row.sha256)
for sec in config()['sections']:
 o=outmap(sec);st=json.loads((o/'RESULT_STATUS.json').read_text());f=review/'06_figures/by_dataset'/sec.split('_')[0]/sec
 ck(sec+' required panels',all((f/(n+'.png')).exists() for n in ['coordinate_basemap','mapping_overview','raw_expression_overview','K32_mean','lower80_mean','Epi_total_mean','I_state_mean','K32_q05','vCAF_mean','state_expression_associations']))
 cap=(f/'mapping_overview_caption.txt').read_text();ck(sec+' final state in caption',st['numerical_status'] in cap and str(st['epochs']) in cap)
wt(pd.DataFrame(checks),R/'DELIVERY_CHECKS.tsv')
result={'n_checks':len(checks),'n_failed':sum(not x['passed'] for x in checks),'zip_bytes':(review/'review.zip').stat().st_size,'zip_sha256':sha(review/'review.zip'),'package_n_files':len(manifest),'checked_epoch':time.time(),'not_inside_zip':'This final external validation avoids a circular archive hash.'};js(R/'DELIVERY_CHECKS.json',result);say('DELIVERY_CHECKS',result);assert all(x['passed'] for x in checks)
