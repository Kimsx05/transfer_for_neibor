from fetch import *
extras=[('GSE149933','GSE149933_family.soft.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149933/soft/GSE149933_family.soft.gz'),('Lawson2020','mValidationAnalysis.Rmd','https://raw.githubusercontent.com/NMikolajewicz/Lawson2020/82662ff8183307ec09439dc001834537ec00bda3/Catherine/mValidationAnalysis.Rmd')]
for d,n,u in extras:
 dst=C/d/n
 if dst.exists():continue
 with tempfile.NamedTemporaryFile(dir=dst.parent,delete=False) as f:p=Path(f.name)
 try:
  sp.run(['/usr/bin/curl','-fsSL','--retry','3','--retry-all-errors','--max-time','180',u,'-o',str(p)],check=True);validate(p,n,d);os.replace(p,dst)
 finally:
  if p.exists():p.unlink()

