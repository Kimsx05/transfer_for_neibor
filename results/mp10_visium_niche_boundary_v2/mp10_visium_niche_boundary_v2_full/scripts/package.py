from pathlib import Path
import hashlib,csv,json,zipfile,os
O=Path(__file__).resolve().parents[1];archive=O/'mp10_visium_niche_boundary_v2_full.zip'
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
skip={archive.name,'file_manifest.tsv','archive_validation.json','portable_test_path.txt'}
files=sorted(p for p in O.rglob('*') if p.is_file() and p.name not in skip and '__pycache__' not in p.parts and p.suffix!='.tmp')
assert json.loads((O/'audit/delivery_validation.json').read_text())['validated']
mf=O/'file_manifest.tsv'
with mf.open('w') as f:
    w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['path','bytes','sha256'])
    for p in files:w.writerow([p.relative_to(O),p.stat().st_size,sha(p)])
temp=archive.with_suffix('.zip.tmp')
with zipfile.ZipFile(temp,'w',zipfile.ZIP_DEFLATED,compresslevel=3,allowZip64=True) as z:
    for p in files+[mf]:z.write(p,p.relative_to(O))
with zipfile.ZipFile(temp) as z:
    assert z.testzip() is None
    names=z.namelist()
    for d,suff in [('spots/','.parquet'),('boundary/','.parquet'),('local/','.parquet'),('review/','.jpg')]:
        matches=[i for i in z.infolist() if i.filename.startswith(d) and i.filename.endswith(suff)]
        assert len(matches)==63 and all(i.file_size>500 for i in matches)
    assert len(names)==len(set(names))==len(files)+1
os.replace(temp,archive)
report=dict(archive=archive.name,sha256=sha(archive),archive_bytes=archive.stat().st_size,uncompressed_members_bytes=sum(p.stat().st_size for p in files+[mf]),members=len(files)+1,crc='all members passed',parquet_sections_each=63,review_JPEGs=63,manifest_excludes_itself_and_archive=True)
r=O/'audit/archive_validation.json';r.write_text(json.dumps(report,indent=2));print(json.dumps(report))

