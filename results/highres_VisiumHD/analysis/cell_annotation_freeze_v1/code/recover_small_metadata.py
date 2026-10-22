# FROZEN_VERSION_GUARD: reproduce in a new version directory; never overwrite v1.
from pathlib import Path as _FreezePath
if (_FreezePath(__file__).resolve().parents[1] / 'metadata/FROZEN.json').exists():
    raise RuntimeError('Frozen v1 is immutable; create v2 for changes.')
from pathlib import Path
import importlib.util,json,zipfile,hashlib
out=Path(__file__).resolve().parents[1];hd=out.parents[1]
spec=importlib.util.spec_from_file_location('validated_range',hd/'code/download_inputs.py');d=importlib.util.module_from_spec(spec);spec.loader.exec_module(d);d.ROOT=out
meta=json.loads((out/'metadata/public_cervilla_zenodo.txt').read_text());entry=next(x for x in meta['files'] if x['key']=='visiumhd_Bladder.zip')
z=zipfile.ZipFile(d.RemoteZip(entry['links']['self'],entry['size']))
selected=['visiumhd_Bladder/probe_set.csv','visiumhd_Bladder/spatial/aligned_fiducials.jpg','visiumhd_Bladder/spatial/detected_tissue_image.jpg']
records=[]
print('Expected selected uncompressed bytes',sum(z.getinfo(n).file_size for n in selected),flush=True)
for n in selected:
 info=z.getinfo(n);dest=out/'inputs/Cervilla_recovered_metadata'/Path(n).name;dest.parent.mkdir(parents=True,exist_ok=True)
 if not dest.exists():
  data=z.read(n) # zipfile validates member CRC
  assert len(data)==info.file_size
  dest.write_bytes(data)
 records.append(dict(member=n,path=str(dest),bytes=dest.stat().st_size,sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),zip_crc32=f'{info.CRC:08x}',member_crc_verified=True,full_archive_md5_verified=False))
 print(n,'OK',flush=True)
(out/'metadata/recovered_small_inputs.json').write_text(json.dumps(records,indent=2))
