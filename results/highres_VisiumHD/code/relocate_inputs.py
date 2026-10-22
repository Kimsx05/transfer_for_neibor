"""Move recovered source assets on the same filesystem, retaining old-path symlinks."""
from pathlib import Path
import os,json
ROOT=Path(__file__).resolve().parents[1]; OLD=ROOT.parent/'highres_ST'
jobs=[(OLD/'highres_spatial_audit_v3_GSE300692/00_input_audit/extracted_source/D1_Johnson',ROOT/'data/GSE300692/D1_Johnson'),(OLD/'VisiumHD_MP10_GSE299573_independent_validation_v1/02_resolution_construction_audit/extracted_8um_spatial/spatial',ROOT/'data/GSE299573/GSM9040999/square_008um/spatial')]
for p in (OLD/'spatial_input_audit_3GEO_v1/audit_GSE318916/_inspection_cache').glob('*.h5'):jobs.append((p,ROOT/'data/GSE318916'/p.name))
for src,dst in jobs:
 if src.is_symlink() and src.resolve()==dst:continue
 if not src.exists():raise FileNotFoundError(src)
 if dst.exists():raise FileExistsError(dst)
 dst.parent.mkdir(parents=True,exist_ok=True)
 before=src.stat();os.rename(src,dst);src.symlink_to(dst,target_is_directory=dst.is_dir())
 assert dst.stat().st_ino==before.st_ino and src.exists()
 with (ROOT/'logs/relocations.jsonl').open('a') as f:f.write(json.dumps({'old':str(src),'new':str(dst),'inode':before.st_ino,'method':'same-filesystem rename + compatibility symlink'})+'\n')
 print(src,'->',dst)
