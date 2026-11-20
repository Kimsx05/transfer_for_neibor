from base6 import *
import importlib.metadata as im,hashlib,platform,subprocess
m=rd(ROOT/'input_manifest.tsv')
extra=[]
def add(s,role,p):
 p=Path(p);extra.append(dict(section=s,role=role,path=str(p.resolve()),readable=p.is_file(),bytes=p.stat().st_size if p.is_file() else None))
for s in sections:
 for f in ['cells_scored.parquet','candidate_ledger.parquet','raw_candidates.parquet','genes.tsv']:add(s,'v5_'+f,V5/s/f)
 for f in ['tables/reference_membership.parquet','models/reference_prepared.json','code/common.py']:add(s,'historical_reference_or_scoring',V/f)
 add(s,'frozen_full_reference',js(V/'models/reference_prepared.json')['reference']);add(s,'patient_mapping',A.parents[1]/'mp10_k32_reference_identifiability_v5/00_audit/provenance/080_sample_patient_mapping.tsv');add(s,'calibrated_H_and_E_model',A/'cell_annotation_freeze_v1/models/keras_cache/models/StarDist2D/2D_versatile_he/2D_versatile_he_extracted/weights_best.h5')
put(pd.concat([m[~m.role.str.startswith(('v5_','historical_reference_or_scoring','frozen_full_reference','patient_mapping','calibrated_H_and_E_model'))],pd.DataFrame(extra)],ignore_index=True),ROOT/'input_manifest.tsv')
vers={'python':platform.python_version(),'proseg':subprocess.check_output([str(ROOT/'.tools/prebuilt/bin/proseg'),'--version'],text=True).strip()}
for pkg in ['scikit-learn','pyarrow','numpy','scipy','h5py','anndata','pillow','shapely','numba','matplotlib']:
 vers[pkg]=im.version(pkg)
put(pd.DataFrame(vers.items(),columns=['package','version']),ROOT/'tables/software.tsv')
