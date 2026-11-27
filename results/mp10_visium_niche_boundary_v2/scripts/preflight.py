from lib import *
for fn in ['v1_per_section_hashes.tsv','v1_input_hashes.tsv','software_inputs.tsv']:
    manifest=read(O/'audit'/fn)
    for r in manifest.itertuples():
        p=Path(r.path)
        if not p.is_file() or sha(p)!=r.sha256:raise ValueError('Input changed or missing: '+str(p))
print('Frozen checksums verified for v1 and scoring inputs.')

