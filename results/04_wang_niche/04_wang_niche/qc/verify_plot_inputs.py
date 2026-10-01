import sys
from pathlib import Path
ROOT=Path('/home/data/t070721/codex_workspace/Bladder Metabolism/spatial_ecology_top20_v2')
sys.path.insert(0,str(ROOT/'code'))
from common import *
OUT=ROOT/'04_wang_niche'
hist=readcsv(ROOT/'00_inputs/histology_manifest.tsv',sep='\t')
provenance=readcsv(ROOT/'00_inputs/input_provenance.tsv',sep='\t')
rows=[]
for _,r in hist.iterrows():
    expected=provenance.loc[provenance.path.eq(r.image_path),'sha256'].dropna().unique()
    assert len(expected)==1,f'No unique frozen image SHA: {r.image_path}'
    current=sha(r.image_path);assert current==expected[0],f'Frozen H&E changed: {r.image_path}'
    rows.append(dict(dataset=r.dataset,section_id=r.section_id,image_path=r.image_path,
        current_sha256=current,upstream_frozen_sha256=expected[0],status='PASS'))
assert len(rows)==63
tsv(OUT/'plot_input_validation.tsv',rows)
savejson(OUT/'plot_input_signature.json',dict(n_images=len(rows),status='PASS',
    plot_input_validation_sha256=sha(OUT/'plot_input_validation.tsv'),
    validator_path=str(Path(__file__).resolve()),validator_sha256=sha(__file__),
    source='Module 00 frozen input_provenance.tsv and histology_manifest.tsv',
    comparison='all actual H&E image bytes remain identical to frozen source SHA; no new image preprocessing'))
print('PLOT_INPUTS_VERIFIED: 63 unchanged accepted H&E images')
