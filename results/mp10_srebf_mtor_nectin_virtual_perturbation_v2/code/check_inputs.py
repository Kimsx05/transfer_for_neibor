from common import *
import argparse
parser=argparse.ArgumentParser();parser.add_argument('--seal',action='store_true');args=parser.parse_args()
record=Path('config/preparation_checkpoint.json')
outputs=['tables/gene_mapping.tsv','tables/coverage.tsv','tables/panel_members.tsv','tables/program_members.tsv','cache/panel.bin','cache/panel_counts.npz','cache/metadata.tsv.gz','qc/input.tsv','cache/base_edges.tsv.gz','config/input_fingerprints.tsv']
for row in rd('config/input_fingerprints.tsv').itertuples():assert sha(row.path)==row.sha256, 'source fingerprint changed: '+row.path
z={f:sha(f) for f in outputs+['code/prepare.py','code/oracle_prior.py','config/main.json']}
if args.seal:
 temp=record.with_suffix('.new');temp.write_text(json.dumps(z,indent=2));os.replace(temp,record)
else:assert json.loads(record.read_text())==z,'prepared features, prior or configuration changed'

