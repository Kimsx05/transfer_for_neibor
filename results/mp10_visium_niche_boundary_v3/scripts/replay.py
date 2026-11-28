# Refit a frozen label configuration in a temporary output directory, from caches only.
from common3 import *
import tempfile,subprocess
stages=['cluster','lock','down','auxiliary','geometry_sensitivity','paired_curves','extra_tests','group_plots']
with tempfile.TemporaryDirectory(prefix='light_replay_',dir=O/'logs') as tmp:
    dest=Path(tmp);env=dict(os.environ,V3_MODEL_OUT=tmp,V3_SEED=str(seed))
    for stage in stages:
        with open(dest/(stage+'.log'),'w') as f:subprocess.run(['bash',str(O/'run.sh'),stage],cwd=O,env=env,check=True,stdout=f,stderr=subprocess.STDOUT)
    for name in ['classification/section_classification.tsv','tables/four_primary_tests.tsv','tables/overall_level_global.tsv','tables/geometry_curve_tests.tsv','tables/paired_curve_summary.tsv.gz']:
        pd.testing.assert_frame_equal(read(O/name),read(dest/name))
    n=len(list((dest/'figures/groups').glob('*.pdf')))
    assert n==len(read(O/'audit/group_figure_index.tsv'))
    (O/'audit/lightweight_replay.json').write_text(json.dumps(dict(ok=True,seed=seed,stages=stages,reproduced_groups=len(read(dest/'classification/section_classification.tsv')),identical_four_primary_tests=True,standalone_group_figures=n,input='immune caches, point-bin/geometry summaries, existing v2 lightweight section tables. No counts, Seurat or H5AD.'),indent=2))

