from prepare import *
from annotate import annotate
import subprocess, argparse, traceback

def execute(folder):
    folder=Path(folder);sig=json.loads((folder/'input_signature.json').read_text())
    cache=folder/'completion.json'
    if cache.exists():
        prior=json.loads(cache.read_text())
        if prior.get('input_signature')==sig['signature'] and all(Path(p).is_file() and sha(p)==h for p,h in prior['output_checksums'].items()):
            print('REUSE_VERIFIED',sig['dataset'],flush=True);return prior
    try:
        with open(folder/'discovery.log','w',buffering=1) as log:
            subprocess.run(['Rscript',str(OUT/'code/discover.R'),str(folder)],stdout=log,stderr=subprocess.STDOUT,check=True)
        status=json.loads((folder/'discovery_status.json').read_text())
        if status['status']=='SKIPPED_NOT_ESTIMABLE': return status
        baseline=OUT/'execution_validation_baseline'/sig['dataset']
        if baseline.is_dir():
            before=readcsv(baseline/'all_candidate_memberships.csv.gz')
            after=readcsv(folder/'all_candidate_memberships.csv.gz')
            assert before.equals(after),'Parallel candidate memberships differ from verified serial baseline'
            old_grid=readcsv(baseline/'silhouette_grid.tsv',sep='\t')
            new_grid=readcsv(folder/'silhouette_grid.tsv',sep='\t')
            err=float(np.max(np.abs(old_grid.mean_silhouette-new_grid.mean_silhouette)))
            assert err==0 and old_grid.n_clusters.equals(new_grid.n_clusters)
            savejson(folder/'execution_parallelism_parity.json',dict(dataset=sig['dataset'],status='PASS',
                n_candidates=len(new_grid),n_spot_memberships_compared=len(after)*(len(after.columns)-1),
                n_membership_differences=0,max_silhouette_difference=err,
                comparison='serial 1-worker versus 4 candidate workers; same graph, seed, features and grid'))
        annotate(folder)
        with open(folder/'plots.log','w',buffering=1) as log:
            subprocess.run(['Rscript',str(OUT/'code/plots.R'),str(folder)],stdout=log,stderr=subprocess.STDOUT,check=True)
        files=[p for p in folder.rglob('*') if p.is_file() and p.name!='completion.json']
        result=dict(dataset=sig['dataset'],status='PASS_WITH_WARNING',input_signature=sig['signature'],
            n_valid=status['n_valid'],n_niches=status['n_niches'],
            output_checksums={str(p):sha(p) for p in files})
        savejson(cache,result);print('DATASET_COMPLETE',sig['dataset'],flush=True);return result
    except Exception as exc:
        savejson(folder/'failure.json',dict(dataset=sig['dataset'],status='BLOCKED',error=str(exc),traceback=traceback.format_exc()))
        raise

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--phase',choices=['prepare','smoke','all'],required=True)
    args=parser.parse_args()
    if args.phase=='prepare': prepare()
    else:
        order=json.loads((OUT/'dataset_execution_order.json').read_text());folders=order['dataset_dirs']
        if args.phase=='smoke':
            result=execute(folders[0]);assert result['status']=='PASS_WITH_WARNING'
            savejson(OUT/'smoke_mechanical_check.json',dict(dataset=order['smoke_dataset'],status='PASS',
                complete_45_candidate_grid=True,checks=['factor-only feature columns','input SHA and real spot keys',
                'positive undirected graph','fixed shared Euclidean silhouette set','assignments retain full universe',
                'asymptotic Wilcoxon matches scipy','PDF and 600 dpi PNG rendered'],completion=result))
        else:
            smoke=json.loads((OUT/'smoke_mechanical_check.json').read_text());assert smoke['status']=='PASS'
            for folder in folders: execute(folder)
