from common import *
from concurrent.futures import ThreadPoolExecutor
import psutil
m=rd('cache/metadata.tsv.gz');v=rd('tables/panel_members.tsv').set_index('gene');v1=V/'tables/matched_controls.tsv'
# Reuse outputs only when the specific script, immutable inputs and model configuration match.
base=stamp([Path(f) for f in ['config/main.json','config/input_fingerprints.tsv','code/prepare.py','code/fit.R','code/knk.py','tables/panel_members.tsv','tables/program_members.tsv','cache/panel.bin','cache/metadata.tsv.gz']])
def run(j,mode):
 guard();d=Path('cache/jobs')/j;d.mkdir(exist_ok=True,parents=True);marker=d/(mode+'.json');key=base+(sha('tables/control_matching.tsv') if mode!='reference' else '')
 if marker.exists():
  z=json.loads(marker.read_text())
  if z['key']==key and all(Path(f).exists() and sha(f)==h for f,h in z['outputs'].items()):return
 t=time.monotonic()
 with open('logs/'+j+'_'+mode+'.log','w') as log:
  cmd=['/usr/bin/time','-v','-o','logs/'+j+'_'+mode+'.usage','Rscript','--vanilla','code/fit.R',j,mode];subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
 if mode=='reference':out=[d/'WT.rds',d/'coverage.tsv.gz']
 else:
  out=[d/'status.tsv',d/'membership.tsv.gz']
  if rd(d/'status.tsv').status.iloc[0]=='complete':out += [d/'heat.tsv.gz',d/'KO_status.tsv',d/'coverage.tsv.gz']
  if j=='clean':out += [d/'manifold.tsv.gz',d/'patient.tsv.gz']
 z=dict(key=key,seconds=time.monotonic()-t,outputs={str(f):sha(f) for f in out});tmp=marker.with_suffix('.tmp');tmp.write_text(json.dumps(z));os.replace(tmp,marker);print(j,mode,'ok',flush=True)
def match():
 clean=rd('cache/jobs/clean/coverage.tsv.gz').set_index('gene');full=rd('cache/jobs/full/coverage.tsv.gz').set_index('gene')
 prog=set(rd('tables/program_members.tsv').gene)|set(E+K);rows=[]
 for k in K:
  typ='sequence_TF' if k in K[:2] else 'non_TF_protein'
  ix=[g for g in clean.index if g not in prog and v.at[g,'control_class']==typ and clean.at[g,'outdegree']>0 and full.at[g,'outdegree']>0]
  def scores(ctr):
   X=ctr.loc[ix+[k],['mean','det','outdegree']].copy();X[['mean','outdegree']]=np.log1p(X[['mean','outdegree']]);return (X-X.loc[k])/X.std()
  z=scores(clean);f=scores(full);dist=z.loc[ix].pow(2).sum(1).pow(.5);ids=dist.sort_values(kind='stable').head(C['controls_per_KO']).index
  for n,g in enumerate(ids,1):
   row=dict(ko=k,control=g,control_class=typ,match_order=n,candidates=len(ix),distance=dist[g],full_distance=np.linalg.norm(f.loc[g]))
   for nm,tab,st in [('clean',clean,z),('full',full,f)]:
    for col in z.columns:row.update({f'{nm}_{col}_SMD':st.at[g,col],f'{nm}_{col}_KO':tab.at[k,col],f'{nm}_{col}_control':tab.at[g,col]})
   rows.append(row)
 save(pd.DataFrame(rows),'tables/control_matching.tsv')
if __name__=='__main__':
 run('clean','reference');run('full','reference');match()
 jobs=['clean','full']+['seed'+str(i) for i in C['other_seeds']]+['drop_'+d for d in sorted(m.dataset.unique())]+['cluster_'+s for s in sorted(m.state.unique())]+['sub%02d'%i for i in range(1,21)]
 with ThreadPoolExecutor(max_workers=C['workers']) as pool:
  fs=[pool.submit(run,j,'KO') for j in jobs]
  for f in fs:f.result()

