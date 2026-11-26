"""Validate the returned form, archive it, and recompute the affected samples in a separate revision."""
from common import *
import hashlib,subprocess,shutil,datetime,shlex
def parse(path,s):
 t=pd.read_csv(path,sep='\t',dtype=str).fillna('');base=pd.read_csv(paths(s)/'human_review.tsv',sep='\t').id
 if not t.id.is_unique or not set(t.id)<=set(base):raise ValueError('duplicate or unrecognized review ID')
 todo=t[t.human_status.isin(['reviewed','human-reviewed'])].copy()
 if not todo.human_reviewer.str.strip().ne('').all():raise ValueError('human reviewer is required for a signed opinion')
 if not todo.human_grade.isin(['supported','possible','unresolved','unsupported']).all():raise ValueError('invalid review grade')
 for k in ['human_reliable_T','human_epi_eligible']:
  if not todo[k].str.lower().isin(['','true','false','1','0']).all():raise ValueError(k)
 return todo
def apply(changes):
 digest=hashlib.sha256()
 for s,t in sorted(changes.items()):digest.update(s.encode());digest.update(t.to_csv(index=False).encode())
 dst=OUT/'revisions'/('review_'+digest.hexdigest()[:12]);dst.mkdir(parents=True,exist_ok=True)
 if (dst/'complete.json').is_file():print('Already completed',dst);return
 base=OUT
 for path in ['audit','logs','reports','niche','figures','hd','xenium','review','config','spatial']:(dst/path).mkdir(exist_ok=True)
 if not (dst/'code').exists():(dst/'code').symlink_to(base/'code',target_is_directory=True)
 cfg=dict(CFG,version=dst.name,human_review_status='partially human-reviewed');(dst/'config/frozen.json').write_text(json.dumps(cfg,indent=2))
 for s in HD+XE:
  v,old=folders(s);sampledir=dst/old.relative_to(base);q=dst/'review'/s
  if s not in changes:
   for path in [old,base/'review'/s,base/'niche'/s,base/'figures'/s]:
    new=dst/path.relative_to(base)
    if path.exists() and not new.exists():new.symlink_to(path,target_is_directory=True)
   continue
  sampledir.mkdir(exist_ok=True);q.mkdir(exist_ok=True);t=changes[s].set_index('id');put(t.reset_index(),q/'signed_changes.tsv')
  i=pd.read_parquet(old/'identity.parquet');prior=i.loc[t.index,['grade','reliable_T','major']].copy()
  i.loc[t.index,'grade']=t.human_grade
  asbool={'true':True,'1':True,'false':False,'0':False}
  for field,col in [('human_reliable_T','reliable_T'),('human_epi_eligible','human_epi_eligible')]:
   val=t[field].str.lower().map(asbool).dropna()
   if len(val):i.loc[val.index,col]=val.astype(bool)
  if (i.grade.isin(['supported','possible'])&~i.reliable_T).any():raise ValueError(s+' Treg component must belong to reliable T; resolve signed parent conflict')
  if s in XE and 'human_epi_eligible' in i and (i.reliable_T&i.human_epi_eligible.eq(True)).any():raise ValueError('T and Epi identities conflict')
  i.loc[t.index,'human_review_status']='human-reviewed';i.loc[t.index,'reviewer']=t.human_reviewer;i.loc[t.index,'review_origin']='signed human form';i.loc[t.index,'review_status']='human-reviewed';i.to_parquet(sampledir/'identity.parquet')
  prior.join(t).to_csv(q/'before_after.tsv',sep='\t');put(pd.DataFrame({'id':i.index,'eligible_override':i.get('human_epi_eligible',pd.Series(index=i.index,dtype='boolean'))}),q/'eligibility.tsv')
  candidate_ids=i.index[i.candidate_union | i.index.isin(t.index)]
  upgraded=i.loc[candidate_ids]
  if s in HD:upgraded=upgraded.join(pd.read_parquet(v/'weights.parquet').add_prefix('weight_'))
  upgraded.to_parquet(q/'candidates.parquet')
  for n in ['atlas.html','window_index.tsv','human_review.tsv','index.html','cards','tiles']:
   if (paths(s)/n).exists() and not (q/n).exists():(q/n).symlink_to(paths(s)/n)
 env=dict(os.environ,MP_V2_OUTPUT=str(dst),MP_V2_BASE=str(base),OUT=str(V))
 commands=[]
 def job(*args):commands.append(shlex.join([str(base/'code/run.sh'),*args]))
 for s in changes:
  if s in XE:job('bash','code/actual_R.sh','code/update_scores.R',s)
  for script in ['score','graph','spatial','reverse']:job('code/'+script+'.py',s)
  job('code/niche.py','mender',s);job('bash','code/actual_R.sh','code/build.R',s)
  # Expression representation and graph depend on the raw units only.
  pp=dst/'niche'/s
  for f in (base/'niche'/s).glob('CellCharter*'):
   if not (pp/f.name).exists():(pp/f.name).symlink_to(f)
  job('code/squid.py',s);job('code/export_tools.py',s);job('bash','code/actual_R.sh','code/spatstat.R',s);job('bash','code/actual_R.sh','code/crawdad.R',s);job('bash','code/actual_R.sh','code/formal_crawdad.R',s);job('code/target_tools.py',s);job('code/kprep.py',s);job('bash','code/actual_R.sh','code/kcross.R',s);job('code/niche_results.py',s);job('code/figures.py',s);job('bash','code/actual_R.sh','code/violin.R',s)
 job('code/check.py');job('code/extras.py');job('code/report.py');job('code/compare.py');job('code/report.py')
 finish=dst/'complete.json';result=json.dumps({'version':dst.name,'samples':list(changes)})
 entry=dst/'recompute.sh'
 header=['#!/bin/bash','set -euo pipefail','cd '+shlex.quote(str(dst))]
 for k in ['MP_V2_OUTPUT','MP_V2_BASE','OUT']:header.append('export '+k+'='+shlex.quote(env[k]))
 entry.write_text('\n'.join(header+commands+["printf '%s\\n' "+shlex.quote(result)+' > '+shlex.quote(str(finish))])+'\n')
 if 'MP_V2_PLAN_FILE' in os.environ:Path(os.environ['MP_V2_PLAN_FILE']).write_text(str(entry)+'\n')
 else:os.execve('/bin/bash',['bash',str(entry)],env)
if __name__=='__main__':
 args=sys.argv[1:];changes={}
 for p in args:
  path=Path(p).resolve();table=pd.read_csv(path,sep='\t',nrows=1);s=next((s for s in HD+XE if set(table.id)<=set(pd.read_csv(paths(s)/'human_review.tsv',sep='\t').id) and (path.parent.name==s or s in path.name)),None)
  if s is None:raise ValueError('Name the form after its sample or keep it in the sample folder')
  t=parse(path,s)
  if len(t):changes[s]=t
 if changes:apply(changes)
 else:print('No signed human decisions. No revision created.')
