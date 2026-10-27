from pathlib import Path
import csv,math,json
BASE=Path('mp10_spatial_communication_v1/spatialcellchat')
def read(p):
 with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(rows,p):
 if not rows:return
 with open(p,'w') as f:
  w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)),delimiter='\t');w.writeheader();w.writerows(rows)
checks=[];coverage=[];runs=[]
for p in sorted(BASE.glob('*/*/COMPLETED.txt')):
 p=p.parent; section=p.parent.name;run=p.name
 stat=read(p/'position_LR_support.tsv');lr={x['interaction']:x for x in read(p/'LR_selected.tsv')};rawflow=sum(float(x['cross_spot_flow']) for x in stat);rawdiag=sum(float(x['same_spot_flow']) for x in stat)
 for level in ['major','fine']:
  e=read(p/f'{level}_edges.tsv');flow=sum(float(x['total_flow']) for x in e);diag=sum(float(x['same_spot_flow_excluded']) for x in e)
  for name,a,b in [('flow_conservation_Q_partition',flow,rawflow),('excluded_diagonal_conservation',diag,rawdiag)]:
   err=abs(a-b);ok=err<=1e-8*max(1,abs(b));assert ok,(section,run,level,name,err)
   checks.append(dict(section=section,run=run,level=level,check=name,value=a,reference=b,abs_error=err,status='PASS'))
 for category in ['Secreted Signaling','ECM-Receptor','Cell-Cell Contact']:
  ss=[x for x in stat if lr[x['interaction']]['annotation']==category]
  coverage.append(dict(section=section,run=run,category=category,n_selected_LR=len(ss),n_with_QC_cross_spot_support=sum(float(x['cross_spot_flow'])>0 for x in ss),cross_spot_links=sum(int(x['cross_spot_links']) for x in ss),cross_spot_flow=sum(float(x['cross_spot_flow']) for x in ss),excluded_same_spot_flow=sum(float(x['same_spot_flow']) for x in ss)))
 for x in read(p/'run_status.tsv'):runs.append(dict(run=run,**x))
write(checks,BASE/'numerical_checks.tsv');write(coverage,BASE/'LR_category_summary.tsv');write(runs,BASE/'run_summary.tsv')
print('complete_runs',len(runs),'checks',len(checks),'all_PASS')
# Same native LR set, same Q and fixed contact radius; first-slice wider diffusible range only.
b=BASE/'GSE171351_A1'
if all((b/r/'COMPLETED.txt').exists() for r in ['main','range3dNN']):
 sens=[]
 for level in ['major','fine']:
  x={(z['source'],z['target']):float(z['total_mean_strength']) for z in read(b/'main'/f'{level}_edges.tsv')}
  y={(z['source'],z['target']):float(z['total_mean_strength']) for z in read(b/'range3dNN'/f'{level}_edges.tsv')}
  nodes=sorted({s for s,t in x if not s.startswith('Epi_MP10_')})
  for node in nodes:
   for dr in ['epi_to_TME','TME_to_epi']:
    key=lambda state:(f'Epi_MP10_{state}',node) if dr=='epi_to_TME' else (node,f'Epi_MP10_{state}')
    lo,hi=x[key('low')],x[key('high')];wl,wh=y[key('low')],y[key('high')]
    sign=lambda v:('High>Low' if v>1e-12 else 'High<Low' if v< -1e-12 else 'zero_difference')
    sens.append(dict(section='GSE171351_A1',level=level,node=node,direction=dr,main_low=lo,main_high=hi,main_delta=hi-lo,wide_low=wl,wide_high=wh,wide_delta=wh-wl,main_direction=sign(hi-lo),wide_direction=sign(wh-wl),direction_agrees=sign(hi-lo)==sign(wh-wl),note='Secreted and ECM 200→300um; Contact105um fixed; same1037LR input'))
 write(sens,BASE/'sensitivity_epi_state_contrasts.tsv')
# GNU time resource records; maximum RSS is the maximum reported process RSS, not summed shared RSS.
import re
resource=[]
def time_record(p):
 s=p.read_text(errors='replace');out={}
 for k,label in [('user_seconds','User time (seconds)'),('system_seconds','System time (seconds)'),('max_process_rss_kib','Maximum resident set size (kbytes)')]:
  a=re.findall(re.escape(label)+r':\s*([\d.]+)',s);out[k]=float(a[-1]) if a else 0.
 a=re.findall(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([^\n]+)',s)
 if a:
  nums=[float(z) for z in a[-1].strip().split(':')];out['wall_seconds']=sum(v*60**i for i,v in enumerate(nums[::-1]))
 else:out['wall_seconds']=0.
 return out
for r in runs:
 section,run=r['section'],r['run'];ps=[Path('mp10_spatial_communication_v1/logs')/f'scc_{section}_{"main" if run=="main" else "range3dNN_parallel"}.log']
 if section=='GSE171351_A1' and run=='main':ps.append(Path('mp10_spatial_communication_v1/logs/scc_GSE171351_A1_aggregation.log'))
 ps=[p if p.exists() else Path(str(p).replace('_range3dNN_parallel.log','_range3dNN.log')) for p in ps];ps=[p for p in ps if p.exists()];data=[time_record(p) for p in ps];cfg={x['parameter']:x['value'] for x in read(BASE/section/run/'parameters.tsv')}
 resource.append(dict(section=section,run=run,wall_seconds=sum(x['wall_seconds'] for x in data),user_seconds=sum(x['user_seconds'] for x in data),system_seconds=sum(x['system_seconds'] for x in data),max_process_rss_kib=max(x['max_process_rss_kib'] for x in data),future_workers=cfg['threads'],BLAS_OMP_threads=1,logs=';'.join(str(p) for p in ps),resource_note='max GNU-time process RSS; shared worker RSS must not be naively summed'))
write(resource,BASE/'resource_summary.tsv')
for r in runs:
 extra=next(x for x in resource if x['section']==r['section'] and x['run']==r['run']);r.update({k:v for k,v in extra.items() if k not in ('section','run')})
write(runs,BASE/'run_summary.tsv')
