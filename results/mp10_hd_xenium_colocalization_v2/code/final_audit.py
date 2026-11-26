from common import *
from report import read,put
import hashlib,importlib.metadata as metadata
from scipy.spatial.distance import jensenshannon
def run():
 root=OUT;out=root/'reports';checks=[];inputs=[];support=[]
 for s in HD+XE:
  hd=s in HD;v,d=folders(s);q=paths(s);r=pd.read_csv(q/'human_review.tsv',sep='\t');idx=pd.read_csv(q/'window_index.tsv',sep='\t');cand=pd.read_parquet(q/'candidates.parquet')
  checks.append((s,'all recalled IDs have a browsable local window',set(cand.index)<=set(idx.id)))
  checks.append((s,'review pairs',all((q/'cards'/(k+ext)).is_file() for k in r.id for ext in ['.png','.pdf'])))
  checks.append((s,'new labels retain pending human signature',r.human_reviewer.isna().all() and r.human_status.eq('pending').all()))
  checks.append((s,'actual reviewed images traced',cand.review_mode.eq('table_and_image').sum()==2))
  n=OUT/'niche'/s;fields=['MENDER.parquet','BuildNicheAssay.tsv.gz','CellCharter_all.parquet','CellCharter_masked.parquet','within_niche_effects.tsv'];checks.extend((s,k,(n/k).exists()) for k in fields)
  for p in [v/'graph_nodes.parquet',v/'window.npz',v/'graph_50.npz',v/'graph_200.npz',v/'major_weights.rds',v/'panel.tsv.gz',v/'conservation.json',v/'RCTD_config.json']+([ROOT/'highres_VisiumHD/analysis/hd_treg_mp10_spatial_all_v1'/s/'T32_members.parquet'] if hd else []):
   if p.exists():inputs.append((s,str(p),p.stat().st_size,p.stat().st_mtime_ns))
  keys=['MP10'] if hd else STATES
  checks.extend((s,k,(d/k).exists()) for k in ['squidpy_50.parquet','CRAWDAD_formal.tsv','spatstat_inhomogeneous_K.tsv','reverse_effects.tsv.gz','fixed_target_curve.tsv']+[f'spatstat_envelope_{k}.tsv' for k in keys])
  # High membership counts and overlap are independently reconstructed from score tables.
  e=pd.read_parquet(d/('centers_16.parquet' if hd else 'centers_cell.parquet'))
  for k in keys:checks.append((s,k+' exact full eligible high',np.array_equal(e[k].gt(e[k].quantile(.95)),e[k+'_high'].to_numpy(bool))))
  if not hd:
   for a,b in itertools.combinations(keys,2):
    h=e[a+'_high'].to_numpy(bool);j=e[b+'_high'].to_numpy(bool);t=read(d/'high_overlap.tsv');row=t[(t.state1==a)&(t.state2==b)].iloc[0];checks.append((s,a+'/'+b+' overlap',row.intersection==(h&j).sum() and row.union==(h|j).sum()))
  cc=read(n/'CellCharter_stability.tsv')
  for name,grp in cc.groupby('genes'):support.append((s,name,len(grp),grp.converged.sum(),int(grp.iterations.max())))
 for p in [V/'shared'/k for k in ['taxonomy.tsv','colors.json','programs.json','common_panel.txt','samples.tsv','frozen_epithelial_truth.parquet','input_manifest.tsv','xenium_inputs.tsv','versions.json']]:
  inputs.append(('shared',str(p),p.stat().st_size,p.stat().st_mtime_ns))
  if p.name in ['taxonomy.tsv','colors.json','programs.json','common_panel.txt']:
   dst=OUT/'config'/p.name
   if not dst.exists():dst.symlink_to(p)
 put(pd.DataFrame(inputs,columns=['sample','path','bytes','mtime_ns']),OUT/'audit/additional_inputs.tsv')
 put(pd.DataFrame(checks,columns=['sample','check','passed']),OUT/'audit/delivery_checks.tsv');put(pd.DataFrame(support,columns=['sample','representation','fits','converged','max_iterations']),out/'CellCharter_convergence.tsv')
 print('verified',len(checks),'checks',sum(v[2] for v in checks),'passing')
if __name__=='__main__':run()
