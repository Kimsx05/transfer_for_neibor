from common import *
from geometry import pixel,clear_paths
from scipy.spatial import cKDTree
from spatial import groups
import importlib.metadata as meta
def run():
 rows=[];bins=[];lost=[];qc=[];dep=[];pool=[];panel=[];rflags=[]
 for s in HD+XE:
  hd=s in HD;v,d=folders(s);z=labels(s);i=pd.read_parquet(d/'identity.parquet');e=i[i.candidate_union];hu=pd.read_csv(paths(s)/'human_review.tsv',sep='\t');ng=pd.read_parquet(v/'graph_nodes.parquet').set_index('id');w=np.load(v/'window.npz');rc=pixel(z[['x_um','y_um']].values,w['affine']);rnd=np.random.default_rng(CFG['seed']);st=e.grade.value_counts();gs=(v/'counts/genes.txt').read_text().splitlines() if hd else (V/'shared/common_panel.txt').read_text().splitlines()
  for k in PAN+CORE:panel.append((s,k,k in gs,i[k].notna().sum(),i[k].gt(0).sum() if k in gs else None))
  assert (i.grade.isin(['supported','possible'])<=i.reliable_T).all() and len(hu.id.unique())==len(hu)
  for a in ['A1_model','A2_FOXP3','A3_combination','A4_secondary']:
   bins.append((s,a,e[a].sum(),(e[a]&e.grade.eq('supported')).sum()))
  for a in [k for k in e if k.startswith('old_') and k!='old_gate_failures']:
   lost.append((s,a,e[a].sum(),(e[a]&e.grade.eq('supported')).sum(),(e[a]&e.A1_model).sum()))
  pool.append(dict(sample=s,units=len(z),candidates=len(e),**{g:int(st.get(g,0)) for g in ['supported','possible','unresolved','unsupported']},reliable_T=i.reliable_T.sum(),old_supported=(z.supported=='Treg').sum(),retained_old=(z.supported.eq('Treg')&i.grade.eq('supported')).sum(),new_supported=(z.supported.ne('Treg')&i.grade.eq('supported')).sum(),human_requested=len(hu),human_pending=int((~i.loc[hu.id,'human_review_status'].eq('human-reviewed')).sum())))
  for res in ([16,32] if hd else ['cell']):
   cs=pd.read_parquet(d/f'states_{res}.parquet');ee=pd.read_parquet(d/f'centers_{res}.parquet');assert cs.index[cs.eligible].equals(ee.index)
   keys=['MP10'] if hd else STATES
   for key in keys:assert np.array_equal(ee[key]>ee[key].quantile(.95),ee[key+'_high'].to_numpy(bool))
   er=pixel(ee[['x_um','y_um']].values,w['affine']);ec=ee.component.to_numpy('i4');nc=ng.component.to_numpy('i4')
   if res==32:parent=pd.read_parquet(d/'parent32.parquet').loc[z.index].region_id32.to_numpy()
   for r in [50,200]:
    g=sparse.load_npz(d/f'epi_{res}_{r}.npz');valid=ee[f'coverage_{r}']>=.5;jj=rnd.choice(g.nnz,min(g.nnz,10000),replace=False);b=g.indices[jj];a=np.searchsorted(g.indptr,jj,side='right')-1;assert clear_paths(a,b,er,rc,w['mask'],ec,nc).all()
    assert np.all(np.linalg.norm(ee[['x_um','y_um']].values[a]-z[['x_um','y_um']].values[b],axis=1)<=r+1e-7)
    # Exhaustive source/component/self checks in chunks of 2048 rows.
    for start in range(0,len(ee),2048):
     sub=g[start:start+2048];u,j=sub.nonzero();u=u+start;assert (ec[u]==nc[j]).all()
     assert ((parent[j] if res==32 else z.index.values[j])!=ee.index.values[u]).all()
    nb=pd.read_parquet(d/f'neighborhood_{res}_{r}.parquet');assert np.array_equal(nb.reliable_T,g@i.reliable_T.to_numpy('i4'));assert (nb.supported<=nb.reliable_T).all() and (nb.broad<=nb.reliable_T).all();assert np.array_equal(nb.valid,valid&(np.diff(g.indptr)>0))
    for k in keys:
     h=ee[k+'_high'].to_numpy(bool)
     for a0,a1 in [('high',h),('low',~h)]:
      nd=nb.loc[nb.valid&a1];dep.append((s,res,r,k,a0,len(nd),len(ee),nd.neighbors.median(),nd.coverage.median(),(nd.known/nd.neighbors).mean(),nd.reliable_T.median(),(nd.reliable_T==0).sum(),(nd.reliable_T>=3).sum(),(nd.reliable_T>=5).sum()))
    rows.append((s,res,r,len(ee),g.nnz,valid.sum(),len(jj),'all source, component, no-self; sampled line visibility passed'))
   ix=np.arange(len(ee));gr=groups(ee,400);groups0=pd.Series(gr).groupby(gr).indices
   for b in groups0.values():ix[b]=rnd.permutation(b)
   h=ee[[k+'_high' for k in keys]].to_numpy(bool);assert np.array_equal(gr,gr[ix]) and len(np.unique(ix))==len(ix)
   for b in groups0.values():assert np.array_equal(h[b].sum(0),h[ix[b]].sum(0))
   qc.append((s,res,'Q95, Epi denominator, joint conditional exchange, IDs passed'))
 put(pd.DataFrame(rows,columns=['sample','res','r','eligible','edges','coverage_ge_half','visibility_sample','checks']),OUT/'audit/graph_checks.tsv')
 put(pd.DataFrame(qc,columns=['sample','res','check']),OUT/'audit/state_checks.tsv')
 put(pd.DataFrame(dep,columns=['sample','res','r','state','group','valid','eligible','all_neighbors_median','coverage_median','known_fraction','T_median','T0','Tge3','Tge5']),OUT/'reports/denominators.tsv')
 put(pd.DataFrame(pool),OUT/'reports/review_summary.tsv');put(pd.DataFrame(bins,columns=['sample','route','called','supported']),OUT/'reports/candidate_routes.tsv');put(pd.DataFrame(lost,columns=['sample','old_gate','all_candidates','new_supported','old_first_Treg']),OUT/'reports/gate_loss.tsv')
 put(pd.DataFrame(panel,columns=['sample','marker','measured','nonmissing_units','detected_units']),OUT/'audit/T_panel.tsv')
if __name__=='__main__':run()
