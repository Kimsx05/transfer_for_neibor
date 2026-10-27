from commot_loader import *
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse,io
from scipy.spatial import cKDTree
import argparse,json,time,resource,gzip,subprocess
p=argparse.ArgumentParser();p.add_argument('section');p.add_argument('--wide',action='store_true');p.add_argument('--parallel-components',action='store_true');args=p.parse_args()
if args.parallel_components:
 from commot_parallel import cot_combine_parallel
 spatial.cot_combine_sparse=cot_combine_parallel
section=args.section;variant='wide3dNN' if args.wide else 'primary2dNN';dest=O/'commot'/section/variant;dest.mkdir(parents=True,exist_ok=True)
if (dest/'COMPLETED.json').exists():print('REUSE',dest);exit()
t0=time.time();inp=O/'input_review/sections'/section
meta=pd.read_csv(inp/'meta.tsv',sep='\t');coords=pd.read_csv(inp/'coords.tsv',sep='\t');genes=pd.read_csv(inp/'genes.tsv',sep='\t')['gene'].astype(str)
assert meta.spot_id.equals(coords.spot_id)
counts=sparse.csr_matrix(io.mmread(inp/'counts.mtx.gz').T,dtype=float);assert counts.shape==(len(meta),len(genes));assert counts.data.min()>=0
libs=np.asarray(counts.sum(1)).ravel();assert (libs>0).all()
x=sparse.diags(1e4/libs)@counts;x=sparse.csr_matrix(x);x.data=np.log1p(x.data)
a=ad.AnnData(x,obs=meta.set_index('spot_id'),var=pd.DataFrame(index=genes.tolist()));a.obsm['spatial']=coords[['x_um','y_um']].to_numpy()
dnn=float(np.median(cKDTree(a.obsm['spatial']).query(a.obsm['spatial'],k=2)[0][:,1]));dnn=round(dnn,8);tol=1e-6;ranges={'Secreted Signaling':(3 if args.wide else 2)*dnn+tol,'ECM-Receptor':2*dnn+tol,'Cell-Cell Contact':1.05*dnn+tol}
db=database.ligand_receptor_database(database='CellChat',species='human',signaling_type=None).copy();db.columns=['ligand','receptor','pathway','category'];db=db[db.category.isin(ranges)]
# Preserve database records and coverage, collapse duplicate ligand-receptor pairs once.
detect=np.asarray((counts>0).sum(0)).ravel()/counts.shape[0];det=dict(zip(genes,detect));gset=set(genes)
db['all_subunits_measured']=[all(g in gset for z in (l,r) for g in z.split('_')) for l,r in zip(db.ligand,db.receptor)]
db['min_subunit_detection']=[min(det.get(g,0) for z in (l,r) for g in z.split('_')) for l,r in zip(db.ligand,db.receptor)]
db['passes_5pct_all_subunits']=db.all_subunits_measured&(db.min_subunit_detection>=.05)
db.to_csv(dest/'gene_lr_coverage.tsv',sep='\t',index=False)
filtered=database.filter_lr_database(db[['ligand','receptor','pathway','category']],a,heteromeric=True,heteromeric_rule='min',heteromeric_delimiter='_',filter_criteria='min_cell_pct',min_cell_pct=.05)
filtered.columns=['ligand','receptor','pathway','category']
# Exact pair is one competitive transport channel; retain all annotations in audit.
filtered.to_csv(dest/'passing_db_records.tsv',sep='\t',index=False)
lr=filtered.groupby(['ligand','receptor'],sort=True).agg(pathway=('pathway',lambda v:'|'.join(sorted(set(v)))),category=('category',lambda v:'|'.join(sorted(set(v))))).reset_index()
lr['distance_um']=lr.category.map(lambda v:min(ranges[c] for c in v.split('|')))
lr.to_csv(dest/'included_lr.tsv',sep='\t',index=False)
# Current official code visits all ligand x receptor combinations, including non-LR pairs.
dis={(l,r):max(ranges.values()) for l in lr.ligand.unique() for r in lr.receptor.unique()}
dis.update({(r.ligand,r.receptor):r.distance_um for r in lr.itertuples()})
params={'method':'COMMOT','version':'0.0.3','commit':subprocess.check_output(['git','-C',str(O/'vendor/COMMOT'),'rev-parse','HEAD'],text=True).strip(),'section':section,'variant':variant,'n_spots':len(meta),'n_measured_genes':len(genes),'n_lr':len(lr),'database':'bundled CellChatDB human','categories':list(ranges),'dNN_um':dnn,'numeric_boundary_tolerance_um':tol,'ranges_um':ranges,'normalization':'log1p(counts / spot library * 10000)','heteromeric':True,'heteromeric_rule':'min','min_subunit_spot_detection':.05,'cot_nitermax':10000,'cot_eps_p':.1,'cot_rho':10,'cot_weights':[.25]*4,'component_scheduling':('4 complete native components in parallel; identical full LR inputs; original weighted sum order' if args.parallel_components else 'native sequential'),'n_component_workers':(4 if args.parallel_components else 1),'smooth':False,'group_permutation':False,'position_diagonal':'retained during native inference; removed only for primary aggregation','primary_metric':'abundance_weighted_mean_cross_spot_strength','solver_convergence_status':'official API does not expose stopping diagnostics; finite/nonnegative output checked'}
(dest/'parameters.json').write_text(json.dumps(params,indent=2));print('START',params,flush=True)
spatial.spatial_communication(a,database_name='CellChat_human_all',df_ligrec=pd.DataFrame(lr[['ligand','receptor','pathway']].to_numpy()),heteromeric=True,heteromeric_rule='min',heteromeric_delimiter='_',dis_thr=dis,cot_nitermax=10000,pathway_sum=False,smooth=False)
print('INFERENCE DONE sec',time.time()-t0,flush=True)
T=a.obsp['commot-CellChat_human_all-total-total'].tocsr();assert np.isfinite(T.data).all() and (T.data>=0).all()
sparse.save_npz(dest/'position_total_with_diagonal.npz',T);diag=T.diagonal().copy();cross=T.copy();cross.setdiag(0);cross.eliminate_zeros();sparse.save_npz(dest/'position_total_cross.npz',cross)
# Freeze all native LR networks without copying measured expression into a second object.
netdir=dest/'position_lr';netdir.mkdir(exist_ok=True);audit=[];sumcheck=sparse.csr_matrix(T.shape)
for i,row in enumerate(lr.itertuples()):
 key='commot-CellChat_human_all-'+row.ligand+'-'+row.receptor
 mat=a.obsp[key].tocsr();sparse.save_npz(netdir/f'{i:04d}.npz',mat);sumcheck+=mat
 audit.append({'lr_index':i,'ligand':row.ligand,'receptor':row.receptor,'category':row.category,'total_flow_with_diagonal':float(mat.sum()),'same_spot_flow':float(mat.diagonal().sum()),'n_links_with_diagonal':mat.nnz,'file':f'position_lr/{i:04d}.npz'})
diff=sumcheck-T;assert diff.nnz==0 or np.max(np.abs(diff.data))<1e-7
pd.DataFrame(audit).to_csv(dest/'position_lr_index.tsv',sep='\t',index=False)
spot=coords.copy();spot['same_spot_flow']=diag;spot['cross_sender']=np.asarray(cross.sum(1)).ravel();spot['cross_receiver']=np.asarray(cross.sum(0)).ravel();spot.to_csv(dest/'spot_communication.tsv',sep='\t',index=False)
for level in ['major','fine']:
 qdf=pd.read_csv(inp/f'Q_{level}.tsv',sep='\t').set_index('spot_id');assert list(qdf.index)==list(meta.spot_id);Q=qdf.to_numpy();assert np.allclose(Q.sum(1),1,atol=1e-8)
 F=Q.T@cross@Q;D=np.outer(Q.sum(0),Q.sum(0))-Q.T@Q;G=np.divide(F,D,out=np.full_like(F,np.nan),where=D>0);same=Q.T@(diag[:,None]*Q)
 B=cross.copy();B.data[:]=1;support=(Q>0).astype(float);K=support.T@B@support
 for name,M in [('total_strength',G),('total_flow',F),('denominator',D),('same_spot_flow',same),('positive_cross_links',K)]:pd.DataFrame(M,index=qdf.columns,columns=qdf.columns).to_csv(dest/f'{level}_{name}.tsv',sep='\t',na_rep='NA',index_label='source')
 rows=[]
 for i,s in enumerate(qdf.columns):
  for j,t in enumerate(qdf.columns):rows.append(dict(section=section,variant=variant,network_level=level,source=s,target=t,strength=G[i,j],total_flow=F[i,j],denominator=D[i,j],same_spot_flow=same[i,j],positive_cross_links=K[i,j],source_effective_abundance=Q[:,i].sum(),target_effective_abundance=Q[:,j].sum(),source_support_spots=int((Q[:,i]>0).sum()),target_support_spots=int((Q[:,j]>0).sum()),computed=bool(D[i,j]>0),group_p_computed=False))
 pd.DataFrame(rows).to_csv(dest/f'{level}_edges.tsv',sep='\t',index=False,na_rep='NA')
 assert np.allclose(F.sum(),cross.sum()) and np.allclose(same.sum(),diag.sum())
result={**params,'elapsed_seconds':time.time()-t0,'peak_rss_GiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,'total_flow_including_same_spot':float(T.sum()),'same_spot_flow':float(diag.sum()),'cross_spot_flow':float(cross.sum()),'status':'complete'}
(dest/'COMPLETED.json').write_text(json.dumps(result,indent=2));print('COMPLETED',result,flush=True)
