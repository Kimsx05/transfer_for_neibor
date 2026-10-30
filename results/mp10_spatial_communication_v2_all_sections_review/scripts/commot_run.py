from commot_loader import *
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.spatial import cKDTree
import argparse,json,time,resource,subprocess,gc
p=argparse.ArgumentParser();p.add_argument('section');p.add_argument('--sequential-components',action='store_true');aargs=p.parse_args();section=aargs.section
out=O/'commot'/section/'primary2dNN';out.mkdir(parents=True,exist_ok=True);inp=O/'input_review/sections'/section;t0=time.time()
def dump(obj,path):
 tmp=Path(str(path)+'.tmp');tmp.write_text(json.dumps(obj,indent=2,ensure_ascii=False));tmp.replace(path)
def save(df,name,index=False):df.to_csv(out/name,sep='\t',index=index,na_rep='NA')
if (out/'COMPLETED.json').exists():print('REUSE',out,flush=True);exit()
meta=pd.read_csv(inp/'meta.tsv',sep='\t');coords=pd.read_csv(inp/'coords.tsv',sep='\t');info=json.loads((inp/'input_summary.json').read_text());assert meta.spot_id.equals(coords.spot_id)
assert (inp/'READY.json').exists()
if not aargs.sequential_components:
 from commot_parallel import cot_combine_parallel
 spatial.cot_combine_sparse=cot_combine_parallel
if (out/'NETWORK_COMPLETED.json').exists():
 params=json.loads((out/'parameters.json').read_text());T=sparse.load_npz(out/'position_total_with_diagonal.npz');print('RESUME AGGREGATION',flush=True)
else:
 genes=pd.read_csv(inp/'genes.tsv',sep='\t')['gene'].astype(str);counts=sparse.csr_matrix(sparse.load_npz(inp/'counts_spot_gene.npz'),dtype=float);assert counts.shape==(len(meta),len(genes));assert np.isfinite(counts.data).all() and (counts.data>=0).all()
 libs=np.asarray(counts.sum(1)).ravel();assert (libs>0).all();x=sparse.csr_matrix(sparse.diags(1e4/libs)@counts);x.data=np.log1p(x.data)
 data=ad.AnnData(x,obs=meta.set_index('spot_id'),var=pd.DataFrame(index=genes.tolist()));data.obsm['spatial']=coords[['x_um','y_um']].to_numpy()
 observed_dnn=float(np.median(cKDTree(data.obsm['spatial']).query(data.obsm['spatial'],k=2)[0][:,1]));ranges={'Secreted Signaling':200.000001,'ECM-Receptor':200.000001,'Cell-Cell Contact':105.000001}
 db=database.ligand_receptor_database(database='CellChat',species='human',signaling_type=None).copy();db.columns=['ligand','receptor','pathway','category'];db=db[db.category.isin(ranges)]
 detect=np.asarray((counts>0).sum(0)).ravel()/counts.shape[0];det=dict(zip(genes,detect));gset=set(genes)
 db['all_subunits_measured']=[all(g in gset for z in (l,r) for g in z.split('_')) for l,r in zip(db.ligand,db.receptor)]
 db['min_subunit_detection']=[min(det.get(g,0) for z in (l,r) for g in z.split('_')) for l,r in zip(db.ligand,db.receptor)]
 db['passes_5pct_all_subunits']=db.all_subunits_measured&(db.min_subunit_detection>=.05);save(db,'gene_lr_coverage.tsv')
 filtered=database.filter_lr_database(db[['ligand','receptor','pathway','category']],data,heteromeric=True,heteromeric_rule='min',heteromeric_delimiter='_',filter_criteria='min_cell_pct',min_cell_pct=.05)
 if len(filtered)==0:raise RuntimeError('No LR survive frozen expression/measured-subunit rules; no alternative selection performed')
 filtered.columns=['ligand','receptor','pathway','category'];save(filtered,'passing_db_records.tsv')
 lr=filtered.groupby(['ligand','receptor'],sort=True).agg(pathway=('pathway',lambda v:'|'.join(sorted(set(v)))),category=('category',lambda v:'|'.join(sorted(set(v))))).reset_index();lr['distance_um']=lr.category.map(lambda v:min(ranges[c] for c in v.split('|')));save(lr,'included_lr.tsv')
 dis={(l,r):200.000001 for l in lr.ligand.unique() for r in lr.receptor.unique()};dis.update({(r.ligand,r.receptor):r.distance_um for r in lr.itertuples()})
 params=dict(method='COMMOT',version='0.0.3',commit='d117445bc07eaa19109c7609e97d9e35b26e99ca',section=section,variant='primary2dNN',n_spots=len(meta),n_measured_genes=len(genes),n_lr=len(lr),database='bundled CellChatDB human',categories=list(ranges),platform_pitch_um=100,observed_dNN_um=observed_dnn,numeric_boundary_tolerance_um=1e-6,ranges_um=ranges,normalization='log1p(counts / all measured gene library * 10000)',heteromeric=True,heteromeric_rule='min',min_subunit_spot_detection=.05,cot_nitermax=10000,cot_eps_p=.1,cot_rho=10,cot_weights=[.25]*4,component_scheduling='native sequential' if aargs.sequential_components else '4 complete native components parallel; same full LR input and original sum order',n_component_workers=1 if aargs.sequential_components else 4,smooth=False,group_permutation=False,position_diagonal='inference retains diagonal; primary aggregation removes it',primary_metric='abundance_weighted_mean_cross_spot_strength',solver_convergence_status='not exposed by official API; finite nonnegative is not proof of mathematical convergence',input_dir=str(inp),input_summary=info,pythonhashseed=os.environ.get('PYTHONHASHSEED','unset'))
 dump(params,out/'parameters.json');print('START',section,'spots',len(meta),'LR',len(lr),flush=True)
 del counts;gc.collect()
 spatial.spatial_communication(data,database_name='CellChat_human_all',df_ligrec=pd.DataFrame(lr[['ligand','receptor','pathway']].to_numpy()),heteromeric=True,heteromeric_rule='min',heteromeric_delimiter='_',dis_thr=dis,cot_nitermax=10000,pathway_sum=False,smooth=False)
 print('INFERENCE DONE seconds',time.time()-t0,flush=True)
 T=data.obsp['commot-CellChat_human_all-total-total'].tocsr();assert np.isfinite(T.data).all() and (T.data>=0).all();sparse.save_npz(out/'position_total_with_diagonal.npz',T)
 netdir=out/'position_lr';netdir.mkdir(exist_ok=True);audit=[];sumcheck=sparse.csr_matrix(T.shape)
 for i,row in enumerate(lr.itertuples()):
  mat=data.obsp['commot-CellChat_human_all-'+row.ligand+'-'+row.receptor].tocsr();sparse.save_npz(netdir/f'{i:04d}.npz',mat);sumcheck+=mat
  off=mat.copy();off.setdiag(0);off.eliminate_zeros()
  audit.append(dict(lr_index=i,ligand=row.ligand,receptor=row.receptor,category=row.category,total_flow_with_diagonal=float(mat.sum()),same_spot_flow=float(mat.diagonal().sum()),cross_spot_flow=float(off.sum()),n_cross_links=off.nnz,n_links_with_diagonal=mat.nnz,file=f'position_lr/{i:04d}.npz'))
 diff=sumcheck-T;assert diff.nnz==0 or np.max(abs(diff.data))<1e-7;save(pd.DataFrame(audit),'position_lr_index.tsv');dump(dict(status='complete',elapsed_seconds=time.time()-t0,n_lr=len(lr)),out/'NETWORK_COMPLETED.json');del data,x,sumcheck;gc.collect()
diag=T.diagonal().copy();cross=T.copy();cross.setdiag(0);cross.eliminate_zeros();sparse.save_npz(out/'position_total_cross.npz',cross)
spot=coords.copy();spot['same_spot_flow']=diag;spot['cross_sender']=np.asarray(cross.sum(1)).ravel();spot['cross_receiver']=np.asarray(cross.sum(0)).ravel();save(spot,'spot_communication.tsv')
support_table=pd.read_csv(inp/'type_support.tsv',sep='\t')
for level in ['major','fine']:
 qdf=pd.read_csv(inp/f'Q_{level}.tsv',sep='\t').set_index('spot_id');assert list(qdf.index)==list(meta.spot_id);Q=qdf.to_numpy();N=qdf.columns;assert np.allclose(Q.sum(1),1,atol=1e-8)
 F=Q.T@cross@Q;D=np.outer(Q.sum(0),Q.sum(0))-Q.T@Q;G=np.divide(F,D,out=np.full_like(F,np.nan),where=D>0);same=Q.T@(diag[:,None]*Q);B=cross.copy();B.data[:]=1;H=(Q>0).astype(float);K=H.T@B@H
 assert np.allclose(F.sum(),cross.sum()) and np.allclose(same.sum(),diag.sum())
 elig=support_table.query('level==@level').set_index('node').type_evaluable.reindex(N).fillna(False).to_numpy(dtype=bool);mask=elig[:,None]&elig[None,:]&(D>0)
 for name,M in [('total_strength',G),('total_flow',F),('denominator',D),('same_spot_flow',same),('positive_cross_links',K)]:
  # Keep raw model aggregates separately; primary inferential tables use existing type eligibility.
  save(pd.DataFrame(M,index=N,columns=N),f'{level}_{name}_unmasked.tsv',True);save(pd.DataFrame(np.where(mask,M,np.nan),index=N,columns=N),f'{level}_{name}.tsv',True)
 rows=[]
 for i,s in enumerate(N):
  for j,t in enumerate(N):
   ok=bool(mask[i,j]);rows.append(dict(section=section,variant='primary2dNN',network_level=level,source=s,target=t,strength=G[i,j] if ok else np.nan,total_flow=F[i,j] if ok else np.nan,denominator=D[i,j],same_spot_flow=same[i,j],positive_cross_links=K[i,j],source_effective_abundance=Q[:,i].sum(),target_effective_abundance=Q[:,j].sum(),source_support_spots=int((Q[:,i]>0).sum()),target_support_spots=int((Q[:,j]>0).sum()),computed=ok,reason_NE='' if ok else 'type_support_or_pair_denominator_insufficient',group_p_computed=False))
 save(pd.DataFrame(rows),f'{level}_edges.tsv')
result={**params,'elapsed_seconds':time.time()-t0,'peak_rss_GiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,'total_flow_including_same_spot':float(T.sum()),'same_spot_flow':float(diag.sum()),'cross_spot_flow':float(cross.sum()),'status':'complete','hl_status':info.get('hl_status')};dump(result,out/'COMPLETED.json');print('COMPLETED',section,flush=True)
