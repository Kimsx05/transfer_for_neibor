from base_data import *
import h5py
from scipy import sparse as sp
from geometry import clear_paths,pixel
checks=[];rng=np.random.default_rng(1729);ap=[]
def check(s,test,value):
 checks.append((s,test,bool(value)));assert value,(s,test)
for s in HD+XE:
 plat='hd' if s in HD else 'xenium';d=paths(s);a=dataset(s);z=pd.read_parquet(d/'analysis_index.parquet')
 check(s,'one original object per row, aligned to counts',a.obs_names.equals(z.index) and z.index.is_unique)
 check(s,'four-level state whitelist',set(z.review_level)<=set(['fine_supported','major_only','mixed_or_suspect','major_unresolved']))
 check(s,'fine labels have parent lineage',z.loc[z.review_level.eq('fine_supported'),'fine_identity'].map(tax.major).eq(z.loc[z.review_level.eq('fine_supported'),'major_identity']).all())
 check(s,'unique qualified nodes receive pooled model predictions',(z.niche_id.ne('unassigned_QC')==z.eligible_expression).all())
 gr=sp.load_npz(d/'niche_graph.npz');gids=pd.read_csv(d/'niche_graph_ids.tsv.gz',header=None)[0]
 check(s,'niche graph IDs include exactly qualified objects in frozen order',np.array_equal(gids,z.index[z.eligible_expression]))
 mat=gr.tocoo();check(s,'symmetric graph and no self edges',(gr!=gr.T).nnz==0 and np.all(mat.row!=mat.col))
 xy=z.loc[gids,['x_um','y_um']].values;check(s,'full graph distance cap',np.linalg.norm(xy[mat.row]-xy[mat.col],axis=1).max()<=30.00001)
 w=np.load(V/plat/s/'window.npz');rc=pixel(xy,w['affine']);co=np.ones(len(xy),np.int32);i=rng.choice(len(mat.data),min(10000,len(mat.data)),replace=False);check(s,'sampled edge interior inside section mask',clear_paths(mat.row[i],mat.col[i],rc,rc,w['mask'],co,co).all())
 if s in HD:
  members=pd.read_parquet(d/'members.parquet');check(s,'every bin once and retained primary assignments identical',members.barcode.is_unique and (members.loc[members.cell_label>0,'final_label']==members.loc[members.cell_label>0,'cell_label']).all())
  check(s,'every nonzero cell label maps to exactly one count row',np.array_equal(np.unique(members.final_label[members.final_label>0]),a.obs.label))
  c=json.loads((d/'conservation.json').read_text());check(s,'gene conservation and assignment totals',c['per_gene_check'] and c['assigned_umi']==a.X.sum() and c['assigned_umi']+c['unassigned_umi']==c['source_umi'])
  ww=pd.read_parquet(d/'tacco_weights.parquet');good=ww.notna().all(1);check(s,'TACCO weights aligned and normalized only for fitted rows',ww.index.equals(z.index) and np.allclose(ww[good].sum(1),1,atol=1e-4))
 else:
  audit=pd.read_parquet(d/'old_gate_audit.parquet');check(s,'original gate reconstruction aligned',audit.index.equals(z.index))
 ap.append(dict(platform=plat,sample=s,original_units=len(z),niche_qualified=int(z.eligible_expression.sum()),fine_supported=int(z.review_level.eq('fine_supported').sum()),major_only=int(z.review_level.eq('major_only').sum()),mixed=int(z.review_level.eq('mixed_or_suspect').sum()),unknown=int(z.review_level.eq('major_unresolved').sum()),index=str(d/'analysis_index.parquet'),counts=str(d/'cells.h5ad') if s in HD else pd.read_csv(V/'shared/xenium_inputs.tsv',sep='\t').query('sample==@s and role=="expression"').path.item(),cell_id='original_cell_id within sample; use composite (sample_id,original_cell_id)',annotation='Codex-reviewed / human pending'))
write(pd.DataFrame(checks,columns=['sample','test','pass']),'delivery_checks.tsv');write(pd.DataFrame(ap),'analysis_objects.tsv')

