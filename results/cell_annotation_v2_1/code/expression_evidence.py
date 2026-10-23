from common import *
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse,stats
from programs import MAJOR,STATE,WEIGHTS
import sys
sys.path.insert(0,str(V2/'code'));from marker_rules import SUBTYPE,ANCHORS
PANELS={**{'major:'+k:v for k,v in MAJOR.items()},**{'subtype:'+k:v for k,v in SUBTYPE.items()},**{'state:'+k:v for k,v in STATE.items()}}
MARKERS=sorted(set(sum(PANELS.values(),[])+sum(ANCHORS.values(),[])))
def load_all():
 xs=[];os=[];zz=[];uu=[];pp=[];syms=None
 for name in ['Cervilla','D1_Johnson']:
  a=read_section(name);g=np.where(a.var.reference_symbol.astype(str).to_numpy()!='',a.var.reference_symbol.astype(str),a.var.gene_symbol.astype(str));genes=pd.Index(pd.unique(g))
  if syms is None:syms=genes
  assert genes.equals(syms)
  t=sparse.csr_matrix((np.ones(a.n_vars,dtype='int32'),(np.arange(a.n_vars),syms.get_indexer(g))),shape=(a.n_vars,len(syms)));xs.append((a.X@t).tocsr());os.append(pd.read_parquet(ROOT/f'metadata/{name}_v2_obs.parquet'));zz.append(a.obsm['X_scVI']);uu.append(a.obsm['X_umap']);pp.append(a.obsm['spatial_um'])
 return sparse.vstack(xs,format='csr'),pd.concat(os),syms,np.vstack(zz),np.vstack(uu),np.vstack(pp)
def normalize(x):
 q=sparse.diags(1e4/np.maximum(np.asarray(x.sum(1)).ravel(),1)).dot(x.astype('float32')).tocsr();q.data=np.log1p(q.data);return q

def group_evidence(x,o,genes,groupfield,output_prefix,universe=None):
 # Raw detection and measured log-count summaries, within the specified expression universe.
 if universe is None:universe=np.arange(len(o))
 q=normalize(x[universe]);ob=o.iloc[universe];det=q.copy();det.data=np.ones_like(det.data);sq=q.copy();sq.data**=2
 totals=np.asarray(q.sum(0)).ravel();total2=np.asarray(sq.sum(0)).ravel();detections=np.asarray(det.sum(0)).ravel();N=len(ob);marker_ix=genes.get_indexer([g for g in MARKERS if g in genes]);rows=[];summary=[]
 for section in ['joint','Cervilla','D1_Johnson']:
  sec=np.ones(N,bool) if section=='joint' else ob.section_id.eq(section).to_numpy();un=np.flatnonzero(sec);n_un=len(un)
  if n_un==0:continue
  total=totals if section=='joint' else np.asarray(q[un].sum(0)).ravel();t2=total2 if section=='joint' else np.asarray(sq[un].sum(0)).ravel();dtotal=detections if section=='joint' else np.asarray(det[un].sum(0)).ravel()
  for cluster,localnames in ob.iloc[un].groupby(groupfield,observed=True).groups.items():
   ix=ob.index.get_indexer(localnames);n=len(ix);nr=n_un-n;s=np.asarray(q[ix].sum(0)).ravel();s2=np.asarray(sq[ix].sum(0)).ravel();d=np.asarray(det[ix].sum(0)).ravel();mean=s/max(n,1);pct=d/max(n,1);rm=(total-s)/max(nr,1);rp=(dtotal-d)/max(nr,1);v=np.maximum((s2-s*s/max(n,1))/max(n-1,1),0);rv=np.maximum((t2-s2-(total-s)**2/max(nr,1))/max(nr-1,1),0);effect=(mean-rm)/np.log(2);se=np.sqrt(v/max(n,1)+rv/max(nr,1)+1e-8);score=(mean-rm)/se;p=2*stats.norm.sf(np.abs(score));rank=np.argsort(p);fdr=np.empty(len(p));fdr[rank]=np.minimum.accumulate((p[rank]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];fdr=np.minimum(fdr,1)
   eligible=(pct>=.03)&(d>=3)&(effect>0);top=np.flatnonzero(eligible)[np.argsort(-score[eligible])[:40]];wanted=np.union1d(top,marker_ix)
   for j in wanted:rows.append(dict(scope=section,cluster=str(cluster),n=n,comparison_n=nr,gene=genes[j],mean_log1p10k=float(mean[j]),detection_fraction=float(pct[j]),other_detection_fraction=float(rp[j]),mean_log2_difference=float(effect[j]),detection_difference=float(pct[j]-rp[j]),descriptive_z=float(score[j]),descriptive_fdr=float(fdr[j]),data_driven_top=bool(j in top)))
   z=ob.iloc[ix];summary.append(dict(scope=section,cluster=str(cluster),n=n,median_counts=float(z.n_counts.median()),median_genes=float(z.n_genes.median()),low_information_fraction=float(z.low_information.mean()),section_composition=json.dumps(z.section_id.astype(str).value_counts().to_dict()),old_major_composition=json.dumps(z.predicted_major.astype(str).value_counts().head(6).to_dict()),v2_labels=json.dumps(z.final_label.astype(str).value_counts().head(6).to_dict()),top_genes=';'.join(genes[top[:12]])))
 pd.DataFrame(rows).to_parquet(ROOT/'metadata'/f'{output_prefix}_markers.parquet',compression='zstd');pd.DataFrame(summary).to_csv(ROOT/'metadata'/f'{output_prefix}_summary.tsv',sep='\t',index=False);return pd.DataFrame(rows),pd.DataFrame(summary)
if __name__=='__main__':
 x,o,genes,z,u,physical=load_all();idx=genes.get_indexer([g for g in MARKERS if g in genes]);present=genes[idx];mx=x[:,idx].tocsr();sparse.save_npz(ROOT/'metadata/observed_marker_counts.npz',mx);js(ROOT/'metadata/observed_marker_genes.json',present.tolist());o.to_parquet(ROOT/'metadata/cell_metadata.parquet',compression='zstd');np.savez_compressed(ROOT/'metadata/existing_coordinates.npz',cell_id=o.index.to_numpy().astype(str),X_scVI=z,X_umap=u,spatial_um=physical)
 normx=normalize(x);m=normx[:,idx].toarray();raw=mx.toarray();scores={}
 for label,gs in PANELS.items():
  jj=present.get_indexer([g for g in gs if g in present]);w=np.array([WEIGHTS.get(g,1) if label.startswith('major:') else 1 for g in present[jj]]);scores[label]=m[:,jj]@w/np.sqrt(max(w.sum(),1));scores[label+'__detected']=(raw[:,jj]>0).sum(1)
 pd.DataFrame(scores,index=o.index).to_parquet(ROOT/'metadata/program_evidence.parquet',compression='zstd');js(ROOT/'metadata/panels.json',dict(panels=PANELS,weights=WEIGHTS,available={k:[g for g in v if g in genes] for k,v in PANELS.items()},missing={k:[g for g in v if g not in genes] for k,v in PANELS.items()},note='Observed all-query genes including those absent from reference; program strengths are descriptive weighted log-expression, not probabilities. COL1A1/downweighted shared programs mitigate widespread ECM background.'))
 group_evidence(x,o,genes,'expression_cluster','global');print('ALL_GLOBAL_EVIDENCE_COMPLETE',o.shape,x.shape,flush=True)
