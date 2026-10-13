from common import *
import gc

def main(out):
 P=out/'02_gene_evidence';d=pd.read_parquet(V4/'04_niche_comparison/anchor_MP10_all_layers.parquet');res=[];cache=[];skip=[]
 for sec,g in d.groupby('section_id',sort=True):
  dest=P/'by_section'/f'{sec}.npz';meta=dest.with_suffix('.parquet')
  if dest.exists() and meta.exists():cache.append(pd.read_parquet(meta));continue
  if g.n_signature_available_model.iloc[0]<20:skip.append(dict(section_id=sec,status='PARTIAL19_EXCLUDED_FROM_FULL20',reason='HMGCS1 absent; not a leave-one-out experiment'));continue
  a,C,s=params(sec);genes=a.var_names.tolist();ix=pd.Index(genes).get_indexer(TOP);ids=a.obs_names.get_indexer(g.barcode);assert (ids>=0).all();w=np.asarray(s['w_sf'],float);mg=np.asarray(s['m_g'],float);amb=np.asarray(s['s_g_gene_add'],float);batch=a.obs['_scvi_batch'].to_numpy(int);arrays={k:[] for k in ['obs_counts','epi_counts','obs_rank','epi_rank','obs_contrib','epi_contrib']};md=g.copy();depth=[];nf=[];typefrac=[]
  for st in range(0,len(g),128):
   ii=ids[st:st+128];Y=sparse.csr_matrix(a.X[ii]).toarray().astype(float);E=allocate(Y,C,w[ii],mg,amb[batch[ii]]);ys,yc,yr=auc_parts(Y,genes);es,ec,er=auc_parts(E,genes);np.testing.assert_allclose(es,g.epi_MP10_model_universe.iloc[st:st+128],rtol=1e-8,atol=1e-10);np.testing.assert_allclose(ys,g.spot_MP10_model_universe.iloc[st:st+128],rtol=1e-8,atol=1e-10)
   for key,val in [('obs_counts',Y[:,ix]),('epi_counts',E[:,ix]),('obs_rank',yr),('epi_rank',er),('obs_contrib',yc),('epi_contrib',ec)]:arrays[key].append(val)
   depth.extend(Y.sum(1));nf.extend((Y>0).sum(1));typefrac.extend(w[ii]/w[ii].sum(1,keepdims=True))
  arrays={k:np.vstack(v) for k,v in arrays.items()};md['model_UMI']=depth;md['model_nFeature']=nf;md['n_model_genes']=len(genes);md['AUC_denominator']=sum(int(np.ceil(.05*len(genes)))-j for j in range(1,21));tf=np.array(typefrac)
  for typ in ['GC B','SPP1 TAM','Neutrophil']:md['center_mean_'+typ]=tf[:,list(C.columns).index(typ)]
  np.savez_compressed(safe(out,dest),**arrays,spot_uid=md.spot_uid.to_numpy(dtype=str));md.to_parquet(safe(out,meta),index=False);cache.append(md);log('gene extraction '+sec);del a;gc.collect()
 md=pd.concat(cache,ignore_index=True);md.to_parquet(safe(out,P/'spatial_technical_metadata.parquet'),index=False)
 for sec,g in md.groupby('section_id'):
  ar=np.load(P/'by_section'/f'{sec}.npz');assert np.array_equal(g.spot_uid,ar['spot_uid']);accept=g.accepted_label.ge(0).to_numpy()
  for n in range(9):
   sel=accept & g.accepted_label.eq(n).to_numpy()
   if not sel.any():continue
   Y=ar['obs_counts'][sel];E=ar['epi_counts'][sel];dep=g.model_UMI.to_numpy()[sel]
   for j,gene in enumerate(TOP):res.append(dict(section_id=sec,dataset_id=g.dataset_id.iloc[0],analysis_unit_uid=g.analysis_unit_uid.iloc[0],niche_name=f'niche{n}',gene=gene,n=int(sel.sum()),raw_detection_rate=(Y[:,j]>0).mean(),raw_mean_counts=Y[:,j].mean(),raw_mean_per10k=(Y[:,j]/dep*1e4).mean(),epi_mean_counts=E[:,j].mean(),epi_median_counts=np.median(E[:,j]),epi_total_counts=E[:,j].sum(),epi_share=E[:,j].sum()/Y[:,j].sum() if Y[:,j].sum()>0 else np.nan,observed_rank_median=np.median(ar['obs_rank'][sel,j]),epi_rank_median=np.median(ar['epi_rank'][sel,j]),observed_cutoff_fraction=(ar['obs_contrib'][sel,j]>0).mean(),epi_cutoff_fraction=(ar['epi_contrib'][sel,j]>0).mean(),observed_auc_mean=ar['obs_contrib'][sel,j].mean(),epi_auc_mean=ar['epi_contrib'][sel,j].mean(),gene_universe_hash=g.model_gene_universe_hash.iloc[0]))
 write_tsv(out,P/'gene_by_section_niche.tsv.gz',res);write_tsv(out,P/'partial_signature_exclusions.tsv',skip);log('gene summaries complete')
if __name__=='__main__':main(OUT)
