from base import *
import anndata as ad,torch,gc,traceback
from scipy.stats import rankdata
from cell2location.models import Cell2location
TOP=pd.read_csv(OLD/'00_preflight/MP10_top20.tsv',sep='\t').gene.tolist()
def gh(g):return hashlib.sha256('\n'.join(map(str,g)).encode()).hexdigest()
def auc(x,genes):
 ix=pd.Index(genes).get_indexer([g for g in TOP if g in genes]);k=int(np.ceil(.05*len(genes)));den=sum(k-r for r in range(1,min(len(ix),k)+1))
 return np.maximum(0,k-rankdata(-x,axis=1,method='average')[:,ix]).sum(1)/den if den else np.full(len(x),np.nan)
def loadmodel(p,a):return Cell2location.load(str(p/'model'),adata=a,accelerator='cpu',device='auto')
def matrix(p,a):
 v=torch.load(p/'model/model.pt',map_location='cpu',weights_only=False);d=v['attr_dict'];ref=d['cell_state_df_'];assert np.array_equal(ref.index.astype(str),a.var_names.astype(str));assert np.array_equal(a.uns['mod']['obs_names'].astype(str),a.obs_names.astype(str));assert np.array_equal(a.uns['mod']['var_names'].astype(str),a.var_names.astype(str));assert np.array_equal(ref.columns.astype(str),a.uns['mod']['factor_names'].astype(str));assert len(ref.columns)==30
 s=a.uns['mod']['post_sample_means'];C=ref.to_numpy(dtype=float);w=np.asarray(s['w_sf'],float);mg=np.asarray(s['m_g'],float);b=np.asarray(s['s_g_gene_add'],float);batch=a.obs['_scvi_batch'].to_numpy(dtype=int);ep=list(ref.columns).index('Epithelial');assert w.shape==(a.n_obs,30);assert d['registry_']['field_registries']['X']['data_registry']['attr_name']=='X';assert d['registry_']['field_registries']['X']['data_registry']['attr_key'] is None
 return C,w,mg,b,batch,ep

def calc(a,t,ids):
 C,w,mg,b,batch,ep=t;Y=sparse.csr_matrix(a.X[ids]).toarray().astype(float);mu=(w[ids]@C.T)*mg+b[batch[ids]];num=(w[ids,ep,None]*C[None,:,ep])*mg
 assert np.isfinite(mu).all() and (mu>=0).all();assert not ((mu<=0)&(Y>0)).any()
 E=np.divide(Y*num,mu,out=np.zeros_like(Y),where=mu>0);A=np.divide(Y*b[batch[ids]],mu,out=np.zeros_like(Y),where=mu>0)
 assert np.isfinite(E).all() and (E>=0).all() and np.all(E<=Y+1e-6)
 return Y,E,A

def main():
 m=pd.read_parquet(V3/'01_features/metadata.parquet');lab=pd.read_parquet(V3/'03_main/SAMPLE_R12_LEIDEN/labels.parquet');assert m.spot_uid.is_unique and lab.spot_uid.is_unique
 for c in lab:
  if c not in m:m=m.merge(lab[['spot_uid',c]],on='spot_uid',validate='one_to_one')
 old=pd.read_parquet(OLD/'03_features/common_anchors.parquet');m=m.merge(old[['spot_uid','nCount','nFeature']],on='spot_uid',validate='one_to_one');m=m.rename(columns={'epithelial_fraction':'epi_fraction_center_mean','nCount':'center_nCount','nFeature':'center_nFeature'});m['niche_id']=m.accepted_label;m['niche_name']=m.accepted_label.map({**{i:f'niche{i}' for i in range(9)},-1:'UNASSIGNED'})
 man=pd.read_csv(ROOT/'cell2location_full_run/v2_30type/01_manifest/01_spatial_manifest_frozen.csv');assets=pd.read_csv(OLD/'00_preflight/MP10_score_assets.tsv',sep='\t').set_index('section');audit=[];models=[];resolved=[]
 for r in man.itertuples():
  p=Path(r.spatial_output_path);raw=ad.read_h5ad(r.raw_input_path,backed='r');rg=raw.var_names.tolist();model=p/'spatial_posterior.h5ad';g=[]
  if model.exists():
   a=ad.read_h5ad(model,backed='r');g=a.var_names.tolist();a.file.close()
  for gene in TOP:audit.append(dict(section_id=r.section_id,dataset_id=r.dataset,gene=gene,raw_present=gene in rg,historical_ranking_present=gene in rg,model_present=gene in g,reason='present' if gene in g else ('absent_raw' if gene not in rg else 'filtered_before_model'),alias_replacement=False))
  raw.file.close();models.append(dict(section_id=r.section_id,dataset_id=r.dataset,n_anchors=int((m.section_id==r.section_id).sum()),raw_path=r.raw_input_path,model_path=str(p/'model/model.pt'),posterior_path=str(model),n_raw_genes=len(rg),n_model_genes=len(g),raw_gene_universe_hash=gh(rg),model_gene_universe_hash=gh(g),n_signature_available_raw=sum(v in rg for v in TOP),n_signature_available_model=sum(v in g for v in TOP),missing_model_genes=','.join(v for v in TOP if v not in g),status='PENDING' if len(g) else 'MODEL_MISSING'))
  resolved.append(dict(section_id=r.section_id,raw_path=r.raw_input_path,frozen_score_path=assets.loc[r.section_id,'path'],model_path=str(p/'model/model.pt'),posterior_path=str(model)))
 tab(OUT/'00_audit_and_freeze/MP10_gene_availability_reaudit.tsv',audit);tab(OUT/'00_audit_and_freeze/resolved_inputs.tsv',resolved)
 mm=pd.DataFrame(models);tab(OUT/'03_epi_attribution/model_and_gene_universe_manifest.tsv',mm)
 pilots=mm[mm.n_anchors>0].sort_values(['n_signature_available_model','section_id'],ascending=[False,True]).groupby('dataset_id',sort=True).head(1);tests=[];passed=set()
 for r in pilots.itertuples():
  try:
   a=ad.read_h5ad(r.posterior_path);p=Path(r.model_path).parent.parent;t=matrix(p,a);sub=m[m.section_id.eq(r.section_id)].sort_values('spot_uid').head(128);ids=a.obs_names.get_indexer(sub.barcode);assert (ids>=0).all()
   from types import SimpleNamespace
   from cell2location.models._cell2location_module import LocationModelLinearDependentWMultiExperimentLocationBackgroundNormLevelGeneAlphaPyroModel as Official
   manager=SimpleNamespace(get_from_registry=lambda key:a.X if key=='X' else a.obs['_scvi_batch'].to_numpy()[:,None]);mod=SimpleNamespace(cell_state_mat=t[0],n_factors=30);official=Official.compute_expected_per_cell_type(mod,a.uns['mod']['post_sample_means'],manager,ind_x=ids);Y,E,A=calc(a,t,ids);oe=official['mu'][t[-1]].toarray();err=float(np.max(np.abs(oe-E)));np.testing.assert_allclose(E,oe,rtol=1e-5,atol=1e-6);total=sum(v.toarray() for v in official['mu'])+A;np.testing.assert_allclose(total,Y,rtol=1e-5,atol=1e-6)
   tests.append(dict(dataset_id=r.dataset_id,section_id=r.section_id,n=len(ids),status='PASS',official_epi_max_abs_error=err,alltypes_ambient_max_abs_error=float(np.max(np.abs(total-Y))),method='installed official compute_expected_per_cell_type; exact saved reference + registry adapter; restored model previously verified'));passed.add(r.dataset_id);del mod,a;gc.collect();log(f'pilot PASS {r.section_id}')
  except Exception as e:tests.append(dict(dataset_id=r.dataset_id,section_id=r.section_id,status='FAIL',reason=repr(e)));traceback.print_exc()
 tab(OUT/'03_epi_attribution/attribution_technical_validation.tsv',tests)
 outputs=[];repro=[];prior_repro=pd.read_csv(OUT/'02_spot_MP10/frozen_reproduction.tsv',sep='\t') if (OUT/'02_spot_MP10/frozen_reproduction.tsv').exists() else pd.DataFrame()
 for ix,r in enumerate(mm.itertuples()):
  out=OUT/'03_epi_attribution/by_section'/f'{r.section_id}.parquet';meta=out.with_suffix('.json');sub=m[m.section_id.eq(r.section_id)].copy()
  if out.exists() and meta.exists() and json.loads(meta.read_text())['status'].startswith('COMPLETE'):outputs.append(pd.read_parquet(out));d=json.loads(meta.read_text());mm.loc[ix,'status']=d['status'];repro.extend(prior_repro[prior_repro.section_id.eq(r.section_id)].to_dict('records') if len(prior_repro) else []);log(f'checkpoint {r.section_id}');continue
  if not len(sub):mm.loc[ix,'status']='NO_COMMON_ANCHORS';continue
  raw=ad.read_h5ad(r.raw_path,backed='r');rg=raw.var_names.tolist();ids=raw.obs_names.get_indexer(sub.barcode);assert (ids>=0).all()
  f=pd.read_csv(assets.loc[r.section_id,'path']);assert not f.barcode.duplicated().any();sub=sub.merge(f[['barcode','EpiMP10_Top20_AUCell']],on='barcode',how='left',validate='one_to_one').rename(columns={'EpiMP10_Top20_AUCell':'spot_MP10_frozen'});assert sub.spot_MP10_frozen.notna().all();sub['spot_MP10_frozen_source']=assets.loc[r.section_id,'path'];sub['frozen_definition_id']='historical_Top20_average_rank_AUC5pct';sub['raw_gene_universe_hash']=r.raw_gene_universe_hash
  check=np.arange(min(12,len(sub)));calcscore=auc(sparse.csr_matrix(raw.X[ids[check]]).toarray(),rg);diff=np.abs(calcscore-sub.spot_MP10_frozen.iloc[check].to_numpy());repro.append(dict(section_id=r.section_id,n=len(check),max_abs_error=float(diff.max()),status='PASS' if diff.max()<1e-10 else 'DRIFT_ORIGINAL_RETAINED'))
  for k in ['spot_MP10_model_universe','epi_MP10_model_universe','spot_MP10_model_partial','epi_MP10_model_partial','signature_epi_count_share_of_observed','signature_nonepi_count_share_of_observed','signature_ambient_count_share_of_observed','allgene_epi_count_share_of_observed','signature_observed_counts','epi_detected_genes','epi_total_counts','template_only_MP10']:sub[k]=np.nan
  sub['n_signature_available_raw']=r.n_signature_available_raw;sub['n_signature_available_model']=r.n_signature_available_model;sub['model_gene_universe_hash']=r.model_gene_universe_hash;sub['model_definition_id']='Top20_average_rank_AUC5pct_matched_model_universe';sub['signature_used_genes_hash']=gh([v for v in TOP if v not in str(r.missing_model_genes).split(',')]);sub['missing_model_genes']=r.missing_model_genes;sub['signature_coverage_status']='FULL20' if r.n_signature_available_model==20 else 'PARTIAL_MODEL_SIGNATURE';sub['epi_score_valid']=False;sub['missing_reason']='';status='FAILED'
  try:
   assert r.dataset_id in passed,'Dataset pilot failed; cannot extend'
   a=ad.read_h5ad(r.posterior_path);t=matrix(Path(r.model_path).parent.parent,a);ids=a.obs_names.get_indexer(sub.barcode);assert (ids>=0).all();genes=a.var_names.tolist();si=pd.Index(genes).get_indexer([v for v in TOP if v in genes]);template=(t[0][:,t[-1]]*np.asarray(t[2]).ravel())[None,:];template_score=float(auc(template,genes)[0])
   for start in range(0,len(sub),128):
    sl=np.arange(start,min(start+128,len(sub)));Y,E,A=calc(a,t,ids[sl]);os=auc(Y,genes);es=auc(E,genes);et=E.sum(1);es[et<=0]=np.nan;sy=Y[:,si].sum(1);ep=E[:,si].sum(1);am=A[:,si].sum(1);none=sy-ep-am;assert (none>=-1e-5).all();div=lambda x,y:np.divide(x,y,out=np.full(len(x),np.nan),where=y>0)
    values={'spot_MP10_model_partial':os,'epi_MP10_model_partial':es,'signature_epi_count_share_of_observed':div(ep,sy),'signature_nonepi_count_share_of_observed':div(np.maximum(none,0),sy),'signature_ambient_count_share_of_observed':div(am,sy),'allgene_epi_count_share_of_observed':div(et,Y.sum(1)),'signature_observed_counts':sy,'epi_detected_genes':(E>0).sum(1),'epi_total_counts':et,'template_only_MP10':template_score}
    if r.n_signature_available_model==20:values.update(spot_MP10_model_universe=os,epi_MP10_model_universe=es)
    for k,v in values.items():sub.loc[sub.index[sl],k]=v
   sub['epi_score_valid']=sub.epi_total_counts.gt(0)&sub.n_signature_available_model.eq(20);sub.loc[sub.epi_total_counts.le(0),'missing_reason']='zero_attributed_epi_counts';sub.loc[sub.n_signature_available_model.lt(20),'missing_reason']='model_missing_signature:'+str(r.missing_model_genes);sub['epi_detected_below_rank_cutoff']=sub.epi_detected_genes.lt(np.ceil(.05*len(genes)));status='COMPLETE_FULL20' if r.n_signature_available_model==20 else 'COMPLETE_PARTIAL_SIGNATURE';del a;gc.collect()
  except Exception as e:sub['missing_reason']=repr(e);traceback.print_exc()
  raw.file.close();sub.to_parquet(out,index=False);js(meta,dict(status=status,n=len(sub),method='posterior means plug-in; not posterior mean of attributed expression',signature=r.missing_model_genes));outputs.append(sub);mm.loc[ix,'status']=status;tab(OUT/'03_epi_attribution/model_and_gene_universe_manifest.tsv',mm);log(f'{ix+1}/63 {r.section_id}: {status}, {len(sub)} anchors')
 allm=pd.concat(outputs,ignore_index=True);assert len(allm)==len(m) and allm.spot_uid.is_unique;overlay=pd.read_csv(OUT/'01_pathology_search/spot_pathology_overlay.tsv',sep='\t');allm=allm.merge(overlay[['spot_uid','pathology_label','pathology_evidence_level','pathology_source_id']],on='spot_uid',validate='one_to_one');allm.to_parquet(OUT/'04_niche_comparison/anchor_MP10_all_layers.parquet',index=False);tab(OUT/'02_spot_MP10/frozen_reproduction.tsv',repro);tab(OUT/'03_epi_attribution/model_and_gene_universe_manifest.tsv',mm);log('ALL LAYERS SAVED')
if __name__=='__main__':main()
