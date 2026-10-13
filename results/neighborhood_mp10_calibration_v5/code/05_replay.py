from common import *
import gc

def main(out):
 P=out/'03_attribution_calibration';base=V4/'05_source_validation';d=pd.read_parquet(V4/'04_niche_comparison/anchor_MP10_all_layers.parquet').set_index('spot_uid');rows=[];genesout=[];aud=[]
 comp=pd.read_csv(base/'composition_only_pseudospots.tsv',sep='\t').set_index('bag_id');state=pd.read_csv(base/'state_change_controls.tsv',sep='\t');state['state_id']=state.blueprint_uid+'::'+state.state;state=state.set_index('state_id')
 for folder in sorted((base/'counts_cache').iterdir()):
  if not folder.is_dir():continue
  axis=json.loads((folder/'axis_manifest.json').read_text());a,C,s=params(folder.name);assert axis['genes']==a.var_names.tolist();genes=axis['genes'];ix=pd.Index(genes).get_indexer(TOP);w=np.asarray(s['w_sf'],float);mg=np.asarray(s['m_g'],float);amb=np.asarray(s['s_g_gene_add'],float);batch=a.obs['_scvi_batch'].to_numpy(int)
  for kind,table in [('composition',comp),('state',state)]:
   keys=axis[kind+'_ids'];p=folder/(kind+'_observed.npz')
   if not keys:continue
   Y=sparse.load_npz(p).tocsr();E=sparse.load_npz(folder/(kind+'_true_epi.npz')).tocsr();assert Y.shape==E.shape==(len(keys),len(genes));assert (Y.data>=0).all() and (E.data>=0).all();assert (Y-E).min()>=0;meta=table.loc[keys];assert meta.index.tolist()==keys;blue=meta.blueprint_uid;ids=a.obs_names.get_indexer(d.loc[blue,'barcode']);assert (ids>=0).all()
   for st in range(0,len(keys),64):
    yy=Y[st:st+64].toarray();ee=E[st:st+64].toarray();ii=ids[st:st+64];eh=allocate(yy,C,w[ii],mg,amb[batch[ii]]);ts,tc,_=auc_parts(ee,genes);hs,hc,_=auc_parts(eh,genes);k=int(np.ceil(.05*len(genes)));tol=1/sum(k-j for j in range(1,21))
    for j in range(len(yy)):
     rr=meta.iloc[st+j];rows.append(dict(bag_id=keys[st+j],blueprint_uid=rr.blueprint_uid,section_id=folder.name,dataset=rr.dataset_id,niche_name=rr.niche_name,branch='FROZEN_WEIGHT_REPLAY',condition='FIXED_COMPOSITION_STATE' if kind=='state' else 'V4_COMPOSITION_STRESS',state=rr.state if kind=='state' else 'fixed',pair_id=rr.blueprint_uid if kind=='state' else keys[st+j],true_score=ts[j],estimated_score=hs[j],score_error=hs[j]-ts[j],depth=yy[j].sum(),true_RNA_fraction=ee[j].sum()/yy[j].sum(),rank_step=tol,background='actual_v4_center',donor='v4_mixed_pool',mismatch='original mg/ambient applied to no-platform/no-ambient v4 simulation; real w not refit'))
     for gi,gene in enumerate(TOP):genesout.append(dict(branch='FROZEN_WEIGHT_REPLAY',bag_id=keys[st+j],background='actual_v4_center',gene=gene,true_counts=ee[j,ix[gi]],estimated_counts=eh[j,ix[gi]],true_auc_contribution=tc[j,gi],estimated_auc_contribution=hc[j,gi]))
   aud.append(dict(section=folder.name,kind=kind,n=len(keys),axis_status='PASS',truth_bounds='PASS'))
  del a;gc.collect();log('replay '+folder.name)
 pd.DataFrame(rows).to_parquet(safe(out,P/'frozen_replay_scores.parquet'),index=False);write_tsv(out,P/'replay_gene_recovery.tsv.gz',genesout);write_tsv(out,P/'replay_axis_audit.tsv',aud);log('all replay complete')
if __name__=='__main__':main(OUT)
