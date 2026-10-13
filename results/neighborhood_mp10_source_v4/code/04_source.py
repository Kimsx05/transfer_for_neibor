from base import *
import anndata as ad
from scipy.stats import rankdata,spearmanr
P=OUT/'05_source_validation';O=OLD/'07_source_validation';TOP=pd.read_csv(OLD/'00_preflight/MP10_top20.tsv',sep='\t').gene.tolist()
def auc(x,genes):
 ix=pd.Index(genes).get_indexer([g for g in TOP if g in genes]);k=int(np.ceil(.05*len(genes)));den=sum(k-r for r in range(1,min(len(ix),k)+1));return np.maximum(0,k-rankdata(-x,method='average')[ix]).sum()/den if x.sum()>0 else np.nan

def main():
 d=pd.read_parquet(OUT/'04_niche_comparison/anchor_MP10_all_layers.parquet');d=d[d.accepted_label.ge(0)&d.epi_score_valid].copy();sc=pd.read_csv(O/'single_cell_selected.tsv.gz',sep='\t');genes=pd.read_csv(O/'selected_genes.tsv',sep='\t').gene.tolist();X=sparse.load_npz(O/'selected_counts.npz');assert len(genes)==X.shape[1] and len(sc)==X.shape[0];mapping=pd.read_csv(ROOT/'subtype/v2/reference_training_v2/11_final_spatial_reference/00_32_to_30_label_mapping.csv');mp=mapping.set_index('reference_celltype_v2').spatial_reference_celltype_v2_final;sc['type30']=sc.reference_celltype_v2.map(mp);assert sc.type30.notna().all();tab(P/'verified_32_to_30_mapping.tsv',mapping)
 usage=pd.read_csv(P/'existing_EpiMP10_usage.tsv',sep='\t');assert usage.cell_id.is_unique;sc=sc.merge(usage,on='cell_id',how='left',validate='one_to_one');sc['signature_detected_genes']=np.asarray((X[:,pd.Index(genes).get_indexer(TOP)]>0).sum(1)).ravel();tab(P/'single_cell_reference_scored.tsv.gz',sc)
 tab(P/'single_cell_type30_summary.tsv',sc.groupby(['dataset','reference_batch','type30']).agg(n=('cell_id','size'),mean_score=('MP10_Top20_AUCell','mean'),median_score=('MP10_Top20_AUCell','median'),mean_signature_detected=('signature_detected_genes','mean'),usage_matched=('MP10_usage','count')).reset_index())
 ep=sc.index[sc.type30.eq('Epithelial')].to_numpy();us=sc.loc[ep].dropna(subset=['MP10_usage']);rho=spearmanr(us.MP10_Top20_AUCell,us.MP10_usage).statistic if len(us)>5 else np.nan
 statevar='MP10_usage' if len(us)>=20 else 'MP10_Top20_AUCell';q=sc.loc[ep,statevar].dropna();lo=q.index[q<=q.quantile(.25)].to_numpy();hi=q.index[q>=q.quantile(.75)].to_numpy();tab(P/'single_cell_usage_check.tsv',[dict(n_epi=len(ep),n_matched_usage=len(us),spearman=rho,state_variable=statevar,usage_definition='Existing fixed-W NNLS projection of original sample-balanced EpiMP basis; same definition, not rebuilt NMF',normal_malignant='Not independently verified; all called Epithelial')])
 # Fixed, non-score selection across dataset x niche and central Epi range.
 b=[]
 for _,g in d.groupby(['dataset_id','accepted_label']):
  g=g.sort_values(['epi_fraction_center_mean','spot_uid']);b.append(g.iloc[np.unique(np.linspace(0,len(g)-1,min(6,len(g))).round().astype(int))])
 blue=pd.concat(b).sort_values(['section_id','spot_uid']).head(500);state_allow=set()
 for _,sg in blue.groupby('niche_name'):
  sg=sg.sort_values(['dataset_id','section_id']);state_allow.update(sg.iloc[np.unique(np.linspace(0,len(sg)-1,min(11,len(sg))).round().astype(int))].spot_uid)
 tab(P/'central_composition_blueprints.tsv',blue[['spot_uid','barcode','dataset_id','section_id','niche_name','epi_fraction_center_mean','center_nCount','model_gene_universe_hash']]);rng=np.random.default_rng(20261013);epids={rep:rng.choice(ep,5,replace=True) for rep in range(5)};rows=[];states=[];members=[];fail=[];state_used=0
 for sec,g in blue.groupby('section_id',sort=True):
  cy=[];ce=[];sy=[];se=[];cid=[];sid=[]
  model=ROOT/'cell2location_full_run/v2_30type/02_sections'/sec;a=ad.read_h5ad(model/'spatial_posterior.h5ad',backed='r');mg=a.var_names.tolist();model_depth=pd.Series(np.asarray(sparse.csr_matrix(a.X[:]).sum(1)).ravel(),index=a.obs_names);a.file.close();gi=pd.Index(genes).get_indexer(mg)
  if (gi<0).any():fail.append(dict(section_id=sec,reason='single_cell_missing_model_genes',genes=','.join(np.array(mg)[gi<0])));continue
  xs=X[:,gi].tocsr();ab=pd.read_csv(model/'abundance_mean.csv.gz',index_col=0);types=list(ab.columns);pools={t:sc.index[sc.type30.eq(t)].to_numpy() for t in types};assert all(len(v)>0 for v in pools.values())
  for row in g.itertuples():
   prop=ab.loc[row.barcode].to_numpy(dtype=float,copy=True);prop/=prop.sum();ei=types.index('Epithelial');target=int(round(model_depth.loc[row.barcode]));assert target>0
   for rep in range(5):
    profiles=[];chosen=[]
    for t in types:
     ids=epids[rep] if t=='Epithelial' else rng.choice(pools[t],3,replace=True);profiles.append(np.asarray(xs[ids].mean(0)).ravel());chosen.append(ids)
    profiles=np.array(profiles);contrib=profiles*prop[:,None];ev=contrib[ei];nv=contrib.sum(0)-ev;total=ev.sum()+nv.sum();joint=rng.multinomial(target,np.r_[ev,nv]/total);E=joint[:len(mg)];N=joint[len(mg):];Y=E+N;bid=f'{row.spot_uid}::r{rep}';scor=auc(Y,mg);es=auc(E,mg)
    cy.append(sparse.csr_matrix(Y[None,:]));ce.append(sparse.csr_matrix(E[None,:]));cid.append(bid)
    rows.append(dict(blueprint_uid=row.spot_uid,bag_id=bid,dataset_id=row.dataset_id,section_id=sec,niche_name=row.niche_name,rep=rep,target_UMI=target,raw_full_universe_UMI=row.center_nCount,depth_definition="observed counts sum in matched model gene universe",attained_UMI=int(Y.sum()),epi_cell_fraction=prop[ei],expected_epi_RNA_fraction=ev.sum()/total,realized_epi_RNA_fraction=E.sum()/target,observed_score=scor,true_epi_score=es,signature_genes=len(TOP),model_gene_universe_hash=row.model_gene_universe_hash,seed=20261013))
    for t,ids,pr in zip(types,chosen,prop):
     members.append(dict(bag_id=bid,celltype=t,cell_fraction=pr,donor_cell_ids=';'.join(sc.loc[ids,'cell_id']),donor_reference_batches=';'.join(sc.loc[ids,'reference_batch'].astype(str)),mean_model_UMI=float(profiles[types.index(t)].sum())))
    if rep==0 and row.spot_uid in state_allow:
     # Exact same non-Epi realized counts and total RNA, Epi state profile replacement only.
     for state,pool in [('low',lo),('high',hi)]:
      ids=rng.choice(pool,5,replace=True);v=np.asarray(xs[ids].sum(0)).ravel();ee=rng.multinomial(int(E.sum()),v/v.sum());yy=ee+N;sy.append(sparse.csr_matrix(yy[None,:]));se.append(sparse.csr_matrix(ee[None,:]));sid.append(row.spot_uid+'::'+state);states.append(dict(blueprint_uid=row.spot_uid,dataset_id=row.dataset_id,section_id=sec,niche_name=row.niche_name,state=state,state_variable=statevar,epi_UMI=int(ee.sum()),total_UMI=int(yy.sum()),observed_score=auc(yy,mg),true_epi_score=auc(ee,mg),donor_cell_ids=';'.join(sc.loc[ids,'cell_id']),model_gene_universe_hash=row.model_gene_universe_hash,scope='paired technical state control; fixed non-Epi counts and total depth, not end-to-end attribution validation'))
     state_used+=1
  cache=P/'counts_cache'/sec;cache.mkdir(parents=True,exist_ok=True)
  for nm,arr in [('composition_observed',cy),('composition_true_epi',ce),('state_observed',sy),('state_true_epi',se)]:
   if arr:sparse.save_npz(cache/(nm+'.npz'),sparse.vstack(arr,format='csr'))
  js(cache/'axis_manifest.json',dict(genes=mg,composition_ids=cid,state_ids=sid))
  log(f'source {sec} cumulative bags={len(rows)}, state_blueprints={state_used}')
 tab(P/'composition_only_pseudospots.tsv',rows);tab(P/'composition_donors_and_cellsize.tsv.gz',members);tab(P/'state_change_controls.tsv',states);tab(P/'simulation_unavailable_sections.tsv',pd.DataFrame(fail,columns=['section_id','reason','genes']))
 r=pd.DataFrame(rows);tab(P/'composition_niche_summary.tsv',r.groupby(['dataset_id','niche_name']).agg(n_bags=('observed_score','size'),n_blueprints=('blueprint_uid','nunique'),median_observed=('observed_score','median'),median_true_epi=('true_epi_score','median')).reset_index())
 s=pd.DataFrame(states);paired=s.pivot(index=['blueprint_uid','dataset_id','section_id','niche_name'],columns='state',values=['observed_score','true_epi_score']);paired.columns=['_'.join(x) for x in paired.columns];paired=paired.reset_index();paired['observed_high_minus_low']=paired.observed_score_high-paired.observed_score_low;paired['true_epi_high_minus_low']=paired.true_epi_score_high-paired.true_epi_score_low;tab(P/'state_control_paired_effects.tsv',paired)
 js(P/'source_validation_manifest.json',dict(n_sc=len(sc),n_types=sc.type30.nunique(),n_blueprints=blue.spot_uid.nunique(),n_bags=len(rows),n_state_blueprints=state_used,n_usage_matched=len(us),state_variable=statevar,seed=20261013,composition='central posterior mean cell abundance, raw model-gene single-cell library-size-weighted counts; exact actual depth',epi_distribution='same five Epi donor sets per repetition across all niches',independence='simulation blueprints/repeats not independent patients; sc not paired to Visium',oracle='not performed; true Epi counts retained for technical scoring check only'))
 log('source controls complete')
if __name__=='__main__':main()
