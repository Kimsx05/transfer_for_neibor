from common import *
import anndata as ad

def main(out):
 P=out/'03_attribution_calibration';rng=np.random.default_rng(123);a,C,s=params('CanDis_2026_Sample1');genes=a.var_names.tolist();types=list(C.columns);assert types[8]=='Epithelial';ref=C.to_numpy(float);libs=ref.sum(0);write_tsv(out,P/'frozen_reference.tsv',C.reset_index());X,eg,em=epi_cache();ix=pd.Index(eg).get_indexer(genes);assert (ix>=0).all();X=X[:,ix].tocsr();fm=pd.read_parquet(out/'01_score_fidelity/cell_level_fidelity.parquet');assert fm.cell_id.tolist()==em.cell_id.tolist();em=fm;old=OLD/'07_source_validation';nm=pd.read_csv(old/'single_cell_selected.tsv.gz',sep='\t');ng=pd.read_csv(old/'selected_genes.tsv',sep='\t').gene.tolist();NX=sparse.load_npz(old/'selected_counts.npz')[:,pd.Index(ng).get_indexer(genes)].tocsr();mp=pd.read_csv(ROOT/'subtype/v2/reference_training_v2/11_final_spatial_reference/00_32_to_30_label_mapping.csv').set_index('reference_celltype_v2').spatial_reference_celltype_v2_final;nm['type30']=nm.reference_celltype_v2.map(mp);profiles=np.stack([np.asarray(NX[nm.type30.eq(t).to_numpy()].mean(0)).ravel() for t in types]);trainlib=profiles.sum(1);assert (trainlib>0).all();scale=trainlib/libs
 # Donors selected only by source coverage and within-source usage availability.
 eligible=[]
 for (dataset,unit),g in em.groupby(['dataset','analysis_unit']):
  q=g.MP10_usage.dropna()
  if len(q)>=100 and q.nunique()>1:eligible.append(dict(dataset=dataset,unit=unit,n=len(q)))
 ed=pd.DataFrame(eligible).sort_values(['dataset','unit']);sel=ed.groupby('dataset').head(1).head(6);write_tsv(out,P/'selected_epi_donors.tsv',sel)
 d=pd.read_parquet(V4/'04_niche_comparison/anchor_MP10_all_layers.parquet');central=[];depth=[]
 for sec,g in d[d.accepted_label.ge(0)&d.epi_score_valid].groupby('section_id'):
  ab=pd.read_csv(ROOT/'cell2location_full_run/v2_30type/02_sections'/sec/'abundance_mean.csv.gz',index_col=0).loc[g.barcode];v=ab.to_numpy(float);v=v/v.sum(1,keepdims=True);central.append(v);depth.extend(g.center_nCount)
 central=np.vstack(central);ne=central.copy();ne[:,8]=0;ne/=ne.sum(1,keepdims=True);base=np.median(ne,axis=0);base/=base.sum();depths=np.quantile(depth,[.25,.5,.75]).round().astype(int);targets=['GC B','SPP1 TAM','Neutrophil'];levels={t:np.quantile(ne[:,types.index(t)][ne[:,types.index(t)]>0],[.1,.9]).tolist() for t in targets};rows=[];Ys=[];Es=[];Ws=[];donors=[]
 def add(meta,Y,E,w):
  assert (E<=Y).all() and np.isfinite(Y).all();bid=f'v5_{len(rows):05d}';rows.append(dict(bag_id=bid,**meta,observed_UMI=float(Y.sum()),epi_UMI=float(E.sum()),cell_equivalent_total=float((w/scale).sum()),true_RNA_fraction=float(E.sum()/Y.sum()),signature_weight_sum=float(w.sum())));Ys.append(sparse.csr_matrix(Y[None,:]));Es.append(sparse.csr_matrix(E[None,:]));Ws.append(w)
 def nonprofile(t,level):
  j=types.index(t);f=base.copy();f[j]=0;f/=f.sum();f*=1-level;f[j]=level;rna=f*trainlib;rna/=rna.sum();return rna
 def draw_non(rna,n):
  quotas=rng.multinomial(n,rna);N=np.zeros(len(genes),int)
  for j in range(30):
   if quotas[j]:N+=rng.multinomial(int(quotas[j]),profiles[j]/trainlib[j])
  return N,quotas
 for donor in sel.itertuples():
  g=em[em.analysis_unit.eq(donor.unit)&em.dataset.eq(donor.dataset)];q=g.MP10_usage;low=g.index[q<=q.quantile(.25)].to_numpy();high=g.index[q>=q.quantile(.75)].to_numpy();allids=g.index.to_numpy()
  for rep in range(2):
   ids={state:rng.choice(pool,min(10,len(pool)),replace=False) for state,pool in [('fixed',allids),('low',low),('high',high)]};ep={state:np.asarray(X[ii].sum(0)).ravel() for state,ii in ids.items()};ep={k:v/v.sum() for k,v in ep.items()}
   for state,ii in ids.items():donors.append(dict(donor=donor.unit,rep=rep,state=state,cell_ids=';'.join(em.loc[ii,'cell_id']),reference_overlap='IN_REFERENCE_OR_UNVERIFIED',MP_basis_overlap='IN_BASIS_OR_UNVERIFIED'))
   for di,D in enumerate(depths):
    Efix=rng.multinomial(int(D//2),ep['fixed']);group=f'{donor.unit}::d{di}::r{rep}'
    for target in targets:
     for li,level in enumerate(levels[target]):
      rna=nonprofile(target,level);N,quotas=draw_non(rna,int(D-Efix.sum()));quota=quotas.astype(float);quota[8]=Efix.sum();w=quota/libs;add(dict(donor=donor.unit,dataset=donor.dataset,rep=rep,depth_band=di,background=target,target_fraction=level,condition='FIXED_STATE_COMPOSITION',state='fixed',pair_id=group+'::'+target+'::'+str(li),fixed_epi_group=group,platform='m_g=1; ambient=0'),N+Efix,Efix,w)
     # Fixed non-Epi and Epi RNA budgets for state pair; state never estimator input.
     rna=nonprofile(target,levels[target][1]);N,quota=draw_non(rna,int(D-D//2));w=quota.astype(float)/libs;w[8]=(D//2)/libs[8]
     for state in ['low','high']:
      E=rng.multinomial(int(D//2),ep[state]);add(dict(donor=donor.unit,dataset=donor.dataset,rep=rep,depth_band=di,background=target,target_fraction=levels[target][1],condition='FIXED_COMPOSITION_STATE',state=state,pair_id=group+'::'+target,fixed_epi_group='',platform='m_g=1; ambient=0'),N+E,E,w)
     for purity in [.25,.7]:
      nE=int(round(D*purity));E=rng.multinomial(nE,ep['fixed']);N,quota=draw_non(rna,int(D-nE));w=quota.astype(float)/libs;w[8]=nE/libs[8];add(dict(donor=donor.unit,dataset=donor.dataset,rep=rep,depth_band=di,background=target,target_fraction=levels[target][1],condition='PURITY_GRADIENT',state='fixed',pair_id=group+'::'+target+'::p'+str(purity),fixed_epi_group='',platform='m_g=1; ambient=0'),N+E,E,w)
 log('real donor panel constructed')
 # Separate model-matched positive control with explicit platform and ambient.
 baseline=[]
 for i in range(24):
  w=rng.dirichlet(np.ones(30))*20;mg=np.asarray(s['m_g']).ravel();ambient=np.asarray(s['s_g_gene_add'])[0];mu=w@ref.T*mg+ambient;muE=w[8]*ref[:,8]*mg;Y=rng.poisson(mu);E=rng.binomial(Y,np.divide(muE,mu,out=np.zeros_like(mu),where=mu>0));eh=allocate(Y[None,:],ref,w[None,:],mg,ambient[None,:])[0];ts,_,_=auc_parts(E[None,:],genes);hs,_,_=auc_parts(eh[None,:],genes);baseline.append(dict(bag_id=i,true_score=ts[0],estimated_score=hs[0],bias=hs[0]-ts[0],scope='MODEL_MATCHED_NUMERICAL_CONTROL_NOT_DONOR_GENERALIZATION'))
 write_tsv(out,P/'model_matched_baseline.tsv',baseline)
 Y=sparse.vstack(Ys,format='csr');E=sparse.vstack(Es,format='csr');W=np.array(Ws);md=pd.DataFrame(rows);assert len(md)<=3000;safe(out,P/'panel');sparse.save_npz(P/'panel_observed.npz',Y);sparse.save_npz(P/'panel_true_epi.npz',E);np.save(safe(out,P/'panel_known_weights.npy'),W);md.to_parquet(safe(out,P/'panel_metadata.parquet'),index=False);write_json(out,P/'panel_axes.json',dict(genes=genes,types=types,bag_ids=md.bag_id.tolist(),reference_section='CanDis_2026_Sample1'));write_tsv(out,P/'simulation_donors.tsv.gz',donors)
 # Verify conditions by original axes, no inner join.
 checks=[]
 for group,g in md[md.condition.eq('FIXED_STATE_COMPOSITION')].groupby('fixed_epi_group'):
  ii=g.index.to_numpy();assert all((E[i]!=E[ii[0]]).nnz==0 for i in ii);checks.append(dict(group=group,check='IDENTICAL_REALIZED_EPI',status='PASS'))
 for pair,g in md[md.condition.eq('FIXED_COMPOSITION_STATE')].groupby('pair_id'):
  ii=g.index.to_numpy();assert len(ii)==2;assert ((Y[ii[0]]-E[ii[0]])!=(Y[ii[1]]-E[ii[1]])).nnz==0;assert E[ii[0]].sum()==E[ii[1]].sum();checks.append(dict(group=pair,check='IDENTICAL_NON_EPI_AND_TOTAL_DEPTH',status='PASS'))
 write_tsv(out,P/'simulation_axis_and_scale_audit.tsv',checks);write_tsv(out,P/'donor_reference_overlap.tsv',[dict(donor=r.unit,dataset=r.dataset,reference_status='IN_REFERENCE_OR_UNVERIFIED',basis_status='IN_BASIS_OR_UNVERIFIED',held_out=False) for r in sel.itertuples()]);write_json(out,P/'calibration_design.json',dict(seed=123,n_panel=len(md),n_donors=len(sel),depths=depths.tolist(),target_levels=levels,platform='real donor panel m_g1 ambient0; baseline original mg and explicit ambient',weights='known expected RNA allocation / fixed reference type library; cell-equivalent numbers use training source mean library, not realized E gene vector',reference_scope='original overlapping reference; no independent donor claim',fit_selection='rep0 only, all donors/backgrounds/depth/state conditions',near_zero='one original normalized rank step; noise summarized across matched technical repetitions'))
 # Known-mixture estimator never reads E for constructing weights.
 score=[]
 for st in range(0,len(md),64):
  yy=Y[st:st+64].toarray();ee=E[st:st+64].toarray();eh=allocate(yy,ref,W[st:st+64],np.ones(len(genes)),np.zeros((1,len(genes))));ts,tc,_=auc_parts(ee,genes);hs,hc,_=auc_parts(eh,genes)
  for j in range(len(yy)):score.append(dict(bag_id=md.bag_id.iloc[st+j],true_score=ts[j],estimated_score=hs[j],score_error=hs[j]-ts[j],top20_count_MAE=float(np.abs(eh[j,pd.Index(genes).get_indexer(TOP)]-ee[j,pd.Index(genes).get_indexer(TOP)]).mean()),branch='KNOWN_MIXTURE_FIXED_SIGNATURE'))
 md.merge(pd.DataFrame(score),on='bag_id',validate='one_to_one').to_parquet(safe(out,P/'known_mixture_scores.parquet'),index=False)
 fitids=np.flatnonzero(md.rep.eq(0));aa=ad.AnnData(Y[fitids].astype('float32'));aa.var_names=genes;aa.obs_names=md.bag_id.iloc[fitids].to_numpy();aa.obs['technical_batch']='synthetic';aa.write_h5ad(safe(out,P/'synthetic_refits/input.h5ad'));log(f'panel complete {len(md)} bags, {len(fitids)} synthetic fit bags')
if __name__=='__main__':main(OUT)
