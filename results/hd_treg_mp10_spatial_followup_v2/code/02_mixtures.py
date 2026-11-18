from fu_common import *
import anndata as ad
F=yaml.safe_load((R/'frozen_config.yaml').read_text());mixcfg=F['mixture']
if (R/'tables/expanded_mixtures.parquet').exists() and os.environ.get('FORCE_MIX')!='1':print('MIXTURE_CACHE');raise SystemExit
previous=pd.read_parquet(R/'tables/expanded_mixtures.parquet',columns=['Treg_main','nFeature'])if (R/'tables/expanded_mixtures.parquet').exists()else None
origin_rng=np.random.default_rng(1729)
def realized(v,prob,pT,pReg):
 # Conditional source attribution preserves the already drawn gene-count matrix exactly.
 t=origin_rng.binomial(v,np.clip(np.divide(pT,prob,out=np.zeros_like(prob),where=prob>0),0,1));reg=origin_rng.binomial(t,np.clip(np.divide(pReg,pT,out=np.zeros_like(prob),where=pT>0),0,1));nt=int(t.sum());nr=int(reg.sum());return dict(realized_T_RNA_UMI=nt,realized_Treg_RNA_UMI=nr,realized_Treg_total_RNA_fraction=nr/int(v.sum()),realized_Treg_T_RNA_fraction=nr/nt if nt else np.nan,truth_definition='component origin sampled conditional on aggregate counts and known mixture profiles')
m=pd.read_parquet(V/'tables/reference_membership.parquet');tr=m.partition.eq('train').to_numpy();assert set(m.loc[tr,'cell_id']).isdisjoint(m.loc[~tr,'cell_id']);assert set(m.loc[tr,'sample']).isdisjoint(m.loc[~tr,'sample'])
p=Path(json.loads((V/'models/reference_prepared.json').read_text())['reference']);a=ad.read_h5ad(p,backed='r');ii=a.obs_names.get_indexer(m.cell_id);assert (ii>=0).all();x=sparse.csr_matrix(a.X[ii,:]);sy=a.var_names.to_numpy();a.file.close();assert len(sy)==17987
selgenes=sorted(set(SIG+CORE+TL));sel=np.array([np.flatnonzero(sy==g)[0]for g in selgenes]);put(m[['cell_id','sample','label','partition','major_celltype']],'reference_membership.tsv.gz')

def scored(mx,meta):
 r=v1.selected_ranks(mx.indptr,mx.indices,mx.data,mx.shape[1],sel);d=pd.DataFrame(meta);raw=mx[:,sel].toarray();ranks=pd.DataFrame(r,columns=selgenes)
 for j,g in enumerate(selgenes):d['raw_'+g]=raw[:,j];d['rank_'+g]=r[:,j]
 d['nFeature']=np.diff(mx.indptr);d['nCount']=np.asarray(mx.sum(1)).ravel();d['T_lineage_detected']=(raw[:,[selgenes.index(g)for g in TL]]>0).sum(1);d['T_lineage_UMI']=raw[:,[selgenes.index(g)for g in TL]].sum(1);d['main_detected']=(raw[:,[selgenes.index(g)for g in SIG]]>0).sum(1);d['main_UMI']=raw[:,[selgenes.index(g)for g in SIG]].sum(1);d['core_detected']=(raw[:,[selgenes.index(g)for g in CORE]]>0).sum(1);d['Treg_main']=v1.ucell(ranks,SIG,1300);d['T_lineage']=v1.ucell(ranks,TL,1300)
 d['T_eligible']=(d.nFeature>=325)&(d.T_lineage_detected>=2)&(d.T_lineage_UMI>=3);d['multigene']=d.main_detected>=3;d['core']=d.core_detected>=1;d['high_old']=d.Treg_main>F['old_threshold'];d['joint_old']=d.T_eligible&d.multigene&d.core&d.high_old;d['score_zero']=d.Treg_main==0;d['zero_reason']=np.where(~d.score_zero,'nonzero',np.where(d.main_detected==0,'no_main_gene_detected','detected_genes_all_average_rank_at_or_beyond1300'));d['within_rank_main']=(ranks[SIG]<1300).sum(1).to_numpy();return d
# Exact replay of v1 RNG consumption, including prior thinning, then original heldout70% mixtures.
rng=np.random.default_rng(123);vv=np.flatnonzero(~tr)
for dep in [500,2000,5000]:
 z=x[vv].copy();tot=np.asarray(z.sum(1)).ravel();prob=np.minimum(1,dep/np.maximum(tot,1));rng.binomial(z.data.astype(int),np.repeat(prob,np.diff(z.indptr)))
pools={l:np.flatnonzero((~tr)&m.label.eq(l).to_numpy())for l in m.label.unique()};rows=[];xs=[]
for comp in mixcfg['competitors']:
 for frac in [0.,.1,.3]:
  for rep in range(100):
   profiles=[]
   for lab in ['Epithelial',comp,'Treg']:
    z=np.asarray(x[rng.choice(pools[lab],min(5,len(pools[lab])),False)].sum(0)).ravel();profiles.append(z/z.sum())
   prob=.7*profiles[0]+(.3-frac)*profiles[1]+frac*profiles[2];prob/=prob.sum()
   for depth in [500,2000,5000]:
    v=rng.multinomial(depth,prob);xs.append(sparse.csr_matrix(v));rows.append(dict(partition='sample_holdout',background='Epi70_original',competitor=comp,Treg_RNA_fraction=frac,Treg_T_RNA_fraction=frac/.3,depth=depth,rep=rep,**realized(v,prob,(.3-frac)*profiles[1]+frac*profiles[2],frac*profiles[2])))
mx=sparse.vstack(xs,format='csr');d=scored(mx,rows);prior=pd.read_parquet(V/'tables/reference_mixture_ranks.parquet');assert np.array_equal(d.nFeature,prior.nFeature)
for g in selgenes:assert np.array_equal(d['rank_'+g],prior[g]),g
d.to_parquet(R/'tables/original70_replayed.parquet',index=False);put(pd.DataFrame([dict(check='original2700_mixtures_exact_full_rank_replay',passed=True,n=len(d),genes=len(selgenes))]),'mixture_replay_validation.tsv');del mx,xs
# Separate partition pools; parent IDs persisted. No derived sample can move partitions.
meta=[];xs=[];parents=[]
for part,mask in [('train',tr),('sample_holdout',~tr)]:
 rng=np.random.default_rng(123 if part=='train'else 456)
 pool={l:np.flatnonzero(mask&m.label.eq(l).to_numpy())for l in m.label.unique()}
 pool['Fibroblast']=np.flatnonzero(mask&m.label.isin(['Matrix_CAF','iCAF','rCAF','vCAF']).to_numpy());pool['Myeloid']=np.flatnonzero(mask&m.label.isin(['C1QC TAM','SPP1 TAM','Classical Mono']).to_numpy())
 for bg in mixcfg['backgrounds']:
  for comp in mixcfg['competitors']:
   for rep in range(mixcfg['replicates']):
    pid=f'{part}_{bg["name"]}_{comp}_{rep}';profiles={};labs=['Epithelial',comp,'Treg']+([bg['other']]if bg['other'] else [])
    for lab in labs:
     inds=rng.choice(pool[lab],min(5,len(pool[lab])),False);z=np.asarray(x[inds].sum(0)).ravel();profiles[lab]=z/z.sum();parents.extend([dict(profile_id=pid,partition=part,role=lab,cell_id=m.cell_id.iloc[j],sample=m['sample'].iloc[j])for j in inds])
    Tmass=1-bg['Epi']-bg['other_fraction']
    for tfrac in mixcfg['Treg_T_RNA_fractions']:
     frac=Tmass*tfrac;prob=bg['Epi']*profiles['Epithelial']+(Tmass-frac)*profiles[comp]+frac*profiles['Treg']
     if bg['other']:prob+=bg['other_fraction']*profiles[bg['other']]
     prob/=prob.sum()
     for dep in mixcfg['depths']:
      z=rng.multinomial(dep,prob);xs.append(sparse.csr_matrix(z));meta.append(dict(profile_id=pid,partition=part,background=bg['name'],Epi_RNA_fraction=bg['Epi'],other_RNA_fraction=bg['other_fraction'],competitor=comp,Treg_RNA_fraction=frac,Treg_T_RNA_fraction=tfrac,T_RNA_fraction=Tmass,depth=dep,rep=rep,**realized(z,prob,(Tmass-frac)*profiles[comp]+frac*profiles['Treg'],frac*profiles['Treg'])))
 print('MIXTURE_PARTITION',part,len(meta),flush=True)
mx=sparse.vstack(xs,format='csr');d=scored(mx,meta)
if previous is not None:assert np.array_equal(previous[['Treg_main','nFeature']].to_numpy(),d[['Treg_main','nFeature']].to_numpy()),'Source attribution must not change scores/counts'
d.to_parquet(R/'tables/expanded_mixtures.parquet',index=False);pd.DataFrame(parents).to_parquet(R/'tables/mixture_parent_membership.parquet',index=False)
put(pd.DataFrame([dict(partition=part,n_reference_cells=int(mask.sum()),n_reference_samples=m.loc[mask,'sample'].nunique(),n_mixtures=int(d.partition.eq(part).sum()),independence='sample split; simulated mixtures conditional on finite sample profiles, not patients')for part,mask in [('train',tr),('sample_holdout',~tr)]]),'mixture_partitions.tsv');print('MIXTURE_DONE',len(d),flush=True)
