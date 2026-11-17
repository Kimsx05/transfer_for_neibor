from common import *
import anndata as ad
from sklearn.metrics import roc_auc_score
from scipy.io import mmwrite
if __name__=='__main__':
 gate()
 if (R/'models/reference_prepared.json').exists():raise SystemExit
 m=pd.read_csv(OLD/'tables/reference_membership.tsv.gz',sep='\t');p=P/'subtype/v2/reference_training_v2/11_final_spatial_reference/reference_v2_final_30type_regression_input.h5ad';a=ad.read_h5ad(p,backed='r');ii=a.obs_names.get_indexer(m.cell_id);assert (ii>=0).all();x=sparse.csr_matrix(a.X[ii,:]);sy=a.var_names.to_numpy();a.file.close();assert np.all(x.data==np.round(x.data));m.to_parquet(R/'tables/reference_membership.parquet',index=False)
 genes=sorted(set(sum(PANELS.values(),[])+literature_genes()+COMM+TOP));avail=[g for g in genes if g in sy];sel=np.array([np.where(sy==g)[0][0]for g in avail]);raw=x[:,sel].toarray();norm=raw*10000/np.maximum(np.asarray(x.sum(1)).ravel()[:,None],1);tr=m.partition.eq('train').to_numpy();treg=m.label.eq('Treg').to_numpy();comp=m.label.str.startswith(('CD4','CD8','Tfh','NK')).to_numpy();rows=[]
 for g in PANELS['Treg_candidates']:
  if g not in avail:continue
  j=avail.index(g);ma=norm[tr&treg,j].mean();compet=[norm[tr&m.label.eq(l).to_numpy(),j].mean()for l in sorted(m.label[comp].unique())];non=norm[tr&~comp&~treg,j].mean();det=(raw[tr&treg,j]>0).mean();maxcomp=max(compet)
  rows.append(dict(gene=g,Treg_detected=det,Treg_mean_per10k=ma,max_competing_T_NK_mean=maxcomp,nonT_mean=non,Treg_vs_max_competitor_ratio=(ma+.01)/(maxcomp+.01),Treg_vs_nonT_ratio=(ma+.01)/(non+.01),selected=bool(det>=.08 and ma>1.2*maxcomp and ma>1.5*non and g not in TOP+COMM)))
 met=pd.DataFrame(rows);put(met,'signature_training_gene_selection.tsv');chosen=met.loc[met.selected,'gene'].tolist();assert len(chosen)>=2,'Cannot establish multigene reference signature'
 excluded=set(TOP+COMM);lit=[g for g in literature_genes()if g not in excluded and g in avail]
 sig={'T_lineage':PANELS['T'],'Treg_main':chosen,'Treg_literature':lit,'Tfh':PANELS['Tfh'],'Activated_T':PANELS['Activated_T']}
 sig={k:[g for g in v if g in avail]for k,v in sig.items()}
 (R/'models/signatures.json').write_text(json.dumps(sig,indent=2));put(pd.DataFrame([dict(signature=k,gene=g,source='reference training candidate selection' if k=='Treg_main' else 'DeSimone Table S4 modified' if k=='Treg_literature' else 'prespecified competitor/lineage panel')for k,v in sig.items()for g in v]),'signature_genes.tsv')
 put(pd.DataFrame([dict(gene=g,reference_present=g in sy,excluded_readout=g in COMM,excluded_MP10=g in TOP,used=g in lit,source='DeSimone full Table S4 309 genes')for g in literature_genes()]),'literature_signature_modifications.tsv')
 classrows=[]
 for lab in sorted(m.label.unique()):
  k=tr&m.label.eq(lab).to_numpy()
  for g in PANELS['Treg_candidates']+PANELS['T']:
   if g in avail:classrows.append(dict(label=lab,gene=g,n=int(k.sum()),mean_per10k=norm[k,avail.index(g)].mean(),detected=(raw[k,avail.index(g)]>0).mean()))
 put(pd.DataFrame(classrows),'signature_reference_expression.tsv')
 # Store only selected full-background ranks, not the reference expression matrix.
 rk=selected_ranks(x.indptr,x.indices,x.data,x.shape[1],sel);rf=pd.DataFrame(rk,columns=avail);rf=pd.concat([m.reset_index(drop=True),rf],axis=1);rf['nFeature']=np.diff(x.indptr);rf['nCount']=np.asarray(x.sum(1)).ravel();rf.to_parquet(R/'tables/reference_full_ranks.parquet',index=False)
 # Held-out samples, expected count thinning; no new selection on validation data.
 vv=np.flatnonzero(~tr);rng=np.random.default_rng(123)
 for dep in [500,2000,5000]:
  z=x[vv].copy();tot=np.asarray(z.sum(1)).ravel();prob=np.minimum(1,dep/np.maximum(tot,1));z.data=rng.binomial(z.data.astype(int),np.repeat(prob,np.diff(z.indptr)));z.eliminate_zeros();rr=selected_ranks(z.indptr,z.indices,z.data,z.shape[1],sel);v=pd.DataFrame(rr,columns=avail);v['label']=m.iloc[vv].label.to_numpy();v['sample']=m.iloc[vv]['sample'].to_numpy();v['cell_id']=m.iloc[vv].cell_id.to_numpy();v['expected_depth']=dep;v['actual_UMI']=np.asarray(z.sum(1)).ravel();v['nFeature']=np.diff(z.indptr);v.to_parquet(R/'tables'/f'reference_holdout_{dep}.parquet',index=False)
 # Synthetic RNA mixtures from held-out profiles. Fractions are RNA, not cell proportions.
 pools={l:np.flatnonzero((~tr)&m.label.eq(l).to_numpy())for l in m.label.unique()};mix=[];labels=[]
 for comp_label in ['CD4 Helper/Conventional','Tfh/CXCL13-like','CD8 Exhausted']:
  for frac in [0.,.1,.3]:
   for rep in range(100):
    profiles=[]
    for l in ['Epithelial',comp_label,'Treg']:
     q=np.asarray(x[rng.choice(pools[l],min(5,len(pools[l])),False)].sum(0)).ravel();profiles.append(q/max(q.sum(),1))
    prob=.7*profiles[0]+(.3-frac)*profiles[1]+frac*profiles[2];prob/=prob.sum()
    for depth in [500,2000,5000]:
     v=rng.multinomial(depth,prob);mix.append(sparse.csr_matrix(v));labels.append(dict(competitor=comp_label,Treg_RNA_fraction=frac,depth=depth,rep=rep))
 mx=sparse.vstack(mix,format='csr');rr=selected_ranks(mx.indptr,mx.indices,mx.data,mx.shape[1],sel);md=pd.concat([pd.DataFrame(labels),pd.DataFrame(rr,columns=avail)],axis=1);md['nFeature']=np.diff(mx.indptr);md.to_parquet(R/'tables/reference_mixture_ranks.parquet',index=False)
 (R/'models/reference_prepared.json').write_text(json.dumps({'reference':str(p),'cells':len(m),'train':int(tr.sum()),'heldout':int((~tr).sum()),'validation_unit':'orig.ident sample, not patient','main_genes':chosen,'literature_original_n':309,'literature_used_n':len(lit),'selection_rule':'Treg detection>=.08; mean normalized >1.2*maximum competing T/NK class and >1.5*pooled nonT; candidates only; no HD associations'},indent=2));print('REFERENCE',chosen,'LITERATURE',len(lit),flush=True)
