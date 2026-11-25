from common import *
import networkx as nx, warnings, gc
import anndata as ad
import celloracle as co
from sklearn.decomposition import PCA
import random
PANEL=rd('tables/panel_members.tsv');META=rd('cache/metadata.tsv.gz')
X=sp.load_npz('cache/panel_counts.npz');base=rd('cache/base_edges.tsv.gz');settings=C['celloracle']
controls=rd('tables/control_matching.tsv');kos=K[:2]+sorted(set(controls[controls.ko.isin(K[:2])].control))
def analyze(job,state):
 seed=456 if job=='seed456' else 123
 rng=np.random.RandomState(seed);np.random.seed(seed);random.seed(seed)
 out=Path('cache/oracle')/job/state;out.mkdir(exist_ok=True,parents=True)
 sig=stamp([Path('code/oracle.py'),Path('cache/base_edges.tsv.gz'),Path('tables/control_matching.tsv'),Path('config/main.json'),Path('config/input_fingerprints.tsv'),Path('tables/panel_members.tsv')]);marker=out/'complete.json'
 if marker.exists():
  j=json.loads(marker.read_text())
  if j.get('signature')==sig and all((out/n).is_file() and sha(out/n)==v for n,v in j['outputs'].items()):return
 def complete(status,why=''):
  files=[f for f in out.glob('*.gz')]+[f for f in out.glob('*.npz')]
  j=dict(status=status,reason=why,cells=len(m),patients=m.patient_id.nunique(),signature=sig,outputs={f.name:sha(f) for f in files})
  q=marker.with_suffix('.new');q.write_text(json.dumps(j));os.replace(q,marker)
 m=META[META.state==state]
 if job.startswith('drop_'):m=m[m.dataset!=job[5:]]
 m=m.iloc[rng.permutation(len(m))]
 if job=='pilot':m=m.head(400)
 if len(m)<80 or m.patient_id.nunique()<4:complete('NE','<80 cells or <4 patients');return
 x=X[m.index];m=m.reset_index(drop=True);aa=x.getnnz(0);b=sum((x[ix].getnnz(0)>0).astype(int) for ix in m.groupby('patient_id').indices.values())
 keep=PANEL.P_clean & ~PANEL.gene.str.startswith('RAW~') & (aa>=10)&(b>=3)&((np.ravel(x.astype(float).power(2).mean(0))-np.ravel(x.mean(0))**2)>0)
 genes=pd.Index(PANEL.gene[keep]);x=x[:,keep.values].copy()
 local=base[base.source.isin(genes)&base.target.isin(genes)]
 tf=local.groupby('target').source.agg(list).to_dict()
 save(pd.DataFrame(dict(gene=PANEL.gene,positive_cells=aa,positive_patients=b,in_oracle=keep)) ,out/'coverage.tsv.gz')
 a=ad.AnnData(X=x.multiply(1e4/m.nCount_RNA.to_numpy()[:,None]).tocsr().astype(np.float64),obs=m.set_index('barcode').copy(),var=pd.DataFrame(index=genes))
 a.X.data=np.log1p(a.X.data);a.layers['raw_count']=x
 a.obs['state']=pd.Categorical(a.obs.state);a.uns['state_colors']=['#195571'];a.obsm['X_pca']=PCA(n_components=2,random_state=seed,svd_solver='randomized').fit_transform(a.X.toarray())
 o=co.Oracle();o.import_anndata_as_normalized_count(a,cluster_column_name='state',embedding_name='X_pca');o.import_TF_data(TFdict=tf)
 o.perform_PCA(n_components=min(settings['PCs'],len(m)-1,len(genes)-1));o.knn_imputation(k=min(settings['neighbors'],len(m)-1),n_pca_dims=o.pcs.shape[1],balanced=False,n_jobs=1)
 links=o.get_links(cluster_name_for_GRN_unit='state',alpha=settings['alpha'],bagging_number=settings['bagging_number'],model_method=settings['model'],n_jobs=1,verbose_level=0)
 links.filter_links(p=settings['p_cutoff'],weight='coef_abs',threshold_number=settings['top_links_per_cluster'])
 o.get_cluster_specific_TFdict_from_Links(links);o.fit_GRN_for_simulation(alpha=settings['alpha'],GRN_unit='cluster',use_cluster_specific_TFdict=True,verbose_level=0)
 coef=o.coef_matrix_per_cluster[state];ii,jj=np.nonzero(coef.values)
 ed=pd.DataFrame(dict(source=coef.index[ii],target=coef.columns[jj],coefficient=coef.values[ii,jj]))
 filt=links.filtered_links[state];assert set(zip(ed.source,ed.target))<=set(zip(filt.source,filt.target))<=set(zip(local.source,local.target))
 save(ed.merge(filt,how='left',on=['source','target'] ),out/'edges.tsv.gz')
 wt=np.asarray(o.adata.layers['imputed_count']);sd=wt.std(0)
 o.simulate_shift(perturb_condition={},n_propagation=3);none_max=float(np.max(np.abs(o.adata.layers['delta_X'])));assert none_max<1e-10
 outputs=[];patient=[];audits=[];qcs=[];gg=nx.from_pandas_edgelist(ed,source='source',target='target',create_using=nx.DiGraph)
 if job=='main':o.calculate_randomized_coef_table(random_seed=settings['randomized_seed'])
 for k in kos:
  valid=k in genes and k in tf_sources(tf) and k in gg and gg.out_degree(k)>0
  variants=[('native',3)]
  if job=='main' and k in K[:2]:variants+=[('native',2),('native',4),('randomized',3)]
  for mode,steps in variants:
   graph=gg
   if mode=='randomized':
    mat=o.coef_matrix_per_cluster_randomized[state];rr,cc=np.nonzero(mat.values);graph=nx.DiGraph();graph.add_edges_from(zip(mat.index[rr],mat.columns[cc]))
   lengths=nx.single_source_shortest_path_length(graph,k,cutoff=steps) if k in graph else {}
   dist=np.array([lengths.get(g,-1) for g in genes]);reach=(dist>=1)
   why='' if valid else ('source absent/ineligible' if k not in genes else 'source absent in mapped prior' if k not in tf_sources(tf) else 'source has no fitted outgoing edge')
   if valid:
    try:o.simulate_shift(perturb_condition={k:0.0},n_propagation=steps,GRN_unit='cluster',use_randomized_GRN=(mode=='randomized'),clip_delta_X=False)
    except ValueError as exc:why=str(exc)
   ok=not bool(why);base_reason=np.where(reach,'','no path within steps') if ok else np.repeat(why,len(genes))
   selfgene=(genes==k);assess=(reach|selfgene)&ok
   frame=pd.DataFrame(dict(gene=genes,ko=k,graph=mode,steps=steps,path_length=np.where(dist>=0,dist,np.nan),reachable=reach,is_KO=selfgene,status=np.where(assess,'EVALUABLE','NE'),reason=base_reason))
   if ok:
    d=np.asarray(o.adata.layers['delta_X']);sim=np.asarray(o.adata.layers['simulated_count']);assert np.allclose(d,sim-wt)
    if (~reach & ~selfgene).any():assert np.max(np.abs(d[:,~reach & ~selfgene]))<1e-8
    frame['native_mean_delta']=d.mean(0);frame['native_mean_abs_delta']=np.abs(d).mean(0)
    for pid,ind in m.groupby('patient_id').indices.items():
     if not assess.any():continue
     sliced=d[ind][:,assess];pp=pd.DataFrame(dict(gene=genes[assess],ko=k,graph=mode,steps=steps,patient_id=pid,cells=len(ind),signed_sum=sliced.sum(0),abs_sum=np.abs(sliced).sum(0),std_sum=np.divide(sliced.sum(0),sd[assess],out=np.full(assess.sum(),np.nan),where=sd[assess]>0)));patient.append(pp)
    qcs.append(dict(ko=k,graph=mode,steps=steps,no_perturbation_max_abs=none_max,below_WT_min=int((sim<wt.min(0)-1e-10).sum()),above_WT_max=int((sim>wt.max(0)+1e-10).sum()),at_zero=int((sim==0).sum()),total=sim.size,clip_delta_X=False))
    if job=='main' and k in K[:2] and mode=='native' and steps==3:
     focus=genes.isin(E+K+list(rd('tables/program_members.tsv').gene.unique()))
     np.savez_compressed(out/(k+'_layers.npz'),barcode=m.barcode.to_numpy(dtype=str),gene=genes[focus].to_numpy(dtype=str),imputed_count=wt[:,focus].astype('float32'),simulated_count=sim[:,focus].astype('float32'),delta_X=d[:,focus].astype('float32'))
   outputs.append(frame);audits.append(dict(ko=k,graph=mode,steps=steps,expression=k in genes,base_regulator=k in tf_sources(tf),outdegree=graph.out_degree(k) if k in graph else 0,status='EVALUABLE' if ok else 'NE',reason=why))
 save(pd.concat(outputs),'%s/all_genes.tsv.gz'%out);save(pd.DataFrame(audits),out/'nodes.tsv.gz');save(pd.DataFrame(qcs,columns=['ko','graph','steps','no_perturbation_max_abs','below_WT_min','above_WT_max','at_zero','total','clip_delta_X']),out/'simulation_qc.tsv.gz')
 cols='gene ko graph steps patient_id cells signed_sum abs_sum std_sum'.split();save(pd.concat(patient) if patient else pd.DataFrame(columns=cols),out/'patient.tsv.gz')
 complete('complete')
def tf_sources(tf):return set(v for a in tf.values() for v in a)
if __name__=='__main__':analyze(sys.argv[1],sys.argv[2])

