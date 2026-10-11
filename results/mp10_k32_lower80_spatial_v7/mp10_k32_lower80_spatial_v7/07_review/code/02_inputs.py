from common import *
import anndata as ad,math,gc
from scipy import sparse
from scipy.stats import rankdata
if (R/'00_manifest/SPATIAL_INPUTS_READY.json').exists():
 for sec in config()['sections']:
  o=outmap(sec);f=json.loads((o/'INPUT.json').read_text());assert sha(o/'counts.h5ad')==f['counts_sha256']
 say('EXISTING SPATIAL INPUTS VERIFIED; frozen display genes preserved');raise SystemExit
cfg=config();md=pd.read_csv(R/'00_manifest/member_label_manifest.tsv.gz',sep='\t',low_memory=False).set_index('cell_id');labels=json.loads((R/'00_manifest/labels.json').read_text());a=ad.read_h5ad(R/'01_reference/TRAIN_counts.h5ad');genes=a.var_names
# TRAIN-only expression effect; full raw UMI is used only for descriptive normalization.
effects=pd.DataFrame({'gene':genes})
for lab,short in [(K,'K32'),(L,'lower80')]:
 x=a.X[a.obs.label_v7.eq(lab).to_numpy()];effects[short+'_mean']=np.asarray(x.mean(0)).ravel();effects[short+'_detection']=np.asarray((x>0).mean(0)).ravel()
effects['log2FC_raw_counts_pc0.1']=np.log2((effects.K32_mean+.1)/(effects.lower80_mean+.1))
src=Path('/home/data/t070721/Metabolistic_scRNA/Multi2/4.3 NMF allgene/output_full_cc_seed5_v2');top=pd.read_csv(src/'07_EpiMP_top50_genes.csv');top20=top.loc[(top.EpiMP=='EpiMP10')&(top['rank']<=20)].sort_values('rank');wt(top20,R/'00_manifest/FROZEN_MP10_TOP20.tsv');ann=pd.read_csv(P/'mp10_k32_validation_v4/02_rna_DEG/full/rest/DE_annotated.tsv.gz',sep='\t',usecols=['gene_id','lipid_pathway_annotations','n_lipid_pathway_annotations']).drop_duplicates('gene_id');effects=effects.merge(ann,left_on='gene',right_on='gene_id',how='left');effects['MP10_Top20']=effects.gene.isin(top20.gene);effects['MP10_Top50']=effects.gene.isin(top.loc[top.EpiMP.eq('EpiMP10'),'gene']);mp=pd.read_csv(P/'mp10_k32_validation_v4/00_audit/top50_signatures.tsv',sep='\t');lookup=mp.groupby('gene').MP.agg(lambda x:';'.join(sorted(set(x))));effects['frozen_14MP_Top50_memberships']=effects.gene.map(lookup);wt(effects,R/'01_reference/TRAIN_expression_effects.tsv.gz')
# Selection fixed before any new spatial prediction exists.
secs=[];raw={};man=pd.read_csv(P/'cell2location_full_run/v2_30type/01_manifest/01_spatial_manifest_frozen.csv').set_index('section_id');available=set(genes)
for sec in cfg['sections']:
 assert sec in man.index
 path=P/'cell2location_preflight_v1/02_visium/anndata'/(sec+'.h5ad');b=ad.read_h5ad(path);assert b.obs_names.is_unique and b.var_names.is_unique and (b.X.data>=0).all() and np.equal(b.X.data,np.floor(b.X.data)).all();assert set(genes)<=set(b.var_names);xy=np.asarray(b.obsm['spatial']);assert xy.shape==(b.n_obs,2) and np.isfinite(xy).all();raw[sec]=b;available&=set(b.var_names)
 priorpath=P/'cell2location_full_run/v2_30type/02_sections'/sec/'model_config.json';prior=json.loads(priorpath.read_text());o=outmap(sec);o.mkdir(parents=True,exist_ok=True);c=ad.AnnData(X=b[:,genes].X.copy(),obs=pd.DataFrame({'sample':sec},index=b.obs_names),var=pd.DataFrame(index=genes));c.write_h5ad(o/'counts.h5ad',compression='gzip');js(o/'INPUT.json',{'section':sec,'dataset':sec.split('_')[0],'input_path':str(path),'input_sha256':sha(path),'counts_sha256':sha(o/'counts.h5ad'),'gene_hash':key(genes),'spot_hash':key(b.obs_names),'prior_source':str(priorpath),'N_cells_per_location':prior['N_cells_per_location'],'detection_alpha':prior['detection_alpha'],'input_provenance':dict(b.uns['input_provenance']),'image_note':'preflight h5ad has coordinates but no embedded image; coordinate-only figures, no unverified H&E overlay'})
 secs.append({'section':sec,'dataset':sec.split('_')[0],'assay_recorded':'fresh-frozen polyA' if sec=='GSE171351_B1' else 'CytAssist FFPE v2','n_spots':b.n_obs,'n_genes_full':b.n_vars,'median_full_UMI':float(np.median(np.asarray(b.X.sum(1)).ravel())),'G0_available':len(genes),'Top20_available':sum(top20.gene.isin(b.var_names)),'N_cells_per_location':prior['N_cells_per_location'],'detection_alpha':prior['detection_alpha'],'counts_path':str(path),'coordinates_valid':True,'image_render':'coordinates only'})
wt(pd.DataFrame(secs),R/'00_manifest/SELECTED_SECTIONS.tsv')
cand=effects[effects.MP10_Top20&effects.gene.isin(available)&effects.n_lipid_pathway_annotations.gt(0)&effects.K32_detection.ge(.1)].sort_values(['log2FC_raw_counts_pc0.1','K32_detection','gene'],ascending=[False,False,True]);display=cand.head(3).copy();display['selection_rule']='Frozen MP10 Top20, measured in both selected slices and G0, inherited lipid pathway annotation, TRAIN K32 detection>=0.1, raw mean log2FC descending';wt(display,R/'00_manifest/DISPLAY_GENES.tsv');assert not any((outmap(x)/'model/model.pt').exists() for x in cfg['sections'])
direction=[]
for sample in a.obs['sample'].unique():
 s=a.obs['sample'].eq(sample).to_numpy();ik=s&a.obs.label_v7.eq(K).to_numpy();il=s&a.obs.label_v7.eq(L).to_numpy();nk=int(ik.sum());nl=int(il.sum())
 for g in top20.gene:
  j=genes.get_indexer([g])[0];ok=nk>=30 and nl>=30 and j>=0;mk=float(a.X[ik,j].mean()) if ok else np.nan;ml=float(a.X[il,j].mean()) if ok else np.nan;direction.append({'sample':sample,'source_group':md.loc[a.obs_names[s][0],'source_group'],'dataset':md.loc[a.obs_names[s][0],'dataset'],'gene':g,'n_K32':nk,'n_lower80':nl,'eligible':ok,'K32_mean':mk,'lower80_mean':ml,'log2FC_pc0.1':np.log2((mk+.1)/(ml+.1)) if ok else np.nan})
wt(pd.DataFrame(direction),R/'01_reference/TRAIN_source_directions.tsv');del a;gc.collect()
# Real spatial scores: reuse exact Top20 only; existing direct Top20 can differ from frozen rank Top20, so verify before reuse.
oldcode=(P/'MP10_spatial_multisection_top20_vs_core19_v1/code/run_benchmark.py').read_text();import ast
oldtop=ast.literal_eval(oldcode.split('TOP20=')[1].split('\n')[0]);scoreaudit=[]
for sec,b in raw.items():
 xy=np.asarray(b.obsm['spatial']);lib=np.asarray(b.X.sum(1)).ravel();d=pd.DataFrame({'spot_id':b.obs_names,'x':xy[:,0],'y':xy[:,1],'nCount_full':lib});previous=P/'MP10_spatial_multisection_top20_vs_core19_v1/02_spot_scores'/(sec.split('_')[0]+'__'+sec+'__spot_scores.csv.gz');reuse=False
 if previous.exists() and set(oldtop)==set(top20.gene):
  pr=pd.read_csv(previous).set_index('spot_id');reuse=set(pr.index)==set(b.obs_names)
  if reuse:d['MP10_Top20_UCell']=pr.loc[b.obs_names,'Top20_A_UCell'].to_numpy();d['MP10_Top20_AUCell']=pr.loc[b.obs_names,'Top20_A_AUCell'].to_numpy()
 if not reuse:
  idx=b.var_names.get_indexer(top20.gene);assert (idx>=0).all();n=len(idx);mr=1500;cut=math.ceil(.05*b.n_vars);den=sum(cut-r for r in range(1,min(n,cut)+1));uu=[];aa=[]
  for i in range(b.n_obs):
   x=b.X.getrow(i).toarray().ravel();r=rankdata(-x,method='average')[idx];uu.append(np.clip(1-(np.minimum(r,mr+1).sum()-n*(n+1)/2)/(n*mr),0,1));aa.append(np.maximum(0,cut-r).sum()/den)
  d['MP10_Top20_UCell']=uu;d['MP10_Top20_AUCell']=aa
 for g in ['EPCAM','KRT8','KRT18']+display.gene.tolist():
  if g in b.var_names:
   x=np.asarray(b[:,[g]].X.toarray()).ravel();d[g+'_counts']=x;d[g+'_log1pCP10k']=np.log1p(x*1e4/np.maximum(lib,1))
  else:d[g+'_counts']=np.nan;d[g+'_log1pCP10k']=np.nan
 wt(d,outmap(sec)/'raw_expression_plot_data.tsv.gz');scoreaudit.append({'section':sec,'reused_exact_spots_and_gene_set':reuse,'frozen_vs_previous_Top20_same_set':set(oldtop)==set(top20.gene),'gene_universe':b.n_vars,'method_source':str(P/'MP10_spatial_multisection_top20_vs_core19_v1/code/run_benchmark.py'),'UCell_maxRank':1500,'ties':'average','AUCell_rank_cutoff':math.ceil(.05*b.n_vars),'matrix':'full measured raw RNA counts; no G0 clipping','note':'Frozen v4/v6 MP10 rank<=20 used unchanged; previous direct Top20 is not substituted'})
wt(pd.DataFrame(scoreaudit),R/'00_manifest/SPATIAL_SCORE_AUDIT.tsv');js(R/'00_manifest/SPATIAL_INPUTS_READY.json',{'sections':cfg['sections'],'display_genes':display.gene.tolist(),'timestamp':time.time()});say('SPATIAL INPUTS READY',display.gene.tolist())
