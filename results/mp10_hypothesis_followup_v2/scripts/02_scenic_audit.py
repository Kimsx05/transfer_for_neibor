import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
import pathlib,json,hashlib,pickle,functools,sys
import numpy as np,pandas as pd,pyarrow as pa
pa.set_cpu_count(1);pa.set_io_thread_count(1)
from pyscenic.utils import modules_from_adjacencies,load_motif_annotations,add_correlation
from pyscenic.transform import module2df,module2features_auc1st_impl,df2regulons
from ctxcore.rnkdb import FeatherRankingDatabase
from ctxcore.recovery import aucs
P=pathlib.Path.cwd();B=P/'results/mp10_hypothesis_followup_v2';O=P/'results/mp10_immune_endothelial_mechanisms_v1/objects/scenic';R=P/'resources/mp10_regulatory_priors';T=B/'tables'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for c in iter(lambda:f.read(8*1024**2),b''):h.update(c)
 return h.hexdigest()
assert (B/'config/M0_FROZEN.json').exists()
meta=json.loads((O/'Epi_regulons.metadata.json').read_text());dbpath=R/'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather';anpath=R/'motifs-v9-nr.hgnc-m0.001-o0.0.tbl'
for f,key in [(O/'Epi_training.tsv.gz','training_sha256'),(O/'Epi_seed123_adjacencies.tsv','adjacency_sha256'),(dbpath,'db_sha256'),(anpath,'motif_sha256')]:assert sha(f)==meta[key],str(f)
key=sha(pathlib.Path(__file__))+meta['input_hash'];cache=B/'objects/SREBF_modules.pkl';cp=B/'objects/SREBF_modules.metadata.json'
if cache.exists() and cp.exists() and json.loads(cp.read_text())['input_hash']==key:
 modules=pickle.load(open(cache,'rb'));print('Matched candidate module cache',flush=True)
else:
 print('Read exact training and full adjacency',flush=True);ex=pd.read_csv(O/'Epi_training.tsv.gz',sep='\t',index_col=0).T;adj=pd.read_csv(O/'Epi_seed123_adjacencies.tsv',sep='\t');assert ex.index.is_unique and ex.columns.is_unique
 tfstats=[]
 for tf in ['SREBF1','SREBF2']:
  x=ex[tf];tfstats.append(dict(TF=tf,n_cells=len(x),detected_cells=int((x>0).sum()),detected_fraction=float((x>0).mean()),mean=float(x.mean()),variance=float(x.var()),TF_predictor=tf in (R/'human_tfs.txt').read_text().split(),raw_edges=int(adj.TF.eq(tf).sum())))
 pd.DataFrame(tfstats).to_csv(T/'scenic_TF_training_audit.tsv',sep='\t',index=False)
 print('Reconstruct full global-quantile and top-regulator candidates',flush=True)
 adj=add_correlation(adj,ex,rho_threshold=.03,mask_dropouts=False);adj[adj.TF.isin(['SREBF1','SREBF2'])].to_csv(T/'scenic_SREBF_edges.tsv.gz',sep='\t',index=False)
 allmods=modules_from_adjacencies(adj,ex,min_genes=20,rho_threshold=.03,rho_mask_dropouts=False);assert len(allmods)==meta['n_modules'],(len(allmods),meta['n_modules']);modules=[m for m in allmods if m.transcription_factor in ['SREBF1','SREBF2']]
 pickle.dump(modules,open(cache,'wb'));cp.write_text(json.dumps(dict(input_hash=key,n_all_modules=len(allmods),n_SREBF_modules=len(modules),training_cells=len(ex),training_genes=ex.shape[1]),indent=2));del ex,adj,allmods
 print('Recovered',len(modules),'SREBF candidate modules',flush=True)
db=FeatherRankingDatabase(str(dbpath),name='hg38_refseq-r80_mc9nr_10kb');ann=load_motif_annotations(str(anpath));annsel=ann[ann.index.get_level_values(0).isin(['SREBF1','SREBF2'])];annsel.to_csv(T/'scenic_SREBF_motif_annotations.tsv',sep='\t')
oldreg=pickle.load(open(O/'Epi_regulons.pkl','rb'));oldgenes={r.transcription_factor:set(r.genes) for r in oldreg if r.transcription_factor in ['SREBF1','SREBF2']};rows=[];members=[];motifs=[];dfs=[]
for module in modules:
 tf=module.transcription_factor;context=';'.join(sorted(module.context));mid=tf+'__'+context+'__'+hashlib.sha256(';'.join(sorted(module.genes)).encode()).hexdigest()[:16];ranks=db.load(module);mapped=set(ranks.columns);n=len(module);coverage=len(mapped)/n;mapping_pass=(n-len(mapped))/n<.20
 av=aucs(ranks,db.total_genes,np.ones(len(mapped)),.05);nes=(av-av.mean())/av.std();dd=pd.DataFrame({'TF':tf,'MotifID':ranks.index,'AUC':av,'NES':nes});dd['NES_pass']=dd.NES>=3;aa=annsel.loc[tf].reset_index();dd=dd.merge(aa,on='MotifID',how='left');dd['annotation_supported']=dd.Annotation.notna();dd['module_id']=mid;dd['mapping_pass']=mapping_pass
 # All motifs retained, including failed NES and absent annotation, to distinguish stages.
 motifs.append(dd);mm=module2df(db,module,ann,module2features_func=functools.partial(module2features_auc1st_impl,rank_threshold=5000,auc_threshold=.05,nes_threshold=3.0,filter_for_annotation=True))
 if len(mm):dfs.append(mm)
 enriched_annotated=int(dd.loc[dd.NES_pass&dd.annotation_supported,'MotifID'].nunique());reason='mapping_coverage_not_above_80_percent' if not mapping_pass else ('no_motif_NES_ge_3' if not dd.NES_pass.any() else 'no_enriched_TF_annotated_motif' if enriched_annotated==0 else 'retained')
 kept=sorted(set(g for targets in mm.get(('Enrichment','TargetGenes'),[]) for g,w in targets))
 rows.append(dict(module_id=mid,TF=tf,context=context,n_original=n,n_mapped=len(mapped),mapping_fraction=coverage,mapping_pass=mapping_pass,n_NES_ge3=int(dd.loc[dd.NES_pass,'MotifID'].nunique()),n_annotated_motifs=int(dd.loc[dd.annotation_supported,'MotifID'].nunique()),max_annotated_NES=dd.loc[dd.annotation_supported,'NES'].max(),n_enriched_annotated=enriched_annotated,n_retained_motifs=len(mm),n_retained_targets=len(kept),retained_targets=';'.join(kept),reason=reason,final_regulon_present=tf in oldgenes,AUCell_status='existing_regulon_scored' if tf in oldgenes else 'not_evaluable_no_regulon'))
 for g in module.genes:members.append(dict(module_id=mid,TF=tf,gene=g,importance=module[g],mapped_exact=g in mapped,TF_added_self=g==tf,retained_target=g in kept))
 print(tf,context,n,len(mapped),reason,enriched_annotated,flush=True)
pd.DataFrame(rows).to_csv(T/'scenic_module_ledger.tsv',sep='\t',index=False);pd.DataFrame(members).to_csv(T/'scenic_module_members.tsv.gz',sep='\t',index=False);pd.concat(motifs).to_csv(T/'scenic_all_motif_diagnostics.tsv.gz',sep='\t',index=False)
if dfs:
 final=df2regulons(pd.concat(dfs));rr=[]
 for r in final:
  assert set(r.genes)==oldgenes.get(r.transcription_factor,set()),'Reconstructed final target mismatch'
  for g,w in r.gene2weight.items():rr.append(dict(TF=r.transcription_factor,regulon=r.name,target=g,weight=w))
 pd.DataFrame(rr).to_csv(T/'scenic_reconstructed_final_targets.tsv',sep='\t',index=False)
(B/'config/scenic_audit_provenance.json').write_text(json.dumps(dict(input_hash=key,original=meta,reference_boundary='missing_fraction >=0.20 excluded (exactly80% also excluded)',reconstructed_final_matches=True,seed456='key-edge check only; not complete ctx/regulon replication',original_files_unchanged=True),indent=2));print('SCENIC STAGE AUDIT COMPLETE',flush=True)
