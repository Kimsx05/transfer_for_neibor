import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
import pathlib,pickle,functools,hashlib,pandas as pd,numpy as np,pyarrow as pa,sqlite3
pa.set_cpu_count(1);pa.set_io_thread_count(1)
from ctxcore.rnkdb import FeatherRankingDatabase
from ctxcore.recovery import aucs
from pyscenic.utils import load_motif_annotations
from pyscenic.transform import module2df,module2features_auc1st_impl
P=pathlib.Path.cwd();B=P/'results/mp10_hypothesis_followup_v2';T=B/'tables';R=P/'resources/mp10_regulatory_priors';mods=pickle.load(open(B/'objects/SREBF_modules.pkl','rb'));ledger=pd.read_csv(T/'scenic_corrected_module_ledger.tsv',sep='\t');mem=pd.read_csv(T/'scenic_corrected_module_members.tsv.gz',sep='\t');ann=load_motif_annotations(str(R/'motifs-v9-nr.hgnc-m0.001-o0.0.tbl'));db=FeatherRankingDatabase(str(R/'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather'),name='hg38_refseq-r80_mc9nr_10kb');rows=[];diag=[]
def annotation_type(s):
 s=str(s).lower();sim=('similar' in s);ort=('ortholog' in s)
 return 'similarity+orthology' if sim and ort else 'similarity' if sim else 'orthology' if ort else 'direct' if 'direct' in s else 'unclassified'
for m in mods:
 context=';'.join(sorted(m.context));mid=m.transcription_factor+'__'+context+'__'+hashlib.sha256(';'.join(sorted(m.genes)).encode()).hexdigest()[:16]
 if not ledger.loc[ledger.module_id.eq(mid),'status'].eq('retained').any():continue
 mp=mem[mem.module_id.eq(mid)].set_index('original_gene').ranking_gene.to_dict();back={v:k for k,v in mp.items()};n=m.copy(gene2weight={mp[g]:w for g,w in m.gene2weight.items()});mm=module2df(db,n,ann,module2features_func=functools.partial(module2features_auc1st_impl,rank_threshold=5000,auc_threshold=.05,nes_threshold=3,filter_for_annotation=True))
 for (tf,motif),row in mm.iterrows():
  pairs=row[('Enrichment','TargetGenes')];text=row[('Enrichment','Annotation')];rows.append(dict(module_id=mid,TF=tf,motif=motif,NES=row[('Enrichment','NES')],AUC=row[('Enrichment','AUC')],annotation=text,annotation_type=annotation_type(text),motif_similarity_q=row[('Enrichment','MotifSimilarityQvalue')],orthologous_identity=row[('Enrichment','OrthologousIdentity')],n_targets=len(pairs),retained_original_targets=';'.join(back[g] for g,w in pairs)))
 df=db.load(n);av=aucs(df,db.total_genes,np.ones(df.shape[1]),.05);ns=(av-av.mean())/av.std();d=pd.DataFrame(dict(module_id=mid,TF=n.transcription_factor,MotifID=df.index,NES=ns,AUC=av));d=d.merge(ann.loc[n.transcription_factor].reset_index(),on='MotifID',how='left');d['annotation_type']=d.Annotation.map(annotation_type);diag.append(d)
pd.DataFrame(rows).to_csv(T/'scenic_corrected_retained_motifs.tsv.gz',sep='\t',index=False);pd.concat(diag).to_csv(T/'scenic_corrected_retained_module_all_motifs.tsv.gz',sep='\t',index=False)
old=pd.read_csv(T/'scenic_SREBF_motif_annotations.tsv',sep='\t');old['annotation_type']=old.Annotation.map(annotation_type);old.to_csv(T/'scenic_motif_annotation_types.tsv',sep='\t',index=False)
mp=pd.read_csv(T/'scenic_gene_mapping_audit.tsv',sep='\t');c=sqlite3.connect('file:/refdir/Rlib/org.Hs.eg.db/extdata/org.Hs.eg.sqlite?mode=ro',uri=True);g=pd.read_sql('select _id,gene_id from genes',c).set_index('_id').gene_id.to_dict();mp['NCBI_Entrez_IDs']=mp.NCBI_internal_ids.fillna('').astype(str).map(lambda s:';'.join(str(g.get(int(float(v)),'')) for v in s.split(';') if v));mp.to_csv(T/'scenic_gene_mapping_audit_with_Entrez.tsv',sep='\t',index=False)
print('RETAINED MOTIF DETAILS',len(rows))
