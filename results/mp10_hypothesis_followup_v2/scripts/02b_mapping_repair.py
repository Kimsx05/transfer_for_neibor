import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
import pathlib,hashlib,json,sqlite3,pickle,functools,re
import pandas as pd,numpy as np,pyarrow as pa
pa.set_cpu_count(1);pa.set_io_thread_count(1)
from ctxcore.rnkdb import FeatherRankingDatabase
from pyscenic.utils import load_motif_annotations
from pyscenic.transform import module2df,module2features_auc1st_impl,df2regulons
P=pathlib.Path.cwd();B=P/'results/mp10_hypothesis_followup_v2';R=P/'resources/mp10_regulatory_priors';O=P/'results/mp10_immune_endothelial_mechanisms_v1';T=B/'tables';dbfile=R/'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather';db=FeatherRankingDatabase(str(dbfile),name='hg38_refseq-r80_mc9nr_10kb');ranking=set(db.genes)
hfile=pathlib.Path('/refdir/Rlib/immunedeconv/extdata/quantiseq/HGNC_genenames_20170418.txt');h=pd.read_csv(hfile,sep='\t',dtype=str).fillna('');sfile=pathlib.Path('/refdir/Rlib/org.Hs.eg.db/extdata/org.Hs.eg.sqlite');c=sqlite3.connect('file:'+str(sfile)+'?mode=ro',uri=True)
info=pd.read_sql('select g._id,g.symbol,t.gene_type from gene_info g left join genetype t on g._id=t._id',c);al=pd.read_sql('select * from alias',c);en=pd.read_sql('select * from ensembl',c);current=info.groupby('symbol')._id.apply(set).to_dict();alias=al.groupby('alias_symbol')._id.apply(set).to_dict();ens=en.groupby('ensembl_id')._id.apply(set).to_dict();itype=info.set_index('_id').gene_type.to_dict();isym=info.set_index('_id').symbol.to_dict();geneids={}
for g in ranking:geneids[g]=current.get(g,alias.get(g,set()))
byid={}
for g,ids in geneids.items():
 if len(ids)==1:byid.setdefault(next(iter(ids)),set()).add(g)
happroved=set(h.loc[h.Status.eq('Approved'),'ApprovedSymbol']);hbyens=h.loc[h.EnsemblGeneID.ne('')].groupby('EnsemblGeneID').ApprovedSymbol.apply(set).to_dict();hpre={};hcert={}
for row in h.itertuples():
 hids=current.get(row.ApprovedSymbol,ens.get(row.EnsemblGeneID,set()))
 for v in [row.ApprovedSymbol]+re.split(r'[|,;]',row.PreviousSymbols)+re.split(r'[|,;]',row.AliasSymbol):
  v=v.strip()
  if v:
   hpre.setdefault(v,set()).add(row.ApprovedSymbol)
   if len(hids)==1:hcert.setdefault(v.upper(),set()).update(hids)
mods=pickle.load(open(B/'objects/SREBF_modules.pkl','rb'));genes=sorted(set(g for m in mods for g in m.genes));rows=[];mapping={}
for g in genes:
 ids=current.get(g,ens.get(g,alias.get(g,set())));typ=';'.join(sorted(str(itype.get(i,'unknown')) for i in ids));sym=';'.join(sorted(isym.get(i,'') for i in ids));candidate=set()
 # Crosswalk requires a unique gene identity AND a ranking name certified by the installed HGNC table (same Ensembl or approved historical name).
 if len(ids)==1:
  possible=byid.get(next(iter(ids)),set());candidate={v for v in possible if hcert.get(v.upper(),set())==ids}
 if g in ranking:out=g;reason='exact_ranking_symbol'
 elif len(ids)>1:out=g;reason='ambiguous_gene_identity'
 elif len(candidate)==1:out=next(iter(candidate));reason='unique_NCBI_identity_and_HGNC_documented_name_to_ranking'
 elif len(candidate)>1:out=g;reason='ambiguous_multiple_ranking_symbols'
 elif not ids:out=g;reason='unresolved_Ensembl_ID' if g.startswith('ENSG') else 'unresolved_clone_or_symbol'
 else:out=g;reason='ranking_not_covered_'+('noncoding_or_pseudogene' if any(k in typ.lower() for k in ['ncrna','pseudo','snrna','snorna']) else 'known_gene')
 mapping[g]=out;rows.append(dict(original_gene=g,current_symbol=sym,NCBI_internal_ids=';'.join(map(str,sorted(ids))),gene_type=typ,ranking_candidates=';'.join(sorted(candidate)),ranking_gene=out,corrected=out!=g,reason=reason,HGNC_source=str(hfile),NCBI_source=str(sfile)))
mapdf=pd.DataFrame(rows);mapdf.to_csv(T/'scenic_gene_mapping_audit.tsv',sep='\t',index=False)
ann=load_motif_annotations(str(R/'motifs-v9-nr.hgnc-m0.001-o0.0.tbl'));ledger=[];dfs=[];mapped_members=[]
for m in mods:
 tf=m.transcription_factor;context=';'.join(sorted(m.context));mid=tf+'__'+context+'__'+hashlib.sha256(';'.join(sorted(m.genes)).encode()).hexdigest()[:16];mp={g:mapping[g] for g in m.genes};rev={}
 for g,v in mp.items():rev.setdefault(v,[]).append(g)
 collision={v for v,gs in rev.items() if len(gs)>1}
 for g,v in list(mp.items()):
  if v in collision:mp[g]=g
 assert len(set(mp.values()))==len(m.genes),'mapping collision remaining'
 renamed=m.copy(gene2weight={mp[g]:w for g,w in m.gene2weight.items()});n=len(m.genes);nm=sum(g in ranking for g in renamed.genes);mm=module2df(db,renamed,ann,module2features_func=functools.partial(module2features_auc1st_impl,rank_threshold=5000,auc_threshold=.05,nes_threshold=3,filter_for_annotation=True));back={v:g for g,v in mp.items()};ret=[]
 if len(mm):
  mm[('Enrichment','TargetGenes')]=mm[('Enrichment','TargetGenes')].map(lambda pairs:[(back[g],w) for g,w in pairs]);dfs.append(mm);ret=sorted(set(g for pp in mm[('Enrichment','TargetGenes')] for g,w in pp))
 ledger.append(dict(module_id=mid,TF=tf,context=context,n_original=n,n_mapped_original=sum(g in ranking for g in m.genes),n_mapped_corrected=nm,mapping_fraction=nm/n,n_names_corrected=sum(g!=v for g,v in mp.items()),n_ambiguous_collision_rejected=len(collision),n_retained_motifs=len(mm),n_retained_targets=len(ret),retained_targets=';'.join(ret),status='retained' if len(mm) else 'mapping_fail' if (n-nm)/n>=.2 else 'motif_annotation_or_NES_fail'))
 for g,v in mp.items():mapped_members.append(dict(module_id=mid,original_gene=g,ranking_gene=v,corrected=g!=v,collision_rejected=mapping[g] in collision))
 print(tf,context,n,nm,len(mm),flush=True)
pd.DataFrame(ledger).to_csv(T/'scenic_corrected_module_ledger.tsv',sep='\t',index=False);pd.DataFrame(mapped_members).to_csv(T/'scenic_corrected_module_members.tsv.gz',sep='\t',index=False)
regs=df2regulons(pd.concat(dfs)) if dfs else [];pickle.dump(regs,open(B/'objects/SREBF_corrected_regulons.pkl','wb'));rows=[]
for r in regs:
 for g,w in r.gene2weight.items():rows.append(dict(TF=r.transcription_factor,regulon=r.name,target=g,weight=w))
pd.DataFrame(rows).to_csv(T/'scenic_corrected_regulon_edges.tsv',sep='\t',index=False)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();meta=dict(correction='HGNC_20170418_approved_names + NCBI_2025Sep24_unique_gene_identity',sources={str(p):sha(p) for p in [hfile,sfile]},reference_unchanged=str(dbfile),original_module_denominators_preserved=True,ambiguous_mappings_rejected=True,original_counts_membership_unchanged=True,new_grn_inference=False,full_HGNC_download='external download unavailable; use installed versioned resources; unresolved names remain unknown',regulons=[dict(name=r.name,TF=r.transcription_factor,n_targets=len(r.genes)) for r in regs],NCBI_web_confirmation_PRAG1='https://www.ncbi.nlm.nih.gov/gene/157285')
(B/'config/scenic_mapping_correction.json').write_text(json.dumps(meta,indent=2));print('MAPPING REPAIR COMPLETE',[(r.name,len(r.genes)) for r in regs],flush=True)
