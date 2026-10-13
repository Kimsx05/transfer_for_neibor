from base import *
import anndata as ad
from PIL import Image
P=OUT/'01_pathology_search/user_review';r=pd.read_csv(P/'ROI_review_to_fill.tsv',sep='\t');paths=set(r.evidence_figure)|set(r.raw_HE)
for rel in paths:
 p=P/rel;dest=p.with_suffix('.jpg')
 with Image.open(p) as im:
  im=im.convert('RGB');im.thumbnail((1600,1600));im.save(dest,quality=85,optimize=True)
for c in ['evidence_figure','raw_HE']:r[c]=r[c].str.replace('.png','.jpg',regex=False)
tab(P/'ROI_review_to_fill.tsv',r);p=P/'review_index.html';p.write_text(p.read_text().replace('.png','.jpg'))
# Portable coordinate evidence for splitting ROIs.
import shutil
for sec in r.section_id.unique():
 src=V3/'05_user_review/evidence'/sec;dst=P/'evidence'/sec;dst.mkdir(exist_ok=True)
 for n in ['image_coordinates.csv','registration.json']:
  if (src/n).exists():shutil.copy2(src/n,dst/n)
rows=[];top=pd.read_csv(OLD/'00_preflight/MP10_top20.tsv',sep='\t').gene.tolist();man=pd.read_csv(OUT/'03_epi_attribution/model_and_gene_universe_manifest.tsv',sep='\t')
for z in man.itertuples():
 for layer,path in [('raw',z.raw_path),('model',z.posterior_path)]:
  a=ad.read_h5ad(path,backed='r');genes=a.var_names.tolist();valid=[v for v in top if v in genes];ix=pd.Index(genes).get_indexer(valid);mat=sparse.csr_matrix(a.X[:,np.sort(ix)]);# backed requires sorted columns
  sums=np.asarray(mat.sum(0)).ravel();sortedgenes=np.asarray(genes)[np.sort(ix)];total=dict(zip(sortedgenes,sums));
  for gene in top:rows.append(dict(section_id=z.section_id,layer=layer,gene=gene,present=gene in genes,total_observed_counts=total.get(gene,np.nan),state='ABSENT' if gene not in genes else ('PRESENT_ALL_ZERO' if total[gene]==0 else 'PRESENT_DETECTED')))
  a.file.close()
 log('gene counts audit '+z.section_id)
tab(OUT/'00_audit_and_freeze/gene_absence_vs_zero.tsv',rows)
tab(OUT/'00_audit_and_freeze/model_checkpoint_hashes.tsv',[dict(section_id=z.section_id,path=z.model_path,sha256=sha(z.model_path)) for z in man.itertuples()])
log('additional audits complete')
