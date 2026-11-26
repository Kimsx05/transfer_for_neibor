from shared import *
import hashlib
src=R/'mp10_pathology_regions_v2/sources/SpaCET/inst/extdata/GeneSets/GOBP.gmt'
hits=[line.rstrip('\n').split('\t') for line in src.open() if line.startswith('CHOLESTEROL_EFFLUX\t')]
assert len(hits)==1
genes=set(filter(None,hits[0][2:]))
sets=read(O/'gene_sets.tsv');sets=sets[sets.module!='Efflux']
rows=[];cov=[];checksum=hashlib.sha256(src.read_bytes()).hexdigest()
for variant,pool in [('full',genes),('minus20',genes-set(TOP))]:
 for gene in sorted(pool):rows.append(dict(gs_name='CHOLESTEROL_EFFLUX',gene_symbol=gene,db_version='SpaCET 1.4.0 bundled GOBP.gmt',gs_collection='GO BP',gs_id='CHOLESTEROL_EFFLUX',module='Efflux',variant=variant,Top20_overlap=gene in TOP))
 for s in S:
  p=O/'.work'/s;idx=pd.Index((p/'genes').read_text().splitlines());ix=np.memmap(p/'i',mode='r',dtype='<i4');det=set(idx[np.bincount(ix,minlength=len(idx))>0])
  cov.append(dict(section_id=s,module='Efflux',variant=variant,n_total=len(pool),n_assayed=len(pool&set(idx)),n_detected=len(pool&det),detected_genes=','.join(sorted(pool&det)),n_background=len(idx),aucMaxRank=int(np.ceil(.05*len(idx))),member_top20=','.join(sorted(pool&set(TOP)))))
save(pd.concat([sets,pd.DataFrame(rows)]),'gene_sets.tsv');b=read(T/'gene_coverage_base.tsv');save(pd.concat([b[b.module!='Efflux'],pd.DataFrame(cov)]),'tables/gene_coverage_base.tsv')
m=read(O/'input_manifest.tsv')
save(pd.concat([m,pd.DataFrame([dict(section_id='shared',role='supplemental efflux GO frozen GMT',path=str(src),size=src.stat().st_size,sha256=checksum)])]).drop_duplicates(['role','section_id','path']),'input_manifest.tsv')

