from common import *
import re,h5py
from collections import defaultdict
hg=rd(prior/'hgnc_complete_set_20261118.txt',low_memory=False).query("status=='Approved'").set_index('symbol')
lookup={}
for col in ['ensembl_gene_id','prev_symbol','alias_symbol']:
 t=defaultdict(set)
 for g,s in hg[col].dropna().items():
  for a in s.split('|'):t[a].add(g)
 lookup[col]=t
def resolve(s):
 if s in hg.index:return s,'approved',s
 for col,t in lookup.items():
  q=re.sub(r'\.\d+$','',s) if col=='ensembl_gene_id' and s.startswith('ENSG') else s
  if q in t:
   hits=t[q]
   if len(hits)==1:return next(iter(hits)),col,next(iter(hits))
   return 'RAW~'+s,'ambiguous','|'.join(sorted(hits))
 return 'RAW~'+s,'unresolved',''
def category(g):
 if g not in hg.index:return 'unclassified'
 h=hg.loc[g]; groups=set(str(h.gene_group).split('|'))
 if h.location=='mitochondria':return 'mitochondrial_genome'
 if h.locus_group=='protein-coding gene':
  if groups&{'L ribosomal proteins','S ribosomal proteins'} or g=='RACK1':return 'cytosolic_ribosome'
  if groups&{'Mitochondrial ribosomal proteins','Large subunit mitochondrial ribosomal proteins','Small subunit mitochondrial ribosomal proteins'}:return 'mitochondrial_ribosome'
 return 'other'
m0=rd(V/'cache/metadata.tsv.gz');m=m0[m0.main].copy().reset_index(drop=True)
assert len(m)==8000 and m.barcode.is_unique and m.patient_id.nunique()==38 and m.sample_uid.nunique()==39
assert m[m.specimen_id.isin(['T4463209a','T4463209b'])].patient_id.nunique()==1
for r in json.loads((V/'config/virtual_perturbation_manifest.json').read_text())['files']:assert sha(V/r['path'])==r['sha256']
fx=rd(V/'cache/features.tsv').feature; X=sp.load_npz(V/'cache/selected_counts.npz')[m0.main.values].tocsr()
assert X.shape==(8000,len(fx)) and X.data.min()>=0 and np.array_equal(X.data,np.floor(X.data))
assert np.array_equal(np.ravel(X.sum(1)),m.nCount_RNA.values) and np.array_equal(X.getnnz(1),m.nFeature_RNA.values)
raw=P/'subtype/v2/reference_training_v2/03_derived_reference/reference_v2_raw_counts_pre_harmonization.h5ad'
with h5py.File(raw) as f:
 allf=f['var/_index/values'].asstr()[:]; assert len(allf)==90011
 fullbar=pd.Index(f['obs/_index/values'].asstr()[:]);assert fullbar.is_unique and (fullbar.get_indexer(m.barcode)>=0).all()
 st=f['X/indptr'][:];assert np.array_equal(np.diff(st)[fullbar.get_indexer(m.barcode)],m.nFeature_RNA)
 # Exact per-cell support and value check, block-local read of sparse integer rows.
 z=f['X'];ix=pd.Index(allf).get_indexer(fx);assert min(ix)>=0
 ri=fullbar.get_indexer(m.barcode)
 for block in np.unique(ri//512):
  start=int(block*512);end=min(start+512,len(fullbar));p=st[start:end+1];s=slice(p[0],p[-1])
  b=sp.csr_matrix((z['data'][s],z['indices'][s],p-p[0]),shape=(end-start,len(allf)))
  hit=np.where(ri//512==block)[0];assert (b[ri[hit]-start][:,ix]!=X[hit]).nnz==0
mapping=pd.DataFrame([ (g,*resolve(g)) for g in allf],columns=['feature_id','gene','mapping','candidate_symbols']).set_index('feature_id')
unique=list(dict.fromkeys(mapping.gene));gi=pd.Index(unique)
code=gi.get_indexer(mapping.loc[fx,'gene'])
project=sp.csr_matrix((np.ones(len(code),dtype=np.int32),(np.arange(len(code)),code)),shape=(len(code),len(gi)))
Y=(X@project).tocsr();assert np.array_equal(np.ravel(Y.sum(1)),m.nCount_RNA)
def statistics(x):
 x=x.astype('float64');avg=np.ravel(x.mean(0))
 return np.asarray(x.getnnz(0)),sum((x[j].getnnz(0)>0).astype(int) for j in m.groupby('patient_id').indices.values()),avg,np.maximum(0,np.ravel(x.power(2).mean(0))-avg**2)
a,b,mu,var=statistics(X)
mapping['raw_positive_cells']=pd.Series(a,index=fx).reindex(mapping.index,fill_value=0)
mapping['raw_positive_patients']=pd.Series(b,index=fx).reindex(mapping.index,fill_value=0)
mapping['raw_variance']=pd.Series(var,index=fx).reindex(mapping.index,fill_value=0)
a,b,mu,var=statistics(Y);coverage=pd.DataFrame(dict(gene=gi,positive_cells=a,positive_patients=b,mean=mu,variance=var,detect=a/len(m))).set_index('gene')
coverage['eligible']=(a>=10)&(b>=3)&(var>0)
coverage['category']=[category(g) for g in gi]
coverage['source_feature_count']=mapping.gene.value_counts()
mapping['category']=mapping.gene.map(coverage.category);mapping['hgnc_id']=mapping.gene.map(hg.hgnc_id);mapping['ensembl_id']=mapping.gene.map(hg.ensembl_gene_id)
save(mapping.reset_index(),'tables/gene_mapping.tsv');rmcats={'cytosolic_ribosome','mitochondrial_ribosome','mitochondrial_genome'}
rlist=pd.DataFrame([(g,category(g)) for g in hg.index if category(g) in rmcats],columns=['gene','category']).merge(coverage.drop(columns='category').reset_index(),on='gene',how='left')
rlist['features']=rlist.gene.map(mapping.reset_index().groupby('gene').feature_id.agg('|'.join));save(rlist,'tables/removed_genes.tsv')
assert not rlist.gene.isin(['RPS6KB1','RPS6KB2','RPS6KA1','RPS6KA2','RPS6KA3','RPS6KA4','RPS6KA5','RPS6KA6']).any()
p=rd(V/'tables/program_members.tsv');top='HMGCS1 MSMO1 INSIG1 FDFT1 SQLE FDPS IDI1 ACSS2 ELOVL5 LDLR MVD SCD HMGCR ACAT2 SLC20A1 ACSL3 ACSL1 ADIPOR2 MXD1 NPC1'.split();assert set(p[p.program=='MP10'].gene)==set(top)
hm=rd(prior/'hallmark.tsv.gz');hm=hm[hm.gs_name=='HALLMARK_MTORC1_SIGNALING'];h=pd.DataFrame([('HALLMARK_MTORC1_SIGNALING',mask,g) for mask in ['Full','M20'] for g in set(hm.gene_symbol) if mask=='Full' or g not in top],columns=p.columns)
p=pd.concat([p,h]);p['frozen_gene']=p.gene;p['gene']=p.gene.map(lambda g:resolve(g)[0]);p=p.drop_duplicates(['program','mask','gene']);save(p,'tables/program_members.tsv')
legacy=rd(V/'tables/matched_controls.tsv').gene.map(lambda g:resolve(g)[0]);old=rd(V/'tables/model_panel.tsv').feature.map(mapping.gene)
tfs=set((prior/'human_tfs.txt').read_text().split());ann=rd(prior/'motifs-v9-nr.hgnc-m0.001-o0.0.tbl',usecols=['gene_name','description'])
lambert=Path('cache/Lambert_TFs_v1.01.txt').read_text().split()
assert len(lambert)>1600 and sum(s.startswith('ENSG') for s in lambert)>1600
strict={resolve(s)[0] for s in lambert if s.startswith('ENSG') and not resolve(s)[0].startswith('RAW~')}
alltf={resolve(g)[0] for g in tfs}|strict
res=set(K+E)|set(p.gene);v=coverage
v['control_class']=['sequence_TF' if g in strict else 'non_TF_protein' if g not in alltf and g in hg.index and hg.at[g,'locus_group']=='protein-coding gene' else 'other' for g in v.index]
eligible=v[v.eligible & ~v.category.isin(rmcats) & ~v.index.isin(res)];extra=set()
for k in K:
 a=eligible[eligible.control_class==('sequence_TF' if k in K[:2] else 'non_TF_protein')];t=a[['mean','detect']].copy();t['mean']=np.log1p(t['mean']);d=((t-[np.log1p(v.at[k,'mean']),v.at[k,'detect']])/t.std()).pow(2).sum(1);extra.update(d.sort_values(kind='stable').head(C['control_candidate_neighbors']).index)
v['v1']=v.index.isin(old);v['reserved']=v.index.isin(res);v['control_supplement']=v.index.isin(extra);v['legacy_control']=v.index.isin(legacy)
v['P_full']=v.eligible & (v.v1|v.reserved|v.control_supplement|v.legacy_control);v['P_clean']=v.P_full & ~v.category.isin(rmcats)
v['status']=np.select([v.variance==0,~v.eligible,~v.P_full,v.P_full&~v.P_clean],['zero variance','low cell/patient detection','not in panel','removed from clean'],default='eligible in clean')
assert v.P_full.sum()<=3000
for g in sorted(res-set(v.index)):v.loc[g,['eligible','P_full','P_clean','status']]=[False,False,False,'no resolved raw feature']
save(v.reset_index(),'tables/coverage.tsv');panel=v[v.P_full==True].sort_index();save(panel.reset_index(),'tables/panel_members.tsv');save(m,'cache/metadata.tsv.gz')
Z=Y[:,gi.get_indexer(panel.index)].tocsr();sp.save_npz('cache/panel_counts.npz',Z)
Z.toarray().astype('<i4').tofile('cache/panel.bin')
save(m.groupby(['dataset','patient_id','specimen_id','state']).size().rename('n_cells').reset_index(),'tables/sample_counts.tsv')
save(pd.DataFrame([dict(n_cells=len(m),n_patients=m.patient_id.nunique(),n_specimens=m.specimen_id.nunique(),datasets=m.dataset.nunique(),raw_features=len(allf),full=int(panel.P_full.sum()),clean=int(panel.P_clean.sum()),cache_matches_h5ad=True,total_umi=int(m.nCount_RNA.sum()))]),'qc/input.tsv')
files=[raw,V/'cache/features.tsv',V/'cache/metadata.tsv.gz',V/'cache/selected_counts.npz',V/'tables/model_panel.tsv',V/'tables/program_members.tsv',V/'tables/matched_controls.tsv',V/'qc/jobs/main_effects.tsv.gz',prior/'hgnc_complete_set_20261118.txt',prior/'hallmark.tsv.gz',prior/'human_tfs.txt',prior/'motifs-v9-nr.hgnc-m0.001-o0.0.tbl',Path('cache/promoter.parquet'),Path('cache/Lambert_TFs_v1.01.txt')]
save(pd.DataFrame([(str(f.resolve()),f.stat().st_size,sha(f)) for f in files],columns=['path','bytes','sha256']),'config/input_fingerprints.tsv')
print(pd.read_csv('qc/input.tsv',sep='\t').to_string(index=False),flush=True)

