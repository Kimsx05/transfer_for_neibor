from common import *
import h5py, scipy.sparse as sp
from collections import defaultdict

def balanced(meta,n,seed):
 rng=np.random.default_rng(seed)
 groups={k:rng.permutation(np.array(sorted(ix)))[:cfg['pool_cap_patient_state']] for k,ix in meta.groupby(['patient_id','state']).groups.items()}
 # deterministic balanced quotas at both levels, with randomized tie order
 patients=rng.permutation(sorted(meta.patient_id.unique()))
 from itertools import cycle
 order={p:cycle([tuple(k) for k in rng.permutation([k for k in groups if k[0]==p])]) for p in patients}
 quotas=defaultdict(int); total=0
 while total<n:
  previous=total
  for p in patients:
   if total==n:break
   for _ in range(sum(1 for k in groups if k[0]==p)):
    k=next(order[p])
    if quotas[k]<len(groups[k]):
     quotas[k]+=1;total+=1;break
  if previous==total:break
 return sorted([cell for k,a in groups.items() for cell in a[:quotas[k]]])

def main():
 t=top.read_text().split()
 assert t=='HMGCS1 MSMO1 INSIG1 FDFT1 SQLE FDPS IDI1 ACSS2 ELOVL5 LDLR MVD SCD HMGCR ACAT2 SLC20A1 ACSL3 ACSL1 ADIPOR2 MXD1 NPC1'.split()
 meta=rd(base/'analysis_cell_metadata.tsv.gz');meta=meta[meta.analysis_major_c2l_final=='Epithelial'].set_index('barcode')
 assert meta.index.is_unique and meta.MP10_cell.notna().all()
 oldmeta=rd(old,usecols=['cell_id','epi_allgene_cluster','annotation']).set_index('cell_id')
 meta=meta.join(oldmeta,validate='one_to_one');meta['state']=meta.epi_allgene_cluster.map(lambda x:'C'+str(int(x)) if pd.notna(x) else 'Unassigned')
 sm=rd(base/'specimen_exposure_eligibility.tsv').set_index('sample_uid')
 meta['specimen_id']=meta.sample_uid.map(sm.specimen_id);meta['patient_id']=meta.sample_uid.map(sm.patient_uid_recorded)
 assert meta.patient_id.eq(meta.patient_uid_recorded).all() and meta.specimen_id.notna().all()
 assert meta[meta.specimen_id.isin(['T4463209a','T4463209b'])].patient_id.nunique()==1
 mem=rd(base/'pseudobulk_membership.tsv.gz')
 for key in ['Epi_all','H05','H10','H50','L50']:
  sub=mem[mem.profile_type==key];meta[key]=meta.index.isin(sub.barcode);assert meta[key].sum()==len(sub)
 assert meta.Epi_all.all() and (meta.H50!=meta.L50).all()
 for name,n,seed in [('pilot',cfg['pilot'],123),('main',cfg['formal'],123),('alt',cfg['formal'],456)]:
  meta[name]=meta.index.isin(balanced(meta[meta.state!='Unassigned'],n,seed))
 alloc=meta.groupby(['dataset','patient_id','state']).agg(eligible_cells=('sample_uid','size'),pilot=('pilot','sum'),main=('main','sum'),alt=('alt','sum')).reset_index()
 for k in ['pilot','main','alt']:alloc[k+'_conditional_probability']=alloc[k]/alloc.eligible_cells
 save(alloc,'tables/allocation.tsv')
 save(meta.reset_index(),'tables/eligible_metadata.tsv.gz')
 chosen=meta.index[meta[['main','pilot','alt']].any(axis=1)]
 def strings(f,k):
  g=f[k];assert not np.any(g['mask'][:]);return pd.Index(g['values'].asstr()[:])
 with h5py.File(raw) as f:
  barcodes=strings(f,'obs/_index');features=strings(f,'var/_index');assert features.is_unique and barcodes.is_unique and meta.index.isin(barcodes).all()
  assert features.equals(strings(f,'var/original_symbol')) and barcodes.equals(strings(f,'obs/barcode'))
  groups=programs();allgenes=set(targets).union(*groups.values());legacy={'NECTIN1':'PVRL1','NECTIN2':'PVRL2','NECTIN3':'PVRL3','NECTIN4':'PVRL4'}
  mapping={g:g if g in features else legacy.get(g,'') for g in sorted(allgenes)}
  mapping={g:v if v in features else '' for g,v in mapping.items()}
  # Prefer exact approved symbol; flag any additional legacy feature rather than silently merging it.
  mapdf=pd.DataFrame([dict(gene=g,raw_feature=v,approved_present=g in features,legacy_present=legacy.get(g) in features) for g,v in mapping.items()])
  save(mapdf,'qc/feature_map.tsv')
  panel=[v for v in mapping.values() if v];idx=features.get_indexer(panel)
  a=f['X'];ptr=a['indptr'][:];assert list(a.attrs['shape'])==[len(barcodes),len(features)]
  assert a.attrs['encoding-type']=='csr_matrix'
  keep=barcodes.isin(chosen); eligible=barcodes.isin(meta.index); parts=[];checks=[]
  for start in range(0,len(barcodes),1000):
   end=min(start+1000,len(barcodes));b=ptr[start:end+1];sl=slice(b[0],b[-1])
   X=sp.csr_matrix((a['data'][sl],a['indices'][sl],b-b[0]),shape=(end-start,len(features)))
   assert np.isfinite(X.data).all() and (X.data>=0).all() and (X.data==np.floor(X.data)).all()
   if keep[start:end].any():parts.append(X[keep[start:end]].astype('int32'))
   mask=eligible[start:end]
   if mask.any():
    sub=X[mask];bar=barcodes[start:end][mask];checks.append(pd.DataFrame(sub[:,idx].toarray(),index=bar,columns=panel));checks[-1]['raw_total']=np.asarray(sub.sum(1)).ravel();checks[-1]['raw_features']=np.asarray((sub>0).sum(1)).ravel()
  pcs=pd.concat(checks);m=meta.loc[pcs.index]
  audit=dict(source=str(raw),shape=list(map(int,a.attrs['shape'])),source_layer='Seurat v5 RNA/counts via pre-harmonization H5AD X',source_layers=rd(base/'source_layer_audit.tsv').to_dict('records'),counts_axis='cells x genes',eligible=len(meta),assigned=int((meta.state!='Unassigned').sum()),ncount_mismatch=int((pcs.raw_total!=m.nCount_RNA).sum()),nfeature_mismatch=int((pcs.raw_features!=m.nFeature_RNA).sum()))
  Path('qc/input.json').write_text(json.dumps(audit,indent=2));print(audit,flush=True)
  coverage=[]
  for g,v in mapping.items():
   for state,ind in [('ALL',m.index),*m.groupby('state').groups.items()]:
    z=pcs.loc[ind,v] if v else None
    coverage.append(dict(gene=g,raw_feature=v,state=state,n_cells=len(ind),n_positive=(z>0).sum() if v else None,detect_frac=(z>0).mean() if v else None,raw_variance=z.var() if v else None,n_positive_patients=m.loc[z[z>0].index,'patient_id'].nunique() if v else None))
  save(pd.DataFrame(coverage),'qc/raw_coverage.tsv')
  X=sp.vstack(parts,format='csr');ids=barcodes[keep];assert len(ids)==len(chosen) and X.shape==(len(ids),len(features))
  cols=np.flatnonzero(np.asarray(X.getnnz(0)).ravel()>0);sp.save_npz('cache/selected_counts.npz',X[:,cols])
  save(pd.DataFrame({'feature':features[cols]}),'cache/features.tsv');save(meta.loc[ids].rename_axis('barcode').reset_index(),'cache/metadata.tsv.gz')
  save(pd.DataFrame([(k,mask,g) for k,members in groups.items() for mask in ['Full','M20'] for g in sorted(members if mask=='Full' else members-set(t))],columns=['program','mask','gene']),'tables/program_members.tsv')
if __name__=='__main__':
 with measured('input_extraction'):main()

