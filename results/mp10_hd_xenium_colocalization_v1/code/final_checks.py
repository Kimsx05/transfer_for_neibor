from base import *
from PIL import Image
from datetime import datetime
import gzip
rows=[]
def chk(s,key,ok,detail):
 rows.append((s,key,bool(ok),str(detail)))
def qc(s):
 hd=s in HD;d=O/('hd' if hd else 'xenium')/s;z=pd.read_parquet(d/'labels.parquet');n=pd.read_parquet(d/'graph_nodes.parquet').set_index('id')
 chk(s,'IDs',z.index.is_unique and z.index.equals(n.index),len(z))
 chk(s,'coverage',len(z)==pd.read_csv(d/'label_counts.tsv',sep='\t').units.sum(),z.flag.value_counts().to_dict())
 if hd:
  w=pd.read_parquet(d/'weights.parquet');chk(s,'weights finite and normalized',np.isfinite(w.to_numpy()).all() and np.allclose(w.sum(axis=1),1),len(w))
  e=pd.read_parquet(d/'epi64_centers.parquet');chk(s,'frozen H05',np.array_equal(e.H05,e.MP10>e.MP10.quantile(.95)),e.H05.sum())
  em=pd.read_parquet(F/s/'Epi64_members.parquet');tm=pd.read_parquet(d/'members.parquet')
  sh=em.merge(tm,on='bin_id',suffixes=('_epi','_bin'),validate='one_to_one')
  sh=sh[sh.region_id_epi.isin(e.region_id)].drop_duplicates(['region_id_epi','region_id_bin'])
  ii=pd.Index(e.region_id).get_indexer(sh.region_id_epi);jj=z.index.get_indexer(sh.region_id_bin)
  for rad in (50,200):
   independent=sparse.load_npz(d/f'epi64_independent_{rad}.npz')
   chk(s,f'no shared expression pairs {rad}',jj.min()>=0 and not np.asarray(independent[ii,jj]).any(),len(ii))
  chk(s,'full mixture weights review',(d/'mixture_weight_review.tsv').exists(),len(w))

 else:
  measured=z.epi_eligible;chk(s,'state population',(z[['Chol_biosynthesis_8','MP10_covered','H05_transfer','REACTOME_CHOLESTEROL_BIOSYNTHESIS','HALLMARK_CHOLESTEROL_HOMEOSTASIS']].notna().sum(axis=1)==5*measured).all(),sum(measured))
  chk(s,'confirmed expression clusters',z.loc[z.nFeature>0,'cluster'].notna().all(),z.cluster.nunique())
  chk(s,'boundary review',(d/'segmentation_ROIs.tsv').exists(),len(z))
 for r in [50,200]:
  if not hd or r==50:
   g=sparse.load_npz(d/f'graph_{r}.npz');chk(s,f'graph {r}',(g!=g.T).nnz==0 and not g.diagonal().any() and g.shape==(len(z),len(z)),g.nnz)
   for k in (['subtype','major'] if hd else ['subtype','major']+(['chol'] if z.chol_group.notna().any() else [])):
    t=pd.read_csv(d/f'{k}_{r}_matrix.tsv.gz',sep='\t');chk(s,f'{k} {r}',t.observed.sum()==t.directed_edges.iloc[0] and np.isclose(t.expected.sum(),t.observed.sum()) and t.loc[t.expected==0,'observed_expected'].isna().all() and (t.permutations==1000).all(),len(t))
  t=pd.read_parquet(d/f'center_composition_{r}.parquet');p=t.filter(regex='^fraction_');chk(s,f'NA for unevaluable {r}',p.loc[~t.evaluable].isna().all().all(),t.evaluable.sum())
 for p in (O/'figures'/s).glob('*.png'):
  with Image.open(p) as im:im.verify()
  chk(s,'figure pair '+p.stem,p.with_suffix('.pdf').stat().st_size>500,p.name)
samples=sys.argv[1:] or HD+XE
for s in samples:qc(s)
val=pd.read_csv(O/'reports/heldout_validation.tsv',sep='\t');chk('shared','holdout both depths',len(val)==6*val.held_out_source.nunique(),len(val));chk('shared','three-state coverage',len(pd.read_csv(O/'reports/xenium_states.tsv',sep='\t'))==5,'five')
put(pd.DataFrame(rows,columns=['sample','check','pass','detail']),O/'reports'/('validation_checks.tsv' if samples==HD+XE else 'xenium_validation_checks.tsv'))
assert all(r[2] for r in rows),[r for r in rows if not r[2]]
print('Verified',len(rows),'checks',flush=True)

