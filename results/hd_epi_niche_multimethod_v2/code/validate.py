from common import *
import subprocess
from load_result import load_niche_metadata,attach_niches
import anndata as ad
from scipy.sparse import csr_matrix
rows=[];a=anchors();T=train_index(a)
for o in (R/'models').iterdir():
 if not o.is_dir() or not (o/'assignments.parquet').exists():continue
 d=pd.read_parquet(o/'assignments.parquet');assert np.array_equal(a.key,d.section+'::'+d.bin_id);assert (d.iloc[T].source=='native_train').all();assert (d.iloc[T].label>=0).all();assert ((d.label<0)==(d.source=='unresolved')).all();rows.append(dict(tag=o.name,unique_IDs=True,all_anchors=len(d),T_preserved=True,source_consistency=True))
pd.DataFrame(rows).to_csv(R/'tables/final_validation.tsv',sep='\t',index=False)
for s in SECTIONS:
 n=pd.read_parquet(V/s/'nodes.parquet');n=n[n.analysis_qc].set_index('bin_id');obj=ad.AnnData(csr_matrix((len(n),1)),obs=n[['x_um','y_um','analysis_qc']].copy());prior=obj.obs.copy();obj=attach_niches(obj,R,s);assert obj.obs[prior.columns].equals(prior);assert obj.obs.loc[~obj.obs.niche_v2_anchor,obj.obs.columns.str.endswith('_label')].isna().all().all()
f=R/'tables/loader_validation.tsv'
if not f.exists():subprocess.run(['Rscript',str(R/'code/validate_loader.R'),str(R)],check=True)
print('VALIDATION_COMPLETE',len(rows),flush=True)
