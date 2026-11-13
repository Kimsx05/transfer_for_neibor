from common import *
import MENDER,anndata as ad
from scipy import sparse
from scipy.spatial import cKDTree
n=pd.read_parquet(V/'Cervilla/nodes.parquet');a=n[n.epi_anchor];center=a.iloc[len(a)//2];core=a[(a.component==center.component)&(abs(a.x_um-center.x_um)<=16)&(abs(a.y_um-center.y_um)<=16)];pool=n[n.analysis_qc&(n.component==center.component)&(abs(n.x_um-center.x_um)<=216)&(abs(n.y_um-center.y_um)<=216)].copy();obs=pd.DataFrame({'ct':pd.Categorical(pool.type_blca.fillna('Unknown'))},index=pool.bin_id);z=ad.AnnData(sparse.csr_matrix((len(pool),1)),obs=obs);z.obsm['spatial']=pool[['x_um','y_um']].to_numpy();m=MENDER.MENDER_single(z,ct_obs='ct',random_seed=123);m.set_MENDER_para(n_scales=8,nn_mode='radius',nn_para=25,count_rep='s',include_self=False);m.run_representation();ix=z.obs_names.get_indexer(core.bin_id);cols=[list(m.ct_unique).index(t) if t in m.ct_unique else -1 for t in TYPES['blca']];c=np.zeros((len(core),len(cols)))
for b in range(8):
 ring=np.zeros_like(c)
 for j,col in enumerate(cols):
  if col>=0:ring[:,j]=z.obsm[f'scale{b}'][ix,col]
 c+=ring
 if b in [1,7]:
  old=pd.read_parquet(V/'Cervilla'/f'r{25*(b+1)}.parquet').set_index('bin_id').loc[core.bin_id];assert np.array_equal(c,old[['blca_C_'+t for t in TYPES['blca']]].to_numpy())
(R/'tables/real_ROI_MENDER_equivalence.json').write_text(json.dumps(dict(section='Cervilla',core_anchors=len(core),full_QC_halo_nodes=len(pool),halo_um=200,rings=8,include_self=False,known_counts_match_frozen_50_200=True,center_x=float(center.x_um),center_y=float(center.y_um)),indent=2));print('REAL_ROI_EQUIVALENCE_PASSED',flush=True)
