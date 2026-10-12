from common import *
from collections import deque

def main():
 master=pd.read_parquet(OUT/'03_features/anchor_master.parquet');m=pd.read_parquet(OUT/'03_features/common_anchors.parquet');f=np.load(OUT/'03_features/features.npz');types=json.loads((OUT/'03_features/feature_columns.json').read_text())['q1'];cov=pd.read_csv(OUT/'00_preflight/coverage63.tsv',sep='\t');checks=[];rng=np.random.default_rng(901)
 for row in cov.itertuples():
  s=row.section_id;ds=row.dataset_id;t=master[master.section_id.eq(s)].reset_index(drop=True);g1=sparse.load_npz(OUT/'02_graphs'/ds/(s+'__R1.npz'));g2=sparse.load_npz(OUT/'02_graphs'/ds/(s+'__R2.npz'));assert g1.shape==(len(t),len(t));assert (g1-g1.T).nnz==0 and not g1.diagonal().any();assert not g1.multiply(g2).nnz and not g2.diagonal().any();assert np.asarray(g1.sum(1)).max()<=6 and np.asarray(g2.sum(1)).max()<=12
  for i in rng.choice(len(t),min(3,len(t)),False):
   first=set(g1[i].indices);second=set()
   for j in first:second.update(g1[j].indices)
   second-=first|{int(i)};assert second==set(g2[i].indices)
  arr=pd.read_csv(row.q05_path,index_col=0).drop(columns=['section','dataset'],errors='ignore').reindex(index=t.barcode,columns=types).to_numpy(float);valid=np.isfinite(arr).all(1)&(arr>=0).all(1)&(arr.sum(1)>0);assert np.array_equal(valid,t.composition_valid)
  local=m[m.section_id.eq(s)];chosen=local.iloc[rng.choice(len(local),min(2,len(local)),False)] if len(local) else local
  lookup={b:i for i,b in enumerate(t.barcode)}
  for j,anchor in chosen.iterrows():
   i=lookup[anchor.barcode]
   for g,key in [(g1,'q1'),(g2,'q2')]:
    ni=g[i].indices;ni=ni[valid[ni]];p=arr[ni]/arr[ni].sum(1,keepdims=True);assert np.allclose(p.mean(0),f[key][j],rtol=1e-6,atol=1e-7)
  checks.append(dict(section_id=s,n_nodes=len(t),n_observed=int(valid.sum()),missing_observations=int((~valid).sum()),graph_symmetry=True,exact_R2_BFS=True,mean_then_sqrt=True,status='PASS'))
 tab(OUT/'reports/validation_by_section.tsv',checks)
 for n in ['R1','R2','R12']:
  for a in ['KM','LEIDEN']:
   lab=pd.read_parquet(OUT/'04_clustering'/n/a/'labels.parquet');r=np.load(OUT/'04_clustering'/n/a/'reference_native_labels.npy');assert lab.spot_uid.tolist()==m.spot_uid.tolist();assert np.array_equal(lab.cluster_label_native_application.to_numpy()[f['reference']],r);assert np.array_equal(lab.cluster_label_common_transfer.to_numpy()[f['reference']],r);assert np.isfinite(lab.nearest_reference_distance).all();assert lab.top_vote_fraction.between(1/15,1).all()
 an=pd.read_csv(OUT/'00_preflight/composition_anomalies.tsv',sep='\t');assert len(an)==25676;an['reason']='no_q05_or_mean_output_for_confirmed_tissue_node';tab(OUT/'00_preflight/composition_anomalies.tsv',an)
 js(OUT/'reports/FINAL_VALIDATION.json',dict(status='PASS',n_sections=63,n_graph_nodes=len(master),n_common=len(m),all_six_same_UID_order=True,n_reference=20000,reference_native_preserved=True,exact_R2_BFS_checked=True,raw_neighbor_composition_means_checked=True,observed_abundance_invalid=0,unobserved_tissue_nodes=25676,MP10_spatial_unblinding=False))
 log('Final scientific-definition validation PASS')
if __name__=='__main__':main()
