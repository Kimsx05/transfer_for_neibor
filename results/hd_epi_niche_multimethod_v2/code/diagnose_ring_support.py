from common import *
a=anchors();d,P,ok=composition('blca',50);old=pd.read_parquet(R/'models/MENDER_adapted_blca_r50_ring_counts/assignments.parquet');new=pd.read_parquet(R/'models/MENDER_adapted_blca_r50_per_ring_normalized/assignments.parquet');known=np.concatenate([np.load(R/'models'/f'{s}_rings.npz')['known'][:,:2,0] for s in SECTIONS]);rows=[]
for lab in [14,5,3,22,0]:
 ix=old.label.to_numpy()==lab
 for nl in sorted(new.loc[ix,'label'].unique()):
  z=ix&(new.label.to_numpy()==nl)
  rows.append(dict(main_label=lab,normalized_label=int(nl),n=int(z.sum()),mean_Epi=float(P[z,0].mean()),pure_Epi_fraction=float((P[z,0]==1).mean()),mean_input_coverage=float(d.loc[z,'input_coverage_component'].mean()),inner_known_mean=float(known[z,0].mean()),outer_known_mean=float(known[z,1].mean()),any_empty_known_ring_fraction=float((known[z]==0).any(1).mean())))
pd.DataFrame(rows).to_csv(R/'tables/ring_support_case_study.tsv',sep='\t',index=False)
