from common import *
import joblib
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score
from matplotlib import colormaps,colors

def main():
 cache=np.load(OUT/'03_features/features.npz');f={k:cache[k] for k in ['q1','q2','reference']};cache.close();ref=f['reference'];m=pd.read_parquet(OUT/'03_features/common_anchors.parquet');rm=m.iloc[ref].reset_index(drop=True);types=json.loads((OUT/'03_features/feature_columns.json').read_text())['q1'];candidate_support=[];candidate_profiles=[];group_stability=[];complete=[];pal=[]
 for name in ['R1','R2','R12']:
  od=OUT/'04_clustering'/name;models,labels,candidates=joblib.load(od/'candidates.joblib')
  for tag,lab in labels.items():
   mm=rm.copy();mm['cluster']=lab
   for c,g in mm.groupby('cluster'):
    known=g[g.source_patient_identity_known & g.patient_uid.notna()];candidate_support.append(dict(input=name,candidate=tag,cluster=c,n_ref=len(g),n_dataset=g.dataset_id.nunique(),n_true_patient=known.patient_uid.nunique(),n_section=g.section_id.nunique(),max_dataset_fraction=g.dataset_id.value_counts(normalize=True).max()))
    for ring,key in [('R1','q1'),('R2','q2')]:
     pp=f[key][ref][lab==c].mean(0)
     for typ,v in zip(types,pp):candidate_profiles.append(dict(input=name,candidate=tag,cluster=c,ring=ring,cell_type=typ,mean_relative_composition=v))
  for alg in ['KM','LEIDEN']:
   path=od/alg;lab=np.load(path/'reference_native_labels.npy');conf=json.loads((path/'working_config.json').read_text());cl=pd.read_parquet(path/'labels.parquet');assert cl.spot_uid.tolist()==m.spot_uid.tolist();k=conf['k'];cmap=colormaps['tab20']
   for c in range(k):pal.append(dict(config=name+'_'+alg,cluster=c,hex=colors.to_hex(cmap(c/(k-1) if k>1 else 0)),note='tab20 repeats colors when K>20; numeric cluster IDs and profile table authoritative; IDs do not align across schemes'))
   for seed in [124,125,126,127,128]:
    z=np.load(path/f'stability_{seed}.npz');sub=z['reference_local_indices'];sl=z['labels'];ct=pd.crosstab(lab[sub],sl);arr=ct.to_numpy();aa,bb=linear_sum_assignment(-arr);match={int(ct.index[a]):int(ct.columns[b]) for a,b in zip(aa,bb)}
    for c in np.unique(lab):
     d=match.get(int(c));a=lab[sub]==c;b=sl==d if d is not None else np.zeros(len(sl),bool);union=(a|b).sum();group_stability.append(dict(config=name+'_'+alg,seed=seed,cluster=int(c),matched_cluster=d,n_original_in_subsample=int(a.sum()),Jaccard=float((a&b).sum()/union) if union else np.nan))
   for ds,ix in m.groupby('dataset_id').groups.items():
    ix=np.asarray(list(ix));for_labels=cl.cluster_label_native_application.to_numpy()[ix];full=m.complete_rings.to_numpy()[ix]
    for c in range(k):complete.append(dict(config=name+'_'+alg,dataset_id=ds,cluster=c,n_all=int((for_labels==c).sum()),n_complete=int(((for_labels==c)&full).sum()),fraction_all=float((for_labels==c).mean()),fraction_complete=float((for_labels[full]==c).mean()) if full.any() else np.nan))
 tab(OUT/'05_comparison/candidate_reference_support.tsv',candidate_support);tab(OUT/'05_comparison/candidate_both_ring_profiles.tsv',candidate_profiles);tab(OUT/'05_comparison/stability_by_cluster.tsv',group_stability);tab(OUT/'05_comparison/complete_ring_label_composition.tsv',complete);tab(OUT/'05_comparison/cluster_palette.tsv',pal)
 # Baseline, pathology and existing frozen posterior audit source inventory.
 paths=[ROOT/'neighborhood_v1/00_inputs_and_audit/PHASE2_INPUT_GATE.md',ROOT/'neighborhood_v1/00_inputs_and_audit/qc/posterior_sampling_manifest.csv',ROOT/'subtype/v2/reference_training_v2/10_final_report/reference_v2_final_30type_freeze_report.md',ROOT/'mp10_destvi_state_pilot_v8/METHODS.md',ROOT/'mp10_destvi_state_pilot_v8/00_manifest/SCORE_CONFIG.json',ROOT/'spatial_ecology_top20_v2/05_07_shared_K12/v1/analysis_contract.json']
 tab(OUT/'00_preflight/provenance_assets.tsv',[dict(path=str(p),exists=p.exists(),sha256=sha(p) if p.exists() else None) for p in paths])
 high=m[m.epithelial_fraction>=.7];tab(OUT/'07_source_validation/high_epi_blind_coverage.tsv',high.groupby(['dataset_id','section_id']).agg(n=('spot_uid','size')).reset_index())
 js(OUT/'07_source_validation/prospective_association_validation.json',dict(frozen_before_spatial_MP10=True,original_top20=True,score='Direct Top20 AUCell',high_epi_threshold=.7,high_epi_n=len(high),labels='reuse six frozen native application labels; no refit',composition_only_design=dict(target='actual center mean c2l30 proportions and six frozen environment labels',n_anchor_per_dataset_max=100,seed=789,n_cells=20,n_technical_repeats=2,source='same-sample scRNA counts; independent source assignment, constant epithelial-state distribution across labels',depth='actual center nCount; thinning only, report unattainable low input depth',genes='actual section original feature universe and original available Top20 set; unavailable features cannot be silently zero-filled',sampling='multinomial center30, without cell duplication within bag; unsupported sample-type requests logged',primary_comparison='patient/source-stratified label contrasts, no spot-independent P values',limitation='Cannot perform patient-matched causal attribution without verified spatial-scRNA linkage'),status='PLANNED_REQUIRED_NOT_EXECUTED',reason='exact per-source feature availability, normal/malignant labels and patient matching need provenance reconciliation; general mixing already run is not this control'))
 log('Candidate, cluster-level stability, full-ring and prospective source diagnostics completed')
if __name__=='__main__':main()
