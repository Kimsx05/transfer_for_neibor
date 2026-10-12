from cluster import *
reps={p.resolve().parent for p in (OUT/'03_main').iterdir()}|{p for p in (OUT/'02_R12_pilot').iterdir() if p.is_dir() and (p/'candidate_parameters.tsv').exists()}
for od in sorted(reps):
 df=pd.read_csv(od/'candidate_parameters.tsv',sep='\t');lv=df[(df.algorithm=='LEIDEN')&df.converged&(df.k>1)&(df.tiny_clusters==0)];ss=lv[(lv.silhouette_embedding>0)&(lv.silhouette_raw_Hellinger>0)];lv=ss if len(ss) else lv;mm=lv[lv.k.between(6,12)];lv=mm if len(mm) else lv
 if not len(lv):lv=df[(df.algorithm=='LEIDEN')&(df.k>1)]
 best=lv.sort_values(['silhouette_embedding','accepted_fraction'],ascending=False).iloc[0];old=json.loads((od/'LEIDEN/working_config.json').read_text())
 if float(old['k_or_resolution'])!=float(best.parameter):
  archive=od/('composition_selection_revision_'+str(int(time.time())));archive.mkdir();(od/'LEIDEN').rename(archive/'LEIDEN');(od/'COMPLETE.json').rename(archive/'COMPLETE.json');log('Reconcile '+od.name);run_rep(old['ring_input'],old['integration_mode'],old['strength'])
# Keep pilot table consistent with final candidate selection; strength selection does not depend on NONE.
pilot=OUT/'02_R12_pilot/pilot_comparison.tsv';df=pd.read_csv(pilot,sep='\t')
for row in range(len(df)):
 rid=df.loc[row,'representation_id'];alg=df.loc[row,'algorithm'];r=json.loads((OUT/'02_R12_pilot'/rid/alg/'working_config.json').read_text())
 for k,v in r.items():
  if k in df.columns:df.loc[row,k]=v
# pandas numeric columns must not receive None; existing NONE NaNs remain NaN.
tab(pilot,df)
notes=[];frozen=json.loads((OUT/'02_R12_pilot/selected_harmony_parameters.json').read_text())
for mode in ['DATASET','SAMPLE']:
 w=df[(df.integration_mode==mode)&(df.strength=='WEAK')].set_index('algorithm');s=df[(df.integration_mode==mode)&(df.strength=='STRONG')].set_index('algorithm');key='dataset_id_same_neighbor' if mode=='DATASET' else 'sample_uid_same_neighbor';mix=float(w[key].mean()-s[key].mean());cost=float(s.raw_neighbor_Hellinger.mean()/w.raw_neighbor_Hellinger.mean()-1);st=bool((s.stability_ARI_median>=w.stability_ARI_median-.02).all());cv=bool((s.accepted_fraction>=w.accepted_fraction-.02).all());co=bool((s.composition_preservation_summary<=w.composition_preservation_summary*1.05).all());choice='STRONG' if mix>.02 and cost<=.05 and st and cv and co and (s.convergence_status=='converged').all() else 'WEAK';assert frozen[mode]['strength']==choice
 notes.append(f'{mode}: selected {choice}; strong mixing improvement={mix:.5f}, raw Hellinger cost={cost:.3%}, stability preserved={st}, coverage preserved={cv}, composition preserved={co}. Conservative WEAK when tradeoffs conflict. No MP10 or UMAP used.')
notes.append('Final finite-grid candidate reconciliation checks positive original-composition separation where supported. This changes some Leiden working granularities but does not change either frozen Harmony strength; the raw neighbor-distance cost alone exceeds the strong-selection bound in both modes. Initial choices remain archived.')
(OUT/'02_R12_pilot/pilot_selection.md').write_text('\n\n'.join(notes)+'\n')
