from cluster import *
import shutil
scope=json.loads((OUT/'00_preflight/RUN_SCOPE.json').read_text());assert scope['stop_after']=='niche' and not scope['mp10_authorized_this_run']
pilots=[]
for mode,strength in [('NONE','NONE'),('DATASET','WEAK'),('DATASET','STRONG'),('SAMPLE','WEAK'),('SAMPLE','STRONG')]:
 od=run_rep('R12',mode,strength)
 for alg in ['KM','LEIDEN']:
  r=json.loads((od/alg/'working_config.json').read_text());integ=pd.read_csv(od/'integration_metrics.tsv',sep='\t');r.update({f'{s.source_key}_same_neighbor':s.same_source_neighbor_fraction for _,s in integ.iterrows()});r['raw_neighbor_Hellinger']=integ.raw_Hellinger_neighbor_mean.iloc[0];pilots.append(r)
tab(OUT/'02_R12_pilot/pilot_comparison.tsv',pilots);df=pd.DataFrame(pilots);selected={};notes=[]
for mode in ['DATASET','SAMPLE']:
 w=df[(df.integration_mode==mode)&(df.strength=='WEAK')].set_index('algorithm');s=df[(df.integration_mode==mode)&(df.strength=='STRONG')].set_index('algorithm');key='dataset_id_same_neighbor' if mode=='DATASET' else 'sample_uid_same_neighbor'
 mixing=float(w[key].mean()-s[key].mean());rawcost=float(s.raw_neighbor_Hellinger.mean()/w.raw_neighbor_Hellinger.mean()-1);stable=bool((s.stability_ARI_median>=w.stability_ARI_median-.02).all());cover=bool((s.accepted_fraction>=w.accepted_fraction-.02).all());comp=bool((s.composition_preservation_summary<=w.composition_preservation_summary*1.05).all());strength='STRONG' if mixing>.02 and rawcost<=.05 and stable and cover and comp and (s.convergence_status=='converged').all() else 'WEAK';selected[mode]=dict(strength=strength,theta=2 if strength=='STRONG' else 1,**{'lambda':1 if strength=='STRONG' else 5});notes.append(f'{mode}: selected {strength}; strong mixing improvement={mixing:.5f}, raw Hellinger cost={rawcost:.3%}, stability preserved={stable}, coverage preserved={cover}, composition preserved={comp}. Conservative WEAK when tradeoffs conflict. No MP10 or UMAP used.')
js(OUT/'02_R12_pilot/selected_harmony_parameters.json',selected);(OUT/'02_R12_pilot/pilot_selection.md').write_text('\n\n'.join(notes)+'\n');log('PILOT SELECTION '+str(selected))
for ring in ['R12','R1','R2']:
 for mode in ['NONE','DATASET','SAMPLE']:
  od=run_rep(ring,mode,'NONE' if mode=='NONE' else selected[mode]['strength'])
  for alg in ['KM','LEIDEN']:
   dest=OUT/'03_main'/f'{mode}_{ring}_{alg}'
   if not dest.exists():dest.symlink_to(od/alg,target_is_directory=True)
rows=[json.loads((p/'working_config.json').read_text()) for p in sorted((OUT/'03_main').iterdir())];tab(OUT/'04_comparison/eighteen_scheme_comparison.tsv',rows);log(f'MAIN COMPLETE {len(rows)}/18')
