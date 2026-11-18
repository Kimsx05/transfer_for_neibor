import pathlib,pandas as pd,json,datetime
B=pathlib.Path('results/mp10_hypothesis_followup_v2');p=B/'tables/model_registry.tsv';d=pd.read_csv(p,sep='\t')
if d.kind.eq('scenic_activity').any():raise SystemExit('Conditional SCENIC registry already present')
rows=[]
for state in ['H05','H10','H50']:
 for mask in ['Full','M20']:
  for setting in ['main','min30','replace_b','common_H05']:
   for tf in ['SREBF1','SREBF2']:
    rows.append(dict(model_id=f'V2S_{len(rows)+1:03d}',kind='scenic_activity',setting=setting,FDR_family=f'recovered_SCENIC_{state}_{mask}_{setting}',x=state,y=tf,state=state,mask=mask,min_cells=30 if setting=='min30' else 50,patient='replacement' if setting=='replace_b' else 'main',primary=state=='H05' and mask=='Full' and setting=='main',planned_n=2,statistical_unit='specimen paired; recorded patient units; a or b'))
r=pd.DataFrame(rows);(B/'tables/model_registry_before_conditional_SCENIC.tsv').write_bytes(p.read_bytes());pd.concat([d,r],ignore_index=True).to_csv(p,sep='\t',index=False);r.to_csv(B/'tables/model_registry_conditional_SCENIC.tsv',sep='\t',index=False)
(B/'config/conditional_SCENIC_analysis.json').write_text(json.dumps(dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),trigger='pre-authorized M1 conditional branch: SREBF2 restored by unique HGNC/NCBI alias mapping, 36 targets',models=len(r),method='paired per-specimen mean AUCell difference; intercept-only OLS HC3 t with df=n-1; raw AUCell units, same original scoring gene order/seed; M20 removes signature genes and rank universe genes only',timing='registered before recovered full-cell AUCell/profile effects available',no_legacy_changes=True),indent=2))
print('Conditional registry',len(r),'models')
