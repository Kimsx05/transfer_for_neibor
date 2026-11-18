import pathlib,pandas as pd,json,datetime
B=pathlib.Path('results/mp10_hypothesis_followup_v2');p=B/'tables/model_registry.tsv';d=pd.read_csv(p,sep='\t')
if d.model_id.str.startswith('V2M_').any():raise SystemExit('mTOR diagnostics already registered')
rows=[]
for st in ['H05','H10','H50']:
 for mask in ['Full','M20']:
  for setting in ['main','min30','quality_RNA','composition_cluster0']:
   rows.append(dict(model_id=f'V2M_{len(rows)+1:03d}',kind='paired_score',setting=setting,FDR_family='mTORC1_paired_diagnostic_'+setting,x=st,y='HALLMARK_MTORC1_SIGNALING',state=st,mask=mask,min_cells=30 if setting=='min30' else 50,patient='main',primary=False,planned_n=6,analysis_role='v2_M2_diagnostic',statistical_unit='specimen paired; recorded patient units; a main'))
r=pd.DataFrame(rows);r.to_csv(B/'tables/mTORC1_diagnostic_model_registry.tsv',sep='\t',index=False);pd.concat([d,r],ignore_index=True).to_csv(p,sep='\t',index=False)
(B/'config/mTOR_diagnostic_addendum.json').write_text(json.dumps(dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),reason='User-prespecified M2 quality/composition and shared-driver checks; exact model enumeration before these additional effects. Main formal lipid/bridge families unchanged.',paired_models=24,EC_masked_GSEA='historical M0073 full gene signed ranking; remove Top20 from all Hallmark sets and ranking universe; all 50 pathways, separate BH50; never a substitute for NE fixed-background score'),indent=2))
