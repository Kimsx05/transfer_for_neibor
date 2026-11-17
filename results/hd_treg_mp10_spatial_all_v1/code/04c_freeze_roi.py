from common import *
assert not (R/'models/association_results_reviewed.json').exists(),'Do not select new ROIs after reviewing association results'
rows=[]
for s in SECTIONS:
 d=pd.read_parquet(R/s/'T32_regions.parquet');v=d[d.measurement_eligible];like=v[v.evidence.eq('Treg_like_multigene')];chosen=[]
 if len(like)<=6:chosen=list(like.index)
 else:
  for qt in [.1,.4,.7,.9]:
   order=(like.Treg_main-like.Treg_main.quantile(qt)).abs().sort_values().index
   i=next((i for i in order if all(np.hypot(d.loc[i,'x_um']-d.loc[j,'x_um'],d.loc[i,'y_um']-d.loc[j,'y_um'])>=256 for j in chosen)),order[0]);chosen.append(i)
 for j,i in enumerate(chosen):
  p=d.loc[i];rows.append(dict(section=s,roi=f'candidate_{j+1}',region_id=p.region_id,x0=p.x_um-128,x1=p.x_um+128,y0=p.y_um-128,y1=p.y_um+128,selection='all_multigene_candidates_when_n_le6'if len(like)<=6 else 'four_signal_quantiles_spatially_separated',evidence=p.evidence))
 # Two low-signature T-matched controls, matched to first/last candidate; no MP10 fields accessed.
 targets=chosen[:1]+chosen[-1:] if chosen else list(v.index[:2]);controls=[];low=v[v.Treg_main<=v.Treg_main.median()]
 for j,target in enumerate(targets):
  if not len(low):continue
  order=(low.T_lineage-d.loc[target,'T_lineage']).abs().sort_values().index;i=next((i for i in order if i not in controls),order[0]);controls.append(i);p=d.loc[i];rows.append(dict(section=s,roi=f'control_{j+1}',region_id=p.region_id,x0=p.x_um-128,x1=p.x_um+128,y0=p.y_um-128,y1=p.y_um+128,selection='matched_T_lineage_lower_Treg',evidence=p.evidence))
put(pd.DataFrame(rows),'roi_manifest.tsv');(R/'models/roi_freeze.json').write_text(json.dumps({'frozen_at':datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),'rules':'all Treg-like multigene regions when <=6, otherwise four score quantiles with256um separation; two matched-T lower-signature controls','MP10_used':False,'association_execution_had_started':(R/'models/association_started.json').exists(),'association_effects_not_read_for_selection':True},indent=2));print('ROI_FROZEN',len(rows))
