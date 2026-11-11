from pathlib import Path
import pandas as pd
r=Path(__file__).resolve().parents[1];notes=pd.read_csv(r/'tables/reviewed_region_notes.tsv',sep='\t',dtype={'cluster':str});choices={'Cervilla':'l02','D1_Johnson':'l08','UC_SCC':'l02','MPBC1':'l02','MPBC2':'l02'};rows=[]
for s,tag in choices.items():
 old=s in ['Cervilla','D1_Johnson']
 if not (r/s/'supplement_metadata.tsv.gz').exists() or (not old and not (r/s/'metadata.tsv.gz').exists()):continue
 m=pd.read_csv((r.parent/'hd_bin8_banksy_v1'/s if old else r/s)/'metadata.tsv.gz',sep='\t',dtype={f'banksy_{tag}_cluster':str});p=r/s/'supplement_metadata.tsv.gz';d=pd.read_csv(p,sep='\t');assert d.bin_id.equals(m.bin_id) and d.bin_id.is_unique
 review_cols=['ext_recommended_branch','ext_recommended_region','ext_region_review_note'];before=d[review_cols].copy() if all(c in d for c in review_cols) else None
 d['ext_recommended_branch']=tag;d['ext_recommended_region']=m[f'banksy_{tag}_region'];d['ext_region_review_note']='Original panel-rule regional summary; not uniform bin identity'
 for _,z in notes[(notes.section==s)&(notes.branch==tag)].iterrows():
  ii=m[f'banksy_{tag}_cluster']==z.cluster;assert ii.any();d.loc[ii,'ext_recommended_region']=z.reviewed_region;d.loc[ii,'ext_region_review_note']=z.review_note
 if before is None or not before.equals(d[review_cols]):
  for fig in (r/s/'figures').glob('recommended_reviewed_regions*'):fig.unlink()
  d.to_csv(p,sep='\t',index=False)
 for label,z in d.groupby('ext_recommended_region'):rows.append(dict(section=s,branch=tag,region=label,n=len(z),review_layer_only=True))
pd.DataFrame(rows).to_csv(r/'tables/recommended_region_counts.tsv',sep='\t',index=False);print('REVIEW_LAYER_ATTACHED; original region and type fields unchanged')
