from common import *
import numpy as np,pandas as pd

o,z,xy=load()
for c in ['low_information','review_strong_competition']:
 assert o[c].dtype==bool, (c,o[c].dtype)
finite=np.isfinite(z).all(1)&np.isfinite(xy).all(1);primary=finite&~o.low_information;conf=o.review_strong_competition.to_numpy();morph=o.morphology_review_flag.astype(str).str.contains('identity_uncertain|ambiguous|possible|irregular',regex=True).to_numpy();disputed=conf|morph
o['interpretation_major']=np.where(disputed,'Disputed',o.working_major);o['interpretation_subtype']=np.where(disputed,'Unknown',o.working_subtype);o['interpretation_display']=np.where(disputed,'Disputed',o.working_display_label);o['interpretation_reason']=np.select([conf&morph,conf,morph],['current_competition_and_morphology','current_strong_competition','morphology_identity_uncertain'],default='inherited_working_label');o['morphology_identity_uncertain']=morph
for name,mask in [('PRIMARY',primary),('EXPANDED_LOWINFO',finite),('NO_CURRENT_CONFLICT',primary&~conf)]:
 o[name+'_included']=mask;o[name+'_reason']=np.select([mask,~finite,name=='PRIMARY' and o.low_information if name=='PRIMARY' else (~finite if name=='EXPANDED_LOWINFO' else o.low_information)],['included','missing_finite_latent_or_coordinates','low_information'],default='current_strong_competition' if name=='NO_CURRENT_CONFLICT' else 'excluded');o[name+'_niche']=pd.Series(pd.NA,index=o.index,dtype='string')
o.to_parquet(ROOT/'cellcharter_assignments.parquet',compression='zstd');print(o.groupby('section_id')[[c for c in o if c.endswith('_included')]+['review_strong_competition','morphology_identity_uncertain']].sum().to_string());js(ROOT/'metadata/input_summary.json',dict(retained=len(o),latent_dim=z.shape[1],finite=int(finite.sum()),primary=int(primary.sum()),no_conflict=int((primary&~conf).sum()),disputed=int(disputed.sum()),morphology=int(morph.sum()),counts='Unmodified integer X in v2 section cells.h5ad',spatial_units='micrometers, spatial_um verified in preceding review',interpretation_not_model_input=True))
