from base import *
from source import kinds
def summarize():
 rows=[];markers=[]
 for s in sections:
  for k in kinds:
   d=pd.read_parquet(R/s/(k+'_evidence.parquet'));assert d.id.is_unique and not (d.supported_tier&~d.model_candidate).any()
   assert (d.FOXP3_raw_screen==d.FOXP3.gt(0)).all()
   for roi in ['all_unique']+rd(R/'roi_manifest.tsv').query('section==@s').roi.tolist():
    t=d if roi=='all_unique'else d[d.roi.str.split('|',regex=False).apply(lambda x:roi in x)]
    rows.append(dict(section=s,roi=roi,branch=k,units=len(t),median_nCount=t.nCount.median(),median_nFeature=t.nFeature.median(),zero_or_no_source=t.nCount.eq(0).sum(),no_source=t.n_supplied_usable.eq(0).sum(),low_information_1_to24=t.nCount.between(1,24).sum(),model_candidates=t.model_candidate.sum(),supported=t.supported_tier.sum(),uncertain_model=(t.model_candidate&~t.supported_tier).sum(),FOXP3_raw=t.FOXP3_raw_screen.sum(),FOXP3_and_T=t.FOXP3_T_screen.sum(),model_with_raw_T=(t.model_candidate&t.T_UMI.gt(0)).sum(),model_with_FOXP3=(t.model_candidate&t.FOXP3_raw_screen).sum(),model_with_FOXP3_and_T=(t.model_candidate&t.FOXP3_T_screen).sum(),screen_union=(t.model_candidate|t.FOXP3_raw_screen).sum()))
   for group,sub in [('all',d),('model_candidate',d[d.model_candidate]),('FOXP3_raw',d[d.FOXP3_raw_screen]),('supported_tier',d[d.supported_tier])]:
    for g in T+core+epi+fib+mye:
     markers.append(dict(section=s,branch=k,group=group,gene=g,n=len(sub),positive=(sub[g]>0).sum(),fraction=(sub[g]>0).mean(),sum_UMI=sub[g].sum(),mean_UMI=sub[g].mean()))
 put(pd.DataFrame(rows),'unit_comparison.tsv');put(pd.DataFrame(markers),'candidate_marker_summary.tsv')
if __name__=='__main__':summarize()
