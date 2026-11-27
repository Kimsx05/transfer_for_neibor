from base_data import *
# Preserve all old rejection reasons on the original Xe panel, with unmeasured anchors distinct from zero.
for s in XE:
 a=dataset(s);z=a.obs;out=pd.DataFrame(index=z.index)
 out['old_major_margin']=z.major_support-z.major_runner;out['old_sub_margin']=z.support-z.runner
 def markers(labels,sets):
  v=np.full(len(z),np.nan);measured=np.zeros(len(z),int)
  for k,ls in sets.items():
   i=np.flatnonzero(labels.eq(k).fillna(False));cols=a.var_names.get_indexer(ls);cols=cols[cols>=0];measured[i]=len(cols)
   if len(i) and len(cols):v[i]=a.X[i][:,cols].getnnz(1)
  return v,measured
 cnt,meas=markers(z.major_candidate,M);anc,num=markers(z.candidate,{k:str(v).split() for k,v in tax.anchors.items()});mark,nm=markers(z.candidate,{k:str(v).split() for k,v in tax.markers.items()})
 out['old_major_marker_detected_raw']=cnt;out['old_fine_anchor_n_measured']=num;out['old_fine_anchor_detected_raw']=anc
 reasons={'unfitted_major':z.major_candidate.isna(),'major_marker_unmeasured':(meas==0)&z.major_candidate.notna(),'major_marker_zero':cnt==0,'major_score_below_0.50':z.major_support<.5,'major_margin_below_0.12':out.old_major_margin<.12,'mixed_competition':z.flag.eq('mixed'),'boundary_degenerate':z.boundary_degenerate,'fine_not_tested':z.candidate.isna(),'fine_marker_zero':mark==0,'fine_anchor_unmeasured':(num==0)&z.candidate.notna(),'fine_anchor_zero':anc==0,'fine_score_below_0.50':z.support<.5,'fine_margin_below_0.15':out.old_sub_margin<.15}
 out['old_gate_reasons']=['|'.join(k for k,r in zip(reasons,row) if r) for row in np.column_stack(list(reasons.values()))]
 out.to_parquet(paths(s)/'old_gate_audit.parquet')

