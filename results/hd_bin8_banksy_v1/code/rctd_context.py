from pathlib import Path
import pandas as pd,sys
r=Path(__file__).resolve().parents[1];s=sys.argv[1];o=r/s;b=r.parent/'hd_bin8_comparison_v1'/s
d=pd.read_csv(o/'metadata.tsv.gz',sep='\t');rt=pd.read_csv(b/'rctd_ROI_results.tsv.gz',sep='\t');z=rt.merge(d,on='bin_id',how='inner',validate='one_to_one',suffixes=('','_attachment'));assert z.qc_keep.all()
rows=[]
for tag in ['l00','l02','l08']:
 col=f'banksy_{tag}_cluster'
 for (roi,cl,group,status),a in z.groupby(['roi',col,'reference_conflict_group','rctd_status']):
  rows.append(dict(branch=tag,roi=roi,cluster=cl,reference_group=group,rctd_status=status,n=len(a),UMI=int(a.nCount_RNA.sum()),first_type_mode=a.first_type.mode().iloc[0] if a.first_type.notna().any() else 'not_fitted',median_weight_1=a.weight_1.median(),median_weight_2=a.weight_2.median()))
pd.DataFrame(rows).to_csv(o/'existing_RCTD_cluster_context.tsv',sep='\t',index=False)
z.groupby(['roi','reference_conflict_group','rctd_status'],dropna=False).agg(n=('bin_id','size'),UMI=('nCount_RNA','sum')).reset_index().to_csv(o/'existing_RCTD_conflict_context.tsv',sep='\t',index=False)
# Pair composition records weights/status as delivered, without relabeling reject or excluded bins.
z.groupby(['roi','first_type','second_type','rctd_status'],dropna=False).agg(n=('bin_id','size'),UMI=('nCount_RNA','sum'),median_weight_1=('weight_1','median'),median_weight_2=('weight_2','median')).reset_index().to_csv(o/'existing_RCTD_pair_context.tsv',sep='\t',index=False)
print(s,'existing RCTD context',len(z),'analysis bins; no new fitting')
