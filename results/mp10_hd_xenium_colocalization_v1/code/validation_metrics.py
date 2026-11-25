from base import *
from sklearn.metrics import roc_auc_score,average_precision_score,precision_recall_curve
d=O/'shared/validation';z=pd.concat([pd.read_csv(p,sep='\t') for p in sorted(d.glob('test_*.tsv.gz'))]);out=[];curves=[]
for (ds,mode),t in z.groupby(['dataset','mode']):
 y=t.H05_truth=='H05'
 for k in ['Chol_biosynthesis_8','MP10_covered','H05_transfer']:
  v=t[k];valid=y.nunique()==2
  out.append((ds,mode,k,len(y),y.sum(),y.mean(),v.corr(t.MP10_cell),v.corr(t.MP10_cell,method='spearman'),roc_auc_score(y,v) if valid else np.nan,average_precision_score(y,v) if valid else np.nan))
  if valid:
   pr,re,th=precision_recall_curve(y,v);curves.append(pd.DataFrame({'source':ds,'mode':mode,'program':k,'threshold':np.r_[th,np.nan],'precision':pr,'recall':re}))
put(pd.DataFrame(out,columns=['held_out_source','mode','score','n','h05','prevalence','Pearson','Spearman','AUROC','average_precision']),O/'reports/heldout_validation.tsv')
put(pd.concat(curves),d/'PR_curves.tsv.gz')

