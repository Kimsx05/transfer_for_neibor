from public_common import *
import openpyxl,gzip,csv
targets=sorted(set(sum(sets.values(),[])+['CHMP2A','SREBF2','SCAP']))
def f(a,pat):return next((C/a).glob(pat))
def meta(a):
 wanted={}
 with gzip.open(f(a,'*series_matrix*'),'rt') as h:
  for r in csv.reader(h,delimiter='\t'):
   if r and r[0] in ('!Sample_title','!Sample_geo_accession'):wanted[r[0][8:]]=r[1:]
 return pd.DataFrame(wanted)
def diff(a,b):
 a,b=np.asarray(a,dtype=float),np.asarray(b,dtype=float)
 av,bv=np.var(a,ddof=1)/len(a),np.var(b,ddof=1)/len(b);delta=np.mean(a)-np.mean(b)
 if av+bv==0:return dict(effect=delta,lo=np.nan,hi=np.nan,p=np.nan)
 df=(av+bv)**2/(av**2/(len(a)-1)+bv**2/(len(b)-1));se=np.sqrt(av+bv);rad=stats.t.ppf(.975,df)*se
 return dict(effect=delta,lo=delta-rad,hi=delta+rad,p=2*stats.t.sf(abs(delta/se),df))
def main():
 a='GSE175862';nk=pd.read_csv(f(a,'*scores*'));nk.columns=['gene','score'];nk['rank_desc']=nk.score.rank(ascending=False);nk['background']=len(nk)
 output(pd.DataFrame({'gene':targets}).merge(nk,on='gene',how='left'),'tables/GSE175862_NK_candidates.tsv');output(meta(a),'tables/GSE175862_samples.tsv')
 a='GSE175863';fn=f(a,'*.xlsx');x=pd.read_excel(fn,sheet_name='TPM Abundance',index_col=0);z=meta(a)
 z['column']=z.title.map({c.split('_',1)[1]:c for c in x});assert z.column.notna().all()
 z['contrast']=z.title.str.extract(r'^(\w+?_(?:Sh|Sg))',expand=False);z['arm']=np.where(z.title.str.contains('NT'),'control','CHMP2A_perturbed')
 output(z,'tables/GSE175863_samples.tsv')
 mat=np.log2(1+x.groupby(level=0).mean());rows=[];ms=[]
 for k,g in z.groupby('contrast'):
  ctrl=g[g.arm.eq('control')].column;ko=g[g.arm.ne('control')].column;assert len(ctrl)==len(ko)==2
  for gene in targets:
   if gene in mat.index:rows.append(dict(contrast=k,gene=gene,n_ctrl=2,n_perturbed=2,**diff(mat.loc[gene,ko],mat.loc[gene,ctrl])))
  for mask in ['Full','M20-masked']:
   for mod,gs in sets.items():
    members=[t for t in gs if t in mat.index and (mask=='Full' or t not in top)]
    if members:
     sc=mat.loc[members].mean();ms.append(dict(contrast=k,module=mod,mask=mask,members=' '.join(members),n_ctrl=2,n_perturbed=2,**diff(sc[ko],sc[ctrl])))
 output(adjust(pd.DataFrame(rows),'contrast'),'tables/GSE175863_genes.tsv');output(adjust(pd.DataFrame(ms),['contrast','mask']),'tables/GSE175863_modules.tsv')
 de=pd.read_excel(fn,sheet_name='DESeq2',usecols=range(7));output(de[de.iloc[:,0].isin(targets)],'tables/GSE175863_author_pooled.tsv')
 w=openpyxl.load_workbook(f('HMGCS1_2026','*.xlsx'),read_only=True,data_only=True)
 output(pd.DataFrame([(s.title,s.max_row,s.max_column) for s in w],columns=['sheet','rows','cols']),'tables/HMGCS1_workbook_index.tsv')
 rr=[]
 for sheet,blocks in {'Fig .1c':[(2,[3,4],'A549'),(9,[10,11],'H1299')],'Fig. 2h':[(4,[5,6],'A549'),(10,[11,12],'H1299')],'Fig. 3a':[(3,[4,5],'A549')],'Fig. 3b':[(4,[5,6],'H1299')],'Fig. 4e':[(3,[4,5],'A549'),(9,[10,11],'H1299')]}.items():
  sh=list(w[sheet].values)
  for header,rs,model in blocks:
   for row in rs:
    for j in range(3):
     name=sh[header-1][2+j*3]
     for t in range(3):
      val=sh[row-1][2+3*j+t];assert isinstance(val,(int,float));rr.append(dict(sheet=sheet,model=model,condition=str(sh[row-1][1]),group=name,source_row=row,source_col=3+3*j+t,value=val))
 rr=pd.DataFrame(rr);output(rr,'plot_data/HMGCS1_source.tsv')
 out=[]
 for (sheet,model,condition),g in rr.groupby(['sheet','model','condition']):
  arms=list(g.group.unique())
  for x1,x0 in [(arms[1],arms[0]),(arms[2],arms[1])]:
   b=g[g.group==x1].value;c=g[g.group==x0].value;out.append(dict(sheet=sheet,model=model,condition=condition,arm=x1,reference=x0,mean_arm=b.mean(),mean_ref=c.mean(),source_values_per_arm=3,**diff(b,c)))
 output(adjust(pd.DataFrame(out),'sheet'),'tables/HMGCS1_source_contrasts.tsv')
if __name__=='__main__':main()

