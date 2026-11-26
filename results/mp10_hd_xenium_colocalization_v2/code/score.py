from common import *
from numba import njit
from scipy.stats import spearmanr
@njit
def rank_score(ptr,ind,val,lookup,n):
 out=np.zeros(len(ptr)-1)
 for i in range(len(out)):
  a,b=ptr[i],ptr[i+1];vs=val[a:b];hist=np.bincount(vs);asc=np.cumsum(hist)
  for j in range(a,b):
   if lookup[ind[j]]:
    k=val[j];rank=(b-a)-asc[k]+(1+hist[k])/2.
    out[i]+=max(905-rank,0.)
 return out/17890.
def matrix(m,ids,z):
 ix=pd.Index(z).get_indexer(m.loc[ids,'region_id']);assert min(ix)>=0
 return sparse.csr_matrix((np.ones(len(ids),np.int32),(ix,np.arange(len(ids)))),shape=(len(z),len(ids)))
def run(s):
 v,d=folders(s);hd=s in HD;z=labels(s);ident=pd.read_parquet(d/'identity.parquet');rngs=CFG['seeds']
 qrows=[];corr=[];overlaps=[]
 def deliver(y,elig,methods,res):
  y['eligible']=elig
  for k in methods:
   a=y.loc[elig,k];h=a>a.quantile(.95);y.loc[elig,k+'_high']=h.to_numpy()
   n=int(h.sum());multi=y.loc[a.index[h],k+'_detected'].ge(2).mean() if k+'_detected' in y else np.nan
   jacs=[];rhos=[]
   for seed in rngs:
    if str(seed) not in y:continue
    b=y.loc[a.index,str(seed)];hh=b>b.quantile(.95);jacs.append((hh&h).sum()/(hh|h).sum());rhos.append(a.corr(b,method='spearman'))
   qc='not_estimable' if n<5 or a.nunique()<3 else 'low_support' if n<20 or (jacs and np.median(jacs)<.5) or (multi<.5) else 'evaluable'
   qrows.append(dict(sample=s,res=res,state=k,n=len(a),high=n,low=len(a)-n,q95=a.quantile(.95),ties=(a==a.quantile(.95)).sum(),zeros=(a==0).sum(),distinct=a.nunique(),high_multigene=multi,high_single_gene_fraction=y.loc[a.index[h],k+'_dominance'].mean() if k+'_dominance' in y else np.nan,QC=qc,jaccard80=','.join(map(str,jacs)),rho80=','.join(map(str,rhos))))
   for var in ['UMI','nFeature','epi_weight','cell_area']+methods:
    if var in y and var!=k:corr.append((s,res,k,var,a.corr(y.loc[elig,var],method='spearman')))
  for a,b in itertools.combinations(methods,2):
   h=y.loc[elig,a+'_high'].astype(bool);j=y.loc[elig,b+'_high'].astype(bool);overlaps.append((s,a,b,(h&j).sum(),(h|j).sum(),(h&j).sum()/(h|j).sum()))
  y.index.name='id';y.to_parquet(d/f'states_{res}.parquet')
 if hd:
  src=json.loads((v/'conservation.json').read_text());x,ids,g=tenx(src['input']);genes=g.name.values
  assert len(genes)==18085 and g.id.is_unique and sum(np.isin(genes,TOP))==20
  m=pd.read_parquet(v/'members.parquet').set_index('bin_id');m16=matrix(m,ids,z.index);x16=m16@x;del x
  legacy=ROOT/'highres_VisiumHD/analysis/hd_treg_mp10_spatial_all_v1'/s
  m32=pd.read_parquet(legacy/'T32_members.parquet').set_index('bin_id');pairs=m.join(m32,lsuffix='16',rsuffix='32');ln=pairs[['region_id16','region_id32']].drop_duplicates().set_index('region_id16')
  assert ln.index.is_unique and set(ln.index)==set(z.index)
  y32=pd.read_parquet(legacy/'T32_measure.parquet',columns=['region_id','component','x_um','y_um','nCount','nFeature','n_bins']).set_index('region_id');g32=matrix(ln.rename(columns={'region_id32':'region_id'}),z.index,y32.index)
  ln.to_parquet(d/'parent32.parquet');weights=ident.epi_weight.copy()
  if 'human_epi_eligible' in ident:
   opinion=ident.human_epi_eligible.dropna();weights.loc[opinion.index]=opinion.astype(float)
  dep=z.nCount.to_numpy();den=g32@(weights.notna()*dep);y32['epi_weight']=(g32@(weights.fillna(0)*dep))/np.where(den>0,den,np.nan);y32['fitted_fraction']=den/y32.nCount
  y16=z[['component','x_um','y_um','nCount','nFeature','n_bins']].copy();y16['epi_weight']=weights;y16['fitted_fraction']=weights.notna().astype(float)
  for res,y,xx in [(16,y16,x16),(32,y32,g32@x16)]:
   xx=xx.tocsr();xx.eliminate_zeros();assert np.array_equal(xx.sum(1).A1,y.nCount) and np.array_equal(np.diff(xx.indptr),y.nFeature)
   y=y.rename(columns={'nCount':'UMI'});y['epi_markers']=xx[:,np.isin(genes,MAJ['Epithelial'])].getnnz(1)
   el=(y.UMI>=100)&(y.nFeature>=50)&(y.epi_markers>=1)&(y.epi_weight>=.50)&(y.fitted_fraction>=.5)
   if res==16 and 'human_epi_eligible' in ident:el=ident.human_epi_eligible.combine_first(el).astype(bool)
   yy=xx[el.values].copy();mpt=np.isin(genes,TOP);mp=yy[:,mpt].toarray();y.loc[el,'MP10_detected']=(mp>0).sum(1);y.loc[el,'MP10_dominance']=mp.max(1)/np.where(mp.sum(1)>0,mp.sum(1),np.nan)
   for key in ['MP10']+list(map(str,rngs)):
    yy2=yy if key=='MP10' else yy.copy()
    if key!='MP10':yy2.data=np.random.default_rng(int(key)).binomial(yy2.data,.8).astype('i4');yy2.eliminate_zeros()
    y.loc[el,key]=rank_score(yy2.indptr,yy2.indices,yy2.data,mpt,len(genes))
   from scipy.stats import rankdata
   chk=np.maximum(905-rankdata(-yy[:16].toarray().astype(float),axis=1,method='average')[:,mpt],0).sum(1)/sum(905-np.arange(1,21));np.testing.assert_allclose(chk,y.loc[el,'MP10'].iloc[:16],atol=1e-12)
   deliver(y,el,['MP10'],res)
 else:
  b=pd.read_parquet(v/'boundary_geometry.parquet').reindex(z.index);el=z.epi_eligible&~b.review_boundary&~ident.reliable_T
  if 'human_epi_eligible' in ident:el=ident.human_epi_eligible.combine_first(el).astype(bool)
  extra=d/'additional_scores.tsv.gz'
  if extra.exists():
   for key,col in pd.read_csv(extra,sep='\t').set_index('id').items():z.loc[col.index,key]=col
  common=z[STATES].notna().all(axis=1);assert (common|~el).all(),'newly eligible Epi needs panel transfer'
  y=z[['x_um','y_um','cell_area','nFeature','nGeneTranscripts']+STATES+[k+'_detected' for k in STATES[:2]]+[k+'_max_gene_fraction' for k in STATES[:2]]].rename(columns={'nGeneTranscripts':'UMI',**{k+'_max_gene_fraction':k+'_dominance' for k in STATES[:2]}}).copy()
  y['component']=pd.read_parquet(v/'graph_nodes.parquet').set_index('id').component;y.loc[~el,STATES]=np.nan;deliver(y,el,STATES,'cell')
 put(pd.DataFrame(qrows),d/'score_QC.tsv');put(pd.DataFrame(corr,columns=['sample','res','state','variable','Spearman']),d/'score_correlations.tsv');put(pd.DataFrame(overlaps,columns=['sample','state1','state2','intersection','union','Jaccard']),d/'high_overlap.tsv')
 print(s,qrows,flush=True)
if __name__=='__main__':
 for s in sys.argv[1:] or HD+XE:run(s)
