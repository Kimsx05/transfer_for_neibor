from common import *
import matplotlib.pyplot as plt

def main(out):
 P=out/'01_score_fidelity';X,genes,m=epi_cache();u=pd.read_csv(P/'original_usage.tsv',sep='\t');assert u.cell_id.is_unique and m.cell_id.is_unique;m=m.merge(u,on='cell_id',how='left',validate='one_to_one');assert len(m)==X.shape[0];ix=pd.Index(genes).get_indexer(TOP);assert (ix>=0).all();k=int(np.ceil(.05*len(genes)));den=sum(k-j for j in range(1,21));C=np.zeros((len(m),20));R=np.zeros_like(C);sig=X[:,ix].toarray();depth=np.asarray(X.sum(1)).ravel();nf=np.diff(X.indptr)
 for i in range(len(m)):
  vals=np.sort(X.data[X.indptr[i]:X.indptr[i+1]]);v=sig[i];gt=len(vals)-np.searchsorted(vals,v,side='right');eq=np.searchsorted(vals,v,side='right')-np.searchsorted(vals,v,side='left');eq[v==0]=len(genes)-len(vals);rank=1+gt+(eq-1)/2;R[i]=rank;C[i]=np.maximum(0,k-rank)/den
 m['Top20_score']=C.sum(1);m.loc[depth==0,'Top20_score']=np.nan;m['nCount']=depth;m['nFeature']=nf;m['n_signature_detected']=(sig>0).sum(1);m['usage_type']='original_stored_usage_raw_fixed14_basis_NNLS';m['analysis_unit']=m['orig.ident'].astype(str);m['identity_scope']='source_sample_not_verified_patient';m['malignant_status']='not_independently_verified';m.to_parquet(safe(out,P/'cell_level_fidelity.parquet'),index=False);np.savez_compressed(safe(out,P/'top20_gene_diagnostics.npz'),counts=sig,rank=R,auc_contribution=C,genes=np.array(TOP),cell_id=m.cell_id.to_numpy(dtype=str));rows=[]
 for (dataset,unit),g in m.groupby(['dataset','analysis_unit']):
  g=g.dropna(subset=['Top20_score','MP10_usage']);valid=len(g)>=30 and g.Top20_score.nunique()>1 and g.MP10_usage.nunique()>1;r=spearmanr(g.Top20_score,g.MP10_usage).statistic if valid else np.nan
  if valid:
   D=np.column_stack([np.ones(len(g)),np.log1p(g.nCount),np.log1p(g.nFeature)]);yr=rankdata(g.Top20_score);ur=rankdata(g.MP10_usage);yr-=D@np.linalg.lstsq(D,yr,rcond=None)[0];ur-=D@np.linalg.lstsq(D,ur,rcond=None)[0];adj=np.corrcoef(yr,ur)[0,1]
  else:adj=np.nan
  rows.append(dict(dataset=dataset,analysis_unit=unit,n=len(g),spearman=r,technical_partial_rank_correlation=adj,score_depth_spearman=spearmanr(g.Top20_score,g.nCount).statistic if valid else np.nan,usage_depth_spearman=spearmanr(g.MP10_usage,g.nCount).statistic if valid else np.nan,status='CALCULABLE' if valid else 'LOW_N_OR_CONSTANT',identity_scope='source_sample_not_verified_patient'))
 write_tsv(out,P/'fidelity_by_patient_dataset.tsv',rows);t=pd.DataFrame(rows);v=m.dropna(subset=['Top20_score','MP10_usage']);write_json(out,P/'fidelity_summary.json',dict(n_counts=len(m),n_usage_matched=len(v),n_genes=len(genes),n_source_units=len(t),overall_spearman=spearmanr(v.Top20_score,v.MP10_usage).statistic,donor_balanced_median_spearman=t.spearman.median(),donor_balanced_mean_spearman=t.spearman.mean(),n_negative_units=int(t.spearman.lt(0).sum()),known_patient_scope='orig.ident sources; no verified patient merge claimed',reference_overlap='check calibration overlap audit',original_stored_usage=True))
 write_tsv(out,P/'usage_provenance.tsv',[dict(path=str(ROOT/'ST/mapping/08_cell_EpiMP_scores.qs'),type='original_stored_usage',column='usage_raw[EpiMP10,]',definition='original14basis fixed-W NNLS, positive-centered RNA per sample',new_projection=False,independence='Top20 belongs to same original basis; association is not independent biological validation')]);write_tsv(out,P/'coverage_and_missing_ids.tsv',[dict(count_cells=len(m),usage_cells=len(u),matched=len(v),counts_without_usage=int(m.MP10_usage.isna().sum()),usage_without_counts=int((~u.cell_id.isin(m.cell_id)).sum()))])
 rng=np.random.default_rng(123);display=[];official=[]
 for _,g in m.groupby('analysis_unit'):
  g=g.sort_values('nCount');display.extend(g.iloc[np.unique(np.linspace(0,len(g)-1,min(100,len(g))).astype(int))].index);official.extend(g.iloc[np.unique(np.linspace(0,len(g)-1,min(3,len(g))).astype(int))].index)
 fig,ax=plt.subplots(figsize=(7,5));z=m.loc[display];ax.scatter(z.MP10_usage,z.Top20_score,s=2,alpha=.25);ax.set(xlabel='Original stored EpiMP10 raw usage',ylabel='Frozen Top20 average-tie score',title='Source-balanced display; all matched cells in summaries');fig.tight_layout()
 for ext in ['png','pdf']:fig.savefig(safe(out,out/'06_figures'/('score_fidelity.'+ext)),transparent=True,dpi=150)
 plt.close(fig)
 for dataset,g in t.groupby('dataset'):
  fig,ax=plt.subplots(figsize=(8,4));ax.scatter(range(len(g)),g.spearman);ax.axhline(0,c='gray');ax.set(xlabel='source unit (not independent cells)',ylabel='within-unit Spearman',title=str(dataset));fig.tight_layout()
  for ext in ['png','pdf']:fig.savefig(safe(out,out/'06_figures'/str(dataset)/('fidelity.'+ext)),transparent=True,dpi=140)
  plt.close(fig)
 # Official comparison subset, fixed source/depth selection before comparison.
 from scipy.io import mmwrite
 official=np.array(official,dtype=int);mmwrite(safe(out,P/'official_input.mtx'),X[official].T);write_tsv(out,P/'official_axes.tsv',m.iloc[official][['cell_id','Top20_score','nCount','nFeature']]);safe(out,P/'official_genes.txt').write_text('\n'.join(genes)+'\n');safe(out,P/'top20.txt').write_text('\n'.join(TOP)+'\n')
 safe(out,P/'score_fidelity_report.md').write_text('Complete epithelial counts scored on one fixed90011-gene background; exact canonical cell ids; original stored raw usage only. Source-unit correlations and technical partial-rank diagnostics are descriptive. No new malignant labels or NMF. Official tie comparison is a separate implementation sensitivity.\n');log('full epithelial fidelity complete')
if __name__=='__main__':main(OUT)
