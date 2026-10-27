from common import *
from scipy import sparse,io
f=P/'neighborhood_mp10_calibration_v5/01_score_fidelity/counts_cache'
cells=(f/'cells.txt').read_text().splitlines();genes=(f/'genes.txt').read_text().splitlines()
ix=np.memmap(f/'indices.bin',mode='r',dtype='<i4');data=np.memmap(f/'data.bin',mode='r',dtype='<f8');ind=np.fromfile(f/'indptr.bin',dtype='<i4')
x=sparse.csr_matrix((data,ix,ind),shape=(len(cells),len(genes)),copy=False)
v=read('score_validation_cells.tsv');new=io.mmread(O/'objects/validation_counts.mtx').tocsr().T.tocsr();old=x[pd.Index(cells).get_indexer(v.barcode)]
diff=new-old;errs=np.asarray(abs(diff).max(axis=1).toarray()).ravel();calc=read('score_recalculation.tsv');v=v.merge(calc,on='barcode');v['count_max_abs_error']=errs;v['score_abs_error']=abs(v.MP10_cell-v.filled_score)
save(v,'score_validation.tsv');assert np.max(errs)==0 and v.score_abs_error.max()<1e-12
m=read('cells_pre_score.tsv.gz');s=read('samples_base.tsv');fill=calc.set_index('barcode').filled_score
sel=m.major.eq('Epithelial')&m.MP10_cell.isna();m.loc[sel,'MP10_cell']=m.loc[sel,'barcode'].map(fill);m.loc[sel&m.MP10_cell.notna(),'score_source']='filled_same_definition'
m['score_unavailable_reason']=np.where(m.major.eq('Epithelial')&m.MP10_cell.isna(),'zero_counts_or_unmeasured','')
e=m[m.major.eq('Epithelial')];g=e.groupby('sample_uid');agg=g.MP10_cell.agg(n_epi_scored='count',MP10_sample_allEpi='mean',MP10_median_allEpi='median');agg['n_epi_total']=g.size()
for src in ['reused','filled_same_definition','unavailable']:agg['n_score_'+src]=e[e.score_source.eq(src)].groupby('sample_uid').size().reindex(agg.index,fill_value=0)
s=s.merge(agg,on='sample_uid',how='left');s['score_coverage']=s.n_epi_scored/s.n_epi_total;s['score_coverage_trusted']=s.score_coverage.eq(1)
s['MP10_eligible']=s.tissue_status.eq('Tumor')&s.n_epi_scored.ge(30)&s.score_coverage_trusted
s['MP10_threshold_allEpi']=s.dataset.map(s[s.MP10_eligible].groupby('dataset').MP10_sample_allEpi.median());s['sample_MP10_group_allEpi']=None
s.loc[s.MP10_eligible,'sample_MP10_group_allEpi']=np.where(s.loc[s.MP10_eligible,'MP10_sample_allEpi']>s.loc[s.MP10_eligible,'MP10_threshold_allEpi'],'high','low')
m['epi_MP10_state_within_specimen']=None;med=m.sample_uid.map(s.set_index('sample_uid').MP10_median_allEpi);z=m.major.eq('Epithelial')&m.tissue_status.eq('Tumor')&m.MP10_cell.notna()
m.loc[z,'epi_MP10_state_within_specimen']=np.where(m.loc[z,'MP10_cell']>med[z],'Epi_MP10_high','Epi_MP10_low')
for state in ['Epi_MP10_low','Epi_MP10_high']:s['n_'+state]=s.sample_uid.map(m[m.epi_MP10_state_within_specimen.eq(state)].groupby('sample_uid').size()).fillna(0).astype(int)
s['epi_state_eligible']=s.score_coverage_trusted&s.n_Epi_MP10_low.ge(30)&s.n_Epi_MP10_high.ge(30)&s.tissue_status.eq('Tumor')
s['score_definition']='frozen average-tie ranking of RNA counts over 90011 genes; K4501; denominator89810'
s['patient_unit_note']='recorded IDs only; all derived/provisional; not verified biological identities'
save(s,'samples_v2.tsv');save(m,'cells_v2.tsv.gz');save(m[m.major.eq('Epithelial')],'epi_MP10_cells_v2.tsv.gz')
old=pd.read_csv(V/'tables/samples.tsv',sep='\t');cols=['sample_uid','n_malignant','n_MP10_scored','MP10_sample','MP10_eligible','MP10_threshold','MP10_group'];old=old[cols].rename(columns={c:'old_'+c for c in cols if c!='sample_uid'})
vs=s.merge(old,on='sample_uid');vs['reentered']=vs.MP10_eligible&~vs.old_MP10_eligible;vs['group_flipped']=vs.MP10_eligible&vs.old_MP10_eligible&vs.sample_MP10_group_allEpi.ne(vs.old_MP10_group);save(vs,'old_vs_new_samples.tsv')
print(s.groupby('dataset').agg(specimens=('sample_uid','size'),eligible=('MP10_eligible','sum'),state=('epi_state_eligible','sum'),filled=('n_score_filled_same_definition','sum')).to_string());print('reentries',vs.reentered.sum(),'flips',vs.group_flipped.sum())
