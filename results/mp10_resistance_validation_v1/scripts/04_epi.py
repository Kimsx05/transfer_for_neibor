exec(open(__file__.replace('04_epi.py','01_prepare.py')).read().split('# Archive contents')[0])
from scipy.stats import spearmanr,rankdata
import time
C=Path(cfg['inputs']['epi_cache']);genes=C.joinpath('genes.txt').read_text().splitlines();cells=C.joinpath('cells.txt').read_text().splitlines();assert len(genes)==cfg['epi']['background'];K=cfg['epi']['aucMaxRank'];assert K==int(np.ceil(.05*len(genes)))
mp=pd.read_csv(B/'00_audit/epi_signature_mapping.tsv',sep='\t');lookup=mp.set_index('gene').input_gene.to_dict();sets=cfg['genes'];un=list(dict.fromkeys(lookup[g] for gs in sets.values() for g in gs if pd.notna(lookup[g])));ix=pd.Index(genes).get_indexer(un);assert (ix>=0).all();rev=np.full(len(genes),-1,dtype=int);rev[ix]=np.arange(len(ix));inds=np.memmap(C/'indices.bin',dtype='<i4',mode='r');dat=np.memmap(C/'data.bin',dtype='<f8',mode='r');ptr=np.memmap(C/'indptr.bin',dtype='<i4',mode='r');assert len(ptr)==len(cells)+1 and ptr[-1]==len(dat)==len(inds)
f=pd.read_parquet(cfg['inputs']['epi_score']).set_index('cell_id').loc[cells].reset_index();assert f.cell_id.tolist()==cells
sample_keys=(f.dataset+'|'+f['orig.ident']);sample_levels=sorted(sample_keys.unique());sample_codes=pd.Categorical(sample_keys,categories=sample_levels).codes;gene_detected=np.zeros((len(sample_levels),len(un)),dtype=int)
scores=np.zeros((len(f),len(sets)));det=np.zeros_like(scores);sigix={s:np.array([un.index(lookup[g]) for g in gs if pd.notna(lookup[g])],dtype=int) for s,gs in sets.items()};coverage=[]
for s,gs in sets.items():coverage.append(dict(signature=s,n_original=len(gs),n_mapped=len(sigix[s]),n_present=len(sigix[s]),missing=';'.join(g for g in gs if pd.isna(lookup[g])),background=len(genes)))
assert len(sigix['MP10'])>=16
started=time.time();print('Streaming sparse counts',len(f),'cells',flush=True)
for i in range(len(f)):
 a,z=ptr[i:i+2];v=np.asarray(dat[a:z]);j=rev[inds[a:z]];vals=np.zeros(len(un));keep=j>=0;vals[j[keep]]=v[keep];sv=np.sort(v);gt=len(sv)-np.searchsorted(sv,vals,side='right');eq=np.searchsorted(sv,vals,side='right')-np.searchsorted(sv,vals,side='left');eq[vals==0]=len(genes)-np.count_nonzero(v);r=1+gt+(eq-1)/2
 gene_detected[sample_codes[i]]+=(vals>0)
 for k,(s,si) in enumerate(sigix.items()):scores[i,k]=np.maximum(0,K-r[si]).sum()/sum(K-jj for jj in range(1,len(si)+1));det[i,k]=np.count_nonzero(vals[si])
 if i%20000==0:print(i,'elapsed',round(time.time()-started,1),flush=True)
rows_gene=[]
for si,sample in enumerate(sample_levels):
 n=int((sample_codes==si).sum())
 for j,gene in enumerate(un):rows_gene.append(dict(sample=sample,dataset=sample.split('|')[0],input_gene=gene,n_cells=n,n_detected=int(gene_detected[si,j]),detection_rate=gene_detected[si,j]/n))
save(pd.DataFrame(rows_gene),'00_audit/epi_gene_detection_by_sample.tsv')
err=np.max(np.abs(scores[:,0]-f.Top20_score));assert err<1e-12,err
out=f[['cell_id','dataset','orig.ident','Type','nCount_RNA','nFeature_RNA']].copy();out['sample']=out.dataset+'|'+out['orig.ident'];pat=pd.read_csv(cfg['inputs']['patient_map'],sep='\t');out=out.merge(pat[['sample','patient','patient_recorded','patient_status','patient_source']],on='sample',how='left',validate='many_to_one');assert out.cell_id.tolist()==cells
for k,s in enumerate(sets):out[s]=f.Top20_score.to_numpy() if s=='MP10' else scores[:,k];out[s+'_n_detected']=det[:,k]
base=Path('/home/data/t070721/Metabolistic_scRNA/Multi2');cm=pd.read_csv(base/'3.3 recognize malignant/rough_cluster_infercnv_score_cell_metadata.csv',usecols=['cell','rough_cluster']);st=pd.read_csv(base/'3.4 malignant cell process/NMF1/rough_cluster_malignant_like_status_summary.csv');cm=cm.merge(st[['rough_cluster','malignant_like_status']],on='rough_cluster',validate='many_to_one');out=out.merge(cm.rename(columns={'cell':'cell_id'}),on='cell_id',how='left',validate='one_to_one');save(out,'01_epi/cell_scores.tsv.gz');save(pd.DataFrame(coverage),'00_audit/epi_signature_coverage.tsv')
rows=[]
for subset,md in [('tumor',out[out.Type.isin(cfg['epi']['tumor'])]),('normal',out[out.Type=='Normal']),('CNV_supported_tumor',out[out.Type.isin(cfg['epi']['tumor'])&out.malignant_like_status.eq('malignant_like_main')])]:
 for sample,g in md.groupby('sample',sort=True):
  for s in sets:
   if s=='MP10':continue
   ok=g[['MP10',s]].notna().all(axis=1);h=g[ok];valid=len(h)>=50 and h.MP10.nunique()>1 and h[s].nunique()>1;rho=spearmanr(h.MP10,h[s]).statistic if valid else np.nan
   base=dict(subset=subset,sample=sample,dataset=g.dataset.iloc[0],patient=g.patient.iloc[0],patient_recorded=g.patient_recorded.iloc[0],patient_status=g.patient_status.iloc[0],signature=s,n_cells=len(g),n_used=len(h),rho=rho,adjustment='unadjusted',MP10_sd=h.MP10.std(),signature_sd=h[s].std(),mean_detection=g[s+'_n_detected'].mean()/len(sigix[s]),reason='OK' if valid else 'n<50_or_constant_or_missing')
   rows.append(base)
   for depth in ['nCount_RNA','nFeature_RNA']:
    row=dict(base,adjustment=depth+'+S.Score+G2M.Score',rho=np.nan,reason='BLOCKED: S.Score and G2M.Score absent; no substitutes calculated');rows.append(row)
eff=pd.DataFrame(rows);save(eff,'01_epi/sample_effects.tsv')
# No verified identities: required patient effects remain explicit NA, not falsely independent specimens.
pr=eff.copy();pr['patient_effect']=np.nan;pr['reason']='BLOCKED: all source identities derived/provisional; independently verified patient crosswalk required';save(pr[['subset','dataset','patient_recorded','patient_status','signature','adjustment','patient_effect','reason']].drop_duplicates(),'01_epi/patient_effects_UNAVAILABLE.tsv')
# Descriptive sample summaries carry no patient-level CI/P.
raw=eff[eff.adjustment=='unadjusted'];summary=raw.groupby(['subset','dataset','signature']).agg(n_samples=('rho','count'),median_sample_rho=('rho','median'),min_rho=('rho','min'),max_rho=('rho','max'),n_positive=('rho',lambda x:(x>0).sum()),n_negative=('rho',lambda x:(x<0).sum())).reset_index();summary['unit']='specimen descriptive; patient independence unverified';summary['P']=np.nan;summary['CI_low']=np.nan;summary['CI_high']=np.nan;save(summary,'01_epi/descriptive_sample_summary.tsv')
qc=dict(n_cells=len(out),n_samples=out['sample'].nunique(),n_datasets=out.dataset.nunique(),MP10_reproduction_max_abs_error=float(err),seconds=time.time()-started,verified_patients=int(out.patient.nunique()),patient_primary_status='BLOCKED',adjustment_status='BLOCKED',n_CNV_supported_tumor=int((out.Type.isin(cfg['epi']['tumor'])&out.malignant_like_status.eq('malignant_like_main')).sum()))
(B/'01_epi/execution_QC.json').write_text(json.dumps(qc,indent=2));print(qc,flush=True)
