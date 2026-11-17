from common import *
from sklearn.metrics import roc_auc_score
import importlib.util
if (R/'frozen_config.yaml').exists():print('FROZEN_ALREADY');raise SystemExit
assert all((R/s/'prepared.json').exists() for s in SECTIONS)
sig=json.loads((R/'models/signatures.json').read_text());hd=[pd.read_csv(R/'tables'/f'{s}_gene_mapping.tsv',sep='\t').symbol.tolist()for s in SECTIONS];common=set.intersection(*[set(v) for v in hd]);sig={k:[g for g in v if g in common]for k,v in sig.items()};assert all(len(v)>=2 for v in sig.values());(R/'models/signatures.json').write_text(json.dumps(sig,indent=2))
qc=[];med=[]
for s in SECTIONS:
 for kind in ['T32','T16','Epi64']:
  d=pd.read_parquet(R/s/f'{kind}_measure.parquet');v=d[d.T_detected>=2] if kind.startswith('T') else d
  qc.append(dict(section=s,branch=kind,n=len(d),T_supported_n=int((d.T_detected>=2).sum()),median_nFeature=float(v.nFeature.median()),median_nCount=float(v.nCount.median()),feature_q10=float(v.nFeature.quantile(.1)),feature_q90=float(v.nFeature.quantile(.9))))
  if kind=='T32':med.append(float(v.nFeature.median()))
# One cross-section setting: minimum section median among T-supported32um regions, rounded to100.
maxrank=int(np.clip(np.floor(min(med)/100)*100,500,3000));minfeature=max(100,maxrank//4)
put(pd.DataFrame(qc),'measurement_calibration.tsv')
rf=pd.read_parquet(R/'tables/reference_full_ranks.parquet');validation=[];scoretables=[]
for dep in ['full',500,2000,5000]:
 d=rf.copy()if dep=='full' else pd.read_parquet(R/'tables'/f'reference_holdout_{dep}.parquet')
 for k,v in sig.items():d[k]=ucell(d,v,maxrank)
 if dep=='full':
  train=d[d.partition.eq('train')];neg=train[train.label.str.startswith(('CD4','CD8','Tfh'))];threshold=float(neg.Treg_main.quantile(.95));vv=d[d.partition.ne('train')]
 else:vv=d
 for label in ['all_nonTreg','T_competitors','CD4 Helper/Conventional','CD4 Naive/Memory-like','Tfh/CXCL13-like','CD8 Exhausted','NK','Epithelial']:
  if label=='all_nonTreg':mask=np.ones(len(vv),bool)
  elif label=='T_competitors':mask=vv.label.eq('Treg')|vv.label.str.startswith(('CD4','CD8','Tfh'))
  else:mask=vv.label.isin(['Treg',label])
  v=vv.loc[mask];truth=v.label.eq('Treg');
  for name in ['Treg_main','Treg_literature']:
   validation.append(dict(depth=str(dep),comparison=label,signature=name,n=len(v),Treg_n=int(truth.sum()),AUROC=roc_auc_score(truth,v[name]) if truth.nunique()==2 else np.nan))
 if dep=='full':
  valid=d[d.partition.ne('train')]
  for g in sig['Treg_main']:
   keep=[k for k in sig['Treg_main']if k!=g];v=ucell(valid,keep,maxrank);validation.append(dict(depth='full_leave_one_gene',comparison='all_nonTreg',signature='minus_'+g,n=len(valid),Treg_n=int(valid.label.eq('Treg').sum()),AUROC=roc_auc_score(valid.label.eq('Treg'),v)))
 d[['cell_id','label','sample']+list(sig)].assign(depth=str(dep)).to_parquet(R/'tables'/f'reference_scores_{dep}.parquet',index=False)
put(pd.DataFrame(validation),'signature_validation.tsv')
mix=pd.read_parquet(R/'tables/reference_mixture_ranks.parquet')
for k,v in sig.items():mix[k]=ucell(mix,v,maxrank)
mix['exceeds_training_T95']=mix.Treg_main>threshold
put(mix.groupby(['competitor','depth','Treg_RNA_fraction']).agg(n=('Treg_main','size'),score_median=('Treg_main','median'),score_q10=('Treg_main',lambda x:x.quantile(.1)),score_q90=('Treg_main',lambda x:x.quantile(.9)),above_T95=('exceeds_training_T95','mean')).reset_index(),'reference_mixture_validation.tsv')
put(mix[['competitor','depth','Treg_RNA_fraction','rep','nFeature']+list(sig)+['exceeds_training_T95']],'reference_mixture_scores.tsv.gz')
# Score calibration only: no paired Treg/MP10 neighborhoods, effects or correlation examined here.
stability=[];roi=[];variogram=[];rng=np.random.default_rng(123)
for s in SECTIONS:
 for kind in ['T32','T16','Epi64']:
  d=pd.read_parquet(R/s/f'{kind}_measure.parquet');rank=d[[c for c in d if c.startswith('rank_')]].rename(columns=lambda c:c[5:])
  for k,v in sig.items():d[k]=ucell(rank,v,maxrank)
  d['MP10']=mp10(rank);d['Treg_main_detected']=(d[['raw_'+g for g in sig['Treg_main']]]>0).sum(1);d['Treg_main_UMI']=d[['raw_'+g for g in sig['Treg_main']]].sum(1);d['Treg_core_detected']=(d[['raw_'+g for g in ['FOXP3','IL2RA','CTLA4','IKZF2']]]>0).sum(1)
  d['measurement_eligible']=(d.nFeature>=minfeature)&(d.T_detected>=2)&(d.T_UMI>=3) if kind.startswith('T') else (d.nFeature>=905)&(d.MP10_detected>=5)&(d.Epi_detected>=2)
  d['evidence']=np.where(~d.measurement_eligible,'insufficient',np.where((d.Treg_main_detected>=3)&(d.Treg_core_detected>=1)&(d.Treg_main>threshold),'Treg_like_multigene','T_lineage_or_activated')) if kind.startswith('T') else np.where(d.measurement_eligible,'Epi_enriched_scoreable','insufficient')
  d['mixed_Epi_Fib']=(d.Epi_detected>=2)&(d.Fib_detected>=2)
  half=pd.read_parquet(R/s/f'{kind}_halfdepth.parquet').set_index('region_id');q=d.set_index('region_id').loc[half.index];hr=half.filter(like='rank_').rename(columns=lambda c:c[5:]);hs=ucell(hr,sig['Treg_main'],maxrank) if kind.startswith('T') else mp10(hr);primary=q.Treg_main if kind.startswith('T')else q.MP10;valid=q.measurement_eligible.to_numpy()
  stability.append(dict(section=s,branch=kind,n_calibration=len(q),eligible_calibration=int(valid.sum()),halfdepth_spearman=float(spearmanr(primary[valid],hs[valid]).statistic) if valid.sum()>5 else np.nan,score_count_spearman=float(spearmanr(d.loc[d.measurement_eligible,'Treg_main' if kind.startswith('T')else 'MP10'],d.loc[d.measurement_eligible,'nCount']).statistic),eligible=int(d.measurement_eligible.sum())))
  pd.DataFrame({'region_id':q.index,'full':primary.to_numpy(),'half':hs,'eligible':valid}).to_parquet(R/s/f'{kind}_depth_validation.parquet',index=False)
  if kind=='Epi64':
   v=d.loc[d.measurement_eligible,'MP10'];quant=v.quantile([.5,.9,.95]);d['L50']=d.measurement_eligible&(d.MP10<=quant.loc[.5]);d['H50']=d.measurement_eligible&(d.MP10>quant.loc[.5]);d['H10']=d.measurement_eligible&(d.MP10>quant.loc[.9]);d['H05']=d.measurement_eligible&(d.MP10>quant.loc[.95])
  # Marginal spatial autocorrelation diagnostic only, independent random pairs per distance band.
  if kind in ['T32','Epi64']:
   v=d[d.measurement_eligible];coords=v[['x_um','y_um']].to_numpy();sc=rankdata(v.Treg_main if kind=='T32'else v.MP10);sc=(sc-sc.mean())/max(sc.std(),1e-9);tree=cKDTree(coords);ix=rng.choice(len(v),min(1200,len(v)),False)
   for radius in [100,200,400,800]:
    products=[]
    for i in ix:
     js=tree.query_ball_point(coords[i],radius);js=[j for j in js if j!=i and v.component.iloc[j]==v.component.iloc[i] and np.linalg.norm(coords[j]-coords[i])>radius/2]
     if js:products.extend((sc[i]*sc[rng.choice(js,min(20,len(js)),False)]).tolist())
    variogram.append(dict(section=s,branch=kind,distance_band_outer_um=radius,n_pairs=len(products),mean_standardized_product=np.mean(products) if products else np.nan))
  if kind=='T32':
   # Three Treg signal quantiles + one T-matched lower-signature control, before MP10 pairing.
   v=d[d.measurement_eligible].copy();chosen=[]
   if len(v):
    for qt in [.25,.5,.9]:
     order=(v.Treg_main-v.Treg_main.quantile(qt)).abs().sort_values().index
     i=next((i for i in order if all(np.hypot(d.loc[i,'x_um']-d.loc[j,'x_um'],d.loc[i,'y_um']-d.loc[j,'y_um'])>=256 for j in chosen)),order[0]);chosen.append(i)
    high=d.loc[chosen[-1]];control=v[v.Treg_main<=v.Treg_main.median()];i=(control.T_lineage-high.T_lineage).abs().idxmin();chosen.append(i)
    for j,i in enumerate(chosen):
     v0=d.loc[i];roi.append(dict(section=s,roi=f'roi_{j+1}',region_id=v0.region_id,x0=v0.x_um-128,x1=v0.x_um+128,y0=v0.y_um-128,y1=v0.y_um+128,selection=['Treg_score_q25','Treg_score_q50','Treg_score_q90','T_lineage_matched_lower_Treg'][j]))
  d.drop(columns=[c for c in d if c.startswith('rank_')]).to_parquet(R/s/f'{kind}_regions.parquet',index=False)
put(pd.DataFrame(stability),'score_stability.tsv');put(pd.DataFrame(roi),'roi_manifest.tsv');vg=pd.DataFrame(variogram);put(vg,'marginal_spatial_correlation.tsv')
# A common block size responds only to marginal dependence at400um, never to cross-score association.
block=800 if (vg.loc[vg.distance_band_outer_um==400,'mean_standardized_product']>.2).any() else 400
cfg=yaml.safe_load((R/'frozen_design.yaml').read_text());cfg.update({'phase':'frozen_before_HD_association','frozen_at':datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),'UCell_version':'2.14.0','UCell_maxRank':maxrank,'T_min_nFeature':minfeature,'T_lineage_detected_min':2,'T_lineage_UMI_min':3,'Treg_like_main_detected_min':3,'Treg_like_core_detected_min':1,'Treg_like_threshold_training_competitor_q95':threshold,'Epi_min_nFeature':905,'Epi_min_MP10_detected':5,'signatures':sig,'block_um':block,'larger_block_um':2*block,'association_min_regions':50,'association_min_blocks':12,'circle_observed_fraction_min':.5,'minimum_independent_Epi_regions':1,'literature_source':'De Simone2016 Table S4 full309 genes; modified exclusions and common measured genes','projection_or_smoothing':False})
(R/'frozen_config.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False));print('FROZEN',maxrank,minfeature,block,threshold,flush=True)
