from common import *
from scipy.stats import rankdata
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
raw='EpiMP10_Top20_AUCell'
support=['H05_support','H05_conditional_support']
keys=['section_id','barcode']
metrics=[]; diagnostics=[]; ablation=[]; locations=[]; regions=[]; non_epi=[]; thresholds=[]; denominator=[]; label_context=[]; conditional_pairs=[]; outside=[]
def corr(x,y,rank=False):
 x=np.asarray(x); y=np.asarray(y)
 if len(x)<3 or np.ptp(x)==0 or np.ptp(y)==0:return np.nan
 if rank:x,y=rankdata(x),rankdata(y)
 return np.corrcoef(x,y)[0,1]
def residual(x,Z):return x-Z@np.linalg.lstsq(Z,x,rcond=None)[0]
def attach(d,p):
 a=read(p);assert len(a)==len(d) and set(zip(*(a[k] for k in keys)))==set(zip(*(d[k] for k in keys)))
 return d.merge(a,on=keys,validate='one_to_one')
for sec in S:
 d=read(O/f'spots/{sec}.input.tsv.gz')
 base=d.copy()
 for branch in ['main','noTop20']:
  d=attach(base,O/f'spots/{sec}.{branch}.tsv.gz')
  e=d.loc[d['query']].copy();e['logdepth']=np.log1p(e.nCount_RNA)
  pair=e[[raw,*support,'Epithelial_fraction','logdepth']].dropna()
  pair_Z=np.column_stack([np.ones(len(pair)),pair[['Epithelial_fraction','logdepth']]])
  for score in support:
   conditional_pairs.append(dict(section_id=sec,branch=branch,score=score,n=len(pair),pearson=corr(pair[raw],pair[score]),spearman=corr(pair[raw],pair[score],True),resid_pearson=corr(residual(pair[raw],pair_Z),residual(pair[score],pair_Z))))
   complete=e[[raw,score,'Epithelial_fraction','logdepth']].replace([np.inf,-np.inf],np.nan).dropna()
   x,y=complete[raw],complete[score]
   Z=np.column_stack([np.ones(len(x)),complete[['Epithelial_fraction','logdepth']]])
   rx,ry=residual(x,Z),residual(y,Z)
   r=dict(section_id=sec,branch=branch,score=score,n_query=len(e),n=len(x),missing=len(e)-len(x),low_den_query=int(e.low_denominator.sum()),low_den_effective=e.loc[x.index,'low_denominator'].sum(),zero_den=(e.Epi_support==0).sum(),pearson=corr(x,y),spearman=corr(x,y,True),resid_pearson=corr(rx,ry),resid_spearman=corr(rx,ry,True))
   metrics.append(r)
   for cov in ['Epithelial_fraction','logdepth']:
    z=complete[cov]
    z1=np.column_stack([np.ones(len(x)),z])
    diagnostics.append(dict(section_id=sec,branch=branch,score=score,covariate=cov,n=len(x),MP10_vs_cov_Pearson=corr(x,z),support_vs_cov_Pearson=corr(y,z),residual_r_adjust_only_this=corr(residual(x,z1),residual(y,z1)),joint_residual_r=r['resid_pearson']))
  if branch=='noTop20':
   for score in support:
    t=e[[*keys,score,raw]].merge(main.loc[main['query'],[*keys,score]],on=keys,suffixes=('_masked','_main'),validate='one_to_one').dropna()
    fm=read(O/f'features/{sec}_main.tsv');fa=read(O/f'features/{sec}_noTop20.tsv')
    assert set(fa.gene)==set(fm.gene)-set(TOP)
    ablation.append(dict(section_id=sec,score=score,n_paired=len(t),main_features=len(fm),masked_features=len(fa),main_MP10_genes=','.join(fm[fm.MP10_member].gene),masked_MP10_genes=int(fa.MP10_member.sum()),main_r=corr(t[raw],t[score+'_main']),masked_r=corr(t[raw],t[score+'_masked']),support_rank_agreement=corr(t[score+'_main'],t[score+'_masked'],True)))
   continue
  main=d.copy()
  for population,g in [('full',d),('frozen_query',e)]:
   for group,b in [('all',g),('Epi_support_lt_0.10',g[g.low_denominator]),('Epi_support_ge_0.10',g[~g.low_denominator])]:
    valid=b.H05_conditional_support.notna()
    denominator.append(dict(section_id=sec,population=population,stratum=group,n=len(b),n_zero=int((b.Epi_support==0).sum()),n_ratio_valid=valid.sum(),n_conditional_ge_0_9=int((b.H05_conditional_support>=.9).sum()),median_raw=b.H05_support.median(),median_den=b.Epi_support.median(),raw_vs_Epi_support_r=corr(b.H05_support,b.Epi_support),MP10_vs_Epi_support_r=corr(b[raw],b.Epi_support)))
   for (dominant,predicted),h in g.groupby(['dominant_c2l_major','predicted_label']):
    label_context.append(dict(section_id=sec,population=population,c2l_dominant=dominant,predicted_label=predicted,n=len(h)))
  m25,m75=e[raw].quantile([.25,.75]);h25,h75,h95=e.H05_support.quantile([.25,.75,.95])
  hi=d['query']&(d[raw]>=m75)&(d.H05_support<=h25);lo=d['query']&(d[raw]<=m25)&(d.H05_support>=h75)
  d['discordance']=np.select([hi,lo],['MP10 high / support low','MP10 low / support high'],default='other')
  d['display_top5_support']=d['query']&(d.H05_support>=h95)
  thresholds.append(dict(section_id=sec,MP10_Q25=m25,MP10_Q75=m75,raw_support_Q25=h25,raw_support_Q75=h75,raw_support_Q95=h95,n_desc_top=d.display_top5_support.sum()))
  d['discordant_component']=pd.Series(pd.NA,index=d.index,dtype='Int64')
  d['top20_UMI']=d[['UMI_'+g for g in TOP]].sum(axis=1,min_count=1)
  for mask,name in [(hi,'MP10 high / support low'),(lo,'MP10 low / support high')]:
   b=d.loc[mask].copy();lookup={(r.array_row,r.array_col):i for i,r in enumerate(b.itertuples())}
   ixs=[];jxs=[]
   for (r,c),i in lookup.items():
    for dr,dc in [(0,2),(1,1),(1,-1)]:
     j=lookup.get((r+dr,c+dc))
     if j is not None and b.tissue_component_id.iloc[i]==b.tissue_component_id.iloc[j]:ixs.append(i);jxs.append(j)
   _,lab=connected_components(coo_matrix((np.ones(len(ixs)),(ixs,jxs)),shape=(len(b),len(b))),directed=False)
   b['component']=lab+1;d.loc[b.index,'discordant_component']=lab+1
   for component,g in b.groupby('component'):
    regions.append(dict(section_id=sec,discordance=name,component=component,n=len(g),x_px_center=g.x.mean(),y_px_center=g.y.mean(),array_row_min=g.array_row.min(),array_row_max=g.array_row.max(),array_col_min=g.array_col.min(),array_col_max=g.array_col.max(),modal_position=g.lesion_location.mode().iat[0] if g.lesion_location.notna().any() else 'unknown',median_epi_fraction=g.Epithelial_fraction.median(),median_top20_UMI=g.top20_UMI.median()))
  for name,g in d[d['query']].groupby('discordance'):
   entry=dict(section_id=sec,discordance=name,n=len(g),median_depth=g.nCount_RNA.median(),median_epi_fraction=g.Epithelial_fraction.median(),median_Epi_support=g.Epi_support.median(),n_inside=(g.lesion_location=='inside').sum(),median_top20_UMI=g.top20_UMI.median())
   for c in d.filter(regex='^UMI_'):
    entry[c+'_median']=g[c].median(); entry[c+'_detected']=(g[c]>0).mean() if g[c].notna().all() else np.nan
   locations.append(entry)
  below=d[(d.Epithelial_fraction<.25)&(d.dominant_c2l_major!='Epithelial')]
  for name,g in [('ALL',below),*below.groupby('dominant_c2l_major')]:
   row=dict(section_id=sec,c2l_dominant=name,n=len(g),median_H05=g.H05_support.median(),Q95_H05=g.H05_support.quantile(.95),max_H05=g.H05_support.max())
   for marker in ['EPCAM','KRT8','COL1A1','PTPRC','PECAM1']:row[marker+'_median_UMI']=g['UMI_'+marker].median()
   outside.append(row)
  non=d[(d.Epithelial_fraction<.1)&(d.dominant_c2l_major!='Epithelial')]
  for name,g in [('ALL',non),*non.groupby('dominant_c2l_major')]:
   entry=dict(section_id=sec,c2l_major=name,n=len(g),median_H05=g.H05_support.median(),Q95_H05=g.H05_support.quantile(.95),max_H05=g.H05_support.max(),median_Epi_support=g.Epi_support.median(),n_low_den=g.low_denominator.sum())
   for gene in ['EPCAM','KRT8','PTPRC','COL1A1','PECAM1']:entry[gene+'_UMI_median']=g['UMI_'+gene].median()
   non_epi.append(entry)
  save(d,f'spots/{sec}.tsv.gz')
save(outside,'below_legacy_gate_diagnostics.tsv');save(conditional_pairs,'common_support_comparison.tsv');save(denominator,'denominator_diagnostics.tsv');save(label_context,'label_context.tsv');save(metrics,'section_metrics.tsv');save(diagnostics,'composition_depth_diagnostics.tsv');save(ablation,'feature_ablation.tsv');save(regions,'discordant_regions.tsv');save(locations,'discordant_expression.tsv');save(non_epi,'non_epi_diagnostics.tsv');save(thresholds,'display_thresholds.tsv')
if all((O/f'pseudo_{d}.tsv.gz').exists() for d in ['medium','low']):
 truth=read(O/'pseudo_truth.tsv')
 pred=pd.concat([read(O/f'pseudo_{d}.tsv.gz') for d in ['medium','low']]).rename(columns={'barcode':'pseudospot_id'})
 out=truth.merge(pred,on='pseudospot_id',validate='one_to_one');assert len(out)==len(truth)==len(pred)
 save(out,'pseudo_spot_predictions.tsv.gz')
 rows=[]
 for depth,d in out.groupby('depth_condition'):
  for src,g in [('all_reused_mixtures',d),*d.groupby('sample_uid')]:
   zero=g[g.n_H05==0];pos=g[g.n_H05>0]
   rows.append(dict(depth=depth,source=src,n=len(g),n_true_zero=len(zero),n_positive=len(pos),zero_support_median=zero.H05_support.median(),zero_support_Q95=zero.H05_support.quantile(.95),zero_support_max=zero.H05_support.max(),positive_support_median=pos.H05_support.median(),fraction_vs_raw_spearman=corr(g.H05_fraction,g.H05_support,True),count_vs_raw_spearman=corr(g.n_H05,g.H05_support,True),median_depth=g.post_shared_UMI.median()))
 save(rows,'pseudo_spot_validation.tsv')
 w=out.pivot(index=['sample_uid','mixture_id'],columns='depth_condition',values='H05_support').reset_index()
 save([dict(source=src,n=len(g),paired_median_low_minus_medium=(g.low-g.medium).median(),paired_rank_agreement=corr(g.medium,g.low,True)) for src,g in w.groupby('sample_uid')],'pseudo_paired_depth.tsv')
else:
 save([{'status':'unavailable','reason':'See logs/pseudo*failure.txt'}],'pseudo_spot_validation.tsv')

