from shared import *
import warnings
warnings.filterwarnings('ignore',category=pd.errors.PerformanceWarning)
mapping=read(O/'annotation_mapping.tsv').query('subtype!="Epithelial"')
fine=list(mapping.subtype); majors=['major_'+x for x in mapping.major.unique()]
sets=read(O/'gene_sets.tsv')
functions=[f'{mod}__{v}' for mod in sets.module.unique() for v in ['full','minus20']]
cov=read(T/'gene_coverage_base.tsv')
measures=['pearson','spearman','resid_pearson','resid_spearman']
full=[];cross=[];sens=[];fc=[];ft=[];groups=[];contr=[];disc=[];niches=[];diag=[];ranges=[];old=[];rings=[];lowtme=[];overlap=[];qcs=[]
ringroot=R/'neighborhood_ring_comparison_v2'
anchor=ringroot/'03_features'
ordered=pd.read_parquet(anchor/'common_anchors.parquet',columns=KEY+['spot_uid'])
ringcols=json.loads((anchor/'feature_columns.json').read_text())['q1']
archive=np.load(anchor/'features.npz')
region=read(ringroot/'00_preflight/coverage63.tsv').set_index('section_id')
def group(s):
 return 'high' if s else 'low'
def bands(x):
 q25,q75=x.quantile([.25,.75])
 if q25==q75:return pd.Series('overlap',index=x.index),q25,q75
 return pd.Series(np.select([x>=q75,x<=q25],['high','low'],default='middle'),index=x.index),q25,q75
def desc(a):
 a=np.asarray(a);a=a[np.isfinite(a)]
 return dict(n=len(a),minimum=np.min(a),q25=np.quantile(a,.25),median=np.median(a),q75=np.quantile(a,.75),maximum=np.max(a),mean=a.mean(),sd=np.std(a),zero_fraction=np.mean(a==0)) if len(a) else dict(n=0)
manifest=read(O/'input_manifest.tsv')
extra=[(anchor/f,'R12 original per-center feature axes','shared') for f in ['features.npz','common_anchors.parquet','feature_columns.json']]+[(ringroot/'00_preflight/coverage63.tsv','original ring q05 provenance','shared')]
for s in S:
 d=read(O/f'spots/{s}.base.tsv.gz')
 for part in ['function','markers','Epi','TME']:d=join(d,read(O/f'spots/{s}.{part}.tsv.gz'))
 for model,cols in [('Epi',[H,N]),('TME',[scol(t) for t in fine])]:
  sc=d[cols];err=(sc.sum(axis=1)-1).abs().max()
  assert sc.notna().all().all() and sc.min().min()>=0 and err<1e-7
  qcs.append(dict(section_id=s,model=model,n=len(d),n_labels=len(cols),missing=int(sc.isna().sum().sum()),min_score=sc.min().min(),max_score=sc.max().max(),max_unit_sum_error=err))
  feat=read(O/f'features/{s}.{model}.tsv')
  for (mod,var),g in sets.groupby(['module','variant']):
   hit=sorted(set(g.gene_symbol)&set(feat.gene))
   overlap.append(dict(section_id=s,model=model,module=mod,variant=var,n_mapping_features=len(feat),mapping_Top20=','.join(feat.loc[feat.Top20,'gene']),n_mapping_intersection=len(hit),mapping_intersection=','.join(hit)))
 for major,g in mapping.groupby('major'):d[scol('major_'+major)]=d[[scol(t) for t in g.subtype]].sum(axis=1,min_count=len(g))
 ratios={}
 for t in fine+majors:
  ratios['mean_all30_fraction__'+t]=d[ccol(t)]/d.c2l_mean_all30.where(d.c2l_mean_all30>0)
  ratios['mean_TME29_fraction__'+t]=d[ccol(t)]/d.c2l_mean_TME29.where(d.c2l_mean_TME29>0)
 d=pd.concat([d,pd.DataFrame(ratios,index=d.index)],axis=1)
 e=d[d['query']].copy();assert len(e)==S[s]
 for x in [H,MP,N]+[scol(t) for t in fine]+[lcol(t) for t in fine]+functions:
  for label,g in [('query',e),('full',d)]:ranges.append(dict(section_id=s,population=label,readout=x,**desc(g[x])))
 for x in [H,MP]:
  label,a,b=bands(e[x]);col='region_'+('H05' if x==H else 'MP10');d.loc[e.index,col]=label;e[col]=label
  for y in functions:
   fc.append(dict(section_id=s,exposure=x,module_score=y,**relate(e,x,y)))
   vals={}
   for g in ['high','low']:
    z=e[label==g];vals[g]=z[y]
    groups.append(dict(section_id=s,exposure=x,module_score=y,group=g,Q25=a,Q75=b,n_group=len(z),**desc(z[y])))
   contr.append(dict(section_id=s,exposure=x,module_score=y,n_high=len(vals['high']),n_low=len(vals['low']),high_minus_low_mean=vals['high'].mean()-vals['low'].mean(),high_minus_low_median=vals['high'].median()-vals['low'].median(),coincident_quartiles=a==b))
 e['discordance']=np.select([(e.region_H05=='high')&(e.region_MP10=='low'),(e.region_H05=='low')&(e.region_MP10=='high')],['H05 high / MP10 low','H05 low / MP10 high'],default='other')
 d.loc[e.index,'discordance_new']=e.discordance
 for group_name,z in e.groupby('discordance'):
  rec=dict(section_id=s,group=group_name,n=len(z),n_inside=(z.lesion_location=='inside').sum(),median_x=z.x.median(),median_y=z.y.median())
  for col in functions+['nCount_RNA','Epithelial_fraction','top20_UMI',H,MP]+list(z.filter(regex='^UMI_')):rec['median_'+col]=z[col].median()
  disc.append(rec)
 for y in ['H05_support','H05_conditional_support',MP]:old.append(dict(section_id=s,reference_readout=y,**relate(e,H,y)))
 for x in [H,MP]+functions:
  for c in ['Epithelial_fraction','nCount_RNA']:
   z=np.log1p(e[c]) if c=='nCount_RNA' else e[c]
   diag.append(dict(section_id=s,readout=x,covariate='log1p_nCount' if c=='nCount_RNA' else c,pearson=corr(e[x],z),spearman=corr(rankdata(e[x]),rankdata(z))))
 full.extend(square(e,s,'A',[H,N]+[scol(t) for t in fine]))
 for t in fine+majors:
  lev='major' if t.startswith('major_') else 'subtype'
  on=[H,MP,scol(t),lcol(t),'Epithelial_fraction','nCount_RNA'];z=e.dropna(subset=on)
  for combo,x,y in [('A',H,scol(t)),('B',MP,lcol(t)),('H05_c2l',H,lcol(t)),('MP10_transfer',MP,scol(t))]:
   cross.append(dict(section_id=s,level=lev,cell_type=t,combination=combo,**relate(z,x,y)))
  for method,v in [('mean_absolute',z[ccol(t)]),('mean_all30_fraction',z[ccol(t)]/z.c2l_mean_all30),('mean_TME29_fraction',z[ccol(t)]/z.c2l_mean_TME29)]:
   for x in [H,MP]:sens.append(dict(section_id=s,level=lev,cell_type=t,exposure=x,TME_units=method,**stat(z[x],v,z.Epithelial_fraction,np.log1p(z.nCount_RNA))))
  for x in functions:
   for method,y in [('transfer',scol(t)),('log1p_mean',lcol(t))]:ft.append(dict(section_id=s,level=lev,cell_type=t,module_score=x,TME_method=method,**relate(z,x,y)))
 for label,z in [('full',d),('query',e)]:
  cutoff=z.c2l_mean_TME29.quantile(.25)
  for strata,b in [('bottom_quartile',z[z.c2l_mean_TME29<=cutoff]),('remaining',z[z.c2l_mean_TME29>cutoff])]:
   a=b[[scol(t) for t in fine]];win=a.idxmax(axis=1)
   for t in fine:lowtme.append(dict(section_id=s,population=label,stratum=strata,Q25_absolute_TME=cutoff,n=len(b),cell_type=t,n_top=(win==scol(t)).sum(),mean_score=b[scol(t)].mean(),median_winner_score=a.max(axis=1).median(),Q95_winner_score=a.max(axis=1).quantile(.95)))
 for k,z in e.dropna(subset=['niche']).groupby('niche'):
  for x in [H,MP]+functions:niches.append(dict(section_id=s,niche=int(k),readout=x,**desc(z[x])))
 # Identical assigned centers for the original q05/all30 composition at center, R1 and R2.
 idx=np.where(ordered.section_id==s)[0];u=ordered.iloc[idx].reset_index(drop=True)
 g=e[e.niche.notna()].merge(u.assign(row_id=idx),on=KEY+['spot_uid'],validate='one_to_one')
 assert len(g)==e.niche.notna().sum()
 qpath=Path(region.loc[s,'q05_path']);extra.append((qpath,'q05 absolute used only for matched-scale context',s))
 q=pd.read_csv(qpath,index_col=0).loc[g.barcode,ringcols];q=q.div(q.sum(axis=1),axis=0)
 mats={'center_q05_all30':q.to_numpy(),'R1_q05_all30':archive['q1'][g.row_id],'R2_q05_all30':archive['q2'][g.row_id]}
 rc=[g[KEY+['niche','ring1_valid_n','ring2_valid_n']].reset_index(drop=True)]
 for name,arr in mats.items():rc.append(pd.DataFrame(arr,columns=[name+'__'+c for c in ringcols]))
 save(pd.concat(rc,axis=1),f'spots/{s}.original_rings.tsv.gz')
 for typ in fine:
  i=ringcols.index(typ)
  for label,arr in mats.items():
   for x in [H,MP,'Synthesis__minus20','Homeostasis__minus20','Fatty_acid__minus20']:
    rings.append(dict(section_id=s,cell_type=typ,exposure=x,scale=label,**stat(g[x],arr[:,i],g.Epithelial_fraction,np.log1p(g.nCount_RNA))))
 save(d,f'spots/{s}.tsv.gz')
for filename,rows in [('A_pairwise.tsv.gz',full),('four_way.tsv',cross),('abundance_sensitivity.tsv.gz',sens),('functional_correlations.tsv',fc),('functional_TME.tsv.gz',ft),('functional_regions.tsv',groups),('functional_contrasts.tsv',contr),('discordant_expression.tsv',disc),('niche_scores.tsv',niches),('diagnostic_covariates.tsv',diag),('dynamic_ranges.tsv',ranges),('old_vs_new.tsv',old),('peripheral_vs_center.tsv.gz',rings),('low_TME.tsv',lowtme),('score_QC.tsv',qcs)]:save(rows,'tables/'+filename)
save(cov.merge(pd.DataFrame(overlap),on=['section_id','module','variant'],validate='one_to_many'),'tables/gene_coverage.tsv')
more=[dict(section_id=s,role=r,path=str(p),size=p.stat().st_size) for p,r,s in extra]
save(pd.concat([manifest,pd.DataFrame(more)]).drop_duplicates(['role','section_id','path']),'input_manifest.tsv')

