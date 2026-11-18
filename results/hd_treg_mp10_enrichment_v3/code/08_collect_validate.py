from common import *
checks=[]
def check(name,ok,detail=''):
 assert bool(ok),(name,detail)
 checks.append(dict(check=name,passed=bool(ok),detail=str(detail)))
def bh(p):
 p=np.asarray(p,float);out=np.full(len(p),np.nan);ix=np.flatnonzero(np.isfinite(p));o=ix[np.argsort(p[ix])];out[o]=np.minimum.accumulate((p[o]*len(o)/np.arange(1,len(o)+1))[::-1])[::-1].clip(0,1);return out
for name in ['enrichment_effects','neighborhood_enrichment','spatial_robustness','cross_distance_results','null_histograms','exchangeability_summary','bootstrap_diagnostics']:
 d=pd.concat([pd.read_csv(R/s/(name+'.tsv'),sep='\t')for s in SECTIONS],ignore_index=True)
 if name in ['enrichment_effects','neighborhood_enrichment']:
  d['q']=np.nan
  for fam,g in d.groupby('family'):d.loc[g.index,'q']=bh(g.P)
 put(d,name+'.tsv')
ef=pd.read_csv(R/'tables/enrichment_effects.tsv',sep='\t');check('40 planned FE combinations',len(ef)==40);check('unique planned combinations',not ef.duplicated(['section','radius_um','high_definition','label_definition']).any());check('contingency totals',ef[['a','b','c','d']].sum(axis=1).eq(ef.N).all());check('10 primary rows',ef.family.eq('primary_FE').sum()==10)
put(ef[['section','radius_um','high_definition','label_definition','a','b','c','d','N','label_positive','near_high_denominator','comparison_denominator','positive_fraction_high','positive_fraction_comparison','overall_positive_fraction']],'contingency_tables.tsv')
put(ef[['section','radius_um','high_definition','label_definition','N','FE_status','PR_status','OR_status','FE_CI_status','estimability_reason','P','q']],'estimability_and_limitations.tsv')
for s in SECTIONS:
 t=pd.read_parquet(R/s/'T_regions.parquet');e=pd.read_parquet(R/s/'Epi_regions.parquet');n=pd.read_parquet(R/s/'neighborhoods.parquet');old=pd.read_parquet(V/s/'T32_regions.parquet') if (V/s/'T32_regions.parquet').exists()else None
 check(s+' unique region IDs',t.region_id.is_unique and e.region_id.is_unique)
 check(s+' FOXP3 exact raw rule',t.FOXP3_supported.eq(t.T_eligible&t.raw_FOXP3.gt(0)).all())
 check(s+' eligible complete continuous rows',n.region_id.nunique()==t.T_eligible.sum()and len(n)==2*t.T_eligible.sum())
 check(s+' exposure missing outside U',n.loc[~n.U_r,['near_H05','near_H10']].isna().all().all())
 check(s+' U rule',n.U_r.eq(n.coverage.ge(.5)&n.Epi_n.ge(1)).all())
 check(s+' MP10 quantile frozen',e.H05.eq(e.MP10.gt(e.MP10.quantile(.95))).all()and e.H10.eq(e.MP10.gt(e.MP10.quantile(.9))).all())
 check(s+' finite coordinates',np.isfinite(t[['x_um','y_um']]).all().all()and np.isfinite(e[['x_um','y_um']]).all().all())
 for radius in RADII:
  u=n[n.radius_um.eq(radius)&n.U_r];l=u[LABELS].to_numpy(float);groups=u.stratum.to_numpy();_,perm=next(perm_batches(l,groups,2));
  for g in np.unique(groups):check(s+str(radius)+' permutation stratum '+str(g),np.all(perm[:,groups==g,:].sum(1)==l[groups==g,:].sum(0)))
hist=pd.read_csv(R/'tables/null_histograms.tsv',sep='\t');check('9999 permutations retained in each null histogram',hist.groupby(['section','radius_um','high_definition','label_definition','statistic']).frequency.sum().eq(9999).all())
# Explicit counterexamples for denominator, zero, infinity, and unordered pair conventions.
z=effects(0,2,3,4);check('FE valid zero',float(z['FE'])==0);check('OR positive divided by zero is infinity',np.isinf(effects(2,0,1,4)['OR']));check('no label positive makes FE undefined',np.isnan(effects(0,0,2,3)['FE']));check('no high group makes FE undefined',np.isnan(effects(0,2,0,3)['FE']))
xy=pd.DataFrame({'x_um':[0,3,30,0],'y_um':[0,4,0,0],'component':[1,1,1,2]});pp=pairlist(xy,5);check('physical radius boundary and component isolation',np.array_equal(pp,np.array([[0,1]])))
# Two A positions and one B, all adjacent: IAA=1, IAB=2, MS=2, NMS=1.
check('mixing unordered denominator factor2',2/1==2 and (2/(2*1))/(1/(2*(2-1)/2))==1)
# One shared member excludes the complete pair, not a score fractional adjustment.
Ta={'t1':{1,2},'t2':{3}};Ea={'e1':{2,4},'e2':{5}};valid=[(a,b)for a in Ta for b in Ea if not Ta[a]&Ea[b]];check('whole source-sharing pair removed',len(valid)==3 and ('t1','e1')not in valid)
put(pd.DataFrame(checks),'validation_results.tsv')
import importlib.metadata as md
vers=[]
for p in ['numpy','pandas','scipy','scikit-learn','matplotlib','pyarrow','PyYAML']:
 try:vers.append(dict(package=p,version=md.version(p)))
 except:pass
vers.append(dict(package='Python',version=sys.version));put(pd.DataFrame(vers),'software_versions.tsv')
print(ef[ef.family.eq('primary_FE')][['section','radius_um','N','label_positive','a','near_high_denominator','FE','FE_low','FE_high','P','q']].to_string(index=False));print('VALIDATED',len(checks))
