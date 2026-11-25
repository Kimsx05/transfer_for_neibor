from common import *
from scipy.stats import spearmanr
import itertools
T=Path('tables');D0=Path('cache/jobs');panel=rd(T/'panel_members.tsv').set_index('gene');gm=rd(T/'gene_mapping.tsv').set_index('feature_id');program=rd(T/'program_members.tsv');ctrl=rd(T/'control_matching.tsv')
old=rd(V/'qc/jobs/main_effects.tsv.gz');old['ko']=old.perturbed_gene.map(gm.gene);old['mapped']=old.gene.map(gm.gene)
old['raw_v1_gene']=old.gene;old['gene']=old.mapped;old['is_KO']=(old.ko==old.gene);old['rank_denominator']=1999;old['rank_fraction']=(old['rank']-1)/1998
focus=set(E+K);jobs=[];tops=[];cats=[];prog=[];extended=[];spec=[];compare=[];rankchange=[];index=[];WT={};roles=K
def programstats(a,j,branch):
 join=program.merge(a[~a.is_KO],on='gene',how='inner')
 join['negative']=(join.directionScore<0).where(join.directionScore.notna());join['positive']=(join.directionScore>0).where(join.directionScore.notna())
 by=['ko','program','mask'];stats=join.groupby(by).agg(n_evaluable=('gene','nunique'),negative_fraction=('negative','mean'),positive_fraction=('positive','mean'),median_rank_fraction=('rank_fraction','median'),median_direction=('directionScore','median'),median_distance=('distance','median')).reset_index()
 sizes=program.groupby(['program','mask']).gene.nunique().to_dict()
 rows=[]
 for k in a.ko.unique():
  for name in program.program.unique():
   for mask in ['Full','M20']:rows.append(dict(ko=k,program=name,mask=mask,n_frozen=sizes.get((name,mask),0)))
 z=pd.DataFrame(rows).merge(stats,on=by,how='left');z.n_evaluable=z.n_evaluable.fillna(0).astype(int);z['n_in_panel']=[program[(program.program==n)&(program['mask']==mm)].gene.isin(set(a.gene)).sum() for n,mm in []] if len(z)==0 else z.apply(lambda q:len(pmembers.get((q.program,q['mask']),set())&set_a),axis=1)
 z['status']=np.where((z.n_evaluable>=3)&(z.n_evaluable*2>=z.n_frozen),'EVALUABLE','NE')
 z.loc[z.status=='NE',['negative_fraction','positive_fraction','median_rank_fraction','median_direction','median_distance']]=np.nan
 z['reason']=np.where(z.status=='NE',np.where(z.n_frozen==0,'empty masked program','coverage below gate'),'');z['analysis']=j;z['branch']=branch
 if j in ('clean','full') and branch=='heat':
  save(join[join.ko.isin(K)][['ko','gene','program','mask','distance','directionScore','direction','rank','rank_fraction','p.adj']].assign(analysis=j),T/('program_members_response_'+j+'.tsv.gz'))
 return z
pmembers={k:set(v.gene) for k,v in program.groupby(['program','mask'])}
def rank_compare(a,b,an,bn,k,other):
 aa=a[a.ko.eq(k)&~a.is_KO].drop_duplicates('gene',keep=False).set_index('gene')
 bb=b[b.ko.eq(other)&~b.is_KO].drop_duplicates('gene',keep=False).set_index('gene')
 ii=aa.index.intersection(bb.index).difference([k,other]);x=aa.loc[ii];y=bb.loc[ii]
 for n in [20,50,100]:
  sx=set(x.nsmallest(n,'rank').index);sy=set(y.nsmallest(n,'rank').index)
  compare.append(dict(left_model=an,right_model=bn,ko=k,other=other,n_common=len(ii),spearman=spearmanr(x['rank'],y['rank']).statistic,top_n=n,top_intersection=len(sx&sy),top_union=len(sx|sy),jaccard=len(sx&sy)/len(sx|sy)))
 if k==other:
  zz=pd.DataFrame(dict(gene=ii,left_rank=x['rank'],right_rank=y['rank'],right_minus_left_fraction=y.rank_fraction-x.rank_fraction,left_q=x['p.adj'],right_q=y['p.adj'])).reset_index(drop=True)
  rankchange.append(zz.assign(left_model=an,right_model=bn,ko=k))
for folder in sorted(D0.iterdir()):
 j=folder.name;status=rd(folder/'status.tsv').iloc[0];m=rd(folder/'membership.tsv.gz');jobs.append(dict(analysis=j,status=status.status,reason=status.reason,n_cells=len(m),n_patients=m.patient_id.nunique(),n_datasets=m.dataset.nunique()))
 if status.status!='complete':continue
 for branch in ['heat']+(['manifold'] if j=='clean' else []):
  f=folder/(branch+'.tsv.gz');a=rd(f);set_a=set(a.gene);index.append(dict(method='scTenifoldKnk',analysis=j,branch=branch,path=str(f),rows=len(a),sha256=sha(f)))
  if j in ('clean','full'):WT[(j,branch)]=a
  prog.append(programstats(a,j,branch));a['category']=a.gene.map(panel.category)
  z=a[~a.is_KO & (a['rank']<=100)].copy();tops.append(z)
  for n in [20,50,100]:
   zz=z[z['rank']<=n].groupby(['ko','category']).size().rename('count').reset_index();zz['fraction']=zz['count']/zz.groupby('ko')['count'].transform('sum');cats.append(zz.assign(analysis=j,branch=branch,top_n=n))
  for k in K:
   q=a[a.ko==k].set_index('gene');cc=list(ctrl[ctrl.ko==k].control);vv=a[a.ko.isin(cc)&~a.is_KO].groupby('gene').agg(control_median_rank=('rank_fraction','median'),control_median_distance=('distance','median'),n_controls=('ko','nunique'))
   yy=q.join(vv);yy['rank_advantage_vs_controls']=yy.control_median_rank-yy.rank_fraction;yy['distance_minus_control_median']=yy.distance-yy.control_median_distance
   scores=a[a.ko.isin(cc)].pivot(index='gene',columns='ko',values='distance');yy['n_controls_lower']=scores.lt(yy.distance,axis=0).sum(1).reindex(yy.index)
   yy=yy.reset_index()
   if j=='clean' and branch=='heat':spec.append(yy[~yy.is_KO])
   zz=yy[yy.gene.isin(focus)&~yy.is_KO].copy();extended.append(zz.assign(n_cells=len(m),n_patients=m.patient_id.nunique()))
# historical baseline has older panel denominators and unmapped/ambiguous features.
oldcat=old.gene.map(panel.category).fillna('unclassified')
for n in [20,50,100]:
 u=old[(~old.is_KO)&(old['rank']<=n)].copy();u['category']=oldcat.loc[u.index];v=u.groupby(['ko','category']).size().rename('count').reset_index();v['fraction']=v['count']/v.groupby('ko')['count'].transform('sum');cats.append(v.assign(analysis='v1',branch='heat',top_n=n))
a=WT['clean','heat'];b=WT['full','heat'];ma=WT['clean','manifold']
for k in K:
 rank_compare(b,a,'P_full','P_clean',k,k);rank_compare(a,ma,'clean_heat','clean_manifold',k,k)
 for ct in ctrl[ctrl.ko==k].control:rank_compare(a,a,'clean_heat','matched_control',k,ct)
 if k in K[:2]:
  rank_compare(old,a,'v1','P_clean',k,k);rank_compare(old,b,'v1','P_full',k,k)
for k,y in itertools.combinations(K,2):rank_compare(a,a,'clean_heat','other_main_KO',k,y)
shared=a[~a.is_KO & a.ko.isin(K)].groupby('gene').agg(n_main_top50=('rank',lambda z:(z<=50).sum()),n_main_ranked=('rank','count'),median_main_rank_fraction=('rank_fraction','median')).reset_index()
cc=a[a.ko.isin(set(ctrl.control)) & ~a.is_KO].groupby('gene')['rank'].apply(lambda z:(z<=50).sum()).rename('n_unique_controls_top50')
shared=shared.merge(cc,on='gene',how='left');shared['shared_main']=shared.n_main_top50>=3
save(shared.sort_values(['n_main_top50','median_main_rank_fraction'],ascending=[False,True]),T/'shared_response.tsv')
save(pd.concat(tops),T/'top100.tsv.gz');save(pd.concat(tops).query('rank<=50'),T/'top50.tsv.gz');save(pd.concat(tops).query('rank<=20'),T/'top20.tsv.gz')
save(pd.concat(cats),T/'rank_category.tsv');save(pd.concat(prog),T/'program_response.tsv');save(pd.concat(spec).merge(shared,on='gene'),T/'control_specificity.tsv')
save(pd.DataFrame(compare),T/'ranking_comparisons.tsv');save(pd.concat(rankchange),T/'paired_rank_changes.tsv.gz')
ex=pd.concat(extended);save(ex,T/'endpoint_model_results.tsv.gz')
# Explicitly retain NE requests for absent sources/targets and skipped clusters.
cross=pd.MultiIndex.from_product([[v['analysis'] for v in jobs],K,sorted(focus)],names=['analysis','ko','gene']).to_frame(index=False);cross=cross[cross.ko!=cross.gene]
ex=cross.merge(ex[ex.branch=='heat'],how='left',on=['analysis','ko','gene']);ex['reason']=np.where(ex.status.isna(),'source/endpoint local coverage, no outgoing source or cluster gate','');ex.status=ex.status.fillna('NE')
def kind(j):
 return 'patient_subsample' if j.startswith('sub') else 'LODO' if j.startswith('drop_') else 'seed' if j.startswith('seed') else 'state' if j.startswith('cluster_') else j
ex['refit_type']=ex.analysis.map(kind);save(ex,T/'robustness.tsv')
base=ex[ex.analysis=='clean'].set_index(['ko','gene']);r=[]
for (k,g,ty),z in ex.groupby(['ko','gene','refit_type']):
 ref=base.loc[(k,g)];ok=z[z.status=='EVALUABLE'];signs=ok.directionScore.dropna()
 r.append(dict(ko=k,endpoint=g,type=ty,n_planned=len(z),n_evaluable=len(ok),direction_same_as_clean=(np.sign(signs)==np.sign(ref.directionScore)).mean() if len(signs) else None,down_fraction=(signs<0).mean(),distance_min=ok.distance.min(),distance_max=ok.distance.max(),rank_fraction_min=ok.rank_fraction.min(),rank_fraction_max=ok.rank_fraction.max(),n_q05=(ok['p.adj']<.05).sum()))
save(pd.DataFrame(r),T/'endpoint_stability.tsv')
prog=pd.concat(prog);pr=prog[(prog.branch=='heat')&prog.ko.isin(K)].copy();pr['refit_type']=pr.analysis.map(kind);bas=pr[pr.analysis=='clean'].set_index(['ko','program','mask']);rr=[]
for (k,n,mask,ty),z in pr.groupby(['ko','program','mask','refit_type']):
 a0=bas.loc[(k,n,mask)];ok=z[z.status=='EVALUABLE']
 rr.append(dict(ko=k,program=n,mask=mask,type=ty,n_configs=len(z),n_evaluable=len(ok),fraction_negative_median=ok.negative_fraction.median(),fraction_negative_min=ok.negative_fraction.min(),fraction_negative_max=ok.negative_fraction.max(),median_rank_fraction=ok.median_rank_fraction.median(),rank_min=ok.median_rank_fraction.min(),rank_max=ok.median_rank_fraction.max(),same_median_sign=(np.sign(ok.median_direction)==np.sign(a0.median_direction)).mean() if len(ok) else None))
save(pd.DataFrame(rr),T/'program_stability.tsv')
p=rd(D0/'clean/patient.tsv.gz');save(p,T/'patient_conditioned.tsv.gz')
# Paired contrasts use the intersection of eligible patients for a given KO and readout.
p=p[p.status=='EXPLORATORY'];idx=['patient_id','ko','gene'];lo=p[p.stratum=='L50'][idx+['distance','n_cells']].rename(columns={'distance':'low','n_cells':'low_cells'})
hi=p[p.stratum.isin(['H05','H10','H50'])][idx+['stratum','distance','n_cells']].rename(columns={'distance':'high','n_cells':'high_cells'})
pairs=hi.merge(lo,on=idx);pairs['paired_high_minus_low']=pairs.high-pairs.low;save(pairs,T/'patient_high_low_pairs.tsv.gz')
save(pairs.groupby(['ko','gene','stratum']).paired_high_minus_low.agg(['count','median','min','max']).reset_index(),T/'patient_paired_summary.tsv')
save(pd.DataFrame(jobs),T/'knk_job_status.tsv');save(pd.DataFrame(index),T/'knk_all_gene_index.tsv')

