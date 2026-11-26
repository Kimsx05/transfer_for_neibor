from common import *
def review(s):
 hd=s in HD;v,d=folders(s);z=labels(s);x,g,ids=expression(s)
 assert np.array_equal(z.index,ids)
 q=paths(s);q.mkdir(parents=True,exist_ok=True)
 genes=list(dict.fromkeys(PAN+CORE+'CD4 CD8A CD8B IL7R GZMK TOX LAG3 NKG7 GNLY KLRD1 EPCAM KRT8 KRT18 LST1 TYROBP'.split()+MAJ['T/NK']))
 a=raw(x,g,ids,genes);pan=(a[PAN]>0).sum(1);comp=(a[CORE[1:]]>0).sum(1);fox=a.FOXP3.gt(0)
 hit1=z.candidate.eq('Treg');hit2=fox;hit3=pan.ge(1)&comp.ge(2)
 e=z[['x_um','y_um','major','major_candidate','major_support','major_second','major_runner','supported','flag','candidate','support','second','runner','nFeature']].copy()
 e['UMI']=z['nCount' if hd else 'transcript_counts'];e['panT_genes']=pan;e['companion_genes']=comp;e['fit_or_transfer']=z.major_candidate.notna()
 if hd:
  w=pd.read_parquet(v/'weights.parquet').reindex(ids);t=w[TAX.index[TAX.major.eq('T/NK')]].sum(1,min_count=1);tw=w.Treg
  rank=w.gt(tw,axis=0).sum(1).add(1).where(tw.notna());e['treg_rank']=rank;e['treg_weight']=tw;e['tnk_weight']=t;e['treg_within_tnk']=tw/t.where(t>0)
  hit4=(tw>=.1)&((rank<=3)|(e.treg_within_tnk>=.30))&((pan+comp+fox)>0)
  ref=hit4|(hit1&(tw>=.2));tref=(t-w.NK)>=.20
  good=pd.Series(True,index=z.index); e['epi_weight']=w.Epithelial
 else:
  boundary=pd.read_parquet(v/'boundary_geometry.parquet').reindex(ids)
  good=~(z.area_outlier|z.boundary_degenerate|boundary.review_boundary)
  e=e.join(z[['cell_area','area_outlier','boundary_degenerate','segmentation_method']]).join(boundary[['review_boundary','polygon_valid']]);hit4=z.second.eq('Treg')
  ref=(hit1&(z.support>=.5))|(hit4&(z.runner>=.25));tref=z.major_candidate.eq('T/NK')
  good &= z.major.isin(['T/NK','Unknown','Mixed'])
 resolved=z.supported.isin(TAX.index[(TAX.major=='T/NK') & (TAX.index!='NK')]);reliable=resolved|(pan>=2)|((pan>=1)&tref)
 if not hd:reliable &= good&~z.supported.eq('NK')
 else:reliable &= ~((z.supported=='NK')&(pan<2))
 cand=hit1|hit2|hit3|hit4
 strong=reliable & ((fox&((comp>=1)|ref))|((comp>=2)&ref)) & good
 possible=((pan>=1)|resolved) & (fox|(comp>=2)|(ref&(comp>=1))) & good & reliable
 status=np.select([~cand,strong,possible,((pan+comp)>0)&((pan>0)|(comp>=2)|ref)],['outside','supported','possible','unresolved'],default='unsupported')
 e['reference_support']=ref;e['reliable_T']=reliable;e['candidate_union']=cand;e['grade']=status
 for k,val in zip(['A1_model','A2_FOXP3','A3_combination','A4_secondary'],[hit1,hit2,hit3,hit4]):e[k]=val
 e['reason']=np.select([~cand,strong&fox,strong,possible,~good,(pan==0)&(comp==0),pan==0,~ref],['outside recall union','T lineage + FOXP3 with companion or reference','T lineage + multiple companions with reference','partial Treg evidence in reliable T','segmentation or alternative major','isolated or absent raw marker evidence','limited T lineage attribution','activation or competing identity'],default='incomplete evidence')
 fails={'unfitted_or_no_fine_transfer':z.candidate.isna(),'parent_not_TNK':z.major.ne('T/NK'),'mixed':z.major.eq('Mixed'),'Treg_not_top1':~hit1,'major_margin':(z.major_support-z.major_runner<.12),'major_below_0.50':z.major_support.lt(.5),'no_measured_parent_marker':(a[MAJ['T/NK']]>0).sum(1).eq(0),'subtype_weight_or_score':z.support.lt(.4 if hd else .5),'subtype_margin':(z.support-z.runner)<(.12 if hd else .15),'missing_anchor':~(fox|(a.IL2RA>0)),'old_combination':~(fox|(comp>=2))}
 for key,val in fails.items():e['old_'+key]=val
 e['old_gate_failures']=['|'.join(np.array(list(fails))[r]) for r in np.array(list(fails.values())).T]
 e['human_review_status']='pending';e['review_origin']='provisional';e['reviewer']='Codex';e['review_mode']='table_only';e['review_status']=np.where(cand,'Codex-reviewed','outside')
 e['evidence_path']='review/'+s+'/candidates.parquet';e=e.join(a)
 assert not ((e.grade.isin(['supported','possible']))&~e.reliable_T).any()
 e.index.name='id';e.to_parquet(d/'identity.parquet')
 ce=e.loc[cand];ce=ce.join(w.add_prefix('weight_')) if hd else ce;ce.to_parquet(q/'candidates.parquet')
 strata=pd.qcut(e.UMI.rank(method='average'),4,labels=False,duplicates='drop');selected=[]
 rng=np.random.default_rng(CFG['seed'])
 for name,mask,n in [('supported',e.grade.eq('supported'),8 if hd else 4),('disputed',e.grade.isin(['possible','unresolved']),10 if hd else 4),('unsupported',e.grade.eq('unsupported'),6 if hd else 2),('T-lineage control',~cand&(pan>0),6 if hd else 2)]:
  pool=e[mask].copy();pool['stratum']=strata[mask].astype(str)+'/'+pool.flag.astype(str)+'/'+pool.A1_model.astype(str)+'/'+pool.A4_secondary.astype(str)
  if name=='T-lineage control' and len(pool):
   targets=strata[e.index.isin([row[0] for row in selected])].dropna().to_numpy()
   for t0 in rng.choice(targets if len(targets) else [0,1,2,3],min(n,len(pool))):
    ix=(strata.loc[pool.index]-t0).abs();i=rng.choice(ix.index[ix==ix.min()]);selected.append((i,name));pool=pool.drop(i)
  else:
   pool['u']=rng.random(len(pool));pool=pool.sort_values('u');pool['cycle']=pool.groupby('stratum').cumcount()
   selected.extend((i,name) for i in pool.sort_values(['cycle','u']).head(n).index)
 b=pd.DataFrame(selected,columns=['id','sampling_stratum']).set_index('id').join(e[['grade','UMI','panT_genes','flag','reason']]);b['human_grade']='';b['human_reliable_T']='';b['human_epi_eligible']='';b['human_reviewer']='';b['comment']='';b['human_status']='pending';
 if (q/'human_review.tsv').exists():
  previous=pd.read_csv(q/'human_review.tsv',sep='\t');assert np.array_equal(previous.id,b.index),'review sampling changed; use a new output revision'
 else:put(b.reset_index(),q/'human_review.tsv')
 out=e.groupby(['grade','major','A1_model','A2_FOXP3','A3_combination','A4_secondary']).size().rename('n').reset_index();put(out,q/'coverage.tsv')
 print(s,ce.grade.value_counts().to_dict(),len(b))
if __name__=='__main__':
 for s in sys.argv[1:] or HD+XE:review(s)
