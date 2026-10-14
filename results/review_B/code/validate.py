from common import *
from association import residualize,correlate
import subprocess,tempfile,argparse,zipfile,re
checks=[]
def ck(name,ok,detail=''):
 checks.append(dict(test=name,status='PASS' if bool(ok) else 'FAIL',detail=detail));assert bool(ok),name

def numeric():
 p=pairs();g=pd.read_csv(OUT/'00_approval/approved_gate_v7.tsv',sep='\t');v6=pd.read_csv(V6/'04_gate/candidate_review_TO_FILL.tsv',sep='\t');ck('approval counts and tiers',len(p)==139 and len(g)==1112 and g.candidate_id.is_unique and p.celltype.nunique()==29 and p.review_tier.value_counts().to_dict()==dict(zip(TIERS,[28,32,39,40])))
 j=g[['candidate_id','scope','niche','celltype','accepted_label']].merge(v6[['candidate_id','scope','niche','celltype','accepted_label']],on='candidate_id',suffixes=('_new','_old'),validate='one_to_one');ck('all 1112 true IDs match v6 keys',len(j)==1112 and all(j[c+'_new'].eq(j[c+'_old']).all() for c in ['scope','niche','celltype','accepted_label']))
 ck('original gate remains unmodified',sha(V6/'04_gate/candidate_review_TO_FILL.tsv')==json.loads((OUT/'00_approval/approval_receipt.json').read_text())['original_v6_gate_sha256'])
 ck('supplemental workbook TSV matched',json.loads((OUT/'00_approval/supplemental_crosscheck.json').read_text())['status']=='PASS')
 a=np.array([[.1,.9],[.2,.8],[.3,.7]]);ck('UID permutation explicit join',np.array_equal(align(['a','b','c'],['b','c','a'],a),a[[1,2,0]]))
 for name,uid,target in [('duplicate UID',['a','a','c'],['a','c']),('wrong UID',['a','b','c'],['x'])]:
  try:align(uid,target,a);ok=False
  except AssertionError:ok=True
  ck(name+' rejected',ok)
 orig=['A','B'];rev=['B','A'];ck('type order restored explicitly',np.array_equal(a[:,::-1][:,pd.Index(rev).get_indexer(orig)],a));ck('wrong column detected',pd.Index(orig).get_indexer(['A','WRONG'])[1]==-1)
 audit=pd.read_csv(OUT/'01_ring123/graph_and_cache_audit.tsv',sep='\t');ck('all63 graph audits',len(audit)==63 and audit.BFS_exact.all() and audit.rings_disjoint.all() and audit.q1_max_error.le(1e-6).all() and audit.q2_max_error.le(1e-6).all())
 m,q,c=load();ck('all ring compositions finite normalized',np.isfinite(q).all() and (q>=0).all() and np.allclose(q.sum(2),1,atol=1e-6));ck('accepted common set identical in all rings',np.array_equal(m.main_eligible,m.accepted_label.ge(0)&m.R1_valid_n.ge(4)&m.R2_valid_n.ge(8)&m.R3_valid_n.ge(12)));ck('closure Epi+TME=1',np.allclose(q[:,:,EPI]+q[:,:,NON].sum(2),1,atol=1e-6));ck('B5 composition retained full20 score excluded',m[m.section_id.eq('GSE238145_B5')].main_eligible.any() and not m[m.section_id.eq('GSE238145_B5')].score_eligible.any())
 # Physical missing nodes persist; stored topology/node axis matches full frozen tissue rows.
 master=pd.read_parquet(V2/'03_features/anchor_master.parquet');ck('missing composition not removed from physical graph',audit.n_tissue_nodes.sum()==len(master) and audit.n_missing_composition.sum()==int((~master.composition_valid).sum()))
 ident=pd.read_csv(OUT/'02_background_composition/background_identity_audit.tsv.gz',sep='\t');ck('background shrinkage identity all eligible cases',ident.max_error.lt(1e-12).all())
 sec=pd.read_csv(OUT/'02_background_composition/section_ring123_composition.tsv.gz',sep='\t');unit=pd.read_csv(OUT/'02_background_composition/unit_ring123_composition.tsv.gz',sep='\t');ck('same section coverage for all three rings',sec.groupby(['analysis_set','section_id','pair_id']).eligible.nunique().eq(1).all());ck('unit aggregation unique across datasets',not unit.duplicated(['analysis_set','scope','summary_type','pair_id','ring','analysis_unit_uid']).any())
 z=sec[sec.analysis_set.eq('MAIN')&sec.eligible].groupby(['pair_id','ring','analysis_unit_uid']).excess.mean();actual=unit[unit.analysis_set.eq('MAIN')&unit.scope.eq('ALL_DATASETS')&unit.summary_type.eq('INFERENTIAL')].set_index(['pair_id','ring','analysis_unit_uid']).excess;ck('overall equal-section mean without dataset splitting',np.allclose(actual.loc[z.index],z,atol=1e-14))
 bg=pd.read_csv(OUT/'03_ring_statistics/ring123_background_tests.tsv',sep='\t');rp=pd.read_csv(OUT/'03_ring_statistics/ring123_pairwise_tests.tsv',sep='\t');mp=pd.read_csv(OUT/'04_MP10_association/MP10_scope_associations.tsv',sep='\t')
 for name,df,keys,n,qcol in [('background',bg,['unit_subset'],417,'background_zero_q'),('ringpair',rp,['unit_subset'],417,'ring_pairwise_q'),('association',mp,['unit_subset','method'],1008,'association_q')]:
  ck(name+' planned families retained',df.groupby(keys+['scope']).size().eq(n).all())
  ck(name+' all scope BH correct',all(np.allclose(bh(d.p,n),d[qcol],equal_nan=True) for _,d in df.groupby(keys+['scope'])))
  ck(name+' auxiliary 8scope BH correct',all(np.allclose(bh(d.p,n*8),d[qcol+'_8scopes'],equal_nan=True) for _,d in df.groupby(keys)))
  ck(name+' auxiliary 7dataset BH correct',all(np.allclose(bh(d.p,n*7),d[qcol+'_7datasets'],equal_nan=True) for _,d in df[df.scope.ne('ALL_DATASETS')].groupby(keys)))
 families=pd.read_csv(OUT/'00_audit/family_manifest.tsv',sep='\t');ck('all80 main and verified family definitions',len(families)==80 and families.groupby(['unit_subset','scope']).size().eq(5).all())
 ck('distinct background and pair q columns','ring_pairwise_q' not in bg and 'background_zero_q' not in rp)
 ck('zero and insufficient statuses',signed([0]*5)['p']==1 and signed([1,0,0,0,0])['status']=='INSUFFICIENT_NONZERO_PAIRS' and signed([1,2])['status']=='INSUFFICIENT_UNITS')
 ck('planned BH not reduced by missing',bh([.001,np.nan],417)[0]==.417 and np.isnan(bh([.001,np.nan],417)[1]))
 # Independent R small-sample and tied/zero numerical parity.
 parity=[]
 cases={'exact_no_ties':[-2,-1,3,4,5,6],'many_ties_zero':[0,0,1,1,-1,2,2,-3,4,4,5,-6],'large_asymptotic':list(range(1,31))+list(range(-1,-31,-1))}
 for name,x in cases.items():
  sx=signed(x);exact=sx['actual_method'].startswith('exact');expr='x<-c('+','.join(map(str,x))+');r<-wilcox.test(x,mu=0,alternative="two.sided",exact='+('TRUE' if exact else 'FALSE')+',correct=TRUE);cat(format(r$p.value,digits=17))';rpv=float(subprocess.check_output(['Rscript','-e',expr],text=True));ck('R signed-rank parity '+name,np.isclose(rpv,sx['p'],atol=1e-12,rtol=1e-10));parity.append(dict(case=name,R_p=rpv,Python_p=sx['p'],method=sx['actual_method']))
 tab('00_audit/R_signed_rank_parity.tsv',parity)
 # Stable ranks/background invariance and independent OLS partial-rank implementation.
 rng=np.random.default_rng(123);X=rng.normal(size=(100,2));Y=rng.normal(size=(100,2));C=rng.normal(size=(100,2));XR=rankdata(X,axis=0);YR=rankdata(Y,axis=0);CR=rankdata(C,axis=0);res,rank,note=residualize(np.column_stack([XR,YR]),CR);D=np.column_stack([np.ones(100),CR]);ols=np.column_stack([XR,YR])-D@np.linalg.lstsq(D,np.column_stack([XR,YR]),rcond=None)[0];ck('QR rank residualization independent OLS parity',np.allclose(correlate(res[:,:2],res[:,2:]),correlate(ols[:,:2],ols[:,2:]),atol=1e-12));ck('subtracting constant background leaves Spearman',np.array_equal(rankdata(X[:,0]),rankdata(X[:,0]-.12345)))
 assoc=pd.read_csv(OUT/'04_MP10_association/MP10_section_associations.tsv.gz',sep='\t');ck('all planned section rows retained',len(assoc)==63*168*3*2*3*5);ck('two score layers share centers and valid sections',assoc.groupby(['analysis_set','section_id','query_id','ring','method']).agg(n=('n_centers','nunique'),valid=('rho',lambda x:x.notna().nunique())).eq(1).all().all());ck('ALL queries deduplicated',assoc[assoc.range.eq('ALL_ELIGIBLE_EPI')].query_id.nunique()==29 and assoc[assoc.range.eq('WITHIN_SELECTED_NICHE')].query_id.nunique()==139)
 au=pd.read_csv(OUT/'04_MP10_association/MP10_unit_associations.tsv.gz',sep='\t');ck('association units globally unique',not au.duplicated(['analysis_set','scope','query_id','ring','method','score_layer','analysis_unit_uid']).any())
 # Recompute actual subset (one fixed query per dataset) from raw frozen arrays.
 records=[]
 for ds in SCOPES[1:]:
  ix=m.index[m.dataset_id.eq(ds)&m.score_eligible].to_numpy()
  if len(ix)==0:continue
  secid=m.loc[ix].groupby('section_id').size().idxmax();ix=m.index[m.section_id.eq(secid)&m.score_eligible].to_numpy();cell=NON[0];X=q[ix,0,cell,None];Y=m.loc[ix,LAYERS].to_numpy();C=np.column_stack([c[ix,EPI],np.log1p(m.loc[ix,'center_nCount'])]);D=np.column_stack([np.ones(len(ix)),rankdata(C,axis=0)]);Z=np.column_stack([rankdata(X,axis=0),rankdata(Y,axis=0)]);R=Z-D@np.linalg.lstsq(D,Z,rcond=None)[0];effect=correlate(R[:,:1],R[:,1:])[0]
  actual=assoc[assoc.analysis_set.eq('MAIN')&assoc.section_id.eq(secid)&assoc.query_id.eq(f'all_c{cell:02d}')&assoc.ring.eq('R1')&assoc.method.eq('ADJUSTED_MAIN')].set_index('score_layer').loc[LAYERS,'rho'].to_numpy();ck('actual rank association '+ds,np.allclose(effect,actual,atol=1e-10,equal_nan=True));records.append(dict(dataset_id=ds,section_id=secid,n=len(ix),max_error=np.nanmax(abs(effect-actual))))
 tab('00_audit/actual_association_recalculation.tsv',records)
 block=pd.read_csv(OUT/'05_sensitivity/spatial_block_selected_anchors.tsv',sep='\t');ck('shared block sample max1 per physical block',not block.duplicated(['section_id','block_row','block_col']).any() and set(block.spot_uid).issubset(set(m.loc[m.score_eligible,'spot_uid'])))
 sm=masks(m,q);ck('TME fixed centers and ring-level normalization',np.array_equal(sm['TME_MATCHED_RAW'],sm['TME_CONDITIONAL']) and (q[sm['TME_CONDITIONAL']][:,:,NON].sum(2)>=.05).all());ck('TME normalized nonEpi sums1',np.allclose((q/q[:,:,NON].sum(2,keepdims=True))[:,:,NON].sum(2),1))
 man=pd.read_csv(OUT/'00_audit/input_manifest.tsv',sep='\t');ck('all inputs unchanged',all(sha(r.path)==r.sha256 for r in man.itertuples()));supp=pd.read_csv(OUT/'00_approval/supplemental_input_manifest.tsv',sep='\t');ck('supplemental sources unchanged',all(sha(r.path)==r.sha256 for r in supp.itertuples()))
 tab('reports/validation_numeric.tsv',checks);log(f'{len(checks)} numeric validation checks passed')

def delivery():
 checks.clear();maps=[]
 for scope in SCOPES:
  d=pd.read_csv(OUT/'06_figures'/scope/'display_name_map.tsv',sep='\t');ck(scope+' complete139 composition and168 MP10 figures',d.kind.eq('composition').sum()==139 and d.kind.eq('MP10').sum()==168 and not d.duplicated(['kind','query_id']).any());ck(scope+' all PNG/PDF exist',all((OUT/p).is_file() for p in d.png) and all((OUT/p).is_file() for p in d.pdf));maps.append(d)
  for p in d.png:
   from PIL import Image
   with Image.open(OUT/p) as im:assert im.info.get('dpi',(0,0))[0]>=299
  ck(scope+' PNG >=300dpi',True)
 mapping=pd.concat(maps);tab('06_figures/display_name_map.tsv',mapping)
 # Figure sample N matches inferential unit statistics (all files are generated from these exact sets).
 u=pd.read_csv(OUT/'02_background_composition/unit_ring123_composition.tsv.gz',sep='\t');u=u[u.analysis_set.eq('MAIN')&u.summary_type.eq('INFERENTIAL')];counts=u.groupby(['scope','pair_id','ring']).size();t=pd.read_csv(OUT/'03_ring_statistics/ring123_background_tests.tsv',sep='\t');t=t[t.unit_subset.eq('ALL_POLICY_UNITS')];ck('all composition plot Ns match tests',all(counts.get((r.scope,r.pair_id,r.ring),0)==r.n_units for r in t.itertuples()))
 ck('visual inspection four scopes',len(pd.read_csv(OUT/'reports/visual_inspection.tsv',sep='\t'))>=4)
 for p in (OUT/'reports').glob('*.zip'):
  with zipfile.ZipFile(p) as z:
   ck('zip CRC '+p.name,z.testzip() is None)
   with tempfile.TemporaryDirectory(prefix='v7_',dir='/tmp') as td:
    z.extractall(td);root=Path(td);missing=[]
    for file in list(root.rglob('*.md'))+list(root.rglob('*.html')):
     if '00_approval' in str(file):continue
     txt=file.read_text();links=re.findall(r'\]\(([^)]+)\)',txt)+re.findall(r'(?:href|src)="([^"]+)"',txt)
     for link in links:
      if link.startswith(('http','mailto','#')):continue
      target=(file.parent/link.split('#')[0]).resolve()
      if not target.exists():missing.append(str(file.relative_to(root))+':'+link)
    ck('portable extracted links '+p.name,not missing,';'.join(missing[:5]))
 tab('reports/validation_delivery.tsv',checks);num=pd.read_csv(OUT/'reports/validation_numeric.tsv',sep='\t');tab('reports/validation_results.tsv',pd.concat([num,pd.DataFrame(checks)],ignore_index=True))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--delivery',action='store_true');a=ap.parse_args();delivery() if a.delivery else numeric()
